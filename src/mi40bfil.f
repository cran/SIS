************************************************************************
*
*     File  mi40bfil fortran.
*
*     m4getb   m4chek   m4id     m4name   m4inst   m4load   m4oldb
*     m4savb   m4dump   m4newb   m4pnch
*     m4rept   m4soln   m4solp   m4stat
*
*     09 Oct 2000: m2xmat called if Report file 91, 92, 93 specified.
*     30 Jan 2004: character*4 xxx  changed to character xxx*4 (etc).
*     14 Feb 2004: misolf implemented and called by m4soln.
*     09 Mar 2004: m4savb now prints biggest x and pi (and where),
*                  because max( pinorm, 1 ) was always annoying!!
*     17 Jun 2004: m4savb saves biggest primal and dual infeasibilities
*                  before and after scaling, so GAMS can use them.
*                  dparm(27:30) saves  Binf1,  Dinf1,  Binf,  Dinf.
*                  iparm(27:30) saves jBinf1, jDinf1, jBinf, jDinf.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4getb( ncycle, istart, m, mbs, n, nb, nn, nname, nscl,
     $                   lcrash, ns,
     $                   ne, nka, a, ha, ka,
     $                   hrtype, hs, kb, ascale, bl, bu,
     $                   pi, xn, y, y2, name1, name2, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            ha(ne), hrtype(mbs), hs(nb)
      integer            ka(nka), kb(mbs), name1(nname), name2(nname)
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb)
      double precision   pi(m), xn(nb), y(m), y2(m), z(nwcore)

*     ------------------------------------------------------------------
*     m4getb is called with ncycle = 0 before the first cycle.
*     A basis file is loaded (if any exists).
*
*     m4getb is called with ncycle > 0 before every cycle.
*     The Jacobian is evaluated and stored in a (unscaled).
*     If gotscl is false, m2scal is called to compute scales.
*     If relevant, the scales are applied to a, bl, bu, xn, pi, fcon.
*     If gotbas is false, m2crsh is called to initialize hs.
*     (kb, y, y2  are used as workspace by m2crsh.)
*
*     In both cases, lcrash is an output parameter to tell m5solv
*     if further calls to crash are needed.
*
*     14 May 1992: pi(1:nncon) assumed to be initialized on entry.
*                  xn passed to m2crsh.
*     04 Jun 1992: lcrash added as an output parameter.
*     09 Jul 1992: istart added as an input parameter.
*                  Used when ncycle = 0 to load a basis file only if
*                  it is a cold start.
*     10 Oct 1993: If nncon > 0, xlam now has length m instead of nncon,
*                  so that m8srch can take a short step in all of pi.
*                  We set pi(nncon+1:m) = 0 and then xlam = pi.
*     08 Oct 1994: pi(iobj) always set to -1.0.  This may be important
*                  if the first major itn takes a step less than 1.
*     18 Nov 1995: pi(iobj) = - minimz  is more like it (for max probs).
*     28 Nov 1996: Jacobian is now initially random.  No longer need
*                  to set it (via old m8ajac).
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m2parm/ dparm(30),iparm(30)
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /m8loc / lfcon ,lfcon2,lfdif ,lfdif2,lfold ,
     $                   lblslk,lbuslk,lxlam ,lrhs  ,
     $                   lgcon ,lgcon2,lxdif ,lxold
      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)
      common    /m8func/ nfcon(4),nfobj(4),nprob,nstat1,nstat2
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave
      logical            gotbas,gotfac,gothes,gotscl
      common    /cycle1/ gotbas,gotfac,gothes,gotscl

      equivalence      ( iparm(1), icrash )
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      lcrash = 0

      if (ncycle .eq. 0) then
*        ---------------------------------------------------------------
*        This call is made before the cycle loop.
*        ---------------------------------------------------------------

*        Initialize  lvldif = 1, nfcon(*) = 0, nfobj(*) = 0.

         lvldif = 1
         call iload1( 4, 0, nfcon, 1 )
         call iload1( 4, 0, nfobj, 1 )

         if (istart .eq. 0) then
*           -----------
*           Cold start.
*           -----------
*           We have to initialize xn(n+1:nb) and pi(nncon+1:m)
*           before the problem is scaled.
*           The basis files initialize all of xn.
*           One day they may load pi for nonlinear problems.

            call dload ( m      , zero, xn(n+1)    , 1 )
            if (nncon .lt. m)
     $      call dload ( m-nncon, zero, pi(nncon+1), 1 )

*           Load a basis file if one exists.

            if      (ioldb .gt. 0) then
               call m4oldb( m, n, nb, ns, hs, bl, bu, xn )
            else if (insrt .gt. 0) then
               call m4inst( m, n, nb, nname, ns,
     &                      hs, bl, bu, xn, name1, name2 )
            else if (iload .gt. 0) then
               call m4load( m, n, nb, nname, ns,
     &                      hs, bl, bu, xn, name1, name2 )
            end if
            nssave = ns
         end if

*        Make sure the nonlinear variables are within bounds.

         do 160 j = 1, nn
            xn(j) = max( xn(j), bl(j) )
            xn(j) = min( xn(j), bu(j) )
  160    continue

      else
*        ---------------------------------------------------------------
*        ncycle > 0.  This call is made every cycle.
*        ---------------------------------------------------------------

*        If nncon > 0, initialize xlam from pi.
*        Do it at the very start, or on later cycles if the present
*        solution is feasible.  (ninf = 0 in both cases.
*        Keep the previous xlam if solution is infeasible.)
*        First change the sign of pi if maximizing.
*
*        If pi seems ridiculously big, assume that it
*        has not been correctly initialized and just use xlam = 0.
                          
         if (nncon .gt. 0  .and.  ninf .eq. 0) then
            if (minimz .lt. 0) then
               call dscal ( m, -one, pi, 1 )
            end if

            toobig = 1.0d+10
            imax   = idamax( m, pi, 1 )
            xlmax  = abs( pi(imax) )
            if (xlmax .lt. toobig) then
               call dcopy ( m, pi, 1, z(lxlam), 1 )
            else
               call dload ( m,  zero, z(lxlam), 1 )
               if (iprint .gt. 0) write(iprint, 1000)
               if (isumm  .gt. 0) write(isumm , 1000)
            end if
         end if

*        08 Oct 1994: Make sure pi(iobj) and xlam(iobj) are correct.
*        18 Nov 1995: Must allow for maximization probs (minimz = -1).

         if (iobj .gt. 0) then
            pi(iobj) = - minimz
            if (nncon .gt. 0) z(lxlam + iobj - 1) = pi(iobj)
         end if

*        Compute scales from  a, bl, bu  (unless we already have them).
*        Then apply them to   a, bl, bu, pi, xn.

         if (lscale .gt. 0) then
            if (.not. gotscl) then
               gotscl = .true.
               call m2scal( m, n, nb, ne, nka, nn, nncon, nnjac,
     $                      hrtype, ha, ka, a, ascale, bl, bu, y, y2 )
            end if

            call m2scla( 1, m, n, nb, ne, nka,
     $                   ha, ka, a, ascale, bl, bu, pi, xn )
         end if

*        ---------------------------------------------------------------
*        If there was no basis file, find an initial basis via Crash.
*        This works best with A scaled (hence much of the complication).
*        xn(1:n) is input.  xn(n+1:nb) is initialized by m2crsh.
*        ---------------------------------------------------------------
         if (.not. gotbas) then
            gotbas = .true.
            if (icrash .gt. 0  .and.  icrash .le. 3) lcrash = icrash

            call m2crsh( lcrash, m, n, nb, nn,
     $                   ne, nka, a, ha, ka,
     $                   kb, hs, hrtype, bl, bu, xn, z, nwcore )
         end if
      end if

*     Exit.

  900 return

 1000 format(' XXX Warning: pi(*) is big (perhaps not initialized).'
     $     / ' XXX Setting  pi(*) = 0  for use as initial lambda.')

      end ! subroutine m4getb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4chek( m, maxs, mbs, n, nb, ns,
     $                   hs, kb, bl, bu, xn )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            kb(mbs)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m4chek  takes hs and xn and checks they contain reasonable values.
*     The entries hs(j) = 2 are used to set ns and nssave and possibly
*     the list of superbasic variables kb(m+1) thru kb(m+ns).
*     Scaling, if any, has taken place by this stage.
*
*     If gotbas and gothes are both true, nssave and the superbasic kb's
*     are assumed to be set.  It must be a Hot start, or ncycle > 1.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave
      logical            gotbas,gotfac,gothes,gotscl
      common    /cycle1/ gotbas,gotfac,gothes,gotscl

      logical            setkb
      parameter        ( zero = 0.0d+0,  tolb = 1.0d-4 )

*     Make sure hs(j) = 0, 1, 2 or 3 only.

      do 5 j = 1, nb
         js  = hs(j)
         if (js .lt. 0) hs(j) = 0
         if (js .ge. 4) hs(j) = js - 4
    5 continue

      setkb  = .not. (gotbas .and. gothes)

*     ------------------------------------------------------------------
*     Make sure the objective is basic and free.
*     Then count the basics and superbasics, making sure they don't
*     exceed m and maxs respectively.  Also, set ns and possibly
*     kb(m+1) thru kb(m+ns) to define the list of superbasics.
*     Mar 1988: Loop 100 now goes backwards to make sure we grab obj.
*     Apr 1992: Backwards seems a bit silly in the documentation.
*               We now go forward through the slacks,
*               then forward through the columns.
*     ------------------------------------------------------------------
   10 nbasic = 0
      ns     = 0
      if (iobj .gt. 0) then
          jobj     =   n + iobj
          hs(jobj) =   3
          bl(jobj) = - plinfy
          bu(jobj) =   plinfy
      end if

*     If too many basics or superbasics, make them nonbasic.
*     Do slacks first to make sure we grab the objective slack.

      j = n

      do 100 jj = 1, nb
         j      = j + 1
         if (j .gt. nb) j = 1
         js     = hs(j)
         if (js .eq. 2) then
            ns     = ns + 1
            if (ns .le. maxs) then
               if ( setkb ) kb(m + ns) = j
            else
               hs(j) = 0
            end if

         else if (js .eq. 3) then
            nbasic = nbasic + 1
            if (nbasic .gt. m) hs(j) = 0
         end if
  100 continue

*     Proceed if the superbasic kbs were reset, or if ns seems to
*     agree with nssave from the previous cycle.
*     Otherwise, give up trying to save the projected Hessian, and
*     reset the superbasic kbs after all.

      if (setkb) then
*        ok
      else if (ns .ne. nssave) then
         setkb  = .true.
         gothes = .false.
         if (iprint .gt. 0) write(iprint, 1000) ns, nssave
         if (isumm  .gt. 0) write(isumm , 1000) ns, nssave
         go to 10
      end if

*     Check the number of basics.

      ns     = min( ns, maxs )
      nssave = ns
      if (nbasic .ne. m ) then
         gothes = .false.
         if (iprint .gt. 0) write(iprint, 1100) nbasic, m
         if (isumm  .gt. 0) write(isumm , 1100) nbasic, m
      end if

*     -----------------------------------------------------------
*     On all cycles, set each nonbasic xn(j) to be exactly on its
*     nearest bound if it is within tolb of that bound.
*     -----------------------------------------------------------
      bplus  = 0.9d+0*plinfy
      do 300 j = 1, nb
         xj    = xn(j)
         if (abs( xj ) .ge.  bplus) xj = zero
         if (hs(j)     .le.  1    ) then
            b1    = bl(j)
            b2    = bu(j)
            xj    = max( xj, b1 )
            xj    = min( xj, b2 )
            if ((xj - b1) .gt. (b2 - xj)) b1 = b2
            if (abs(xj - b1)  .le.  tolb) xj = b1
            hs(j) = 0
            if (xj .gt. bl(j)) hs(j) = 1
         end if
         xn(j) = xj
  300 continue

      return

 1000 format(/ ' WARNING:', i6, ' superbasics in hs(*);',
     $         ' previously ns =', i6, '.  Hessian not saved')
 1100 format(/ ' WARNING:', i7, ' basics specified;',
     $         ' preferably should have been', i7)

      end ! subroutine m4chek

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )

      integer            name1(nname), name2(nname)

*     ------------------------------------------------------------------
*     m4id   returns a name id1-id2 for the j-th variable.
*     If nname = nb, the name is already in name1, name2.
*     Otherwise nname = 1. Some generic column or row name is cooked up.
*
*     26 Oct 1996: Use . instead of blank in generic names.
*     ------------------------------------------------------------------

      character          gname*8

      if (nname .gt. 1) then
         id1 = name1(j)
         id2 = name2(j)
      else
         if (j .le. n) then
            write(gname, '(a1,i7)') 'x', j
         else
            write(gname, '(a1,i7)') 'r', j - n
         end if

*        This gives names like 'x   1234' or 'r     99'.
*        Change them to be     'x...1234' or 'r.....99'.

         do 10 k = 2, 7
            if (gname(k:k) .eq. ' ') then
                gname(k:k)   =  '.'
            else
                go to 20
            end if
   10    continue
   20    read (gname, '(2a4)') id1, id2
      end if

      end ! subroutine m4id

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4name( m, n, nb, nname, name1, name2, id1, id2,
     $                   ncard, notfnd, maxmsg, j1, j2, jmark, jfound )

      implicit           double precision (a-h,o-z)
      integer            name1(nname), name2(nname)

*     ------------------------------------------------------------------
*     m4name searches for id1-id2 in arrays name1-2(j), j = j1, j2.
*     jmark  will probably speed the search on the next entry.
*     Used by subroutines m3mpsc, m4inst, m4load.
*
*     Left-justified alphanumeric data is being tested for a match.
*     On Burroughs B6700-type systems, one could replace .eq. by .is.
*
*     Parameter    In/Out  Meaning
*     ----------------------------
*     m            I       Number of rows.
*     n            I       Number of columns.
*     nb           I       Length of list of names to be searched.
*                          Might be n or n+m.
*     nname        I       Length of name1(*), name2(*).
*                          nname < nb means there are no real names.
*     name1, name2 I       List of names (left and right halves).
*     id1  , id2   I       Name being searched for.
*     ncard        I       A number to identify this search.
*     notfnd       I O     Counts number of failed searches.
*     maxmsg       I       Limits number of error messages.
*     j1           I       Marks the first name to be searched.
*     j2           I       Marks the last  name to be searched.
*     jmark        I O     Says where to start the search.
*                          Changed to jfound if search is successful.
*                          Changed to j1     otherwise.
*     jfound         O     Position of name (id1, id2) in list of names.
*     29 Apr 2001: Allow for missing names.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm

      do 50 j = jmark, j2
         call m4id  ( j, m, n, nb, nname, name1, name2, id1j, id2j )
         if (id1 .eq. id1j  .and.  id2 .eq. id2j) go to 100
   50 continue

      do 60 j = j1, jmark
         call m4id  ( j, m, n, nb, nname, name1, name2, id1j, id2j )
         if (id1 .eq. id1j  .and.  id2 .eq. id2j) go to 100
   60 continue

*     Not found.

      jfound = 0
      jmark  = j1
      notfnd = notfnd + 1
      if (notfnd .le. maxmsg) then
         if (iprint .gt. 0) write(iprint, 1000) ncard, id1, id2
      end if
      return

*     Got it.

  100 jfound = j
      jmark  = j
      return

 1000 format(' XXX  Line', i6, '  --  name not found:', 8x, 2a4)

      end ! subroutine m4name

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4inst( m, n, nb, nname, ns,
     $                   hs, bl, bu, xn, name1, name2 )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            name1(nname), name2(nname)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     This impression of INSERT reads a file produced by  m4pnch.
*     It is intended to read files similar to those produced by
*     standard MPS systems.  It recognizes SB as an additional key.
*     Also, values are extracted from columns 25--36.
*
*     17 May 1992: John Stone (Ketron) mentioned trouble if rows and
*                  columns have the same name.  The quick fix is to
*                  search column names from the beginning always,
*                  rather than from position jmark.  Just set jmark = 1.
*     29 Apr 2001: Allow for missing names.
*     30 Apr 2001: Initialize hs(*) before reading file.  (Why not long ago??)
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj

      integer            id(5)
      character          key*4
      character          lll*4       , lul*4       , lxl*4       ,
     $                   lxu*4       , lsb*4       , lend*4
      data               lll /' LL '/, lul /' UL '/, lxl /' XL '/,
     $                   lxu /' XU '/, lsb /' SB '/, lend/'ENDA'/

      bplus  = 0.9d+0*plinfy
      if (iprint .gt. 0) write(iprint, 1999) insrt
      if (isumm  .gt. 0) write(isumm , 1999) insrt
                         read (insrt , 1000) id
      if (iprint .gt. 0) write(iprint, 2000) id
      l1    = n + 1

*     Make logicals basic.
*     30 Apr 2001: Also make columns nonbasic (ignoring given hs(*)).
*                  Might as well assume the basis info is essentially complete.

      call iload1( n, 0, hs    , 1 )
      call iload1( m, 3, hs(l1), 1 )
      ignord = 0
      nbs    = 0
      ns     = 0
      notfnd = 0
      ncard  = 0
*     jmark  = 1
      lmark  = l1
      ndum   = n + 100000

*     Read names until ENDATA

      do 300 nloop = 1, ndum
         read(insrt, 1020) key, name1c, name2c, name1r, name2r, xj
         if (key .eq. lend) go to 310

*        Look for name1.  It may be a column or a row,
*        since a superbasic variable could be either.
*        17 May 1992: Set jmark = 1 so columns are searched first.
*        This avoids trouble when columns and rows have the same name.

         ncard  = nloop
         jmark  = 1
         call m4name( m, n, nb, nname, name1, name2, name1c, name2c,
     $                ncard, notfnd, mer, 1, nb, jmark, j )
         if (   j  .le. 0) go to 300
         if (hs(j) .gt. 1) go to 290

         if (key .eq. lxl  .or.  key .eq. lxu) then
*           ------------------------------------------------------------
*           XL, XU (exchange card) -- make col j basic, row l nonbasic.
*           ------------------------------------------------------------
                  
*           Look for name2.  It has to be a row.
         
            call m4name( m, n, nb, nname, name1, name2, name1r, name2r,
     $                   ncard, notfnd, mer, l1, nb, lmark, l )
            if (l  .le.  0  ) go to 300
            if (l  .eq. jobj) go to 290
            if (hs(l) .ne. 3) go to 290

            nbs    = nbs + 1
            hs(j)  = 3
            if (key .eq. lxl) then
               hs(l) = 0
               if (bl(l) .gt. -bplus) xn(l) = bl(l)
            else
               hs(l) = 1
               if (bu(l) .lt.  bplus) xn(l) = bu(l)
            end if

*        ---------------------------------------------------------------
*        else LL, UL, SB  --  only  j  and  xj  are relevant.
*        ---------------------------------------------------------------
         else if (key .eq. lll) then
            hs(j) = 0
         else if (key .eq. lul) then
            hs(j) = 1
         else if (key .eq. lsb) then
            hs(j) = 2
            ns    = ns + 1
         else
            go to 290
         end if

*        Save xj.

         if (abs(xj) .lt. bplus) xn(j) = xj
         go to 300

*        Card ignored.

  290    ignord = ignord + 1
         if (iprint .gt. 0  .and.  ignord .le. mer) then
            write(iprint, 2010) ncard, key, name1c,name2c,name1r,name2r
         end if
  300 continue

  310 ignord = ignord + notfnd
      if (iprint .gt. 0) write(iprint, 2050) ncard, ignord, nbs, ns
      if (isumm  .gt. 0) write(isumm , 2050) ncard, ignord, nbs, ns
      if (insrt  .ne. iread) rewind insrt
      return

 1000 format(14x, 2a4, 2x, 3a4)
 1005 format(2a4)
 1020 format(3a4, 2x, 2a4, 2x, e12.5)
 1999 format(/ ' INSERT file to be input from file', i4)
 2000 format(/ ' NAME', 10x, 2a4, 2x, 3a4)
 2010 format(' XXX  Line', i6, '  ignored:', 8x, 3a4, 2x, 2a4)
 2050 format(/ ' No. of lines read      ', i6, '  Lines ignored', i6
     $       / ' No. of basics specified', i6, '  Superbasics  ', i6)

      end ! subroutine m4inst

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4load( m, n, nb, nname, ns,
     $                   hs, bl, bu, xn, name1, name2 )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            name1(nname), name2(nname)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m4load  inputs a load file, which may contain a full or partial
*     list of row and column names and their states and values.
*     Valid keys are   BS, LL, UL, SB.
*
*     07 Oct 1994: Fixed columns are made nonbasic, even if they are
*                  specified to be BS or SB.  This may lead to a
*                  singular initial basis (to be patched up by the
*                  first LU factorize), but it could be more helpful
*                  later if the problem is infeasible.  Previously,
*                  Yan Wang, ECNZ, Hamilton, New Zealand found that
*                  some fixed variables were still basic (at infeasible
*                  values) in the final infeasible basis.  In terms of
*                  finding the cause of infeasibility, it is more
*                  intuitive to let non-fixed variables be infeasible.
*                  Note: We let fixed slacks be basic, since they are
*                  often needed to give a nonsingular basis.
*     29 Apr 2001: Allow for missing names.
*     30 Apr 2001: Initialize hs(*) before reading file.  (Why not long ago??)
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj

      integer            id(5)
      character          key*4
      character          lbs*4       , lll*4       , lul*4       ,
     $                   lsb*4       , lend*4
      data               lbs /' BS '/, lll /' LL '/, lul /' UL '/,
     $                   lsb /' SB '/, lend/'ENDA'/

      bplus  = 0.9d+0*plinfy
      if (iprint .gt. 0) write(iprint, 1999) iload
      if (isumm  .gt. 0) write(isumm , 1999) iload
                         read (iload , 1000) id
      if (iprint .gt. 0) write(iprint, 2000) id

*     30 Apr 2001: Initialize everything to be nonbasic (ignoring given hs(*)).
*                  Might as well assume the basis info is essentially complete.

      call iload1( nb, 0, hs, 1 )

      ignord = 0
      nbs    = 0
      ns     = 0
      notfnd = 0
      ncard  = 0
      jmark  = 1
      ndum   = n + 100000

*     Read names until ENDATA is found.

      do 300 nloop = 1, ndum
         read (iload, 1020) key, id1, id2, xj
         if (key .eq. lend) go to 310

         ncard  = nloop
         call m4name( m, n, nb, nname, name1, name2, id1, id2,
     $                ncard, notfnd, mer, 1, nb, jmark, j )
         if (j .le. 0) go to 300

*        The name id1-id2 belongs to the j-th variable.

         if (hs(j) .gt. 1) go to 290
         if (j   .eq.jobj) go to  90

*        07 Oct 1994:  Make fixed columns nonbasic.

         if (j .le. n) then
            if (bl(j) .eq. bu(j)) go to 100
         end if

*        Test for BS, LL, UL, SB.

         if (key .eq. lbs) go to  90
         if (key .eq. lll) go to 100
         if (key .eq. lul) go to 150
         if (key .eq. lsb) go to 200
         go to 290

*        Make basic.

   90    nbs    = nbs + 1
         hs(j)  = 3
         go to 250

*        LO or UP.

  100    hs(j)  = 0
         go to 250

  150    hs(j)  = 1
         go to 250

*        Make superbasic.

  200    ns     = ns + 1
         hs(j)  = 2

*        Save  x  values.

  250    if (abs(xj) .lt. bplus) xn(j) = xj
         go to 300

*        Card ignored.

  290    ignord = ignord + 1
         if (ignord .le. mer) then
            if (iprint .gt. 0) write(iprint, 2010) ncard, key, id1, id2
         end if
  300 continue

  310 ignord = ignord + notfnd
      if (iprint .gt. 0) write(iprint, 2050) ncard, ignord, nbs, ns
      if (isumm  .gt. 0) write(isumm , 2050) ncard, ignord, nbs, ns

*     Make sure the linear objective is basic.

      if (iobj  .gt. 0  .and.  hs(jobj) .ne. 3) then
         hs(jobj) = 3

*        Swap obj with last basic variable.

         do 850 j = nb, 1, -1
            if (hs(j) .eq. 3) go to 860
  850    continue

  860    hs(j)  = 0
      end if

      if (iload .ne. iread) rewind iload
      return

 1000 format(14x, 2a4, 2x, 3a4)
 1005 format(2a4)
 1020 format(3a4, 12x, e12.5)
 1999 format(/ ' LOAD file to be input from file', i4)
 2000 format(/ ' NAME', 10x, 2a4, 2x, 3a4)
 2010 format(' XXX  Line', i6, '  ignored:', 8x, 3a4, 2x, 2a4)
 2050 format(/ ' No. of lines read      ', i6, '  Lines ignored', i6
     $       / ' No. of basics specified', i6, '  Superbasics  ', i6)

      end ! subroutine m4load

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4oldb( m, n, nb, ns,
     $                   hs, bl, bu, xn )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m4oldb  inputs a compact basis file from file  ioldb.
*
*     07 Oct 1994: Fixed columns are made nonbasic, even if they are
*                  specified to be BS or SB.  See m4load.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m5log1/ idebug,ierr,lprint

      character          id*80

      bplus  = 0.9d+0*plinfy
      if (iprint .gt. 0) write(iprint, 1999) ioldb
      if (isumm  .gt. 0) write(isumm , 1999) ioldb
         read (ioldb , 1000) id
      if (iprint .gt. 0) then
         write(iprint, 2000) id
      end if

         read (ioldb , 1005) id(1:52), newm, newn, ns
      if (iprint .gt. 0) then
         write(iprint, 2005) id(1:52), newm, newn, ns
      end if

      if (newm .ne. m  .or.  newn .ne. n) go to 900
      read (ioldb , 1010) hs

*     Set values for nonbasic variables.

      do 200 j = 1, nb
         js    = hs(j)
         if (js .le. 1) then
            if (js .eq. 0) xj = bl(j)
            if (js .eq. 1) xj = bu(j)
            if (abs(xj) .lt. bplus) xn(j) = xj
         end if
  200 continue

*     Load superbasics (and other numerical values).
*     02 May 2001: This version leaves hs(j) unaltered.
*                  (Earlier manuals implied that entries here
*                  would be made superbasic regardless.)

      ns     = 0
      ndummy = m + n + 10000

      do 300 idummy = 1, ndummy
         read(ioldb, 1020, end=310) j, xj
         if (j .le.  0) go to  310
         if (j .le. nb) then
            xn(j)  = xj
            if (hs(j) .eq. 2) ns = ns + 1
         end if
  300 continue

  310 if (ns .gt. 0) then
         if (iprint .gt. 0) write(iprint, 2010) ns
         if (isumm  .gt. 0) write(isumm , 2010) ns
      end if

*     07 Oct 1994:  Make sure fixed columns are nonbasic.

      nfix   = 0
      do 400 j = 1, n
         if (hs(j) .gt. 1) then
            if (bl(j) .eq. bu(j)) then
               nfix   = nfix + 1
               hs(j)  = 0
            end if
         end if
  400 continue

      if (nfix .gt. 0) then
         if (iprint .gt. 0) write(iprint, 2040) nfix
         if (isumm  .gt. 0) write(isumm , 2040) nfix
      end if

      go to 990

*     Error exits.

  900 call m1page( 1 )
      if (iprint .gt. 0) write(iprint, 3000)
      if (isumm  .gt. 0) write(isumm , 3000)
      ierr   = 30

  990 if (ioldb .ne. iread) rewind ioldb
      return

 1000 format(a80)
 1005 format(a52, 2x, i7, 3x, i7, 4x, i5)
 1010 format(80i1)
 1020 format(i8, e24.14)
 1999 format(/ ' OLD BASIS file to be input from file', i4)
 2000 format(1x, a80)
 2005 format(1x, a52, 'M=', i7, ' N=', i7, ' SB=', i5)
 2010 format(' No. of superbasics loaded:      ', i7)
 2040 format(' No. of fixed cols made nonbasic:', i7)
 3000 format(/ ' EXIT -- the basis file dimensions do not match',
     $   ' this problem')

      end ! subroutine m4oldb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4savb( mode, m, mbs, n, nb, nn, nname, nscl,
     $                   msoln, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, kb, ascale, bl, bu,
     $                   name1, name2,
     $                   pi, rc, xn, y, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            ha(ne), hs(nb)
      integer            ka(nka), kb(mbs), name1(nname), name2(nname)
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb)
      double precision   pi(m), rc(nb), xn(nb), y(m), z(nwcore)

*     ------------------------------------------------------------------
*     m4savb  saves basis files  and/or  prints the solution.
*
*     If mode = 1, the problem is first unscaled, then from 0 to 4 files
*     are saved (PUNCH file, DUMP file, SOLUTION file, REPORT file,
*     in that order).
*     A new BASIS file, if any, will already have been saved by m5solv.
*     A call with mode = 1 must precede a call with mode = 2.
*
*     If mode = 2, the solution is printed under the control of msoln
*     (which is set by the Solution keyword in the SPECS file).
*
*     18 Nov 1991: Scaled pinorm saved for use in m4soln.
*     25 Nov 1991: rc added as parameter to return reduced costs.
*     31 Jan 1991: Call m2rcA  to get the reduced costs.
*     18 Dec 1992: Maximum primal and dual infeasibilities computed
*                  and printed here.
*     14 Feb 2004: The nonlinear slack values are now copied into the
*                  correct part of xn.  m4soln prints them from there.
*                  This may result in slightly different behavior
*                  when xn is examined or used for restarts.
*                  However, it seems reasonable to set xn this way.
*
*                  m4soln now deals with the UNSCALED solution,
*                  so we no longer have to save the scaled pinorm.
*     09 Mar 2004: Print max elements of x and pi (and which ones).
*                  This is more helpful than piNorm >= 1.
*     17 Jun 2004: Save biggest primal and dual infeasibilities
*                  before and after scaling, so GAMS can use them.
*                  dparm(27:30) saves  Binf1,  Dinf1,  Binf,  Dinf.
*                  iparm(27:30) saves jBinf1, jDinf1, jBinf, jDinf.
*     ------------------------------------------------------------------

      logical            alone, AMPL, GAMS, MINT, page1, page2
      common    /m1env / alone, AMPL, GAMS, MINT, page1, page2
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m2parm/ dparm(30),iparm(30)
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      common    /m5step/ featol, tolx0,tolinc,kdegen,ndegen,
     $                   itnfix, nfix(2)
      common    /m5tols/ toldj(3),tolx,tolpiv,tolrow,rowerr,xnorm
      common    /m7len / fobj  ,fobj2 ,nnobj ,nnobj0
      common    /m7loc / lgobj ,lgobj2
      common    /m7tols/ xtol(2),ftol(2),gtol(2),pinorm,rgnorm,tolrg
      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /m8loc / lfcon ,lfcon2,lfdif ,lfdif2,lfold ,
     $                   lblslk,lbuslk,lxlam ,lrhs  ,
     $                   lgcon ,lgcon2,lxdif ,lxold
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave

      character          istate*12
      logical            feasbl, prnt
      save               pnorm1, pnorm2
      parameter        ( one = 1.0d+0 )

      feasbl = ninf .eq. 0
      k      = ierr + 1
      call m4stat( k, istate )

      if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        mode = 1.
*        Compute rc and unscale everything.
*        Then save basis files.
*        ---------------------------------------------------------------

*        Compute reduced costs rc(*) for all columns and rows.
*        Find the maximum bound and dual infeasibilities.

         call m2rcA ( feasbl, featol, minimz,
     $                m, n, nb, nnobj, nnobj0,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, z(lgobj), pi, rc, xn )
         call m2Binf( nb, bl, bu, xn, Binf, jBinf )
         call m2Dinf( nb, jobj, bl, bu, rc, xn, Dinf, jDinf )

         Binf1     =  Binf
         Dinf1     =  Dinf
         jBinf1    = jBinf
         jDinf1    = jDinf
         dparm(27) =  Binf
         dparm(28) =  Dinf
         iparm(27) = jBinf
         iparm(28) = jDinf

         maxx   = idamax( n, xn, 1 )
         maxp   = idamax( m, pi, 1 )
         xnorm  = abs( xn(maxx) )
         pimax  = abs( pi(maxp) )
         pinorm = max( pimax, one )

         maxx1  = maxx
         maxp1  = maxp
         xnorm1 = xnorm
         pimax1 = pimax
         pnorm1 = pinorm

*        Unscale a, bl, bu, pi, xn, rc, fcon, gobj and xnorm, pinorm.
*        (Previously, m4soln used the scaled pinorm, but no more.)

         if (lscale .gt. 0) then
            call m2scla( 2, m, n, nb, ne, nka,
     $                   ha, ka, a, ascale, bl, bu, pi, xn )

            call dddiv ( nb, ascale, 1, rc, 1 )

            if (lscale .eq. 2) then
               if (nncon .gt. 0)
     $         call ddscl ( nncon, ascale(n+1), 1, z(lfcon), 1 )

               if (nnobj .gt. 0)
     $         call dddiv ( nnobj, ascale     , 1, z(lgobj), 1 )
            end if

*           09 Mar 2004: Previously we called m5setp to redefine pinorm.
*                        y was not used.
*                        Now use max elements.

            maxx   = idamax( n, xn, 1 )
            maxp   = idamax( m, pi, 1 )
            xnorm  = abs( xn(maxx) )
            pimax  = abs( pi(maxp) )
            pinorm = max( pimax, one )

            call m5setp( 3, m, y, pi, z, nwcore )
            call m2Binf( nb, bl, bu, xn, Binf, jBinf )
            call m2Dinf( nb, jobj, bl, bu, rc, xn, Dinf, jDinf )
         end if

         pnorm2 = pinorm
         dparm(29) =  Binf
         dparm(30) =  Dinf
         iparm(29) = jBinf
         iparm(30) = jDinf

*        ---------------------------------------------------------------
*        Print various scaled and unscaled norms.
*        ---------------------------------------------------------------
         if (lscale .gt. 0) then
            if (iprint .gt. 0)
     &         write(iprint, 1010) maxx1, xNorm1, maxp1, pimax1
            if (isumm  .gt. 0)
     &         write(isumm , 1010) maxx1, xNorm1, maxp1, pimax1
         end if
            if (iprint .gt. 0)
     &         write(iprint, 1020) maxx , xNorm , maxp , pimax
            if (isumm  .gt. 0)
     &         write(isumm , 1020) maxx , xNorm , maxp , pimax
         if (lscale .gt. 0) then
            if (iprint .gt. 0) write(iprint, 1030) jBinf1, Binf1 ,
     $                                             jDinf1, Dinf1
            if (isumm  .gt. 0) write(isumm , 1030) jBinf1, Binf1 ,
     $                                             jDinf1, Dinf1
         end if
            if (iprint .gt. 0) write(iprint, 1040) jBinf , Binf  ,
     $                                             jDinf , Dinf
            if (isumm  .gt. 0) write(isumm , 1040) jBinf , Binf  ,
     $                                             jDinf , Dinf

*        Change the sign of pi and rc if feasible and maximizing.

         if (ninf .eq. 0  .and.  minimz .lt. 0) then
            call dscal ( m , -one, pi, 1 )
            call dscal ( nb, -one, rc, 1 )
         end if

*        Compute nonlinear constraint infeasibilities (violations).
*        14 Feb 2004: Copy them into xn(n+1:n+nncon-1)

         if (nncon .gt. 0) then
            call m8Cinf( n, nb, nncon, Cinf, iCinf,
     $                   ne, nka, a, ha, ka,
     $                   bl, bu, z(lfcon), xn, y, z, nwcore )
            call dcopy ( nncon, y, 1, xn(n+1), 1 )
            vimax  = Cinf
            maxvi  = iCinf
            virel  = vimax / (one + xnorm)
            if (iprint .gt. 0) write(iprint, 1080) vimax
            if (isumm  .gt. 0) write(isumm , 1080) vimax
         end if

*        ---------------------------------------------------------------
*        Output PUNCH, DUMP, SOLUTION and/or REPORT files.
*        ---------------------------------------------------------------
         if (ipnch .gt. 0) then
            call m4pnch( ipnch, m, n, nb, nname,
     $                   hs, bl, bu, xn, name1, name2 )
         end if

         if (idump .gt. 0) then
            call m4dump( idump, m, n, nb, nname,
     $                   hs, bl, bu, xn, name1, name2 )
         end if

         ! 14 Feb 2004: No longer worry about scaled solution.
         ! pinorm = pnorm1  
         if (isoln .gt. 0) then
            call m4soln( .true., m, n, nb, nname, nscl,
     $                   nn, nnobj, nnobj0, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, ascale, bl, bu,
     $                   z(lgobj), pi, rc, xn, y,
     $                   name1, name2, istate, z, nwcore )
         end if

         ! 09 Oct 2000: Export A, B or (B S) to Report file 91, 92 or 93.

         if (ireprt .ge. 91  .and.  ireprt .le. 93) then
            call m2xmat( ireprt, n, nb, ne, nka, 
     &                   a, ha, ka, hs )

         else if (ireprt .gt. 0) then
            call m4rept( .true., m, n, nb, nname, nscl,
     $                   nn, nnobj, nnobj0, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, ascale, bl, bu,
     $                   z(lgobj), pi, rc, xn, y,
     $                   name1, name2, istate, z, nwcore )
         end if

         ! 14 Feb 2004: No longer worry about scaled solution.
         ! pinorm = pnorm2

      else
*        ---------------------------------------------------------------
*        mode = 2.    Print solution if requested.
*
*        msoln  = 0   means   no
*               = 1   means   if optimal, infeasible or unbounded
*               = 2   means   yes
*               = 3   means   if error condition
*        ---------------------------------------------------------------
         prnt   = iprint .gt. 0  .and.  msoln .gt. 0
         if ((msoln .eq. 1  .and.  ierr .gt. 2)  .or.
     $       (msoln .eq. 3  .and.  ierr .le. 2)) prnt = .false.
         if ( prnt ) then

            ! 14 Feb 2004: No longer worry about scaled solution.
            ! pinorm = pnorm1
            call m4soln( .false., m, n, nb, nname, nscl,
     $                   nn, nnobj, nnobj0, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, ascale, bl, bu,
     $                   z(lgobj), pi, rc, xn, y,
     $                   name1, name2, istate, z, nwcore )
            ! 14 Feb 2004: No longer worry about scaled solution.
            ! pinorm = pnorm2
            if (isumm  .gt. 0) write(isumm, 1200) iprint

         else if (.not. (GAMS .or. AMPL)) then
            if (isumm  .gt. 0) write(isumm, 1300)
         end if
      end if

      return

 1010 format(  ' Max x       (scaled)', i9, 1p, e8.1,
     &     2x, ' Max pi      (scaled)', i9,     e8.1)
 1020 format(  ' Max x               ', i9, 1p, e8.1,
     &     2x, ' Max pi              ', i9,     e8.1)
 1030 format(  ' Max Prim inf(scaled)', i9, 1p, e8.1,
     $     2x, ' Max Dual inf(scaled)', i9,     e8.1)
 1040 format(  ' Max Primal infeas   ', i9, 1p, e8.1,
     $     2x, ' Max Dual infeas     ', i9,     e8.1)
 1080 format(  ' Nonlinear constraint violn', 1p, e11.1)
 1100 format(2a4)
 1200 format(/ ' Solution printed on file', i4)
 1300 format(/ ' Solution not printed')

      end ! subroutine m4savb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4dump( idump, m, n, nb, nname,
     &                   hs, bl, bu, xn, name1, name2 )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            name1(nname), name2(nname)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m4dump outputs basis names in a format compatible with m4load.
*     This file is normally easier to modify than a punch file.
*     29 Apr 2001: Allow for missing names.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax

      character          key(4)*4
      data               key   /' LL ', ' UL ', ' SB ', ' BS '/

      write(idump, 2000) name

      do 500 j = 1, nb
         call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
         k     = hs(j) + 1
         write(idump, 2100) key(k), id1, id2, xn(j)
  500 continue

      write(idump , 2200)
      if (iprint .gt. 0) write(iprint, 3000) idump
      if (isumm  .gt. 0) write(isumm , 3000) idump
      if (idump .ne. iprint) rewind idump
      return

 2000 format('NAME', 10x, 2a4, 2x, '   DUMP/LOAD')
 2100 format(3a4, 12x, 1p, e12.5)
 2200 format('ENDATA')
 3000 format(/ ' DUMP file saved on file', i4)

      end ! subroutine m4dump

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4newb( mode, inewb, m, n, nb, nn, ns, ms, nscl, fsub,
     $                   kb, hs, ascale, bl, bu, x, xn, istate )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            kb(ms)
      double precision   ascale(nscl), bl(nb), bu(nb)
      double precision   x(ms), xn(nb)
      character          istate*12

*     ------------------------------------------------------------------
*     m4newb  saves a compact basis on file inewb.  Called from m5solv.
*     If mode = 1, the save is a periodic one due to the save frequency.
*     If mode = 2, m5solv has just finished the current problem.
*                 (Hence, mode 2 calls occur at the end of every cycle.)
*
*     02 May 2001: Following hs(*), superbasics used to be output first,
*                  then certain other values.  For simplicity,
*                  now output all wanted values in natural order.
*                  For LPs, this is just superbasics (typically none).
*                  For NLPs, include basics and nonlinear variables.
*                  Also output hs(j) at end of line to be helpful.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5lp1 / itn,itnlim,nphs,kmodlu,kmodpi

      logical            scaled

      scaled = lscale .gt. 0
      obj    = sinf
      if (ninf .eq. 0) obj = minimz * fsub

*     Output header cards and the state vector.

      write(inewb, 1000) name, itn, istate, ninf, obj
      write(inewb, 1005) mobj, mrhs, mrng, mbnd, m, n, ns
      write(inewb, 1010) hs

*     02 May 2001: Output all relevant variables in natural order.

      if (nn .le. 0) then       ! LPs
         if (ns .gt. 0) then    ! Output superbasics only
            do j = 1, nb
               js     = hs(j)
               if (js .eq. 2) then
                  xj     = xn(j)
                  if (scaled) xj = xj * ascale(j)
                  write(inewb, 1020) j, xj, js
               end if
            end do
         end if

      else                 ! NLPs
         do j = 1, nb
            js     = hs(j)
            if (j. le. nn  .or.  js .ge. 2) then
               xj     = xn(j)
               if (scaled) xj = xj * ascale(j)
               write(inewb, 1020) j, xj, js
            end if
         end do
      end if

*     Terminate the list with a zero.

      j     = 0
      write(inewb, 1020) j
      if (inewb .ne. iprint) rewind inewb
      if (iprint .gt. 0) write(iprint, 1030) inewb, itn
      if (isumm  .gt. 0) write(isumm , 1030) inewb, itn
      return

 1000 format(2a4, '  ITN', i8, 4x, a12, '  NINF', i7,
     $       '      OBJ', 1p, e21.12)
 1005 format('OBJ=', 2a4, ' RHS=', 2a4, ' RNG=', 2a4, ' BND=', 2a4,
     $       ' M=', i7,  ' N=', i7, ' SB=', i5)
 1010 format(80i1)
 1020 format(i8, 1p, e24.14, i3)
 1030 format(/ ' NEW BASIS file saved on file', i4, '    itn =', i7)

      end ! subroutine m4newb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4pnch( ipnch, m, n, nb, nname,
     &                   hs, bl, bu, xn, name1, name2 )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      integer            name1(nname), name2(nname)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m4pnch  outputs a PUNCH file (list of basis names, states and
*     values) in a format that is compatible with MPS/360.
*     29 Apr 2001: Allow for missing names.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax

      parameter        ( zero = 0.0d+0 )
      character          key(5)*4, ibl*4
      data               key    /' LL ', ' UL ', ' SB ', ' XL ', ' XU '/
      data               ibl    /'    '/

      write(ipnch, 2000) name
      irow   = n

      do 500  j = 1, n
         call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )

         k      = hs(j)
         if (k .eq. 3) then

*           Basics -- find the next row that isn't basic.

  300       irow   = irow + 1
            if (irow .le. nb) then
               k      = hs(irow)
               if (k .eq. 3) go to 300

               call m4id  ( irow, m, n, nb, nname,
     &                      name1, name2, id1r, id2r )
               if (k .eq. 2) k = 0
               write(ipnch, 2100) key(k+4), id1, id2, id1r, id2r, xn(j)
            end if
         else

*           Skip nonbasic variables with zero lower bounds.

            if (k .le. 1) then
               if (bl(j) .eq. zero  .and.  xn(j) .eq. zero) go to 500
            end if
            write(ipnch, 2100) key(k+1), id1, id2, ibl, ibl, xn(j)
         end if
  500 continue

*     Output superbasic slacks.

      do 700 j = n + 1, nb
         if (hs(j) .eq. 2) then
            call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
            write(ipnch, 2100) key(3), id1, id2, ibl, ibl, xn(j)
         end if
  700 continue

      write(ipnch , 2200)
      if (iprint .gt. 0) write(iprint, 3000) ipnch
      if (isumm  .gt. 0) write(isumm , 3000) ipnch
      if (ipnch .ne. iprint) rewind ipnch
      return

 2000 format('NAME', 10x, 2a4, 2x, 'PUNCH/INSERT')
 2100 format(3a4, 2x, 2a4, 2x, 1p, e12.5)
 2200 format('ENDATA')
 3000 format(/ ' PUNCH file saved on file', i4)

      end ! subroutine m4pnch

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4rept( ondisk, m, n, nb, nname, nscl,
     $                   nn, nnobj, nnobj0, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, ascale, bl, bu, gobj, pi, rc, xn, y,
     $                   name1, name2, istate, z, nwcore )

      implicit           double precision (a-h,o-z)
      logical            ondisk
      integer            ha(ne), hs(nb)
      integer            ka(nka), name1(nname), name2(nname)
      character          istate*12
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb)
      double precision   gobj(nnobj0), pi(m), rc(nb), xn(nb), y(m),
     $                   z(nwcore)

*     ------------------------------------------------------------------
*     m4rept  has the same parameter list as m4soln, the routine that
*     prints the solution.  It will be called if the SPECS file
*     specifies  REPORT file  n  for some positive value of  n.
*
*     pi contains the unscaled dual solution.
*     xn contains the unscaled primal solution.  There are n + m = nb
*        values (n structural variables and m slacks, in that order).
*     y  contains the true slack values for nonlinear constraints
*        in its first nncon components (computed by m8viol).
*
*     This version of m4rept does nothing.    Added for PILOT, Oct 1985.
*     31 Oct 1991: Name changed from "report" to "m4rept".
*                  Parameters altered to allow for MPS or generic names.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
*     ------------------------------------------------------------------

      if (iprint .gt. 0) write(iprint, 1000)
      if (isumm  .gt. 0) write(isumm , 1000)
      return

 1000 format(/ ' XXX Report file requested.  m4rept does nothing.')

      end ! subroutine m4rept

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4soln( ondisk, m, n, nb, nname, nscl,
     $                   nn, nnobj, nnobj0, ns,
     $                   ne, nka, a, ha, ka,
     $                   hs, ascale, bl, bu, gobj, pi, rc, xn, y,
     $                   name1, name2, istate, z, nwcore )

      implicit           double precision (a-h,o-z)
      logical            ondisk
      integer            ha(ne), hs(nb)
      integer            ka(nka), name1(nname), name2(nname)
      character          istate*12
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb)
      double precision   gobj(nnobj0), pi(m), rc(nb), xn(nb), y(m),
     $                   z(nwcore)

*     ------------------------------------------------------------------
*     m4soln  is the standard output routine for printing the solution.
*
*     On entry,
*     pi    contains the dual solution.
*     xn    contains the primal solution.  There are n + m = nb values
*           (n structural variables and m slacks, in that order).
*     rc    contains reduced costs for all variables:
*           rc(1:n)    =  gobj + c - A'pi  for structurals,
*           rc(n+1:nb) = -pi  for slacks.
*           pi and rc corresond to the Phase-1 objective
*           if the solution is infeasible.
*     y     is not used.
*
*     All quantities a, bl, bu, pi, rc, xn, fcon, gobj are unscaled
*     and adjusted in sign if maximizing.  (fcon is not used here.)  
*
*     If ondisk is true, the solution is output to the solution file.
*     Otherwise, it is output to the printer.
*
*     31 Jan 1991: rc is now an input parameter.  It is set in m4savb.
*     14 Feb 2004: Solution flags are now set by misolf.
*                  They are defined by the UNSCALED solution!
*                  ascale(*) is no longer referenced.
*                  Also y(*) is no longer used.
*                  Previously it input true nonlinear slack values,
*                  but now m4savb puts them into the right part of xn.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5lp1 / itn,itnlim,nphs,kmodlu,kmodpi
      common    /m5tols/ toldj(3),tolx,tolpiv,tolrow,rowerr,xnorm
      common    /m7tols/ xtol(2),ftol(2),gtol(2),pinorm,rgnorm,tolrg
      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /cycle2/ objtru,suminf,numinf

      logical            feasbl, infsbl, maximz !!, scaled
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------

      bplus  = 0.9d+0*plinfy
   !! scale  = one
      feasbl = ninf   .eq. 0
      infsbl = .not. feasbl
      maximz = minimz .lt. 0
   !! scaled = lscale .gt. 0
      lpr    = iprint
      if (ondisk) lpr = isoln

      call m1page( 1 )
      if (infsbl) write(lpr, 1000) name, ninf, sinf
      if (feasbl) write(lpr, 1002) name, objtru
      write(lpr, 1004) istate, itn, ns
      write(lpr, 1005) mobj, minmax, mrhs, mrng, mbnd
      write(lpr, 1010)
**    tolfea = 0.1 * tolx
**    tolopt = 0.1 * toldj(3) * pinorm
** 05 Oct 1991: Might as well flag according to the tolerances
*               actually used.
      tolfea = tolx
      tolopt = toldj(3) * pinorm

*     ------------------------------------------------------------------
*     Output the ROWS section.
*     ------------------------------------------------------------------
      do 300 iloop = 1, m
         i      = iloop
         j      = n + i
      !! if (scaled) scale = ascale(j)
         js     = hs(j)
         b1     = bl(j)
         b2     = bu(j)
         xj     = xn(j)
         py     = pi(i)
         dj     = rc(j)

*        Define row and slack activities.

         ! 14 Feb 2004: Nonlinear slacks are now put into xn(n+1...)
         !              by m4savb.
      !! if (i .le. nncon) xj = y(i)

         row    = - xj
         d1     =   b1 - xj
         d2     =   xj - b2
         slk    = - d1
         if (abs( d1  )  .gt.  abs( d2 )) slk =   d2
         if (abs( slk )  .ge.  bplus    ) slk = - row

      !! SCALING NOW IGNORED.
      !! d1     = d1 / scale
      !! d2     = d2 / scale
      !! djtest = dj * scale

      !! NEXT LINES NOW DONE IN misolf.
      !! if (feasbl) then         
      !!    if (   maximz   ) djtest =  - djtest
      !! end if

         ! Change slack bounds into row bounds.

         b1     = - b2
         b2     = - bl(j)

         call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
         call misolf( m, n, nb, j, jkey, jstate,
     &                hs, bl, bu, rc, xn )
         call m4solp( ondisk, bplus, tolfea, tolopt,
     $                jkey, jstate,
     $                j, id1, id2, row, slk, b1, b2, py, i )
  300 continue

*     ------------------------------------------------------------------
*     Output the COLUMNS section.
*     ------------------------------------------------------------------
      call m1page( 1 )
      write(lpr, 1020)

      do 400 jloop = 1, n
         j      = jloop
      !! if (scaled) scale = ascale(j)
         js     = hs(j)
         b1     = bl(j)
         b2     = bu(j)
         xj     = xn(j)
         cj     = zero
         dj     = rc(j)

         do 320 k = ka(j), ka(j+1) - 1
            ir    = ha(k)
            if (ir .eq. iobj) cj = a(k)
  320    continue

         d1     =   (b1 - xj) !! / scale
         d2     =   (xj - b2) !! / scale
      !! djtest = - dj * scale
         if (feasbl) then
            if (j .le. nnobj) cj     =    cj + gobj(j)
         !! if (   maximz   ) djtest =  - djtest  !! Now done in misolf
         end if

         call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
         call misolf( m, n, nb, j, jkey, jstate,
     &                hs, bl, bu, rc, xn )
         call m4solp( ondisk, bplus, tolfea, tolopt,
     $                jkey, jstate,
     $                j, id1, id2, xj, cj, b1, b2, dj, m+j )
  400 continue

      if (ondisk) then
         if (isoln .ne. iprint) rewind isoln
         if (iprint .gt. 0) write(iprint, 1400) isoln
         if (isumm  .gt. 0) write(isumm , 1400) isoln
      end if
      return

 1000 format(' NAME', 11x, 2a4, 13x,
     $   ' INFEASIBILITIES', i7, 1p, e16.4)
 1002 format(' NAME', 11x, 2a4, 13x,
     $   ' OBJECTIVE VALUE', 1p, e23.10)
 1004 format(/ ' STATUS', 9x, a12, 9x,
     $   ' ITERATION', i7, '    SUPERBASICS', i7)
 1005 format(/
     $   ' OBJECTIVE', 6x, 2a4, ' (', a3, ')' /
     $   ' RHS      ', 6x, 2a4 /
     $   ' RANGES   ', 6x, 2a4 /
     $   ' BOUNDS   ', 6x, 2a4)
 1010 format(/ ' SECTION 1 - ROWS' //
     $   '  NUMBER  ...ROW.. STATE  ...ACTIVITY...  SLACK ACTIVITY',
     $   '  ..LOWER LIMIT.  ..UPPER LIMIT.  .DUAL ACTIVITY    ..I' /)
 1020 format(  ' SECTION 2 - COLUMNS' //
     $   '  NUMBER  .COLUMN. STATE  ...ACTIVITY...  .OBJ GRADIENT.',
     $   '  ..LOWER LIMIT.  ..UPPER LIMIT.  REDUCED GRADNT    M+J' /)
 1400 format(/ ' SOLUTION file saved on file', i4)

      end ! subroutine m4soln

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4solp( ondisk, bplus, tolfea, tolopt,
     $                   jkey, jstate,
     $                   j, id1, id2, xj, cj, b1, b2, dj, k )

      implicit           double precision (a-h,o-z)
      logical            ondisk

*     ------------------------------------------------------------------
*     m4solp  prints one line of the Solution file.
*
*     The following conditions are marked by key:
*
*        D  degenerate basic or superbasic variable.
*        I  infeasible basic or superbasic variable.
*        A  alternative optimum      (degenerate nonbasic dual).
*        N  nonoptimal superbasic or nonbasic (infeasible dual).
*
*     Prior to 14 Feb 2004,
*     tests for these conditions were performed on scaled quantities
*     d1, d2, djtest,
*     since the correct indication was then more likely to be given.
*     On badly scaled problems, the unscaled solution could then appear
*     to be flagged incorrectly, but it would be just an "illusion".
*
*     19 Mar 1994: Use variable format to print big reals in e format.
*     14 Feb 2004: Now use jkey and jstate from misolf.
*                  d1, d2, djtest are no longer used in here.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt

      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
      parameter        ( big  = 1.0d+9                )

      character          line*111
      character          ckey  (0:4)*1
      character          cstate(0:5)*4
      character          lzero*16, lone*16, lmone*16, none*16
      character          e*10
      character          form*82, form1*82

      data               ckey   /' ', 'A', 'D', 'I', 'N'/

      data               cstate /' LL ', ' UL ', 'SBS ',
     $                           ' BS' , ' EQ' , ' FR '/

      data               lzero  /'          .     '/
      data               lone   /'         1.0    '/
      data               lmone  /'        -1.0    '/
      data               none   /'           None '/
      data               e      /' 1p,e16.6,'/
      data               form1  /'(i8, 2x, 2a4, 1x, a1, 1x,a3, 0p,f16.5,
     $ 0p,f16.5, 0p,f16.5, 0p,f16.5, 0p,f16.5, i7)'/
*     ------------------------------------------------------------------

      ! Select format for printing.

      if ( ondisk ) then
         write(line, 1000) j, id1, id2, ckey(jkey), cstate(jstate),
     &                     xj, cj, b1, b2, dj, k
      else
         form = form1
         if (abs( xj ) .ge. big) form(29:38) = e
         if (abs( cj ) .ge. big) form(39:48) = e
         if (abs( b1 ) .ge. big) form(49:58) = e
         if (abs( b2 ) .ge. big) form(59:68) = e
         if (abs( dj ) .ge. big) form(69:78) = e
         write(line, form) j, id1, id2, ckey(jkey), cstate(jstate),
     &                     xj, cj, b1, b2, dj, k
      end if

*     Test for 0.0, 1.0 and -1.0

      if (xj .eq. zero) line(25:40) = lzero
      if (xj .eq.  one) line(25:40) = lone
      if (xj .eq. -one) line(25:40) = lmone
      if (cj .eq. zero) line(41:56) = lzero
      if (cj .eq.  one) line(41:56) = lone
      if (cj .eq. -one) line(41:56) = lmone
      if (b1 .eq. zero) line(57:72) = lzero
      if (b1 .eq.  one) line(57:72) = lone
      if (b1 .eq. -one) line(57:72) = lmone
      if (b2 .eq. zero) line(73:88) = lzero
      if (b2 .eq.  one) line(73:88) = lone
      if (b2 .eq. -one) line(73:88) = lmone
      if (dj .eq. zero) line(89:104)= lzero
      if (dj .eq.  one) line(89:104)= lone
      if (dj .eq. -one) line(89:104)= lmone

      if ( ondisk ) then
         write(isoln , 2000) line
      else
         if (b1 .lt. -bplus) line(57:72) = none
         if (b2 .gt.  bplus) line(73:88) = none
         write(iprint, 2000) line
      end if
      return

 1000 format(i8, 2x, 2a4, 1x, a1, 1x, a3, 1p, 5e16.6, i7)
 2000 format(a)

      end ! subroutine m4solp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m4stat( k, istate )

      integer            k
      character          istate*12

*     ------------------------------------------------------------------
*     m4stat loads istate with words describing the current state.
*     ------------------------------------------------------------------

      character          c(0:5)*12
      data               c     /'PROCEEDING  ',
     $                          'OPTIMAL SOLN',
     $                          'INFEASIBLE  ',
     $                          'UNBOUNDED   ',
     $                          'EXCESS ITNS ',
     $                          'ERROR CONDN '/

      j      = min( k, 5 )
      istate = c(j)

      end ! subroutine m4stat
