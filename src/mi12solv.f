************************************************************************
*
*     File  mi12solv fortran
*
*     misolv
*
*     12 Jul 2000: mi12solv.f separates misolv from calls in
*                  minos3 and minoss.
*     03 Jun 2003: Call m2amat BEFORE m8augl to avoid random Jacobian.
*     02 Feb 2004: (At GAMS) Parameters added to misolv for minose.
*     18 Aug 2005: misolv: Initialize vimax, virel, maxvi.
************************************************************************

      subroutine misolv( mimode, start,
     $                   mxx, nxx, nbxx, nexx, nkax, nnamex,
     $                   iobjxx, objadd,
     $                   a, ha, ka, bl, bu, name1, name2,
     $                   hs, xn, pi, rc, 
     $                   inform, ns, z, nwcore,
     $                   ifuser, m1user, majitn, minitn )

      implicit           double precision (a-h,o-z)
      character*(*)      start
      integer            ha(nexx), hs(nbxx)
      integer            ka(nkax), name1(nnamex), name2(nnamex)
      double precision   a(nexx) , bl(nbxx), bu(nbxx)
      double precision   xn(nbxx), pi(mxx) , rc(nbxx), z(nwcore)
      integer            ifuser, nstat, majitn, minitn
      external           m1user

*     ------------------------------------------------------------------
*     misolv solves the current problem.
*
*     On entry,
*     the SPECS file has been read,
*     all data items have been loaded (including a, ha, ka, ...),
*     and workspace has been allocated within z.
*
*     mimode  =  1 if the call is from minos3 (stand-alone MINOS).
*            ge  2 if the call is from minoss.
*
*     On exit,
*     inform  =  0 if an optimal solution was found,
*             =  1 if the problem was infeasible,
*             =  2 if the problem was unbounded,
*             =  3 if the Iteration limit was exceeded,
*            ge  4 if iterations were terminated by some other
*                  error condition (see the MINOS user's guide).
*
*     01 Oct 1991: minoss, mispec and misolv implemented.
*                  minos1, minos2, minos3 reorganized
*                  to facilitate calling MINOS as a subroutine.
*     25 Nov 1991: nname and rc added as parameters of matmod.
*     10 Apr 1992: objadd added as input parameter.
*     20 Apr 1992: Parameter list revised.  nname, name1, name2 added.
*     27 Jun 1992: Cold, Warm, Hot start implemented.
*     09 Jul 1992: ns initialized here only for Cold starts, just to
*                  help debugging.  m4chek always sets it later.
*     28 Nov 1996: Jacobian now initialized to random numbers in a(*).
*                  Gradient checking now done in m5solv when the
*                  linear constraints are satisfied for the first time.
*     23 Dec 1996: 1. Finally debugged previous mods.
*                  2. Arrays fcon, gcon etc now passed to m5solv.
*     18 Aug 2005: Initialize vimax, virel, maxvi here.
*     ------------------------------------------------------------------
*
*     All common blocks are listed here for reference.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File mcommon fortran.


*  Machine-dependent items.

      logical            alone, AMPL, GAMS, MINT, page1, page2
      common    /m1env / alone, AMPL, GAMS, MINT, page1, page2
      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m1savz/ nbytes,newz
      parameter        ( ntime = 5 )
      common    /m1tim / tlast(ntime), tsum(ntime), numt(ntime), ltime
      common    /m1word/ nwordr,nwordi,nwordh


*  Files, maps, parameters.

      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m2len / mrows,mcols,melms
      common    /m2lu1 / minlu,maxlu,lena,nbelem,ip,iq,lenc,lenr,
     $                   locc,locr,iploc,iqloc,lua,indc,indr
      common    /m2lu2 / factol(5),lamin,nsing1,nsing2
      common    /m2lu3 / lenl,lenu,ncp,lrow,lcol
      common    /m2lu4 / parmlu(30),luparm(30)
      common    /m2mapa/ ne    ,nka   ,la    ,lha   ,lka
      common    /m2mapz/ maxw  ,maxz
      common    /m2parm/ dparm(30),iparm(30)


*  Problem size, MPS names, Scale options.

      common    /m3len / m     ,n     ,nb    ,nscl
      common    /m3loc / lascal,lbl   ,lbu   ,lbbl  ,lbbu  ,
     $                   lhrtyp,lhs   ,lkb
      common    /m3mps1/ lname1,lname2,lkeynm,nname
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3mps5/ aelem(2), id(6), iblank
      common    /m3scal/ sclobj,scltol,lscale


*  LP items.

      common    /m5len / maxr  ,maxs  ,mbs   ,nn    ,nn0   ,nr    ,nx
      common    /m5loc / lpi   ,lpi2  ,lw    ,lw2   ,
     $                   lx    ,lx2   ,ly    ,ly2   ,
     $                   lgsub ,lgsub2,lgrd  ,lgrd2 ,
     $                   lr    ,lrg   ,lrg2  ,lxn
      common    /m5freq/ kchk,kinv,ksav,klog,ksumm,i1freq,i2freq,msoln
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      common    /m5log2/ jq1,jq2,jr1,jr2,lines1,lines2
      common    /m5log3/ djq,theta,pivot,cond,nonopt,jp,jq,modr1,modr2
      logical            prnt0 ,prnt1 ,summ0 ,summ1 ,newhed
      common    /m5log4/ prnt0 ,prnt1 ,summ0 ,summ1 ,newhed
      common    /m5lp1 / itn,itnlim,nphs,kmodlu,kmodpi
      common    /m5lp2 / invrq,invitn,invmod
      common    /m5prc / nparpr,nmulpr,kprc,newsb
      common    /m5step/ featol, tolx0,tolinc,kdegen,ndegen,
     $                   itnfix, nfix(2)
      common    /m5tols/ toldj(3),tolx,tolpiv,tolrow,rowerr,xnorm


*  Nonlinear objective.

      logical            conv,restrt
      common    /m7len / fobj  ,fobj2 ,nnobj ,nnobj0
      common    /m7loc / lgobj ,lgobj2
      common    /m7cg1 / cgbeta,itncg,msgcg,modcg,restrt
      common    /m7cg2 / lcg1,lcg2,lcg3,lcg4,modtcg,nitncg,nsubsp
      common    /m7conv/ etash,etarg,lvltol,nfail,conv(4)
      common    /m7phes/ rgmin1,rgnrm1,rgnrm2,jz1,jz2,labz,nfullz,mfullz
      common    /m7tols/ xtol(2),ftol(2),gtol(2),pinorm,rgnorm,tolrg


*  Nonlinear constraints.

      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /m8loc / lfcon ,lfcon2,lfdif ,lfdif2,lfold ,
     $                   lblslk,lbuslk,lxlam ,lrhs  ,
     $                   lgcon ,lgcon2,lxdif ,lxold
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /m8al2 / radius,rhsmod,modpen,modrhs
      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)
      common    /m8func/ nfcon(4),nfobj(4),nprob,nstat1,nstat2
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave
      common    /m8veri/ jverif(4),lverif(2)


*  Miscellaneous.

      logical            gotbas,gotfac,gothes,gotscl
      common    /cycle1/ gotbas,gotfac,gothes,gotscl
      common    /cycle2/ objtru,suminf,numinf
      common    /cyclcm/ cnvtol,jnew,materr,maxcy,nephnt,nphant,nprint

***** end of file mcommon fortran.
*     ------------------------------------------------------------------

      character*1        ch1
      logical            finish, nlncon, nlnobj

*     For minoss we have to copy m, n, etc into common.

      m      = mxx
      n      = nxx
      nb     = nbxx
      ne     = nexx
      nka    = nkax
      nname  = nnamex
      iobj   = iobjxx

*     Initialize a few things.

      nlncon = nncon  .gt. 0
      nlnobj = nnobj  .gt. 0
      ierr   = 0
      lenl   = 0
      lenu   = 0
      ncycle = 0
      ninf   = 0
      njac0  = max( njac, 1 )
      nstat1 = 1
      nstat2 = 1
      sclobj = 1.0d+0
      jobj   = 0
      if (iobj .gt. 0) jobj = n + iobj
      majitn = 0
      minitn = 0
      vimax  = 0.0d+0
      virel  = 0.0d+0
      maxvi  = 0

*     ------------------------------------------------------------------
*     Decode 'start'.
*     ------------------------------------------------------------------
      gotbas = .false.
      gotfac = .false.
      gothes = .false.
      gotscl = .false.
      ch1    = start(1:1)

      if      (ch1 .eq. 'C'  .or.  ch1 .eq. 'c'  .or.
     $         ch1 .eq. 'B'  .or.  ch1 .eq. 'b') then

*        Cold start  or  Basis file.

         istart = 0
         gotbas = (ioldb + insrt + iload) .gt. 0
         ns     = 0

      else if (ch1 .eq. 'W'  .or.  ch1 .eq. 'w') then

*        Warm start.

         istart = 1
         gotbas = .true.

      else if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') then

*        Hot start.
*        'Hot' is the same as 'Hot FHS'.
*        Look for 'Hot F', 'Hot FH', etc.

         istart = 2
         gotbas = .true.
         nchar  = len( start )
         if (nchar .le. 4) then
            gotfac = .true.
            gothes = .true.
            gotscl = .true.
         else
            do 100 j = 5, nchar
               ch1 = start(j:j)
               if (ch1 .eq. 'F'  .or.  ch1 .eq. 'f') gotfac = .true.
               if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') gothes = .true.
               if (ch1 .eq. 'S'  .or.  ch1 .eq. 's') gotscl = .true.
  100       continue
         end if

         if (iprint .gt. 0) then
            write(iprint, 1020) gotbas, gotfac, gothes, gotscl
         end if
         
      else
         istart = 0
         if (iprint .gt. 0) write(iprint, 1030) start
         if (isumm  .gt. 0) write(isumm , 1030) start
      end if

      nssave = ns

*     ------------------------------------------------------------------
*     1. Fiddle with partial price parameter to avoid foolish values.
*        We reduce nparpr if both the row and column section sizes
*        would be smaller than minprc (= 10 say).
*     2. Change Scale option 1 to 0 if all variables are nonlinear.
*     ------------------------------------------------------------------
      minprc = 10
      npr1   = n / nparpr
      npr2   = m / nparpr
      if (max( npr1, npr2 ) .lt. minprc) then
         maxmn  = max( m, n )
         nparpr = maxmn / min( maxmn, minprc )
         npr1   = n / nparpr
         npr2   = m / nparpr
      end if

      if (lscale .eq. 1  .and.  nn .eq. n) lscale = 0

      if (iprint .gt. 0) write(iprint, 1100) lscale, nparpr, npr1, npr2
      if (isumm  .gt. 0) write(isumm , 1110) lscale, nparpr

*     ------------------------------------------------------------------
*     Set the vector of row types and print the matrix statistics.
*     03 Jun 2003: Do this BEFORE m8augl(-2) to avoid random Jacobian.
*                  (Thanks to Michael Friedlander.)
*     ------------------------------------------------------------------
      call m2amat( 1, m, n, nb,
     $             ne, nka, a, ha, ka,
     $             bl, bu, z(lhrtyp) )

*     ------------------------------------------------------------------
*     Save constant Jacobian elements from A in gcon2,
*     then load random Jacobian elements into A.
*     ------------------------------------------------------------------
      if (nlncon) then
         call m8augl( -2, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )
      end if

*     ------------------------------------------------------------------
*     Input a basis file if one exists, thereby defining hs and xn.
*     (Otherwise, m2crsh will be called later to define hs.)
*     At this stage, ncycle = 0.
*     ------------------------------------------------------------------
      if (iprint .gt. 0) then
         write(iprint, 1200)
         if (istart .eq. 0) then
            if (.not. gotbas) write(iprint, 1210)
         else
            write(iprint, 1220) start
         end if
      end if

      call m4getb( ncycle, istart, m, mbs, n, nb, nn, nname, nscl,
     $             lcrash, ns,
     $             ne, nka, a, ha, ka,
     $             z(lhrtyp), hs, z(lkb), z(lascal), bl, bu,
     $             pi, xn, z(ly), z(ly2), name1, name2,
     $             z, nwcore )
      if (ierr .ne. 0) go to 900


*     ------------------------------------------------------------------
*                             CYCLE  PROCEDURE
*
*     The following notes are relevant if Cycle limit = 2 or more.
*
*  1. Scaling and/or Crash are controlled on each cycle by the following
*     logical variables:
*
*     If gotscl is true, scales are retained from the previous cycle.
*                        Otherwise, scales are recomputed (if lscale>0).
*
*     If gotbas is true, the basis is retained.  Otherwise, Crash is
*                        called.
*
*
*  2. When m5solv is called, Flying Starts are controlled by the
*     following logical variables:
*
*     If gotfac is true, an LU factorization of the basis is assumed
*                        to be present.  (Ignored if there are any
*                        nonlinear constraints.)
*
*     If gothes is true, z(lr) is assumed to contain a useful
*                        reduced-Hessian approximation.
*
*
*  3. For the next cycle,
*        m4getb sets gotscl and gotbas to be true, and
*        m5solv sets gotfac and gothes to current values (usually true).
*     These values will often be appropriate.  However, the expert user
*     of matmod must set some or all of the logical variables to .false.
*     if the problem data or state have been significantly altered.
*
*     For example, if the Jacobian was used by the scaling routine
*     (Scale option 2) and if the Jacobian could be rather different
*     from its value at the start of the previous cycle, it may be
*     advisable to request new scales by setting gotscl = .false.
*
*     Similarly, if matmod alters some matrix elements in columns that
*     are currently basic, one should set gotfac = .false. to force
*     refactorization.  In particular, if the linear objective row c is
*     altered, gotfac should be set to .false., since c is part of the
*     LU factors.
*     ------------------------------------------------------------------

      finish = .false.
      jnew   = n - nphant
      materr = 0
      nprntd = 0
      nsolvd = ncycle

*     If Cycle limit is more than 1, call matmod with ncycle = 0 in case
*     the user wants a chance to set things up before any solves.

      if (maxcy .gt. 1) then
         if (iprint .gt. 0) write(iprint, 3000) ncycle
         if (isumm  .gt. 0) write(isumm , 3000) ncycle
      
         call matmod( ncycle, nprob, finish,
     $                m, n, nb, ne, nka, ns, nscl, nname,
     $                a, ha, ka, bl, bu,
     $                z(lascal), hs, name1, name2,
     $                xn, pi, rc, z, nwcore )
         if (finish) go to 800
      end if

*     ==================================================================
*     Start of the Cycle loop.
*     ==================================================================
      do 600 kcycle = 1, maxcy
         ncycle = kcycle
         nsolvd = kcycle
         ierr   = 0
         if (.not. gotbas) then
            gotfac = .false.
            gothes = .false.
         end if
         if (ncycle .ge. 2) then
            call m1page( 1 )
            if (iprint .gt. 0) then
               write(iprint, 2000) ncycle
               write(iprint, 1020) gotbas, gotfac, gothes, gotscl
            end if
            if (isumm  .gt. 0) then
               write(isumm , 2000) ncycle
            end if
         end if
      
*        Make sure the Jacobian variables are inside their bounds.
      
         if (nlncon) then
            call m8augl( 2, m, n, nb, ns, inform,
     $                   ne, nka, a, ha, ka,
     $                   hs, bl, bu, xn, z, nwcore )
         end if
      
*        For the first cycle, the row types have been set by m2amat.
*        Reset them for later cycles in case m2scal or m2crsh are
*        called.
      
         if (ncycle .ge. 2) then
            call m2amat( 2, m, n, nb,
     $                   ne, nka, a, ha, ka,
     $                   bl, bu, z(lhrtyp) )
         end if
      
*        ---------------------------------------------------------------
*        Compute scales from a, bl, bu (except if gotscl is true).
*        Scale a, bl, bu, xn, pi and fcon.
*        Initialize xlam from pi.
*        Call CRASH if a basis file was not supplied
*        (or if gotbas is false).
*        ---------------------------------------------------------------
         call m4getb( ncycle, istart, m, mbs, n, nb, nn, nname, nscl,
     $                lcrash, ns,
     $                ne, nka, a, ha, ka,
     $                z(lhrtyp), hs, z(lkb), z(lascal), bl, bu,
     $                pi, xn, z(ly), z(ly2), name1, name2,
     $                z, nwcore )
         if (ierr .ne. 0) go to 900
      
*        1. Set ns to match hs(*).
*        2. Set kb(m+1) thru kb(m+ns) to define the initial set of
*           superbasics, except if a Hot start
*           (gotbas and gothes are both true).
*        3. Check that nonbasic xn are within bounds.
      
         call m4chek( m, maxs, mbs, n, nb, ns,
     $                hs, z(lkb), bl, bu, xn )
      
*        ---------------------------------------------------------------
*        Solve the current problem.
*        Bail out if there is a fatal error.
*        ---------------------------------------------------------------
         call m1page( 1 )
         if (iprint .gt. 0) write(iprint, 2100)
      
         call m1time( 2,0 )
         call m5solv( lcrash, ns, objadd,
     $                m, maxr, maxs, mbs, 
     $                n, nb, nn, nr, nscl, nx,
     $                nncon, nnjac, njac, nnobj,
     $                nn0, nncon0, njac0, nnobj0, 
     $                ne, nka, a, ha, ka,
     $                z(lhrtyp), hs, z(lkb), z(lascal), bl, bu,
     $                z(lbbl), z(lbbu), fsub, z(lgsub),
     $                z(lgrd), z(lgrd2),
     $                pi, z(lr), rc, z(lrg), z(lrg2),
     $                z(lx), xn, z(ly), z(ly2),
     $                z(lfcon), z(lfcon2), z(lfold),
     $                z(lgcon), z(lgcon2), z(lgobj), z(lgobj2),
     $                z(lxlam), z(lrhs  ), z(lxdif), z(lxold),
     $                z, nwcore,
     $                ifuser, m1user )
         call m1time(-2,0 )

         majitn = majitn + majits
         minitn = minitn + minits
      
         if (ierr .ge. 30                   ) go to 900
         if (ierr .ge. 20  .and.  itn .eq. 0) go to 900
      
*        ---------------------------------------------------------------
*        Unscale, compute nonlinear constraint violations,
*        save basis files and prepare to print the solution.
*        Clock 3 is "Output time".
*        ---------------------------------------------------------------
         call m1time( 3,0 )
         call m4savb( 1, m, mbs, n, nb, nn, nname, nscl,
     $                msoln, ns,
     $                ne, nka, a, ha, ka,
     $                hs, z(lkb), z(lascal), bl, bu,
     $                name1, name2,
     $                pi, rc, xn, z(ly), z, nwcore )
      
*        In some Cycling applications, it may be desirable to suppress
*        the printing of intermediate solutions.  Otherwise if mode = 2,
*        m4savb prints the solution under the control of msoln
*        (which is set by the  Solution  keyword in the SPECS file).
*        The printed solution may or may not be wanted, as follows:
*     
*        msoln  = 0   means      No
*               = 1   means      If optimal, infeasible or unbounded
*               = 2   means      Yes
*               = 3   means      If error condition
*     
*        This call normally prints the solution when there is no
*        Cycling, because the default values are  maxcy = nprint = 1.
      
         if (ncycle .gt. maxcy - nprint) then
            nprntd = nsolvd
            call m4savb( 2, m, mbs, n, nb, nn, nname, nscl,
     $                msoln, ns,
     $                ne, nka, a, ha, ka,
     $                hs, z(lkb), z(lascal), bl, bu,
     $                name1, name2,
     $                pi, rc, xn, z(ly), z, nwcore )
         end if
         call m1time(-3,0 )
      
*        ---------------------------------------------------------------
*        Call the functions one last time with  nstate .ge. 2.
*        We have to disable scaling.
*        mode = 0  tells the functions that gradients are not required.
*        03 Nov 2000: Do this only if ninf = 0!
*        ---------------------------------------------------------------
         if (ierr .eq. 6) go to 800
         if (ninf .eq. 0) then
            lssave = lscale
            lscale = 0
            nstat1 = 2 + ierr
            nstat2 = nstat1
            mode   = 0
            if (nlncon) then
               call m6fcon( mode, nncon, nnjac, njac,
     $                      z(lfcon), z(lgcon2),
     $                      ne, nka, ha, ka, xn, z, nwcore )
            end if
            if (nlnobj) then
               call m6fobj( mode, nnobj,
     $                      fobj, z(lgobj2), xn, z, nwcore )
            end if
            lscale = lssave
            if (mode .lt. 0) go to 800
            nstat1 = 0
            nstat2 = 0
         end if
         
*        Terminate Cycles if m5solv gave a serious error.

         if (ierr .ge. 20) go to 800

*        ---------------------------------------------------------------
*        Let the user modify the problem for the next Cycle.
*        ---------------------------------------------------------------
         if (ncycle .lt. maxcy) then
            if (iprint .gt.  0) write(iprint, 3000) ncycle
            if (isumm  .gt.  0) write(isumm , 3000) ncycle
      
            call matmod( ncycle, nprob, finish,
     $                   m, n, nb, ne, nka, ns, nscl, nname,
     $                   a, ha, ka, bl, bu,
     $                   z(lascal), hs, name1, name2,
     $                   xn, pi, rc, z, nwcore )
            if (finish) go to 800
         end if
  600 continue
*     ==================================================================
*     End of the Cycle loop.
*     ==================================================================


*     Print the final solution if it has not already been printed.

  800 if (nprntd .ne. nsolvd) then
         call m1time( 3,0 )
         call m4savb( 2, m, mbs, n, nb, nn, nname, nscl,
     $                msoln, ns,
     $                ne, nka, a, ha, ka,
     $                hs, z(lkb), z(lascal), bl, bu,
     $                name1, name2,
     $                pi, rc, xn, z(ly), z, nwcore )
         call m1time(-3,0 )
      end if
*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------

  900 inform = ierr
      return

 1020 format(/ ' gotbas =', l2, 4x, ' gotfac =', l2, 4x,
     $         ' gothes =', l2, 4x, ' gotscl =', l2)
 1030 format(/ ' XXX Start parameter not recognized:  ', a)
 1100 format(/ ' Scale option', i3, ',      Partial price', i8
     $       / ' Partial price section size (A) ', i12
     $       / ' Partial price section size (I) ', i12)
 1110 format(/ ' Scale option', i3, ',    Partial price', i4)
 1200 format(/ ' Initial basis' / ' -------------')
 1210 format(  ' No basis file supplied')
 1220 format(  ' start = ', a)
 2000 format(/ ' Start of Cycle', i5,
     $       / ' -------------------')
 2100 format(  ' Iterations' / ' ----------')
 3000 format(/ ' matmod called with ncycle =', i5)

*     end of misolv
      end
