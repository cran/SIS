************************************************************************
*
*     File  mi26bfac fortran.
*
*     m2bmap   m2belm   m2newB   m2bsol   m2sing
*
* 12 Jul 2000: mi26bfac.f contains most basis factor routines
*              (except m2bfac and LUSOL).
* 17 Nov 2001: m2sing uses LUSOL lprint.
* 27 Nov 2002: In m2bsol when big elements of U occur,
*              TPP switches to TRP with Lmax =  10.0
*              (previously  to TCP with Lmax = 100.0).
*              Apparently it never switches back to TPP.
*              Also, BSfac is done with current TPP, TRP or TCP
*              (not always TRP, which would probably be better)
*              but with tight Ltol that does revert to previous value.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2bmap( mode, m, n, ne, minz, maxz, nguess )

*     ------------------------------------------------------------------
*     m2bmap sets up the core allocation for the basis factors.
*     It is called by m2core.
*     
*     Normally the storage is for B = LU.
*     For nonlinear problems, we may also need to factorize (B S)' = LU,
*     where S has at most maxs columns.
*
*     29 Oct 1993: Generalized to allow room for (B S)^T.
*     08 Mar 2002: MORE storage needed (minA) for BS factors.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m1word/ nwordr,nwordi,nwordh
      common    /m2lu1 / minlu,maxlu,lena,nbelem,ip,iq,lenc,lenr,
     $                   locc,locr,iploc,iqloc,lua,indc,indr
      common    /m5len / maxr  ,maxs  ,mbs   ,nn    ,nn0   ,nr    ,nx


*     Allocate arrays for an  ms x m  matrix.
*     We need ms bigger than m only for nonlinear problems.
*     08 Mar 2002: Just in case, might as well use m + maxs always.

!!    if (nn .eq. 0) then
!!       ms     = m
!!    else
!!       ms     = m + maxs
!!    end if
      ms     = m + maxs

      minlu  = minz
      maxlu  = maxz
      mh     = (ms - 1)/nwordh + 1
      mi     = (ms - 1)/nwordi + 1
      nh     = (m  - 1)/nwordh + 1
      ni     = (m  - 1)/nwordi + 1
      ip     = minlu
      iq     = ip     + mh
      lenc   = iq     + nh
      lenr   = lenc   + nh
      locc   = lenr   + mh
      locr   = locc   + ni
      iploc  = locr   + mi
      iqloc  = iploc  + nh
      lua    = iqloc  + mh
      lena   = (maxlu - lua - 1)*nwordh/(nwordh + 2)
      indc   = lua    + lena
      indr   = indc   + (lena - 1)/nwordh + 1

*     Estimate the number of nonzeros in the basis factorization.
*     necola = estimate of nonzeros per column of  a.
*     We guess that the density of the basis factorization is
*     2 times as great, and then allow 1 more such lot for elbow room.
*     18 sep 1989: Tony and Alex change m to min( m, n ) below.
*     14 Aug 2000: Double the estimate.  Memory is much cheaper now!!
*     08 Mar 2002: Allow for tiny m but big m + maxs.
*     30 Mar 2003: Erwin warns of overflow with tiny n but big m and ms.

      necola = max( (ne/n), 10 )
!!!   minA   = 6 * min( m, n ) * necolA  ! Too little for BSfac
!!!!  minA   = 6 *      ms     * necolA  ! ms is biggest LU dimension
      minA   = 6 * min( ms,n ) * necolA  ! ms is biggest LU dimension
      minA   = minA + m + n + 10000      ! So tiny problems have plenty
      nguess = lua + minA + 2*minA/nwordh
      if (mode .ge. 3) then
         if (iprint .gt. 0) write(iprint, 1000) lena
      end if
      return

 1000 format(/ ' Nonzeros allowed for in LU factors', i9)

*     end of m2bmap
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2belm( modeLU, ms, m, n, nz,
     $                   ne, nka, a, ha, ka, kb,
     $                   alu, indc, indr, ip, lena )

      implicit           double precision (a-h,o-z)
      character*2        modeLU
      double precision   a(ne), alu(lena)
      integer            ha(ne), indc(lena), indr(lena), ip(ms)
      integer            ka(nka), kb(ms)

*     ------------------------------------------------------------------
*     m2belm  extracts the basis elements from the constraint matrix,
*     ready for use by the LU factorization routines.
*
*     If modeLU = 'B ', we extract B (the normal case).
*     If modeLU = 'BS', we extract (B S) and transpose it:  (B S)'.
*
*     29 Oct 1993: modeLU options implemented.  
*                  nz is returned to m2bfac (as nbelem),
*                  so that nbelem is defined in m2bsol for 'BS' mode.
*     23 Apr 1994: iobj needed to ensure that the objective slack
*                  remains basic during a 'BS' factorize.
*                  This is only true if slack rows are kept.
*     09 Apr 1996: Went back to excluding slack rows in 'BS' factorize.
*                  iobj is no longer needed.
*     ------------------------------------------------------------------

      parameter        ( one = 1.0d+0 )

      if (modeLU .eq. 'B ') then
*        ---------------------------------------------------------------
*        Normal case.
*        ---------------------------------------------------------------
         nz     = 0
         do 200 k = 1, m
            j     = kb(k)
            if (j .eq. 0) go to 200
            if (j .le. n) then
               do 150  i   = ka(j), ka(j+1)-1
                  ir       = ha(i)
                  nz       = nz + 1
                  alu(nz)  = a(i)
                  indc(nz) = ir
                  indr(nz) = k
  150          continue
            else

*              Treat slacks specially.

               nz       = nz + 1
               alu(nz)  = one
               indc(nz) = j - n
               indr(nz) = k
            end if
  200    continue

      else if (modeLU .eq. 'BS') then
*        ---------------------------------------------------------------
*        Extract (B S)'.
*        ip is needed for workspace.
*        ip(i) = 0 except for rows containing a basic slack.
*        We can ignore all of these rows except for the slack itself.
*        01 Mar 1994: Try keeping them in anyway.
*        23 Apr 1994: Row iobj must be treated specially to ensure
*                     that its slack gets into the basis.
*                     We should really do the same for all free rows.
*        09 Apr 1996: Went back to excluding slack rows.
*                     iobj is no longer needed.
*        ---------------------------------------------------------------
         call hload ( m, 0, ip, 1 )
         do 300 k = 1, ms
            j     = kb(k)
            if (j .gt. n) ip(j-n) = 1
  300    continue

         nz     = 0
         do 400 k = 1, ms
            j     = kb(k)
            if (j .le. n) then
               do 350  i   = ka(j), ka(j+1)-1
                  ir       = ha(i)
                  if (ip(ir) .eq. 0) then
                     nz       = nz + 1
                     alu(nz)  = a(i)
                     indc(nz) = k
                     indr(nz) = ir
                  end if
  350          continue
            else

*              Treat slacks specially.

               nz       = nz + 1
               alu(nz)  = one
               indc(nz) = k
               indr(nz) = j - n
            end if
  400    continue
      end if

*     end of m2belm
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2newB( ms, m, nb, hs, ip, kb, kbsold, locr, nswap )

      implicit           double precision (a-h,o-z)
      integer            hs(nb), ip(ms)
      integer            kb(ms), kbsold(ms), locr(ms)

*     ------------------------------------------------------------------
*     m2newB  permutes kb(*) to reflect the permutation (B S)P,
*     where P is in ip(*).  It updates hs(*) accordingly.
*     kbsold(*) and locr(*) are needed for workspace.
*
*     30 Oct 1993: First version.
*     04 Nov 1993: kbsold, nswap used to save old R if there's no
*                  change in the set of superbasics.
*     ------------------------------------------------------------------
      nswap = 0
      m1    = m  + 1
      ns    = ms - m
      call icopy ( ms, kb    , 1, locr      , 1 )
      call icopy ( ns, kb(m1), 1, kbsold(m1), 1 )

      do 100 k = 1, ms
         i        = ip(k)
         j        = locr(i)
         kb(k)    = j
         if (k .le. m) then
            hs(j) = 3
         else
            if (hs(j) .ne. 2) nswap = nswap + 1
            hs(j) = 2
         end if
  100 continue

*     Restore the old S ordering if S contains the same variables.

      if (nswap .eq. 0) then
         call icopy ( ns, kbsold(m1), 1, kb(m1), 1 )
      end if

*     end of m2newB
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2bsol( mode, m, w, y, z, nwcore )

      implicit           double precision (a-h,o-z)
      double precision   w(m), y(m), z(nwcore)

*     ------------------------------------------------------------------
*     m2bsol  calls up the relevant basis-factorization routines.
*
*     mode
*     ----
*      0    Factorize current basis from scratch, so that  B = L*U.
*      1    Solve  L*w = w(input).  y is not touched.
*      2    Solve  L*w = w(input)  and solve  B*y = w(input).
*      3    Solve  B(transpose)*y = w.  Note that w is destroyed.
*      4    Update the LU factors when the jp-th column is replaced
*           by a vector v.  jp is in common block /m5log3/.
*           On input, w must satisfy L*w = v.  w will be destroyed.
*      5    Solve  L(transpose)*w = w(input).  y is not touched.
*      8    Factorize transpose of (B S), so that  (B') = L*U,
*                                                  (S')
*           without saving L and U.  Get a new partition of (B S).
*
*     The following tolerances are used...
*
*     luparm(3) = maxcol   lu1fac: maximum number of columns
*                          searched allowed in a Markowitz-type
*                          search for the next pivot element.
*     luparm(6) = 0, 1, 2  for TPP, TRP, TCP.
*     luparm(8) = keepLU   lu1fac: keepLU = 1 means keep L and U,
*                                           0 means discard them.
*     parmlu(1) = Lmax1  = maximum multiplier allowed in  L  during
*                          refactorization.
*     parmlu(2) = Lmax2  = maximum multiplier allowed during updates.
*     parmlu(3) = small  = minimum element kept in  B  or in
*                          transformed matrix during elimination.
*     parmlu(4) = utol1  = abs tol for flagging small diagonals of  U.
*     parmlu(5) = utol2  = rel tol for flagging small diagonals of  U.
*     parmlu(6) = uspace = factor allowing waste space in row/col lists.
*     parmlu(7) = dens1  = the density at which the Markowitz strategy
*                          should search maxcol columns and no rows.
*     parmlu(8) = dens2  = the density at which the Markowitz strategy
*                          should search only 1 column.
*                          (In one version of lu1fac, the remaining
*                          matrix is treated as dense if there is
*                          sufficient storage.)
*
*     25 Nov 1991: parmlu(1,2,4,5,8) are now defined by the SPECS file
*                  via  LU Factorization tol
*                       LU Update        tol
*                       LU Singularity   tol
*                       LU Singularity   tol
*                       LU Density       tol
*                  respectively.
*     12 Jun 1992: Decided lu1fac's switch to dense LU was giving
*                  trouble.  Went back to treating all of LU as sparse.
*     26 Oct 1993: mode 8 implemented.
*     27 Feb 1994: Test for excessive growth in U.  Reduce LU Factor Tol
*                  if it is not already near one.
*                  Prompted by problem "maros-r7" in netlib/lp/data.
*                  LU Factor tol = 100.0 gives Umax = 1.0e+10 or worse.
*                                      Many singularities, then failure.
*                  LU Factor tol =  10.0 gives a clean run.
*                  inform = ierr = 2 tells m2bfac to try again.
*     21 Aug 2000: CUTE Problem LISWET3 revealed that the BS Factorize
*                  must use LU Factor Tol < 2.0 for stability, not = 2.0 !!!
*                  Changed it to 1.9.
*     13 Sep 2000: CUTE Problem SVANBERG still has trouble even with 1.1 !!!
*                  LU Complete Pivoting invoked if big growth in U.
*     22 Sep 2000: Print of LU statistics moved into lu1fac.
*     25 Oct 2000: Fixed bug: "growth" was not initialized.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2lu1 / minlu,maxlu,lena,nbelem,ip,iq,lenc,lenr,
     $                   locc,locr,iploc,iqloc,lua,indc,indr
      common    /m2lu2 / factol(5),lamin,nsing1,nsing2
      common    /m2lu3 / lenl,lenu,ncp,lrow,lcol
      common    /m2lu4 / parmlu(30),luparm(30)
      common    /m5freq/ kchk,kinv,ksav,klog,ksumm,i1freq,i2freq,msoln
      common    /m5loc / lpi   ,lpi2  ,lw    ,lw2   ,
     $                   lx    ,lx2   ,ly    ,ly2   ,
     $                   lgsub ,lgsub2,lgrd  ,lgrd2 ,
     $                   lr    ,lrg   ,lrg2  ,lxn
      common    /m5log1/ idebug,ierr,lprint
      common    /m5log3/ djq,theta,pivot,cond,nonopt,jp,jq,modr1,modr2
      common    /m5lp2 / invrq,invitn,invmod
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave

      logical            TPP

      if (mode .eq. 0) then
*        ---------------------------------------------------------------
*        mode = 0.    Factorize the basis.
*        Note that luparm(6) and parmlu(1,2,4,5,8) are defined
*        by the SPECS file options.
*        ---------------------------------------------------------------
         luparm(1) = iprint                     ! nout
         kprint    = mod(lprint,10)             ! lprint
         luparm(2) = 0                          ! errors only
         if (iprint .le.  0) luparm(2) = -1     ! no output
         if (kprint .gt.  0) luparm(2) = 10     ! errors + statistics
         if (idebug .eq. 50) luparm(2) = 50     ! + debug info
         luparm(3) = 5                          ! maxcol
         if (i1freq .gt.  0) luparm(3) = i1freq ! can set in SPECS file
!!       luparm(6) = 0 or 1 or 2                ! TPP, TRP or TCP
         luparm(8) = 1                          ! 1 = keepLU

         parmlu(3) = eps0                       ! small (drop tolerance)
!!       parmlu(4)                              ! utol1
!!       parmlu(5)                              ! utol2
         parmlu(6) = 3.0d+0                     ! uspace
         parmlu(7) = 0.3d+0                     ! dens1
!!       parmlu(8) = 0.3d+0                     ! dens2

         call lu1fac( m,  m   , nbelem  , lena   , luparm , parmlu,
     $                z(lua  ), z(indc ), z(indr), z(ip)  , z(iq) ,
     $                z(lenc ), z(lenr ), z(locc), z(locr),
     $                z(iploc), z(iqloc), z(ly  ), z(ly2 ), w, inform )

         ierr   = inform
         TPP    = luparm(6) .eq. 0
         lamin  = luparm(13)
         lenl   = luparm(23)
         lenu   = luparm(24)
         lrow   = luparm(25)
         ncp    = luparm(26)
         growth = parmlu(16)

*        Test for excessive growth in U.
*        Reduce LU Factor tol and LU Update tol if necessary.
*        (Default values are 100.0 and 10.0)

         if (inform .eq. 0  .and.  growth .ge. 1.0d+8) then
            elmax1  = parmlu(1)
            elmax2  = parmlu(2)

            if ( TPP ) then    ! Switch from TPP to TRP (no longer TCP).
                               ! Lmax = 10.0 is ok with TRP.
                luparm(6) = 1
                elmax1    = min( elmax1, 10.0d+0 )
                parmlu(1) = elmax1
                inform    = 2
                if (iprint .gt. 0) write(iprint, 1005) elmax1
                if (isumm  .gt. 0) write(isumm , 1005) elmax1
                
            else if (elmax1 .ge. 2.0d+0) then
                elmax1    = sqrt( elmax1 )
                parmlu(1) = elmax1
                inform    = 2
                if (iprint .gt. 0) write(iprint, 1010) elmax1
                if (isumm  .gt. 0) write(isumm , 1010) elmax1
            end if

            if (elmax2 .gt. elmax1) then ! Keep Update tol <= Factor tol
                parmlu(2) = elmax1
                if (iprint .gt. 0) write(iprint, 1020) elmax1
                if (isumm  .gt. 0) write(isumm , 1020) elmax1
            end if
         end if            

      else if (mode .le. 2) then
*        ---------------------------------------------------------------
*        mode = 1 or 2.    Solve   L*w = w(input).
*        When LU*y = w is being solved in MINOS, norm(w) will sometimes
*        be small (e.g. after periodic refactorization).  Hence for
*        mode 2 we scale parmlu(3) to alter what lu6sol thinks is small.
*        ---------------------------------------------------------------
         small  = parmlu(3)
         if (mode .eq. 2) parmlu(3) = small * dnormi( m, w, 1 )

         call lu6sol( 1, m, m, w, y, lena, luparm, parmlu,
     $                z(lua ), z(indc), z(indr), z(ip), z(iq),
     $                z(lenc), z(lenr), z(locc), z(locr), inform )

         parmlu(3) = small

         if (mode .eq. 2) then
*           ------------------------------------------------------------
*           mode = 2.    Solve  U*y = w.
*           ------------------------------------------------------------
            call lu6sol( 3, m, m, w, y, lena, luparm, parmlu,
     $                   z(lua ), z(indc), z(indr), z(ip), z(iq),
     $                   z(lenc), z(lenr), z(locc), z(locr), inform )
         end if

      else if (mode .eq. 3) then
*        ---------------------------------------------------------------
*        mode = 3.    Solve  B(transpose)*y = w.
*        ---------------------------------------------------------------
         call lu6sol( 6, m, m, y, w, lena, luparm, parmlu,
     $                z(lua ), z(indc), z(indr), z(ip), z(iq),
     $                z(lenc), z(lenr), z(locc), z(locr), inform )

      else if (mode .eq. 4) then
*        ---------------------------------------------------------------
*        mode = 4.    Update the LU factors of  B  after basis change.
*        ---------------------------------------------------------------
         invmod = invmod + 1
         call lu8rpc( 1, 2, m, m, jp, w, w,
     $                lena, luparm, parmlu,
     $                z(lua ), z(indc), z(indr), z(ip), z(iq),
     $                z(lenc), z(lenr), z(locc), z(locr),
     $                inform, diag, wnorm )
         if (inform .ne. 0) invrq = 7
         lenl   = luparm(23)
         lenu   = luparm(24)
         lrow   = luparm(25)
         ncp    = luparm(26)

      else if (mode .eq. 8) then
*        ---------------------------------------------------------------
*        mode = 8.    Factorize (B S)' = LU without keeping L and U.
*        ---------------------------------------------------------------
         luparm(1) = iprint                     ! nout
         kprint    = mod(lprint,10)             ! lprint
         luparm(2) = 0                          ! errors only
         if (iprint .le.  0) luparm(2) = -1     ! no output
         if (kprint .gt.  0) luparm(2) = 10     ! errors + statistics
         if (idebug .eq. 50) luparm(2) = 50     ! + debug info
         luparm(3) = 5                          ! maxcol
         luparm(8) = 0                          ! Don't keepLU

         if (iprint .gt. 0  .and.  kprint .gt. 0) then
            write(iprint, 1800)
         end if

*        Save tolfac (the existing LU Factor tol) and set it to a small
*        value for this LU, to give a good (B S) partitioning.

         tolfac    = parmlu(1)
         parmlu(1) = 1.9d+0              ! MUST BE LESS THAN 2.0!!!!
         parmlu(3) = eps0
         parmlu(6) = 3.0d+0
         parmlu(7) = 0.3d+0
         ns        = nssave
         ms        = m + ns

         call lu1fac( ms, m   , nbelem  , lena   , luparm , parmlu,
     $                z(lua  ), z(indc ), z(indr), z(ip)  , z(iq) ,
     $                z(lenc ), z(lenr ), z(locc), z(locr),
     $                z(iploc), z(iqloc), z(ly  ), z(ly2 ), w, inform )

         ierr   = inform
         lamin  = luparm(13)
         parmlu(1) = tolfac

         if (iprint .gt. 0  .and.  kprint .gt. 0) then
            write(iprint, 1800)
         end if
      end if

      return

 1005 format(/ ' LU Rook Pivoting invoked.  LU Factor tol =', f10.2)
 1010 format(/ ' LU Factor tol reduced to', f10.2)
 1020 format(/ ' LU Update tol reduced to', f10.2)
 1800 format(' ')

      end ! subroutine m2bsol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2sing( lprint, m, n, nb,
     $                   w, ip, iq, bl, bu, hs, kb, xn )

      implicit           double precision (a-h,o-z)
      integer            ip(m), iq(m), hs(nb)
      integer            kb(m)
      double precision   bl(nb), bu(nb), w(m), xn(nb)

*     -----------------------------------------------------------------
*     m2sing  is called if the LU factorization of the basis appears
*     to be singular.   If  w(j)  is not positive, the  j-th  basic
*     variable  kb(j)  is replaced by the appropriate slack.
*     If any kb(j) = 0, only a partial basis was supplied.
*
*     08 Apr 1992: Now generate internal values for hs(j) if necessary
*                  (-1 or 4) to be compatible with m5hs.
*     17 Nov 2001: Output only if lprint >= 10.
*     -----------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm

      parameter        ( zero = 0.0d+0,  nprint = 5 )

      nsing  = 0
      do 100 k  = 1, m
         j      = iq(k)
         if (w(j) .gt. zero) go to 100
         j      = kb(j)
         if (j    .gt.  0  ) then

*           Make variable  j  nonbasic (and feasible).
*           hs(j) = -1 means xn(j) is strictly between its bounds.

            if      (xn(j) .le. bl(j)) then
               xn(j) =  bl(j)
               hs(j) =  0
            else if (xn(j) .ge. bu(j)) then
               xn(j) =  bu(j)
               hs(j) =  1
            else
               hs(j) = -1
            end if

            if (bl(j) .eq. bu(j)) hs(j) = 4
         end if

*        Make the appropriate slack basic.

         i       = ip(k)
         hs(n+i) = 3
         nsing   = nsing + 1
         if (lprint .ge. 10  .and.  nsing .le. nprint) then
            if (iprint .gt. 0) write(iprint, 1000) j, i
            if (isumm  .gt. 0) write(isumm , 1000) j, i
         end if
  100 continue

      if (lprint .ge. 10  .and.  nsing .gt. nprint) then
         if (iprint .gt. 0) write(iprint, 1100) nsing
         if (isumm  .gt. 0) write(isumm , 1100) nsing
      end if
      return

 1000 format(' Column', i7, '  replaced by slack', i7)
 1100 format(' and so on.  Total slacks inserted =', i6)

      end ! subroutine m2sing
