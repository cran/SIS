************************************************************************
*
*     File  mi25bfac fortran.
*
*     m2bfac
*
*     12 Jul 2000: mi25bfac.f separates m2bfac from most factor routines.
*     17 Nov 2001: m2sing uses LUSOL lprint.
*     03 Apr 2003: That caused a bug.  If m2sing was needed,
*                  the Print level lprint got overwritten.
*     28 May 2004: Don't let the very first LU be badly conditioned.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2bfac( modeLU, gotfac, nfac, nswap, 
     $                   m, mbs, n, nb, nr, nn, ns,
     $                   lcrash, fsub, objadd,
     $                   ne, nka, a, ha, ka,
     $                   kb, hs, bl, bu, bbl, bbu,
     $                   r, w, x, xn, y, y2, z, nwcore )

      implicit           double precision (a-h,o-z)
      character*2        modeLU
      logical            gotfac
      integer            ha(ne), hs(nb)
      integer            ka(nka),kb(mbs)
      double precision   a(ne),  bl(nb), bu(nb), bbl(mbs), bbu(mbs),
     $                   r(nr),  w(m),   x(mbs), xn(nb),
     $                   y(mbs), y2(m),  z(nwcore)

*     ------------------------------------------------------------------
*     m2bfac  computes a factorization of the current basis, B.
*
*     If modeLU = 'B ', the usual B = LU is computed.
*     If modeLU = 'BS', there are some superbasics and we want to
*                       choose a good basis from the columns of (B S).
*                       We first factorize (B S)' to obtain a new B.
*                       Then B = LU is computed as usual.
*     If modeLU = 'BT', we should TRY 'B ' first and go back to 'BS'
*                       only if B seems ill-conditioned.
*
*     gotfac  must be false the first time m2bfac is called.
*     It might be true after the first cycle.
*
*     If lcrash = 3, linear inequality rows (LG rows) are to be
*     treated as free rows.
*
*     04 May 1992: invmod and invitn not reset if gotfac is true.
*                  This will help m5solv keep track of LU properly
*                  during later cycles.
*     04 Jun 1992: lcrash added.
*     29 Oct 1993: modeLU options 'B ' and 'BS' implemented.
*                  nswap returns the number of (B S) changes.
*     28 Feb 1994: Retry with reduced LU Factor tol
*                  if m2bsol says there was large growth in U.
*     20 Mar 1994: 'BT' option implemented to save R more often each
*                  major iteration.
*     09 Apr 1996: kobj now used to keep track of x(iobj) in B.
*     14 Jul 1998: (m*m) overflow fixed in forming bdnsty (Art Cohen).
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2lu1 / minlu,maxlu,lena,nbelem,ip,iq,lenc,lenr,
     $                   locc,locr,iploc,iqloc,lua,indc,indr
      common    /m2lu4 / parmlu(30),luparm(30)
      common    /m2mapz/ maxw,maxz
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      common    /m5lp1 / itn,itnlim,nphs,kmodlu,kmodpi
      common    /m5lp2 / invrq,invitn,invmod
      common    /m5prc / nparpr,nmulpr,kprc,newsb

      character*10       cobj
      logical            BS, BT, modtol, prnt
      integer            nBfac
      double precision   Umin
      save               nBfac, Umin
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------

*     Initialize Umin and nBfac on first entry.
*     nBfac  counts consecutive B factorizations (reset if BS is done).
*     Umin   is the smallest diagonal of U after last BS factor.

      if (nfac .eq. 0) then
         Umin  = zero
         nBfac = 0
      end if

      nfac   = nfac  + 1
      nBfac  = nBfac + 1
      nswap  = 0
      ntry   = 0
      maxtry = 10
      ms     = m + ns

      obj    = sinf
      if (ninf .eq. 0) obj = minimz * fsub  +  objadd
      prnt   = iprint .gt. 0  .and.  mod(lprint,10) .gt. 0
      if (prnt) write(iprint, 1000) nfac, invrq, itn, ninf

      if (gotfac  .and.  invrq .eq. 0) go to 500

*     ------------------------------------------------------------------
*     Set local logicals to select the required type of LU.
*     We come back to 100 if a BT factorize looks doubtful.
*     If BT was requested but we haven't done BS yet,
*     might as well do BS now.
*     ------------------------------------------------------------------
      BT     =  modeLU .eq. 'BT'  .and.  ns   .gt. 0
      BS     = (modeLU .eq. 'BS'  .and.  ns   .gt. 0   )  .or.
     $         (       BT         .and.  Umin .eq. zero)

  100 if ( BS ) then
*        ---------------------------------------------------------------
*        Repartition (B S) to get a better B.
*        ---------------------------------------------------------------
         BT     = .false.
         nBfac  = 1

*        Load the basics into kb, since kb(1:m) isn't defined
*        on the first and second major iteration.

         k      = 0
         do 110 j  = 1, nb
            if (hs(j) .eq. 3) then
               k     = k + 1
               kb(k) = j
            end if
  110    continue

         if (k .eq. m) then
*           ------------------------------------------------------------
*           We have the right number of basics.
*           1. Extract the elements of (B S).
*           2. Factorize (B S)'.
*           3. Apply the resulting row permutation to the cols of (B S).
*           4. If S changed, the Hessian won't be any good, so reset it.
*           ------------------------------------------------------------
            call m2belm( 'BS', ms, m, n, nbelem,
     $                   ne, nka, a, ha, ka, kb,
     $                   z(lua), z(indc), z(indr), z(ip), lena )

            call m2bsol( 8, m, w, y, z, nwcore )

            call m2newB( ms, m, nb, hs, z(ip), kb, y, z(locr), nswap )
            if (nswap .gt. 0) then
               r(1) = zero
            end if
         end if
      end if

*     ------------------------------------------------------------------
*     Normal B = LU.
*     Load the basic variables into kb(1:m), slacks first.
*     Set kobj to tell us where the linear objective is.
*     ------------------------------------------------------------------
  200 ntry   = ntry + 1
      invrq  = 0
      invitn = 0
      invmod = 0
      ierr   = 0
      kobj   = 0
      k      = 0

      do 220 j  = n+1, nb
         if (hs(j) .eq. 3) then
            k     = k + 1
            kb(k) = j
            if (j .eq. jobj) kobj = k
         end if
  220 continue

      nslack = k
      nonlin = 0
      do 240 j = 1, n
         if (hs(j) .eq. 3) then
            k  = k + 1
            if (k  .le. m) then
               kb(k) = j
               if (j .le. nn) nonlin = nonlin + 1
            else
               hs(j) = 0
            end if
         end if
  240 continue

      nbasic = k

      if (nbasic .lt. m) then
*        --------------------------------------------------
*        Not enough basics.
*        Set the remaining kb(k) = 0 for m2belm and m2sing.
*        --------------------------------------------------
         do 250 k = nbasic + 1, m
            kb(k) = 0
  250    continue
*--      call iload ( m-nbasic, 0, kb(nbasic+1), 1 )

      else if (nbasic .gt. m) then
*        --------------------------------------------------
*        Too many basics.
*        This is best treated as a fatal error, since
*        m4getb, etc, aim to keep nbasic .le. m.
*        Something must have gone wrong unintentionally.
*        --------------------------------------------------
         go to 930
      end if

*     -----------------------------------------------------------------
*     Load the basis matrix into the LU arrays.
*     -----------------------------------------------------------------
      minlen = nbelem*5/4
      if (minlen .gt. lena) go to 940

      call m2belm( 'B ', m, m, n, nbelem,
     $             ne, nka, a, ha, ka, kb,
     $             z(lua), z(indc), z(indr), z(ip), lena )

      lin    = max( nbasic - nslack - nonlin, 0 )
      if (ninf .gt. 0) then
         cobj   = 'Sum infeas'
      else
         cobj   = 'Objective '
      end if
      if (prnt) write(iprint, 1010) nonlin, lin, nslack, cobj, obj

*     -----------------------------------------------------------------
*     Now factorize B.
*     modtol says if this is the first iteration after Crash.
*     If so, we use big singularity tols to prevent getting an
*     unnecessarily ill-conditioned starting basis.
*     -----------------------------------------------------------------
      modtol = itn .eq. 0  .and.  lcrash .gt. 0

      if (modtol) then
         utol1     = parmlu(4)
         utol2     = parmlu(5)
         parmlu(4) = max( utol1, eps3 )
         parmlu(5) = max( utol2, eps3 )
      end if

      call m2bsol( 0, m, w, y, z, nwcore )

      if (modtol) then
         parmlu(4) = utol1
         parmlu(5) = utol2
      end if

      nsing  = luparm(11)
      minlen = max( minlen, luparm(13) )
      DUmax  = parmlu(13)
      DUmin  = parmlu(14)
      condU  = DUmax / max( DUmin, 1.0d-20 )

      ! 28 May 2004: Don't let the very first LU be badly conditioned.

      if (nfac .eq. 1  .and.  ns .gt. 0  .and.  .not. BS) then
         if (condU .ge. 1.0d+5) then
            BS = .true.
            go to 100
         end if
      end if

      if (ierr .ge. 7) go to 940
      if (ierr .ge. 3) go to 950
      if (ierr .eq. 2) then
*        --------------------------------------------------------------
*        m2bsol says there was large growth in U.
*        LU Factor tol has been reduced.  Try again.
*        --------------------------------------------------------------
         ierr   = 0
         ntry   = 0
         go to 200
      end if

      ierr   = 0

      if ( BS ) then
*        --------------------------------------------------------------
*        We did a BS factorize this time.  Save the smallest diag of U.
*        --------------------------------------------------------------
         Umin   = dumin

      else if ( BT ) then
*        --------------------------------------------------------------
*        (We come here only once.)
*        See if we should have done a BS factorize after all.
*        In this version we do it if any of the following hold:
*           1. dumin (the smallest diagonal of U) is noticeably smaller
*              than Umin (its value at the last BS factor).
*           2. dumin is pretty small anyway.
*           3. B was singular.
*        nBfac  makes BS increasingly likely the longer we
*        keep doing B and not BS.
*        --------------------------------------------------------------
         BT     = .false.
         Utol   = Umin * 0.1d+0 * nBfac
         BS     = dumin .le. Utol   .or.
     $            dumin .le. eps2   .or.
     $            nsing .gt. 0
         if ( BS ) go to 100
      end if

      if (nsing .gt. 0) then
         if (ntry .gt. maxtry) go to 960
*        --------------------------------------------------------------
*        The basis appears to be singular.
*        Suspect columns are indicated by non-positive components of w.
*        Replace them by the relevant slacks and try again.
*        03 Apr 2003: luparm(2) is lprint for LUSOL.
*                     Next line mustn't change the Print level lprint.
*        --------------------------------------------------------------
!!!      lprint = luparm(2)  ! Bug
         call m2sing( luparm(2), m, n, nb,
     $                w, z(ip), z(iq), bl, bu, hs, kb, xn )

*        See if any superbasics slacks were made basic.

         if (ns .gt. 0) then
            ns0    = ns
            do 410 jq = ns0, 1, -1
               j      = kb(m+jq)
               if (hs(j) .eq. 3) then
                  call m6rdel( m, 1, nr, ns, ms,
     $                         kb, bbl, bbu, x, r, x, x, jq, .false. )
                  ns    = ns - 1
                  ms    = m  + ns
               end if
  410       continue
            if (ns .lt. ns0) r(1) = zero
         end if
         go to 200
      end if

*     ------------------------------------------------------------------
*     Compute the basic variables and check that  A * xn = 0.
*     If gotfac was true, ntry = 0.
*     ------------------------------------------------------------------
  500 gotfac = .false.
      call m5setx( 1, m, n, nb, ms, kb,
     $             ne, nka, a, ha, ka,
     $             bl, bu, x, xn, y, y2, z, nwcore )
      if (ierr .gt. 0  .and.  ntry .eq. 0) go to 200
      if (ierr .gt. 0) go to 980

*     Load basic and superbasic bounds into bbl, bbu.

      do 550 k  = 1, ms
         j      = kb(k)
         bbl(k) = bl(j)
         bbu(k) = bu(j)
  550 continue

*     For Crash option 3, linear LG rows should appear to be free.

      if (lcrash .eq. 3) then
         do 600 k = 1, ms
            j     = kb(k)
            if (j .gt. n) then
               if (bl(j) .lt. bu(j)) then
                  bbl(k) = - plinfy
                  bbu(k) = + plinfy
               end if
            end if
  600    continue
      end if

*     Normal exit.

      if (idebug .eq. 100) then
         if (iprint .gt. 0) write(iprint, 2000) (kb(k), x(k), k = 1, ms)
      end if
      return

*     -------------------------------------------------
*     Error exits.
*     m1page( ) decides whether to write on a new page.
*     m1page(2) also writes '=1' if GAMS.
*     -------------------------------------------------

*     Wrong number of basics.

  930 ierr   = 32
      call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1300) nbasic
      if (isumm  .gt. 0) write(isumm , 1300) nbasic
      return

*     Not enough memory.

  940 ierr   = 20
      more   = maxz + 3*(minlen - lena)
      call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1400) maxz, more
      if (isumm  .gt. 0) write(isumm , 1400) maxz, more
      return

*     Error in the LU package.

  950 ierr   = 21
      call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1500)
      if (isumm  .gt. 0) write(isumm , 1500)
      return

*     The basis is structurally singular even after the third try.
*     Time to give up.

  960 ierr   = 22
      call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1600) ntry
      if (isumm  .gt. 0) write(isumm , 1600) ntry
      return

*     Fatal row error.

  980 ierr   = 10
      return


 1000 format(/ ' Factor', i7, '  Demand', i7, '  Itn', i11,
     &         '  Infeas', i8)
 1010 format(  ' Nonlin', i7, '  Linear', i7, '  Slacks', i8,
     &         '  ', a, 1p, e20.8)
 1300 format(/ ' EXIT -- system error.  Too many basic variables:',
     $         i8)
 1400 format(/ ' EXIT -- not enough storage for the basis factors'
     $      // ' Words available =', i8
     $      // ' Increase Workspace (total) to at least', i8)
 1500 format(/ ' EXIT -- error in basis package')
 1600 format(/ ' EXIT -- the basis is singular after', i4,
     $         '  factorization attempts')
 2000 format(/ ' BS and SB values:' // (5(i7, g17.8)))

*     end of m2bfac
      end
