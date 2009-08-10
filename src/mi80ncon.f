************************************************************************
*
*     File  mi80ncon fortran.
*
*     m8augl   m8chkj   m8chkm   m8prtj   m8sclj
*     m8Cinf   m8Dinf   m8rand   m8rc     m8setj   m8srch
*
*     12 Jul 2000: mi81ncon.f now contains m8aug1.
*     23 Jun 2001: m8sclj split into m8sclj and m7sclg.
*     02 Feb 2004: (At GAMS) ifuser, m1user, ierror added to m8setj.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8augl( mode, m, n, nb, ns, inform,
     $                   ne, nka, a, ha, ka,
     $                   hs, bl, bu, xn, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            ha(ne), hs(nb)
      integer            ka(nka)
      double precision   a(ne), bl(nb), bu(nb), xn(nb), z(nwcore)

*     ------------------------------------------------------------------
*     m8augl  is a front-end for m8aug1.  It is called by  misolv  and
*     m5solv  at various stages of the augmented Lagrangian algorithm.
*     11 Oct 1991: a, ha, ka etc added as parameters.
*     05 Mar 1992: mode 0 now makes nonlinear rows free for first major.
*     11 Sep 1992: hs added as parameter.
*     29 Nov 1996: fcon, fold, xold passed to m8aug1.
*     08 Sep 1997: Pass nncon0, njac0 to m8aug1 rather than
*                       nncon , njac
*                  because m8augl( 2 ) is now called by m5solv
*                  with nncon = 0.
*     ------------------------------------------------------------------

      common    /m3loc / lascal,lbl   ,lbu   ,lbbl  ,lbbu  ,
     $                   lhrtyp,lhs   ,lkb
      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /m8loc / lfcon ,lfcon2,lfdif ,lfdif2,lfold ,
     $                   lblslk,lbuslk,lxlam ,lrhs  ,
     $                   lgcon ,lgcon2,lxdif ,lxold

      ms     = m + ns
      njac0  = max( njac, 1 )
      call m8aug1( mode, ms, nncon0, nnjac, njac0, n, nb, inform,
     $             ne, nka, a, ha, ka,
     $             hs, z(lkb), bl, bu, z(lbbl), z(lbbu),
     $             z(lblslk), z(lbuslk),
     $             z(lgcon ), z(lgcon2), xn )

      end ! subroutine m8augl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8chkj( m, n, njac, nx,
     $                   ne, nka, ha, ka,
     $                   bl, bu, f, f2, g, g2,
     $                   x, y, y2, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka)
      double precision   bl(n), bu(n), f(m), f2(m), g(njac), g2(njac),
     $                   x(n), y(n), y2(nx), z(nwcore)

*     ------------------------------------------------------------------
*     m8chkj  verifies the Jacobian  gcon  using
*     finite differences on the constraint functions  f.
*     The various verify levels are described in  m7chkg.
*
*     27 Nov 1996: m8chkj is now called from m5solv after the
*                  linear constraints are satisfied for the first time.
*                  Tests added to allow for scaling.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5log1/ idebug,ierr,lprint
      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)
      common    /m8veri/ jverif(4),lverif(2)

      parameter        ( zero = 0.0d+0,  one = 1.0d+0,  ok = 0.1d+0 )

      logical            cheap, first
      character*4        key
      character*4        lbad         , lgood
      data               lbad /'bad?'/, lgood /'ok  '/

      lvl    = lverif(1)
      if (lvl .lt. 0) return
      j1     = max( jverif(3), 1 )
      j2     = min( jverif(4), n )
      cheap  = lvl .le. 1  .or.  j1 .gt. j2

*     Evaluate f and g at the base point  x.
*     We have to do it with scaling disabled.

      lssave = lscale
      lscale = 0
      if (lderiv .le. 1) call m6dmmy( njac, g )
      call m6fcon( 2, m, n, njac, f, g,
     $             ne, nka, ha, ka,
     $             x, z, nwcore )

      if (ierr      .ne. 0) go to 900
      if (knowng(2) .eq. 0) go to 900

      if (iprint .gt. 0) then
         if (cheap) then
            write(iprint, 1100)
         else
            write(iprint, 1000)
         end if
      end if

*     --------------------------
*     Cheap test.
*     --------------------------

*     Generate a direction in which to perturb  x.

      yj    = one/n
      do 10 j  =  1, n
         y(j)  =   yj
         y2(j) =   yj
         yj    = - yj * 0.99999d+0
   10 continue

*     Define a difference interval.
*     If needed, alter y to ensure that it will be a feasible direction.
*     If this gives zero, go back to original y and forget feasibility.

      dx     = difint(1) * (one  +  dnormi( n, x, 1 ))
      call m7chkd( n, bl, bu, x, dx, y, nfeas )
      if (nfeas .eq. 0) call dcopy ( n, y2, 1, y, 1 )

*     ------------------------------------------------------------------
*     We must not perturb x(j) if the j-th column of the Jacobian
*     contains any unknown elements.
*     ------------------------------------------------------------------
      nder   = 0
      l      = 0
      do 30 j = 1, n
         k1   = ka(j)
         k2   = ka(j+1) - 1

         do 20 k = k1, k2
            ir   = ha(k)
            if (ir .gt. m) go to 25
            l    = l + 1
            if (g(l) .eq. gdummy) y(j) = zero
   20    continue

   25    if (y(j) .ne. zero) nder = nder + 1
   30 continue

      if (nder .eq. 0) go to 900

*     Set f2 = constraint values at a short step along  y.

      do 40 j  = 1, n
         y2(j) = x(j) + dx*y(j)
   40 continue

      call m6fcon( 0, m, n, njac, f2, g2,
     $             ne, nka, ha, ka,
     $             y2, z, nwcore )
      if (ierr .ne. 0) go to 900

*     Set   y2  =  (f2 - f)/dx  -  Jacobian*y.  This should be small.
*     At the same time, find the first Jacobian element in column j1.

      do 60 i  = 1, m
         y2(i) = (f2(i) - f(i)) / dx
   60 continue

      l      = 0
      lsave  = 0
      do 100 j = 1, n
         yj    = y(j)
         k1    = ka(j)
         k2    = ka(j+1) - 1
         do 80 k = k1, k2
            ir   = ha(k)
            if (ir .gt. m) go to 100
            l    = l + 1
            if (j .lt. j1) lsave = l
            y2(ir) = y2(ir) - g(l)*yj
   80    continue
  100 continue

      imax   = idamax( m,y2,1 )
      gmax   = (f2(imax) - f(imax)) / dx
      emax   = abs( y2(imax) ) / (one + abs( gmax ))
      if (emax .le. ok) then
         if (iprint .gt. 0) write(iprint, 1400)
      else
         if (iprint .gt. 0) write(iprint, 1500)
         if (isumm  .gt. 0) write(isumm , 1500)
      end if
      if (iprint .gt. 0) write(iprint, 1600) emax, imax
      if (cheap) go to 900

*     ----------------------------------------------------
*     Proceed with the verification of columns j1 thru j2.
*     ----------------------------------------------------
      if (iprint .gt. 0) write(iprint, 2000)
      l      = lsave
      nwrong = 0
      ngood  = 0
      emax   = - one
      jmax   = 0

      do 200 j = j1, j2

*        See if there are any known gradients in this column.

         k1    = ka(j)
         k2    = ka(j+1) - 1
         lsave = l

         do 120 k = k1, k2
            ir    = ha(k)
            if (ir .gt. m) go to 200
            l     = l + 1
            if (g(l) .ne. gdummy) go to 140
  120    continue
         go to 200

*        Found one.

  140    xj     = x(j)
         dx     = difint(1) * (one + abs( xj ))
         if (bl(j) .lt. bu(j)  .and.  xj .ge. bu(j)) dx = -dx
         x(j)   = xj + dx
         call m6fcon( 0, m, n, njac, f2, g2,
     $                ne, nka, ha, ka,
     $                x, z, nwcore )
         if (ierr .ne. 0) go to 900

*        Check nonzeros in the jth column of the Jacobian.
*        Don't bother printing a line if it looks like an exact zero.

         l     = lsave
         first = .true.

         do 160 k = k1, k2
            ir    = ha(k)
            if (ir .gt. m) go to 180
            l     = l + 1
            gi    = g(l)
            if (gi .eq. gdummy) go to 160
            gdiff = (f2(ir) - f(ir))/dx
            err   = abs( gdiff - gi ) / (one + abs( gi ))

            if (emax .lt. err) then
               emax  = err
               imax  = ir
               jmax  = j
            end if

            key   = lgood
            if (err .gt. eps5) key = lbad
            if (key .eq. lbad) nwrong = nwrong + 1
            if (key .eq.lgood) ngood  = ngood  + 1
            if (abs( gi ) + err  .le.  eps0) go to 160
            if (iprint .gt. 0) then
               if (first) then
                  write(iprint, 2100) j,xj,dx,l,ir,gi,gdiff,key
               else
                  write(iprint, 2200)         l,ir,gi,gdiff,key
               end if
            end if
            first = .false.
  160    continue

  180    x(j)  = xj
  200 continue

      if (iprint .gt. 0) then
         if (nwrong .eq. 0) then
            write(iprint, 2500) ngood ,j1,j2
         else
            write(iprint, 2600) nwrong,j1,j2
         end if
         write(iprint, 2700) emax,imax,jmax
      end if

      if (emax .ge. one) then

*        Bad gradients in  funcon.

         ierr   = 8
         call m1envt( 1 )
         if (iprint .gt. 0) write(iprint, 3800)
         if (isumm  .gt. 0) write(isumm , 3800)
      end if

*     Exit.

  900 lscale = lssave
      return

 1000 format(/ ' Verification of constraint gradients',
     $   ' returned by subroutine funcon.')
 1100 format(/ ' Cheap test on funcon...')
 1400 format(  ' The Jacobian seems to be OK.')
 1500 format(  ' XXX  The Jacobian seems to be incorrect.')
 1600 format(  ' The largest discrepancy was', 1p, e12.2,
     $   '  in constraint', i6)
 2000 format(// ' Column       x(j)        dx(j)', 3x,
     $   ' Element no.    Row    Jacobian value    Difference approxn')
 2100 format(/ i7, 1p, e16.8, e10.2, 2i10, 2e18.8, 2x, a4)
 2200 format(           33x, 2i10, 1pe18.8, e18.8, 2x, a4)
 2500 format(/ i7, '  Jacobian elements in cols ', i6, '  thru', i6,
     $         '  seem to be OK.')
 2600 format(/ ' XXX  There seem to be', i6,
     $   '  incorrect Jacobian elements in cols', i6, '  thru', i6)
 2700 format(/ ' XXX  The largest relative error was', 1p, e12.2,
     $   '   in row', i6, ',  column', i6 /)
 3800 format(// ' EXIT -- subroutine funcon appears to be',
     $   ' giving incorrect gradients')

      end ! subroutine m8chkj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8chkm( m, n, nb, dmaxlm, hs, pi )

      implicit           double precision (a-h,o-z)
      integer            hs(nb)
      double precision   pi(m)

*     ------------------------------------------------------------------
*     m8chkm  checks the multipliers for the next major iteration.
*     pi(*) contains an estimate of lambda(*).
*     Make sure nonbasic elements of pi(*) have the right sign.
*
*     22 May 1992: If the subproblem is optimal and a slack is
*                  superbasic, the pi(i) should again be zero.
*     02 Mar 1994: That seems to be screwing things up!
*     23 Dec 1996: This is now a subroutine so m8setj can call it twice.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      do 100 i = 1, m
         j     = n + i
         js    = hs(j)
         py    = pi(i)
         if (js .eq. 0  .and.  py .gt. zero) py = zero
         if (js .eq. 1  .and.  py .lt. zero) py = zero
***---   if (js .eq. 2  .and.     optsub   ) py = zero
         py    = min( py,    dmaxlm  )
         py    = max( py, (- dmaxlm) )
         pi(i) = py
  100 continue

      end ! subroutine m8chkm

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8prtj( n, nb, nncon, nnjac, lprint, majits, nscl,
     $                   ne, nka, a, ha, ka, hs,
     $                   ascale, bl, bu, fcon, xlam, xn )

      implicit           double precision (a-h,o-z)
      integer            ha(ne), hs(nnjac)
      integer            ka(nka)
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb),
     $                   fcon(nncon), xlam(nncon), xn(nb)

*     ------------------------------------------------------------------
*     m8prtj  prints x, lambda, f(x) and/or the Jacobian
*     at the start of each major iteration.
*
*     08 Apr 1992: Internal values of hs(*) allow more values for key.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m3scal/ sclobj,scltol,lscale

      logical            prtx, prtl, prtf, prtj, scaled
      character*2        key(-1:4)
      data               key/'FR', 'LO', 'UP', 'SB', 'BS', 'FX'/

      if (iprint .le. 0) return

*     Unscale everything if necessary.

      scaled = lscale .ge. 2

      if (scaled) then
         call m2scla( 2, nncon, n, nb, ne, nka,
     $                ha, ka, a, ascale, bl, bu, xlam, xn )
         call ddscl ( nncon, ascale(n+1), 1, fcon, 1 )
      end if

      l      = lprint/10
      prtx   = mod( l,10 ) .gt. 0
      l      = l/10
      prtl   = mod( l,10 ) .gt. 0  .and.  majits .gt. 1
      l      = l/10
      prtf   = mod( l,10 ) .gt. 0
      l      = l/10
      prtj   = mod( l,10 ) .gt. 0

      if ( prtx ) write(iprint, 1100) (xn(j), j=1,nnjac)
      if ( prtl ) write(iprint, 1200) xlam
      if ( prtf ) write(iprint, 1300) fcon
      if ( prtj ) then
         write(iprint, 1400)
         do 100 j = 1, nnjac
            k1    = ka(j)
            k2    = ka(j+1) - 1
            do 50 k = k1, k2
               ir   = ha(k)
               if (ir .gt. nncon) go to 60
   50       continue
            k    = k2 + 1
   60       k2   = k  - 1
            l    = hs(j)
            write(iprint, 1410) j,xn(j),key(l), (ha(k),a(k), k=k1,k2)
  100    continue
      end if
   
*     Rescale if necessary.
         
      if (scaled) then
         call m2scla( 1, nncon, n, nb, ne, nka,
     $                ha, ka, a, ascale, bl, bu, xlam, xn )
         call dddiv ( nncon, ascale(n+1), 1, fcon, 1 )
      end if
      return

 1100 format(/ ' Jacobian variables'
     $       / ' ------------------'   / 1p, (5e16.7))
 1200 format(/ ' Multiplier estimates'
     $       / ' --------------------' / 1p, (5e16.7))
 1300 format(/ ' Constraint functions'
     $       / ' --------------------' / 1p, (5e16.7))
 1400 format(/ ' x  and  Jacobian' / ' ----------------')
 1410 format(i6, 1p, e13.5, 1x, a2, 4(i9, e13.5)
     $       / (22x, 4(i9, e13.5)))

      end ! subroutine m8prtj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8sclj( nncon, nn, ng, n, nb, ne, nka,
     $                   ascale, ha, ka, g )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka)
      double precision   ascale(nb), g(ng)

*     ------------------------------------------------------------------
*     m8sclj  scales the Jacobian.
*     nn, ng, g  are  nnjac, njac, gcon  respectively.
*
*     m8sclj is called by m6fcon only if modefg = 2.
*     Hence, it is used to scale known gradient elements (if any),
*     but is not called when missing gradients are being estimated
*     by m6dcon.
*
*     23 Jun 2001: m7sclg replaces m8sclj( 1, ...).
*     ------------------------------------------------------------------

      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)

      if (knowng(2) .eq. 0) return

*     ------------------------------------------------------------------
*     Scale known Jacobian elements.
*     ------------------------------------------------------------------
      l     = 0

      do 300 j  = 1, nn
         cscale = ascale(j)

         do 250 k = ka(j), ka(j+1) - 1
            ir    = ha(k)
            if (ir .gt. nncon) go to 300
            l     = l + 1
            grad  = g(l)
            if (grad .ne. gdummy) then
               g(l)  = grad * cscale / ascale(n + ir)
            end if
  250    continue
  300 continue

      end ! subroutine m8sclj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8Cinf( n, nb, nncon, Cinf, iCinf,
     $                   ne, nka, a, ha, ka,
     $                   bl, bu, fcon, xn, y, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka)
      double precision   a(ne), bl(nb), bu(nb), fcon(nncon),
     $                   xn(nb), y(nncon), z(nwcore)

*     ------------------------------------------------------------------
*     m8Cinf  finds the maximum nonlinear constraint infeasibility.
*     y(*)    is output as the vector of violations.
*     Cinf    is the biggest infeasibility.
*     iCinf   is the corresponding row.
*
*     01 Apr 1994: First version, replaces m8viol.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

*     Set  y  =  - fcon - (linear A)*xn,   excluding slacks.

      do 520 i = 1, nncon
         y(i)  = - fcon(i)
  520 continue

      call m2aprd( 8, xn, n, y, nncon,
     $             ne, nka, a, ha, ka,
     $             z, nwcore )

*     See how much  y  violates the bounds on the nonlinear slacks.

      Cinf   = zero
      iCinf  = 1

      do 550 i = 1, nncon
         j     = n + i
         slack = y(i)
         viol  = max( zero, bl(j) - slack, slack - bu(j) )
         if (Cinf  .lt. viol) then
             Cinf    =  viol
             iCinf   =  i
         end if
  550 continue

      end ! subroutine m8Cinf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8Dinf( nb, bl, bu, rc, xn, Dinf, jDinf )

      implicit           double precision (a-h,o-z)
      double precision   bl(nb), bu(nb), rc(nb), xn(nb)

*     ------------------------------------------------------------------
*     m8Dinf  computes the maximum dual infeasibility, measured
*     as the maximum departure from complementarity.
*     For vanilla variables xj with reduced cost zj,
*     the dual infeasibility is   xj * zj.
*     More generally, if xj has a lower bound lj, the dual infeas
*     would be   (xj - lj) * zj,   but to allow for infinite bounds,
*     we use min( xj - lj, one ) * zj.
*
*     Note that (xj - lj) * zj varies smoothly as xj varies near lj,
*     whereas the usual "largest zj of appropriate sign" jumps from
*     zero to |zj| as xj increases from lj.
*
*     m8Dinf  is called by m8setj and m8srch.
*     
*     On exit,
*      Dinf  is the maximum dual infeasibility (departure from
*            complementarity).
*     jDinf  is the corresponding variable.
*
*     03 Apr 1994: First version, derived from m2Dinf and
*                  snopt routine s8infs.
*     17 Nov 1995: m8rc fixed to treat linear obj correctly.
*                  No longer need to treat jobj specially here.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

*      if (jobj .gt. 0) then
*         blobj    = bl(jobj) 
*         bl(jobj) = bu(jobj)
*      end if

      jDinf = 0
      Dinf  = zero
      
      do 500 j = 1, nb
         if (bl(j) .lt. bu(j)) then
            dj     = rc(j)
            if      (dj .gt. zero) then
               dj  =   dj * min( xn(j) - bl(j), one )
            else if (dj .lt. zero) then
               dj  = - dj * min( bu(j) - xn(j), one )
            end if
            
            if (Dinf .lt. dj) then
                Dinf   =  dj
                jDinf  =  j
            end if
         end if
  500 continue

*     if (jobj .gt. 0) then
*        bl(jobj) = blobj
*     end if

      end ! subroutine m8Dinf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8rand( ia, ib, ic, n, x, incx )

      implicit           double precision (a-h,o-z)
      double precision   x(*)

*     ------------------------------------------------------------------
*     m8rand fills a vector x with uniformly distributed random numbers
*     in the interval (0, 1) using a method due to Wichman and Hill.
*
*     ia, ib and ic should be set to integer values between 1 and 30000
*     before the first entry.
*
*     Integer arithmetic up to 30323 is required.
*
*
*     28 Nov 1996: Originally copied from Wichman and Hill 19-Jan-1987.
*                  This version derived from SNOPT routine ddrand.
*                  ia, ib, ic now parameters to avoid new common block.
*     ------------------------------------------------------------------

      if (n .lt. 1) return

      do 100, ix = 1, 1+(n-1)*incx, incx
         ia     = 171*mod(ia, 177) -  2*(ia/177)
         ib     = 172*mod(ib, 176) - 35*(ib/176)
         ic     = 170*mod(ic, 178) - 63*(ic/178)
         
         if (ia .lt. 0) ia = ia + 30269
         if (ib .lt. 0) ib = ib + 30307
         if (ic .lt. 0) ic = ic + 30323
         
         x(ix)  = mod( real(ia)/30269.0 +
     $                 real(ib)/30307.0 + real(ic)/30323.0, 1.0)
  100 continue

      end ! subroutine m8rand

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8rc  ( jobj, minimz, sclobj,
     $                   m, n, nb, nnobj, nnobj0, nncon, nnjac, njac,
     $                   ne, nka, a, ha, ka,
     $                   gobj, gcon, pi, rc )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka)
      double precision   a(ne)
      double precision   gobj(nnobj0), gcon(njac),
     $                   pi(m), rc(nb)

*     ------------------------------------------------------------------
*     m8rc   computes reduced costs rc = gobj - (A I)'pi
*     for all columns of (A I),
*     using  gcon  as the top left-hand corner of A.
*     gcon, gobj and pi are assumed to exist.
*
*     m8rc   is called by m8setj and m8srch.
*
*     28 Sep 1993: First version, derived from m4rc (now m2rcA).
*     17 Nov 1995: jobj needed as extra parameter for problems
*                  with a linear objective.  (Found by PEG in SNOPT.)
*     13 Mar 1996: sclobj also needed for scaled problems.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      l     = 0

      do 200 j = 1, nnjac
         dj    = zero
         do 150 k = ka(j), ka(j+1) - 1
            i     = ha(k)
            if (i .le. nncon) then
               l  = l + 1
               dj = dj  +  pi(i) * gcon(l)
            else
               dj = dj  +  pi(i) * a(k)
            end if
  150    continue
         rc(j) = - dj
  200 continue

      do 300 j = nnjac+1, n
         dj    = zero
         do 250 k = ka(j), ka(j+1) - 1
            i     = ha(k)
            dj    = dj  +  pi(i) * a(k)
  250    continue
         rc(j) = - dj
  300 continue

      do 320 i = 1, m
         rc(n+i) = - pi(i)
  320 continue

*     Include the nonlinear objective gradient.

      sgnobj = minimz
      if (nnobj .gt. 0) then
         call daxpy ( nnobj, sgnobj, gobj, 1, rc, 1 )
      end if

*     Include the linear objective.

      if (jobj .gt. 0) rc(jobj) = rc(jobj)  -  sgnobj * sclobj

      end ! subroutine m8rc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8setj( lcstat, m, n, nb, nscl,
     $                   nn, nncon, nnjac, njac, nnobj, nnobj0,
     $                   lcrash, ns, nswap0, objadd, sameSB,
     $                   nconvg, nomove, nreduc, optsub, convgd,
     $                   ne, nka, a, ha, ka,
     $                   hrtype, hs, kb,
     $                   ascale, bl, bu,
     $                   fcon, fcon2, fold,
     $                   gcon, gcon2, gobj, gobj2,
     $                   xlam, rhs  , xdif, xold,
     $                   pi, rc, xn, y, y2, z, nwcore,
     $                   ifuser, m1user, ierror )

      implicit           double precision (a-h,o-z)
      logical            convgd, optsub, sameSB
      integer            ha(ne), hrtype(m), hs(nb)
      integer            ka(nka), kb(m)
      double precision   a(ne), ascale(nscl), bl(nb), bu(nb)
      double precision   fcon (nncon), fcon2(nncon), fold(nncon),
     $                   gcon (njac ), gcon2(njac ),
     $                   gobj(nnobj0), gobj2(nnobj0),
     $                   xlam (m)    , rhs  (nncon),
     $                   xdif (nb)   , xold (nb)   ,
     $                   pi(m), rc(nb), xn(nb), y(m), y2(m), z(nwcore)
      integer            ifuser, ierror
      external           m1user

*     ------------------------------------------------------------------
*     m8setj  sets up the linearly constrained subproblem for the
*     next major iteration of the augmented Lagrangian algorithm.
*
*     Input parameters:
*
*     lcstat  gives the state of the linearly constrained subproblem
*             just solved.
*             lcstat = 0  It terminated normally.
*                    = 1  It was infeasible.
*                    = 2  It was unbounded.
*                    = 3  It terminated early (too many minor itns).

*     lcrash  If lcrash = 5, the nonlinear constraints have been
*             relaxed so far.  On major itn 2, Crash is called
*             to find a basis including them.
*
*     sameSB  is true if ns remained the same (typically 0) throughout
*             the major iteration.  This is intended to detect square
*             systems of nonlinear equations.  m8srch doesn't need to
*             worry about the dual infeasibility increasing.
*             We allow ns to be any fixed value, since sometimes there
*             are extraneous guys sitting there doing nothing.
*
*     nconvg  counts no. of times optimality test has been satisfied.
*     nreduc  counts no. of times  penpar  has been reduced.
*     optsub  says if the subproblem just solved was optimal.
*
*     Output parameters:
*
*     convgd  says if it is time to stop.
*
*     Local variables:
*
*     Prinf   is the Primal infeasibility (nonlinear constraint violn).
*     Duinf   is the Dual   infeasibility.
*     Prbest  is the smallest Primal infeasibility so far.
*     Dubest  is the smallest Dual   infeasibility so far.
*     PDbest  is the smallest sum of Prinf + Duinf.
*     best    is true if the new point is the best so far.
*     better  is true (from m8srch) if a step of 1.0 was taken and
*             the new point seems better than the previous one.
*
*     06 Mar 1992: The first major iteration now relaxes nonlinear rows
*                  and terminates as soon as the linear constraints are
*                  satisfied.  xold is not defined until Major 2.
*     08 Apr 1992: Internal values of hs(*) simplify setting py2.
*     22 May 1992: penpar reduced for more majors before setting to zero.
*                  Stops Steinke2 (and maybe OTIS?) from suddenly
*                  getting a huge x and lambda.
*     04 Jun 1992: lcrash added.
*     03 Sep 1992: Penalty parameter is now a multiple of the
*                  default value 100/nncon.
*     13 Aug 1993: Implemented key (the "Print level 0 key").
*                  If Print level = 0, there is only one line of output.
*                  If Print level > 0, the same line is included with
*                  key = '>' at the end of it.  A "grep" or "search"
*                  can extract such lines from the full listing.
*     03 Oct 1993: Implemented m8srch to backtrack 0 or 1 times.
*                  Changed dmaxlm from e+8 to e+12 to help Square Root
*                  problems in ProbsB.
*                  (05 Nov 1993: This gave trouble with HS 109.
*                   Settle for e+9.)
*     10 Oct 1993: xlam now has length m, to allow m8srch to take a
*                  short step in all of pi.  (This affects the reduced
*                  costs found by m8rc and m8Dinf, and hence alters
*                  Duinf0 at the next major iteration.)  xlam still has
*                  length nncon elsewhere, e.g. in the augmented
*                  Lagrangian objective for the subproblems.
*     02 Mar 1994: Reduce Penalty parameter if m8srch gives better=true.
*     04 Mar 1994: But not if PDbest says we're not at the best point
*                  so far.
*     15 Mar 1994: Define penmax and penmin from initial Penalty param.
*                  No longer reduce penpar to zero. (User may still set
*                  it to zero initially, but it might be raised later.)
*     25 Mar 1994: xold and xdif now have length n, so we can take a
*                  step in all of x.  We can then decrease the
*                  constraint violation with a sufficiently short step.
*                  m8srch modified to backtrack more thoroughly.
*     11 Sep 1994: Don't reduce penalty if step < 1.0.
*     11 Nov 1994: Bug -- must save local variables Prinf and Duinf.
*     24 Jan 1995: xold, xdif now have length nb.
*                  (m8Dinf uses all of xn, so n was not enough.)
*                  Even if step < 1.0, m8srch resets xn to the full step
*                  to initiate the next major.
*     13 Mar 1996: m8rc needs sclobj.
*     29 Nov 1996: The Jacobian is now random at the start of major1.
*                  We no longer print anything involving the functions.
*     02 Feb 2004: (At GAMS) ifuser, m1user, ierror added permanently.
*     15 Jul 2004: Print relative infeasibilities Prrel, Durel
*                  rather than absolute ones (so we can see how close
*                  we are to terminating).
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2lu3 / lenl,lenu,ncp,lrow,lcol
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      logical            prnt0 ,prnt1 ,summ0 ,summ1 ,newhed
      common    /m5log4/ prnt0 ,prnt1 ,summ0 ,summ1 ,newhed
      common    /m5lp1 / itn,itnlim,nphs,kmodlu,kmodpi
      common    /m5step/ featol, tolx0,tolinc,kdegen,ndegen,
     $                   itnfix, nfix(2)
      common    /m5tols/ toldj(3),tolx,tolpiv,tolrow,rowerr,xnorm
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /m8al2 / radius,rhsmod,modpen,modrhs
      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)
      common    /m8func/ nfcon(4),nfobj(4),nprob,nstat1,nstat2
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave
      common    /cycle2/ objtru,suminf,numinf

      parameter        ( zero   = 0.0d+0,  one    =  1.0d+0,
     $                   dmaxlm = 1.0d+12 )

      logical            better, best  , conv1 , conv2 , convt , convu,
     $                   grdcon, grdobj, major1, major2, Phead , Shead

      double precision   Prbest, Dubest, PDbest, penmax, penmin
      save               Prbest, Dubest, PDbest, penmax, penmin

      double precision   Prinf , Duinf
      save               Prinf , Duinf

      integer            ntight
      save               ntight

* EK oct 2002:
      double precision   dumr

      character*10       line
      character*1        lstate(4), label
      character*2        Pkey, Skey

      data               line   /'----------'/
      data               lstate /' ','I','U','T'/

      rhsmod = zero
      modrhs = 0
      dlmax  = zero
      dlrel  = zero
      dxmax  = zero
      dxrel  = zero
      step   = one
      grdcon = lderiv .ge. 2
      grdobj = nnobj  .eq. 0  .or.  lderiv .eq. 1  .or.  lderiv .eq. 3
      major1 = majits .eq. 1
      major2 = majits .eq. 2
      ierror = 0

      if (major1) then
*        ---------------------------------------------------------------
*        Beginning of first major iteration.
*        misolv has already called m8augl( 2,... ) to make sure that
*        xn is inside its bounds.
*        ---------------------------------------------------------------

*        14 Sep 2000: Initialize gobj in case optimality conditions
*                     are checked before the functions are evaluated.

         if (nnobj .gt. 0) then
            call dload ( nnobj, zero, gobj, 1 )
         end if

*        Set limits for penalty parameter, based on initial value.

         pen0   = penpar
         if (penpar .le. zero) pen0 = one
         penmax = pen0 * 1.0d+1
         penmin = pen0 * 1.0d-4

*        Initialize fold, xold and a few other other things.

         call dload ( nncon, zero,    fold, 1 )
         call dcopy ( nb   , xn  , 1, xold, 1 )
         Prbest = plinfy
         Dubest = plinfy
         PDbest = plinfy

*        Save the nonlinear constraint bounds (unscaled),
*        then relax them.
*        (The first major will just satisfy the linear constraints.)
*        Then set internal values for hs.

         call m8augl(-1, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )
         call m8augl( 0, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )
         call m5hs  ( 'Internal', nb, bl, bu, hs, xn )

*        Set rhs = nonlinear rows of (A I)*xn.  (m2aprd gives y = -A*x.)
*        This stops nonlinear rows interfering on a Warm or Hot start.

         call dload ( m, zero, y, 1 )
         call m2aprd( 5, xn, n, y, m,
     $                ne, nka, a, ha, ka,
     $                z, nwcore )
         call dscal ( nncon, (- one), y, 1 )
         call daxpy ( nncon, (+ one), xn(n+1), 1, y, 1 )
         call dcopy ( nncon, y, 1, rhs, 1 )
         go to 500
      end if

*     ------------------------------------------------------------------
*     Output heading for detailed log.
*     ------------------------------------------------------------------
      if (prnt1) write(iprint, 1010) (line, j=1,12)
      if (summ1) write(isumm , 1010) (line, j=1,7)

      if (major2) then
*        ---------------------------------------------------------------
*        Beginning of Major 2.
*        ---------------------------------------------------------------

*        Check xlam as initialized in m4getb (unscaled).

         call dcopy ( m, xlam, 1, pi, 1 )
         call m8chkm( m, n, nb, dmaxlm, hs, pi )

*        ---------------------------------------------------------------
*        Copy the Jacobian into A (unscaled).
*        ---------------------------------------------------------------
         if (.not. grdcon) then
            lssave = lscale
            lscale = 0
            call m6dcon( nncon, nnjac, njac,
     $                   ne, nka, ha, ka,
     $                   fcon, fcon2, gcon, gcon2,
     $                   xn, y, z, nwcore )
            lscale = lssave
            if (ierr .ne. 0) go to 900
         end if

         call m8augl( 1, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )

*        ---------------------------------------------------------------
*        Scale Jacobian if requested.
*        ---------------------------------------------------------------
         if (lscale .eq. 2) then
            call m2amat( 2, m, n, nb,
     $                   ne, nka, a, ha, ka,
     $                   bl, bu, hrtype )
            call m2scal( m, n, nb, ne, nka, nn, nncon, nnjac,
     $                   hrtype, ha, ka, a, ascale, bl, bu, y, y2 )
            call m2scla( 1, m, n, nb, ne, nka,
     $                   ha, ka, a, ascale, bl, bu, pi, xn )

*           Evaluate functions again so they will be scaled.

            if (.not. grdcon) then
               call m6dmmy( njac, gcon )
            end if
            call m6fcon( 2, nncon, nnjac, njac, fcon, gcon,
     $                   ne, nka, ha, ka,
     $                   xn, z, nwcore )
            if (ierr .ne. 0) go to 900
            if (.not. grdcon) then
               call m6dcon( nncon, nnjac, njac,
     $                      ne, nka, ha, ka,
     $                      fcon, fcon2, gcon, gcon2,
     $                      xn, y, z, nwcore )
            end if
            if (ierr .ne. 0) go to 900

            if (nnobj .gt. 0) then
               if (.not. grdobj) then
                   call m6dmmy( nnobj, gobj )
               end if
               call m6fobj( 2, nnobj, fobj, gobj, xn, z, nwcore )
               if (ierr .ne. 0) go to 900
               if (.not. grdobj) then
                  call m6dobj( nnobj, fobj, gobj, gobj2, xn, z, nwcore )
               end if
               if (ierr .ne. 0) go to 900
            end if
         end if

*        ---------------------------------------------------------------
*        Initialize xlam and xold (possibly scaled).  fold is set below.
*        ---------------------------------------------------------------
         call dcopy ( m , pi, 1, xlam, 1 )
         call dcopy ( nb, xn, 1, xold, 1 )

*        ---------------------------------------------------------------
*        Save bounds on nonlinear rows (now possibly scaled).
*        ---------------------------------------------------------------
         call m8augl(-1, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )

*        ---------------------------------------------------------------
*        Compute reduced costs rc(*) for all columns and rows.
*        Find the maximum primal and dual infeasibilities.
*        ---------------------------------------------------------------
         call m8rc  ( jobj, minimz, sclobj,
     $                m, n, nb, nnobj, nnobj0, nncon, nnjac, njac,
     $                ne, nka, a, ha, ka,
     $                gobj, gcon, pi, rc )

         call m8Cinf( n, nb, nncon, Cinf, iCinf,
     $                ne, nka, a, ha, ka,
     $                bl, bu, fcon, xn, y, z, nwcore )

         call m8Dinf( nb, bl, bu, rc, xn, Dinf, jDinf )

         step   = zero
         Prinf  = Cinf
         Duinf  = Dinf
         Prbest = min( Prbest, Prinf )
         Dubest = min( Dubest, Duinf )
         PDbest = min( PDbest, Prinf + Duinf )
         ntight = 1
*        ---------------------------------------------------------------
*        End of Major 2.
*        ---------------------------------------------------------------
      else
*        ---------------------------------------------------------------
*        Major 3 and later.
*        ---------------------------------------------------------------

*        Check sign of multipliers.

         call m8chkm( m, n, nb, dmaxlm, hs, pi )

*        Restore the bounds on the nonlinear constraints.

         call m8augl( 6, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )

*        Make sure all variables are inside their bounds.
*        Then set internal values for hs.

         call m8augl( 2, m, n, nb, ns, inform,
     $                ne, nka, a, ha, ka,
     $                hs, bl, bu, xn, z, nwcore )
         call m5hs  ( 'Internal', nb, bl, bu, hs, xn )

*        Set xdif = xn - xold.
*        Note that xdif and xold allow for ALL variables,
*        in case "step" turns out not to be one.

         call dcopy ( nb,         xn  , 1, xdif, 1 )
         call daxpy ( nb, (-one), xold, 1, xdif, 1 )

*        Find the relative change in the Jacobian variables.

         imax   = idamax( nnjac, xdif, 1 )
         dxmax  = abs   ( xdif(imax) )
         dxnorm = dnormi( nnjac, xdif, 1 )
         xonorm = dnormi( nnjac, xold, 1 )
         dxrel  = dxnorm / (one + xonorm)

*        Set y2 = pi - xlam for ALL rows.

         call dcopy ( m,         pi  , 1, y2, 1 )
         call daxpy ( m, (-one), xlam, 1, y2, 1 )

*        Find the relative change in  xlam  for just NONLINEAR rows.
*        If xlam is very small, it is probably still zero (never set).

         imax   = idamax( nncon, y2  , 1 )
         dlmax  = abs   ( y2(imax) )
         dlnorm = dnormi( nncon, y2  , 1 )
         xlnorm = dnormi( nncon, xlam, 1 )
         dlrel  = dlnorm / (one + xlnorm)
         if (xlnorm .le. eps2) dlrel = one

*        ---------------------------------------------------------------
*        Determine a step to be used as follows:
*           xn  =  xold + step*xdif,    xlam  =  xlam + step*y2.
*        m8srch evaluates the functions at the chosen steplength.
*        y(1:nncon) returns the nonlinear constraint violations.
*        ---------------------------------------------------------------
         call m8srch( nlag  , ns    , ntight,
     $                objadd, sameSB,
     $                dxrel , dlrel,
     $                Prinf , Duinf ,
     $                step  , better,
     $                m, n, nb, nnobj, nnobj0,
     $                nncon, nnjac, njac,
     $                ne, nka, a, ha, ka,
     $                bl, bu,
     $                fcon, fcon2, gcon, gcon2, gobj,
     $                xlam, pi   , rc  , 
     $                xdif, xn   , xold,
     $                y   , y2   , z   , nwcore )

         call m5hs  ( 'Internal', nb, bl, bu, hs, xn )
         if (ierr .ne. 0) go to 900

         Prbest = min( Prbest, Prinf )
         Dubest = min( Dubest, Duinf )
         PD     = Prinf + Duinf
         PDbest = min( PDbest, PD )
         best   = PDbest .eq. PD

         if (dxrel .le. eps0) then
            if ( prnt1 ) write(iprint, 2100)
            if ( summ1 ) write(isumm , 2100)
            nomove = nomove + 1
         else
            nomove = 0
         end if
      end if

*     ------------------------------------------------------------------
*     All major iterations (except first).
*     ------------------------------------------------------------------

*     Save fcon as fold.  xold and xlam are already ok.
*     Copy the Jacobian into  A  to form the next linearization.

      call dcopy ( nncon, fcon, 1, fold, 1 )
      call m8augl( 1, m, n, nb, ns, inform,
     $             ne, nka, a, ha, ka,
     $             hs, bl, bu, xn, z, nwcore )

*     ------------------------------------------------------------------
*     If Print level = 0, print one line.
*     If Print level > 0, print same line flagged with ">".
*     ------------------------------------------------------------------
      label  = lstate(lcstat + 1)
      lu     = lenl + lenu

      xnorm  = dnormi( nb, xold, 1 )
      pinorm = dnormi( m , xlam, 1 )
      Prrel  = Prinf / (one + xnorm )
      Durel  = Duinf / (one + pinorm)
      vimax  = Prinf
      virel  = Prrel

*     Output heading for terse log.

      if ( prnt0 ) then
         Phead  = mod( majits, 20 ) .eq. 2
     $            .or.  (newhed  .and.  majits .ge. 2)
         if ( Phead ) then
            write(iprint, 1100)
            newhed = .false.
         end if
      end if

      if ( summ0 ) then
         Shead  = mod( majits, 10 ) .eq. 2
         if ( Shead ) write(isumm , 1110)
      end if

*     Output heading for detailed log.

      Pkey   = ' '
      Skey   = ' '
      if ( prnt1 ) then
         Pkey = '>'
         if ( major2 ) then
            write(iprint, 1100) Pkey
         else
            write(iprint, 1100)
         end if
      end if
      if ( summ1 ) then
         Skey = '>'
         if ( major2 ) then
            write(isumm , 1110) Skey
         else
            write(isumm , 1110)
         end if
      end if

      if (iprint .gt. 0)
     $write(iprint, 1200) majits - 1, minits, label, itn, ninf,
     $                    step, objtru, Prrel, Durel, nssave,
     $                    nfcon(1), lu, penpar, nswap0, Pkey

      if (isumm  .gt. 0)
     $write(isumm , 1210) majits - 1, minits, label,
     $                    step, objtru, Prrel, Durel, nssave,
     $                    nfcon(1),     penpar, nswap0, Skey

      if ( prnt1 ) then
         write(iprint, 1000)
         if (.not. major1)
     $   write(iprint, 2400) dxmax, dxrel, dlmax, dlrel
         write(iprint, 2600) vimax, virel
      end if

      if ( summ1 ) then
         write(isumm , 1000)
         if (.not. major1)
     $   write(isumm , 2410) dxmax, dlmax
      end if


*     ------------------------------------------------------------------
*     User formatted output
*     EK oct 2002
*     02 Feb 2004: Added permanently to m8setj.
*                  m1user may return
*                     ierror = 0  All is well
*                            = 1  Resource limit: CPU time exceeded
*                            = 2  User interrupt: control-c.
*                  These are handled by m5solv.
*     ------------------------------------------------------------------
      if (ifuser .gt. 0) then
         call m1user( 2, ierror, majits-1, minits, dumr, dumr, label,
     $                step, objtru, Prinf, Duinf, nssave,
     $                nfcon(1), penpar, nswap0 )
      end if


*     ------------------------------------------------------------------
*     Compute  rhs  =  fcon  -  Jacobian*xn.
*     This is minus the required  rhs  for the new subproblem.
*     ------------------------------------------------------------------
      if (lprint .gt. 1) then
         call m8prtj( n, nb, nncon, nnjac, lprint, majits, nscl,
     $             ne, nka, a, ha, ka, hs,
     $             ascale, bl, bu, fcon, xlam, xn )
      end if

      call dcopy ( nncon, fcon, 1, rhs, 1 )
      call m2aprd( 7, xn, nnjac, rhs, nncon,
     $             ne, nka, a, ha, ka,
     $             z, nwcore )
      call dscal ( nncon, (- one), rhs, 1 )

*     ------------------------------------------------------------------
*     Do Crash on nonlinear rows at start of Major 2
*     unless a basis is known already.
*     ------------------------------------------------------------------
      if (major2  .and.  lcrash .eq. 5) then

*        m8Cinf has set y(*) to be the exact slacks on nonlinear rows.
*        Copy these into xn(n+i) and make sure they are feasible.
*        m2crsh uses them to decide which slacks to grab for the basis.

         do 450 i = 1, nncon
            j     = n + i
            xn(j) = max(  y(i), bl(j) )
            xn(j) = min( xn(j), bu(j) )
  450    continue

*        Set row types in hrtype.
*        m2crsh uses kb and hrtype as workspace.
*        It may alter xn(n+i) for nonlinear slacks.
*        Set internal values for hs again to match xn.

         call m2amat( 2, m, n, nb,
     $                ne, nka, a, ha, ka,
     $                bl, bu, hrtype )
         call m2crsh( lcrash, m     , n     , nb    , nn,
     $                ne    , nka   , a     , ha    , ka,
     $                kb    , hs    , hrtype, bl , bu,
     $                xn    , z     , nwcore )
         lcrash = 0
         call m5hs  ( 'Internal', nb, bl, bu, hs, xn )
      end if

*     ------------------------------------------------------------------
*     Test for optimality.
*     For safety, we used to require the convergence test to be
*     satisfied for 2 subproblems in a row.
*     If Prinf and Duinf are pretty small, 1 in a row is enough.
*
*        conv1: infeas < 0.1            (set Completion Full)
*        conv2: infeas < Radius of conv (start reducing Penalty param)
*        convt: infeas < tight tols     (1 optsub is enough)
*        convu: infeas < user  tols     (2 optsubs for safety)
*     02 Mar 1994: radius = 0.01 delays reduction of Penpar.
*                  radius = 0.1  set as new default.
*     20 Mar 1994: Define conv2 to match 1982 Minos/Augmented paper.
*                  This is happy with small dlrel, not Durel.
*     20 Mar 1994b:Also wait until    small dxrel.
*     22 Dec 1995: Don't allow convergence if lcstat = 1 or 2
*                  (subproblem was infeasible or unbounded).
*     ------------------------------------------------------------------
  500 if (major1  .or.  major2
     $            .or.  lcstat .eq. 1
     $            .or.  lcstat .eq. 2) then
         convgd = .false.
         nconvg = 0
      else
*-       convgd = dxrel .le. 0.1     .and.  dlrel .le. 0.1  .and.
*-   $            virel .le. rowtol  .and.  optsub          .and.
*-   $            dfrel .le. 0.001

         conv1  = Prrel .le. 0.1     .and.  Durel .le. 0.1
*--      conv2  = Prrel .le. radius  .and.  Durel .le. radius
         conv2  = Prrel .le. radius  .and.  dlrel .le. radius
     $                               .and.  dxrel .le. radius
         convt  = Prrel .le. eps2    .and.  Durel .le. eps2
         convu  = Prrel .le. rowtol  .and.  Durel .le. toldj(3)

         if (convt) then
            mxconv = 1
         else
            mxconv = 2
         end if         

         if (convu) then
            nconvg = nconvg + 1
         else
            nconvg = 0
         end if

         convgd = nconvg .ge. mxconv

*        ========================================================
*        Change from Partial Completion to Full Completion if
*        the iterations seem to be converging.
*        (This may be begging the question!  More theory needed.)
*        ========================================================
         if (ncom .eq. 1) then
*           Relax -- it's already Completion Full.
         else if (optsub .and. conv1) then      
            ncom = 1
            if (iprint .gt. 0) write(iprint, 2800)
            if (isumm  .gt. 0) write(isumm , 2800)
         end if
      
*        ========================================================
*        Reduce the penalty parameter if possible.
*        To help a small "Radius of convergence" delay reduction,
*        conv2 always has to be true.
*
*        15 Mar 1994: Don't reduce penalty all the way to zero.
*                     nreduc is no longer needed.
*        11 Sep 1994: Don't reduce penalty if step < 1.0.
*        ========================================================
         if (optsub  .and.  penpar .gt. penmin) then
            if (conv2  .and.  (better .or. best)
     $                 .and.  (step   .eq. one )) then
               nreduc = nreduc + 1
               fac    = 10.0d+0
               penpar = max( penpar / fac, penmin )

*              The next line caused trouble on Steinke2
*              (suddenly a huge search direction).
*              For safety we should always reduce penpar gradually.
*------>       if (nconvg .ge. 1  .or.  nreduc .ge. 5) penpar = zero

*-             if (penpar .le. 1.0d-10) then
*-                penpar = zero
*-             end if
            end if
         end if      

*        ==================================================
*        Otherwise, perhaps increase the penalty parameter.
*        ratio = (infeas at step 1) / (infeas at step 0)
*        used to be given by m8srch.
*        ==================================================
         if (majits .ge. 5  .and.  dlrel  .ge. 1.0d+5) then
*---        if (ratio  .ge. 1.0d+3) then
            if (penpar .lt. penmax) then
               nreduc = 0
               oldpen = penpar
               if (penpar .le. zero) penpar = one

*              04 Sep 1992: Raising the penalty parameter (below)
*                           seems to cause more trouble than it fixes.
*                           Give it up for now.
*              04 Nov 1993: Try again using "big ratio" instead of
*                           "big dlrel" to get into here.
*              04 Dec 1996: Try bigger "dlrel".

               fac    = 2.0d+0
               penpar = min( penpar * fac, penmax )
               if (penpar .gt. oldpen) then
                  if (iprint .gt. 0) write(iprint, 2900) penpar
                  if (isumm  .gt. 0) write(isumm , 2900) penpar
               end if
            end if
         end if
      end if

*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
  900 return

 1000 format(' ')
 1010 format(/ 1x, 13a10)
 1100 format(/ ' Major minor   total ninf',
     $         ' step     objective     Feasible Optimal  nsb',
     $         '   ncon     LU penalty BSwap', a)
 1110 format(/ ' Major minor ',
     $         ' step     objective  Feasible Optimal  nsb',
     $         '   ncon penalty BSwap', a)
 1200 format(1p, 2i6, a1, i7, i5,
     $       e8.1, e16.8, 2e8.1, i5,
     $       i7, i7, e8.1, i6, a)
 1210 format(1p, 2i6, a1,
     $       e8.1, e13.5, 2e8.1, i5,
     $       i7,     e8.1, i6, a)
 2100 format(' Jacobian variables have not changed')
 2400 format(' Maximum change in Jacobian vars =', 1p, e12.2,
     $   4x,  ' ( =', e11.2, ' normalized)'
     $     / ' Maximum change in multipliers   =',     e12.2,
     $   4x,  ' ( =', e11.2, ' normalized)' )
 2410 format(' Change in Jacobn vars =', 1p, e12.2
     $     / ' Change in multipliers =',     e12.2)
 2600 format(' Maximum constraint violation    =', 1p, e12.2,
     $   4x, ' ( =', e11.2, ' normalized)' )
 2800 format(' Completion Full    now requested')
 2900 format(' Penalty parameter increased to', 1p, e12.2)

      end ! subroutine m8setj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8srch( nlag  , ns    , ntight,
     $                   objadd, sameSB,
     $                   dxrel , dlrel ,
     $                   Prinf , Duinf ,
     $                   step  , better,
     $                   m, n, nb, nnobj, nnobj0,
     $                   nncon, nnjac, njac,
     $                   ne, nka, a, ha, ka,
     $                   bl, bu,
     $                   fcon, fcon2, gcon, gcon2, gobj,
     $                   xlam, pi   , rc  , 
     $                   xdif, xn   , xold,
     $                   y   , y2   , z   , nwcore )

      implicit           double precision (a-h,o-z)
      logical            better, sameSB
      integer            ha(ne)
      integer            ka(nka)
      double precision   a(ne), bl(nb), bu(nb)
      double precision   fcon(nncon), fcon2(nncon),
     $                   gcon(njac ), gcon2(njac ), gobj(nnobj0),
     $                   xlam(m)    , pi(m)       , rc(nb)      , 
     $                   xdif(nb)   , xn(nb)      , xold(nb)    ,
     $                   y(m)       , y2(m)       , z(nwcore)

*     ------------------------------------------------------------------
*     m8srch  determines a steplength between major iterations.
*     It finds a step to be used as follows:
*        xn  =  xold + step*xdif,    xlam  =  xlam + step*y2.
*     The preferred steplength is step = 1.0.
*
*     On entry,
*     Prinf, Duinf, xold, xlam
*                  have values from the previous major.
*     xn, pi       have new values.
*     dxrel, dlrel are the relative change in the Jacobian variables
*                  and lambda (after a unit step).
*     optsub       says if the major iteration terminated optimal.
*     sameSB       says if ns was constant during the major.
*
*     On exit,
*     better = .true.  if a unit step is taken and both primal and
*                  dual infeasibilities are reduced or at least not
*                  significantly bigger than before.
*
*     Note that the current working set is retained to start the
*     next major iteration.  If the final step is less than 1.0,
*     the functions will be evaluated there, but we then have a choice
*     on what to do with xn.  Define
*         xfull  = xold + 1.0 *xdif = current xn,
*         xshort = xold + step*xdif.
*     We used to set xn = xshort (with some components between their
*     bounds).  Now it seems better to set xn = xfull to match the
*     working set.  Reason: If step is very small, the next major
*     will do almost the same as the previous one, terminating close
*     to xfull.  Hence, we might as well start there.
*
*     02 Oct 1993 -- Jan 1995: Tried to choose step by monitoring
*                  Prinf and Duinf, the primal and dual infeasibilities.
*                  These seem jump around too much on some problems.
*
*     24 Jan 1995: xold, xdif now have length nb.
*                  (m8Dinf uses all of xn, so n was not enough.)
*
*     27 Feb 1995: Reverted to limiting the change in
*                  x(1:nnjac) and pi(1:nncon), as in 1982 paper.
*                  Still not sure if xn should be xfull or xshort
*                  to initiate the next major.  (See notes above.)
*
*     01 Mar 1995: Realised that m8setj sets up the new rhs using xn!
*                  We should keep xn(1:nnjac) at xshort always.
*
*     04 Feb 1998: Initialize flin = objadd (not zero) to get correct
*                  printed objective value if there is no linear obj.
*
*     30 Jan 2004: (At GAMS) Take shorter step if user asks for it.
*     ------------------------------------------------------------------

      common    /m2parm/ dparm(30),iparm(30)
      common    /m3scal/ sclobj,scltol,lscale
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      common    /m8diff/ difint(2),gdummy,lderiv,lvldif,knowng(2)
      common    /cycle2/ objtru,suminf,numinf

      logical            square
      parameter        ( zero =  0.0d+0,  one = 1.0d+0 )

*     ------------------------------------------------------------------
*     On entry, the objective and its gradient are current.
*     So are the constraint functions and Jacobian, EXCEPT if
*     Lagrangian = No  (nlag = 0).
*     ------------------------------------------------------------------
      if (nlag .eq. 0) then
         if (lderiv .lt. 2)
     $      call m6dmmy( njac, gcon )

         call m6fcon( 2, nncon, nnjac, njac, fcon, gcon,
     $                ne, nka, ha, ka,
     $                xn, z, nwcore )

         if (lderiv .lt. 2)
     $      call m6dcon( nncon, nnjac, njac,
     $                   ne, nka, ha, ka,
     $                   fcon, fcon2, gcon, gcon2, xn, y, z, nwcore )
         if (ierr .ne. 0) go to 900
      end if

*     ------------------------------------------------------------------
*     Save previous primal and dual infeasibilities.
*     (Not used now.)
*     Define  square  to try and detect the nonlinear equations case.
*     ------------------------------------------------------------------
      Prinf0 = Prinf
      Duinf0 = Duinf
      square = sameSB  .and.  ns .eq. 0
      damp   = dparm(4)

*     ------------------------------------------------------------------
*     Choose a step according to dxrel and dlrel, the relative change
*     in xn(1:nnjac) and pi(1:nncon), as in 1982 paper.
*     In the nonlinear equations case, we can ignore dlrel.
*     ------------------------------------------------------------------
      step   = one
      dampx  = damp
      dampl  = damp
      tiny   = 1.0d-15
      step   =  min( step, dampx / (dxrel + tiny) )
      if (square) then
*        Relax
      else
         step = min( step, dampl / (dlrel + tiny) )
      end if

   50 if (step .ge. 0.99d+0) then
*        ---------------------------------------------------------------
*        Take a unit step.
*        ---------------------------------------------------------------
         step = one

      else
*        ---------------------------------------------------------------
*        Take a step < 1.0.
*        ---------------------------------------------------------------
         call dcopy ( nb,       xold, 1, xn, 1 )
         call daxpy ( nb, step, xdif, 1, xn, 1 )
         call dcopy ( m ,       xlam, 1, pi, 1 )
         call daxpy ( m , step, y2  , 1, pi, 1 )

*        ===============================================================
*        Evaluate the constraints and Jacobian.
*        ===============================================================
         if (lderiv .lt. 2)
     $      call m6dmmy( njac, gcon )

         call m6fcon( 2, nncon, nnjac, njac, fcon, gcon,
     $                ne, nka, ha, ka,
     $                xn, z, nwcore )
         if (ierr .ne. 0) go to 900

         if (lderiv .lt. 2)
     $      call m6dcon( nncon, nnjac, njac,
     $                   ne, nka, ha, ka,
     $                   fcon, fcon2, gcon, gcon2, xn, y, z, nwcore )
         if (ierr .ne. 0) go to 900

*        ===============================================================
*        Evaluate the nonlinear objective and gradient (if any).
*        Use rc as workspace for gobj2.
*        ===============================================================
         fobj   = zero
         if (nnobj .gt. 0) then
            call m6fobj( 2, nnobj, fobj, gobj, xn, z, nwcore )
            if (ierr .ne. 0) go to 900

            if (lderiv .ne. 1  .and.  lderiv .ne. 3) then
               call m6dobj( nnobj, fobj, gobj, rc, xn, z, nwcore )
               if (ierr .ne. 0) go to 900
            end if
         end if

         flin   = objadd
         if (jobj .gt. 0) flin = - xn(jobj) * sclobj  +  objadd
         objtru = flin + fobj
      end if

*     ------------------------------------------------------------------
*     Find the infeasibilities at the current step.
*     The linear constraints and bounds are assumed to be satisfied.
*     The primal infeasibility is therefore the maximum violation of
*                              the nonlinear constraints.
*     The dual   infeasibility is the departure from complementarity.
*     ------------------------------------------------------------------
      call m8rc  ( jobj, minimz, sclobj,
     $             m, n, nb, nnobj, nnobj0, nncon, nnjac, njac,
     $             ne, nka, a, ha, ka,
     $             gobj, gcon, pi, rc )

      call m8Cinf( n, nb, nncon, Cinf, iCinf,
     $             ne, nka, a, ha, ka,
     $             bl, bu, fcon, xn, y, z, nwcore )

      call m8Dinf( nb, bl, bu, rc, xn, Dinf, jDinf )
      Prinf  = Cinf
      Duinf  = Dinf

*     ------------------------------------------------------------------
*     Set xold = xn and xlam = pi.
*     Then decide what to do about xn(nnjac+1:nb).
*     nshort is a local variable to test strategies.
*     If nshort = nb   , we leave all of xn at xshort.
*     If nshort = nnjac, we leave only xn(1:nnjac) at the short step,
*     and move the remaining xn(*) back to step 1.0 for the next major.
*     This is tricky because we need to set xold = xn at the same time.
*     ------------------------------------------------------------------
      call dcopy ( m, pi, 1, xlam, 1 )

      if (step .eq. one) then
         call dcopy ( nb, xn, 1, xold, 1 )
      else
         nshort = nnjac
*---     nshort = nb

         call dcopy ( nshort, xn, 1, xold, 1 )

         do 500 j = nshort + 1, nb
            xoldj   = xold(j)
            xold(j) = xn(j)
            xn(j)   = xoldj + xdif(j)
  500    continue
      end if

*     ------------------------------------------------------------------
*     Decide whether the effect of damp should be tightened or relaxed.
*     >>> ntight isn't being used yet.
*     ------------------------------------------------------------------
      if (Prinf .le. Prinf0  .and.  Duinf .le. Duinf0) then
         ntight = ntight - 1
      else
         ntight = ntight + 1
      end if

      ntight = max( ntight, 1  )
      ntight = min( ntight, 10 )

*     ------------------------------------------------------------------
*     Decide if step 1.0 looks good enough to reduce the Penalty param.
*     If the infeasibilities are already as small as tolabs, it doesn't
*     matter if they go up slightly.  Similarly if they increase by a
*     small relative amount, tolrel.
*     (This is a hangover from the previous search strategy.
*     Will probably change it later.)
*     ------------------------------------------------------------------
      tolabs = 1.0d-1
      tolrel = 1.1d+0
      better = step  .eq.  one                         .and.
     $         Prinf .le. (Prinf0 + tolabs) * tolrel   .and.  
     $         Duinf .le. (Duinf0 + tolabs) * tolrel
      return

*     30 Jan 2004: (At GAMS) Take shorter step if user asks for it.

  900 if (ierr .eq. -1) then
         if (step .ge. 1.0d-5) then   ! <-- Smallest step allowed.
            ierr   = 0
            step   = 0.1d+0 * step
            go to 50
         end if
      end if

      end ! subroutine m8srch
