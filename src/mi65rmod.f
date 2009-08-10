************************************************************************
*
*     File  rmod fortran.
*
*     m6bfgs   m6bswp   m6radd   m6rcnd   m6rdel
*     m6rmod   m6rset   m6rsol   m6swap
*
* 27 Jun 2003: Eliminated dnrm2  from m6bfgs.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6bfgs( maxR, nS, lenR, R, g, g2, p, v,
     &     step, told, tolz, inform )

      implicit
     &     none
      integer
     &     inform, lenR, maxR, nS
      double precision
     &     step, told, tolz, g(nS), g2(nS), p(nS), R(lenR), v(nS)

*     ==================================================================
*     m6bfgs  applies the BFGS update to the upper-triangular matrix R,
*     which holds the Cholesky factor of the quasi-Newton approximation
*     of the reduced Hessian.
*
*     R contains a triangle of size nR = min( nS, maxR ).
*     If nS .gt. maxR, R also contains a diagonal of size nS - maxR.
*
*     p       holds the search direction.  It is overwritten.
*     v       must satisfy  R' v = g.      It is overwritten.
*
*     On exit,
*     inform = 0  if no update was performed,
*            = 1  if the update was successful,
*            = 2  if it was nearly singular.
*
*            1983: First version of m6bfgs.
*     12 Jun 2001: First version of SNOPT routine s6Rqn based on m6bfgs.
*     16 Jul 2001: Version of s6Rqn used for new m6bfgs.
*     23 Jul 2001: s6Rqn  modified to become new m6bfgs.
*     ==================================================================
      integer
     &     j, j1, l, lastv, nR
      double precision
     &     d, ddot, delta1, delta2, gtp, gtp2, vlast, vnz
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      inform = 0
      nR     = min( nS, maxR )
      gtp    = ddot  ( nS, g , 1, p, 1 )
      gtp2   = ddot  ( nS, g2, 1, p, 1 )
      if (gtp2 .le. 0.91d+0*gtp) return

      delta1 = one / sqrt( abs( gtp ) )
      delta2 = one / sqrt( step*(gtp2 - gtp) )

*     Normalize  v  and change its sign.

      call dscal ( nS, (- delta1), v, 1 )

*     Avoid cancellation error in forming the new vector  p.

      if ( abs( delta1/delta2 - one ) .ge. 0.5d+0) then
         do j = 1, nS
            p(j) = delta2*( g2(j) - g(j) )  +  delta1*g(j)
         end do
      else
         d = delta1 - delta2
         do j = 1, nS
            p(j) = delta2*g2(j)  +  d*g(j)
         end do
      end if

*     Apply the update in the form  R + v p',  where  v  is held
*     in  v  and  vlast.

      if (nS .gt. maxR) then
       ! vlast = dnrm2 ( nS-maxR, v(maxR+1), 1 )
         vlast = ddot  ( nS-maxR, v(maxR+1), 1, v(maxR+1), 1 )
         vlast = sqrt( vlast )
      else
         vlast = zero 
      end if

*     ---------------------------------------------
*     Find the last nonzero in v (including vlast).
*     ---------------------------------------------
      vnz   = vlast
      lastv = nR + 1

*+    while (lastv .gt. 1  .and.  vnz .le. tolz) do
  100 if    (lastv .gt. 1  .and.  vnz .le. tolz) then
         lastv = lastv - 1
         vnz   = abs(v(lastv)) 
         go to 100
*+    end while
      end if

*     Triangularize   R  +  v p'.

      call m6rmod( maxR, nR, lenR, R, v, p, lastv, vlast,
     &             told, tolz, inform )

*     Deal with surplus diagonals of  R.

      if (nS .gt. maxR) then
         j1   = maxR + 1
         l    = maxR*j1/2
         do j = j1, nS
            l = l + 1
            R(l) = sqrt(R(l)**2 + (p(j) - delta1*g(j))**2
     &                                  -(delta1*g(j))**2)
         end do
      end if

      end ! subroutine m6bfgs

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6bswp( maxr, n, nr, r, v, w, lastv,
     $                   told, tolz, inform )

      implicit           double precision (a-h,o-z)
      double precision   r(nr), v(n), w(n)

*     ------------------------------------------------------------------
*     m6bswp  modifies the upper-triangular matrix R
*     to account for a basis exchange in which the lastv-th
*     superbasic variable becomes basic.  R is changed to
*     R + vw', which is triangularized by r1mod,
*     where v is the lastv-th column of R, and w is input.
*
*     n       is the current size of R.  It is always less than maxr.
*
*     17 Apr 1994: Converted to row-wise storage.
*     18 Sep 2001: vlast = zero needed for calling 21 Jul 2001 version
*                  of m6rmod (= s6Rmod from SNOPT).
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------

*     Set v = lastv-th column of R and find its norm.

      lr     = lastv
      incr   = maxr - 1

      do 100 i = 1, lastv
         v(i)  = r(lr) 
         lr    = lr   + incr
         incr  = incr - 1
  100 continue

      vnorm  = dasum ( lastv, v, 1 )
      vlast  = zero
      call m6rmod( maxr, n, nr, r, v, w, lastv, vlast,
     $             (vnorm*told), (vnorm*tolz), inform )

      end ! subroutine m6bswp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6radd( maxr, nr, ns, r )

      implicit           double precision (a-h,o-z)
      double precision   r(nr)

*     ------------------------------------------------------------------
*     m6radd  adds column ns to the upper triangular matrix R.
*     In this version (Sep 1984) it is just a unit vector.
*
*        Jan 1983: Modified to add a diagonal only, if ns > maxr.
*     17 Apr 1994: Converted to row-wise storage.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      if (ns .le. maxr) then
         lr     = ns
         incr   = maxr
         do 100 k = 1, ns - 1
            r(lr) = zero
            incr  = incr - 1
            lr    = lr   + incr
  100    continue
      else
         lr     = maxr*(maxr + 1)/2  +  (ns - maxr)
      end if

      r(lr) = one

      end ! subroutine m6radd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6rcnd( maxr, nr, ns, r, dmax, dmin, cond )

      implicit           double precision (a-h,o-z)
      double precision   r(nr)

*     ------------------------------------------------------------------
*     m6rcnd  finds the largest and smallest diagonals of the
*     upper triangular matrix R, and returns the square of their ratio.
*     This is a lower bound on the condition number of R'R.
*
*     17 Apr 1994: Converted to row-wise storage.
*     ------------------------------------------------------------------

      ncolr  = min( maxr, ns )
      dmax   = abs( r(1) )
      dmin   = dmax
      lr     = 1
      incr   = maxr

      do 100 j = 2, ncolr
         lr    = lr   + incr
         incr  = incr - 1
         d     = abs( r(lr) )
         dmax  = max( dmax, d )
         dmin  = min( dmin, d )
  100 continue

      cond   = (dmax / dmin)**2

      end ! subroutine m6rcnd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6rdel( m, maxr, nr, ns, ms,
     $                   kb, bbl, bbu, grd, r, rg, x, jq, rset )

      implicit           double precision (a-h,o-z)
      integer            kb(ms)
      double precision   bbl(ms), bbu(ms), grd(ms)
      double precision   r(nr), rg(ns), x(ms)
      logical            rset

*     ------------------------------------------------------------------
*     m6rdel  deletes the jq-th superbasic variable from R
*     and from various arrays kb, bbl, bbu, grd, rg, x.
*     The dimension of R decreases from ns to ns-1.
*
*     17 Apr 1994: Converted to row-wise storage.
*     30 Jul 1994: Bottom part of R moved north-west AFTER the sweep
*                  of rotations has eliminated the jq-th row.
*                  This usually means double handling of that part of R,
*                  but it's the only way to skip identity rotations.
*                  (Who knows if there are any.)
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy

      parameter        ( zero = 0.0d+0 )

      if (jq .eq. ns) return
      if (.not. rset) go to 500

*     ------------------------------------------------------------------
*     Delete the jq-th column of R from the top rows.
*     For the first jq-1 rows, elements R(*,jq+1:ncolr) of each row
*     are shifted 1 place to the left.
*     ------------------------------------------------------------------
      ncolr  = min( maxr, ns )
      lr     = jq
      incr   = maxr
      nmove  = ncolr - jq

      do 40 i = 1, jq - 1
         do 20 k = lr, lr + nmove - 1
            r(k) = r(k+1)
   20    continue
         incr = incr - 1
         lr   = lr   + incr
   40 continue

*     ------------------------------------------------------------------
*     Triangularize the remaining rows of R,
*     using a partial forward sweep of rotations.
*
*     x x x x x x     becomes   x x x x x x  
*       x x x x x                 x x x x x  
*         . - - - -                 . 0 0 0 0  
*           x x x x                   + + + +
*             x x x                     + + +
*               x x                       + +
*                 x                         +
*         |                         |
*        jq                        jq
*
*     The . is not touched because it is later overwritten.
*     ls marks the - being eliminated.
*     lr marks the start of the next + + + row.
*     ------------------------------------------------------------------
      lrsav  = lr
      insav  = incr
      ls     = lr
      tolz   = eps0

      do 100 j = jq + 1, ncolr
         ls    = ls   + 1
         lr    = lr   + incr
         incr  = incr - 1
         b     = r(ls)
         if (abs( b ) .le. tolz) go to 100
         a     = r(lr)
         diag  = sqrt( a**2 + b**2 )
         r(lr) = diag

         if (j .lt. ncolr) then
            cs    = a / diag
            sn    = b / diag
            jr    = lr
            js    = ls

            do  50 k = j + 1, ncolr
               jr    = jr + 1
               js    = js + 1
               rk    = r(jr)
               sk    = r(js)
               r(jr) = cs * rk  +  sn * sk
               r(js) = sn * rk  -  cs * sk
   50       continue
         end if
  100 continue

*     ------------------------------------------------------------------
*     Shift the + + + triangle up and left.
*     lr marks the start of each + + + row being moved.
*     ls marks the start of its final position.
*     ------------------------------------------------------------------
      lr     = lrsav
      incr   = insav
      nmove  = ncolr - jq

      do 200 j = jq + 1, ncolr
         ls    = lr
         lr    = lr + incr
         call dcopy ( nmove, r(lr), 1, r(ls), 1 )
         incr  = incr  - 1
         nmove = nmove - 1
  200 continue

*     ------------------------------------------------------------------
*     Deal with surplus diagonals of R.
*     ------------------------------------------------------------------
      if (ns .gt. maxr) then
         if (jq .le. maxr) then

*           Clear out the last column of R.

            lr     = maxr
            incr   = maxr
            do 250 k = 1, maxr
               r(lr) = zero
               incr  = incr - 1
               lr    = lr   + incr
  250       continue
         end if

*        Shift surplus diagonals of R to the left.

         k      = max( maxr, jq )
         lr     = maxr*(maxr + 1)/2  +  (k - maxr)
         nmove  = ns - k
         do 300 k = lr, lr + nmove - 1
            r(k)  = r(k+1)
  300    continue
      end if

*     ------------------------------------------------------------------
*     Shift all the arrays one place to the left.
*     ------------------------------------------------------------------
  500 do 550 j  = jq, ns - 1
         rg(j)  = rg(j+1)
         k      = m + j
         kb(k)  = kb(k+1)
         bbl(k) = bbl(k+1)
         bbu(k) = bbu(k+1)
         grd(k) = grd(k+1)
         x(k)   = x(k+1)
  550 continue

      end ! subroutine m6rdel

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6rmod( maxR, nS, lenR, R, u, v, lastnz, ulast,
     &                   told, tolz, inform )

      implicit
     &     none
      integer
     &     inform, lenR, lastnz, maxR, nS
      double precision
     &     told, tolz, ulast, R(lenR), u(nS), v(nS)

*     ==================================================================
*     m6rmod  modifies the (nS+1) x nS upper-triangular matrix R so that
*     Q*(R + u v') is upper triangular,  where  Q  is orthogonal.
*     The arrays v and u hold v and the first nS elements of u 
*     respectively.  The (nS+1)th component of u is held in ulast.
*     The new R overwrites the old.
*
*     Q is the product of two sweeps of plane rotations (not stored).
*     These affect the (lastnz)th row of R, which is temporarily held
*     in the array u.  Thus, u is overwritten.  v is not altered.
*
*     ulast  holds  u(nS+1).   It is overwritten. 
*
*     lastnz points to the last nonzero of u.  The value lastnz = nS+1
*            would always be ok, but sometimes it is known to be less
*            than nS+1, so Q reduces to two partial sweeps of rotations.
*
*     told   is a tolerance on the lastv-th diagonal of R.
*     tolz   is a tolerance for negligible elements in  u.
*
*     On exit,
*     inform = 1  if the diagonal of R is larger than told,
*            = 2  if not (the diagonal is not modified).
*
*     06 Sep 1991: First version of SNOPT routine s6Rmod
*                  based on Minos routine m6rmod.
*     11 Sep 1994: s6Rmod modified to update a principal submatrix of R.
*     03 Dec 2000: s6Rmod converted to row-wise storage.
*     17 Jun 2001: Version of s6Rmod used for new m6rmod. 
*     21 Jul 2001: SNOPT routine s6Rmod becomes new m6rmod.
*                  ulast is a new output parameter.
*     ==================================================================
      integer
     &     i, incR, j, l, lastR, ldiag, lm1, nmove
      double precision
     &     root, cs, s, sn, t, u2
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------

      if (lastnz .le. nS) then
         ulast = u(lastnz)
      end if

*     Copy the (lastnz)th row of R into the end of u.

      lm1    = lastnz - 1
      lastR  = lm1*maxR  +  (3 - lastnz)*lastnz/2
      nmove  = nS - lm1
      if (nmove .gt. 0) call dcopy ( nmove, R(lastR), 1, u(lastnz), 1 )

*     ------------------------------------------------------------------
*     Reduce u to a multiple of e(lastnz) using a partial backward sweep
*     of rotations.  This fills in the (lastnz)th row of R (held in u).
*     ------------------------------------------------------------------
      if (lastnz .gt. 1) then
         u2     = ulast**2
         ldiag  = lastR
         incR   = maxR - lm1

         do i = lm1, 1, -1
            incR  = incR  + 1
            ldiag = ldiag - incR
            s     = u(i)
            u(i)  = zero
            if (abs(s) .gt. tolz) then
               u2    = s**2 + u2
               root  = sqrt(u2)
               cs    = ulast/root
               sn    = s    /root
               ulast = root
               l     = ldiag

               do j = i, nS
                  s    = u(j)
                  t    = R(l)
                  u(j) = cs*s + sn*t
                  R(l) = sn*s - cs*t
                  l    = l + 1
               end do
            end if
         end do
      end if

      call daxpy ( nS, ulast, v, 1, u, 1 )      ! Set u = u  +  ulast*v.

*     ------------------------------------------------------------------
*     Eliminate the front of the (lastnz)th row of R (held in u) using a
*     partial forward sweep of rotations.
*     ------------------------------------------------------------------
      if (lastnz .gt. 1) then
         ldiag  = 1
         incR   = maxR

         do i = 1, lm1
            t     = u(i)
            if (abs(t) .gt. tolz) then
               s        = R(ldiag)
               root     = sqrt(s**2 + t**2)
               cs       = s/root
               sn       = t/root
               R(ldiag) = root
               l        = ldiag

               do j = i+1, nS
                  l    = l + 1
                  s    = R(l)
                  t    = u(j)
                  R(l) = cs*s + sn*t
                  u(j) = sn*s - cs*t
               end do
            end if
            ldiag = ldiag + incR
            incR  = incR  - 1
         end do
      end if

*     Insert the new (lastnz)th row of  R.

      if (nmove .gt. 0) then
         call dcopy ( nmove, u(lastnz), 1, R(lastR), 1 )

*        ---------------------------------------------------------------
*        Test for (unlikely) singularity.
*        ---------------------------------------------------------------
         inform = 1
         if (abs( R(lastR) ) .le. told) then
            inform   = 2
*           r(lastr) = told
         end if
      end if

      end ! subroutine m6rmod

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6rset( maxr, nr, ns, r, w, cond )

      implicit           double precision (a-h,o-z)
      double precision   r(nr), w(nr)

*     ------------------------------------------------------------------
*     m6rset  alters R, the upper-triangular factor of the
*     approximation to the reduced hessian.
*
*     If  r(1) = zero, R does not exist.
*     In this case, r is initialized to the identity matrix.
*
*     Otherwise, R already exists and we attempt to make it better
*     conditioned by scaling its columns by the square roots of its
*     current diagonals.
*
*     09 Jun 1994: Converted to row-wise storage.
*                  Parameter w added for workspace.
*     14 Jul 2001: Bug report from Philip Gill.
*                  Last element of diagonal part was not being reset.
*                  Fixed by updating ldiag at end of loop 120, not start.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm

      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      cond   = one
      ncolr  = min( maxr, ns )
      if (ncolr .eq.  0  ) return

      dmax   = abs( r(1) )
      if (dmax  .eq. zero) then

*        Set R = the identity.

         ldiag  = 1
         incr   = maxr
         nzero  = ncolr - 1

         do 100 k = 1, ncolr - 1
            r(ldiag) = one
            call dload ( nzero, zero, r(ldiag+1), 1 )
            ldiag    = ldiag + incr
            incr     = incr  - 1
            nzero    = nzero - 1
  100    continue

         r(ldiag) = one

         do 120 k = maxr + 1, ns
            ldiag    = ldiag + 1
            r(ldiag) = one
  120    continue
      else
*        ---------------------------------------------------------------
*        Scale the columns of R.
*        ---------------------------------------------------------------

*        Find dmin and dmax, and set w = set of scale factors.

         dmin   = dmax
         ldiag  = 1
         incr   = maxr

         do 200 k = 1, ncolr
            diag  = abs( r(ldiag) )
            dmax  = max( dmax, diag )
            dmin  = min( dmin, diag )
            w(k)  = one / sqrt( diag )
            ldiag = ldiag + incr
            incr  = incr  - 1
  200    continue

*        Apply column scales to each row of R.

         ldiag  = 1
         incr   = maxr

         do 250 k = 1, ncolr
            l     = ldiag
            do 240 j = k, ncolr          
               r(l)  =  r(l) * w(j)
               l     =  l + 1
  240       continue
            ldiag    = ldiag + incr
            incr     = incr  - 1
  250    continue

         cond   = dmax / dmin
         if (iprint .gt. 0) write(iprint, 2000) cond
      end if
      return

 2000 format(' Hessian modified.  CondR =', 1p, e8.1)

      end ! subroutine m6rset

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6rsol( mode, maxr, nr, ns, r, y )

      implicit           double precision (a-h,o-z)
      double precision   r(nr), y(ns)

*     ------------------------------------------------------------------
*     m6rsol  solves Rx = y or R'x = y, where R is an
*     upper-triangular matrix of dimension ns stored by rows in r(nr).
*     The solution x overwrites y.
*
*     29 Jul 1994: Converted to row-wise storage.
*     ------------------------------------------------------------------

      ncolr  = min( maxr, ns )

      if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        mode = 1  --  solve Ry = y.
*        ---------------------------------------------------------------
         lr       = (ncolr - 1)*maxr + (3 - ncolr)*ncolr/2
         y(ncolr) = y(ncolr) / r(lr)
         incr     = maxr + 1 - ncolr
         numr     = 0

         do 100 j  = ncolr-1, 1, -1
            numr   = numr + 1
            incr   = incr + 1
            lr     = lr   - incr
            sum    = y(j) - ddot  ( numr, r(lr+1), 1, y(j+1), 1 )
            y(j)   = sum / r(lr)
  100    continue
      else
*        ---------------------------------------------------------------
*        mode = 2  --  solve R'y = y.
*        ---------------------------------------------------------------
         lr    = 1
         incr  = maxr
         numr  = ncolr - 1

         do 200 i = 1, ncolr-1
            y(i)  = y(i) / r(lr)
            call daxpy ( numr, (- y(i)), r(lr+1), 1, y(i+1), 1 )
            lr    = lr   + incr
            incr  = incr - 1
            numr  = numr - 1
  200    continue

         y(ncolr) = y(ncolr) / r(lr)
      end if

*     Deal with surplus diagonals of R.

      if (ns .gt. maxr) then
         l     = maxr*(maxr + 1)/2
         do 320 j = maxr + 1, ns
            l     = l + 1
            y(j)  = y(j) / r(l)
  320    continue
      end if

      end ! subroutine m6rsol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m6swap( m, maxr, nr, ns, ms,
     $                   kb, bbl, bbu, grd, r, rg, x )

      implicit           double precision (a-h,o-z)
      integer            kb(ms)
      double precision   bbl(ms), bbu(ms), grd(ms)
      double precision   r(nr), rg(ns), x(ms)

*     ------------------------------------------------------------------
*     m6swap  (superbasic swap)  finds the largest reduced gradient in
*     the range  rg(maxr+1), ..., rg(ns)  and swaps it into position
*     maxr + 1  (so we know where it is).
*     ------------------------------------------------------------------

      k1     = maxr + 1
      if (ns .gt. k1) then
         nz2    = ns - maxr
         k2     = maxr + idamax( nz2, rg(k1), 1 )
         if (k2 .gt. k1) then
            j      = m + k1
            k      = m + k2
            lastr  = maxr*k1/2
            ldiag1 = lastr + 1
            ldiag2 = lastr + (k2 - maxr)

            rdiag1    = r(ldiag1)
            rg1       = rg(k1)
            j1        = kb(j)
            bl1       = bbl(j)
            bu1       = bbu(j)
            grd1      = grd(j)
            x1        = x(j)

            r(ldiag1) = r(ldiag2)
            rg(k1)    = rg(k2)
            kb(j)     = kb(k)
            bbl(j)    = bbl(k)
            bbu(j)    = bbu(k)
            grd(j)    = grd(k)
            x(j)      = x(k)

            r(ldiag2) = rdiag1
            rg(k2)    = rg1
            kb(k)     = j1
            bbl(k)    = bl1
            bbu(k)    = bu1
            grd(k)    = grd1
            x(k)      = x1
         end if
      end if

      end ! subroutine m6swap
