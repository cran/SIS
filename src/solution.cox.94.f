      subroutine coxsolution ( k, n, nb, ne,
     &                         hs, xn, zsmall, lenz, inform )

      implicit           double precision (a-h,o-z)
      integer            k, n, nb, ne, lenz, inform
      integer*4          hs(nb)
      double precision   xn(nb), zsmall(lenz)
*     ------------------------------------------------------------------
      integer            nwcore
      parameter          (nwcore = 10000000)
*     ------------------------------------------------------------------
      integer            m, p, fk, nname, nncon, nnobj, nnjac, iobj,
     &                   ka(n+1), name1, name2,
     &                   ns, mincor, ninf,
     &                   iprint, isumm, ispecs, i, ii, j, nobs
      integer*4          ha(ne)
      double precision   objadd, a(ne), bl(nb), bu(nb), pi(2*k+1),
     &                   rc(nb), sinf, obj, z(nwcore)
      character*8        names(5)
*     zero, one, infinity-----------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
*     ------------------------------------------------------------------

      iprint = 0   ! The MINOS PRINT   file.
      isumm  = 0   ! The MINOS SUMMARY file.
      ispecs = 0   ! The MINOS SPECS   file.

      call mistart( iprint, isumm, ispecs )  ! Initialize MINOS and open
*     ------------------------------------------------------------------
*       1 + 1  + 1 + 1  + 1 + 1 + 1 + 3 + (p+1)*nobs+nobs+nobs+nobs+nobs
*     nobs lam lam2 mthd  k  lp  prec b   x, y data   d   rep  eta  wsum
*     ------------------------------------------------------------------
      call miopti( 'Workspace (user) ', lenz, 0, 0, inform )
      call miopti( 'LOG FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'PRINT LEVEL ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FILE ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'SUPERBASICS LIMIT ', n+30, 0, 0, inform )
      call miopti( 'PROBLEM NUMBER ', 2, 0, 0, inform )
      call mioptr( 'FUNCTION PRECISION ', zsmall(7), 0, 0, inform )
*     ------------------------------------------------------------------
*     Now set parameters for moniss
*     ------------------------------------------------------------------
      do i = 1, lenz
         z(i) = zsmall(i)
      end do
      p = n - k
      fk = p - k
      m = 2*k + 1
      nname = 1
      nncon = 0
      nnobj = n
      nnjac = 0
      iobj = 0
      objadd = 0.0d+0
      if (fk .gt. 0) then
         do i = 1, fk
            ha(i) = 1
            a(i) = zero
         end do
      end if
      if (k .gt. 0) then
         do i = 1, k
            ii = fk+(i-1)*3+1
            ha(ii) = 1
            ha(ii+1) = i*2
            ha(ii+2) = i*2+1
            a(ii) = zero
            a(ii+1) = -one
            a(ii+2) = one
         end do
         do i = 1, k
            ii = fk+(k+i-1)*3+1
            ha(ii) = 1
            ha(ii+1) = i*2
            ha(ii+2) = i*2+1
            a(ii) = zero
            a(ii+1) = -one
            a(ii+2) = -one
         end do
      end if
      if (fk .gt. 0) then
         do i = 1, fk
            ka(i) = i
         end do
      end if
      if (k .gt. 0) then
         j = fk
         ii = fk-2
         do i = 1, 2*k
            j = j+1
            ii = ii+3
            ka(j) = ii
         end do
      end if
      ka(n+1) = ne+1
      do i = 1, p
         bl(i) = bminus
         bu(i) = bplus
      end do
      do i = (p+1), nb
         bl(i) = zero
         bu(i) = bplus
      end do
      ns = 0
      do i = 1, n
         if (hs(i) .eq. 2) then
            ns = ns + 1
         end if
      end do

      call minoss( 'Warm', m, n, nb, ne, nname,
     $             nncon, nnobj, nnjac,
     $             iobj, objadd, names,
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc,
     $             inform, mincor, ns, ninf, sinf, obj,
     $             z, nwcore )

      nobs = int(z(1))
      ii = 10 + (p+3)*nobs
      do i = 1, 2*nobs
         zsmall(ii+i) = z(ii+i)
      end do
      zsmall(6) = z(6)

      end ! subroutine coxsolution

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine subfunobjcox( n, x, f, g, z, nwcore,
     $                         nobs, lam, lam2, method, kk, p )

      integer            n, nwcore, nobs, kk, p, method
      double precision   x(n), g(n), f, z(nwcore), lam, lam2

      integer            i, j, k, ii, ii0, ii1, ii2, ii3, ii4, ii5wu
      double precision   d(nobs), rept(nobs),
     &                   eta(nobs), mu(nobs),
     &                   wx(p), dwx(p), wsum, dwsum, xk(p),
     &                   ddot, norm1, norm2
      double precision   nobsdble
*     -------------------------------------------------------------
      double precision   zero,        half,        one
      parameter         (zero=0.0d+0, half=5.0d-1, one=1.0d+0)
*     -------------------------------------------------------------
*     ------------------------------------------------------------------
*       1 + 1  + 1 + 1  + 1 + 1 + 1 + 3 + (p+1)*nobs+nobs+nobs+nobs+nobs+ k
*     nobs lam lam2 mthd  k  lp  prec b   x, y data   d   rep  eta  wsum  wt
*     ------------------------------------------------------------------

      ii = 10
      ii0 = ii + p*nobs
      ii1 = ii0 + nobs
      ii2 = ii1 + nobs
      ii3 = ii2 + nobs
      ii4 = ii3 + nobs
      ii5wu = ii4 + nobs
      nobsdble = dble(nobs)
      do 200 i = 1, nobs
         d(i) = z(ii1 + i)
         rept(i) = z(ii2 + i)
         do 100 j = 1, p
            xk(j) = z(ii + (j-1)*nobs + i)
 100     continue
         mu(i) = ddot(p, x, 1, xk, 1)
         eta(i) = exp(mu(i))
         z(ii3+i) = eta(i)
 200  continue
      f = zero
      do 300 j = 1, p
         g(j) = zero
 300  continue
      ii1 = 0
      do 600 i = 1, nobs
         if (d(i) .eq. one) then
            if (method .eq. 1) then
               if (rept(i) .ne. zero .and. ii1 .eq. 0) then
                  ii1 = int(rept(i)) - 1
                  do 320 j = 1, p
                     wx(j) = zero
 320              continue
                  wsum = zero
                  do 360 k = 1, (i + ii1)
                     do 340 j = 1, p
                        xk(j) = z(ii + (j-1)*nobs + k)
                        wx(j) = wx(j) + xk(j)*eta(k)
 340                 continue
                     wsum = wsum + eta(k)
 360              continue
                  do 380 j = 1, p
                     wx(j) = wx(j)/wsum
 380              continue
               else if (ii1 .gt. 0) then
                  ii1 = ii1 - 1
               end if
            else if (method .eq. 2) then
               if (rept(i) .ne. zero .and. ii1 .eq. 0) then
                  ii0 = int(rept(i))
                  ii1 = ii0 - 1
                  do 410 j = 1, p
                     wx(j) = zero
                     dwx(j) = zero
 410              continue
                  wsum = zero
                  dwsum = zero
                  do 430 k = 1, (i + ii1)
                     do 420 j = 1, p
                        xk(j) = z(ii + (j-1)*nobs + k)
                        wx(j) = wx(j) + xk(j)*eta(k)
                        if (ii1 .gt. 0 .and. k .ge. i) then
                           dwx(j) = dwx(j) + xk(j)*eta(k)
                        end if
 420                 continue
                     wsum = wsum + eta(k)
                     if (ii1 .gt. 0 .and. k .ge. i) then
                        dwsum = dwsum + eta(k)
                     end if
 430              continue
                  do 440 j = 1, p
                     wx(j) = wx(j)/wsum
 440              continue
               else if (ii1 .gt. 0) then
                  do 450 j = 1, p
                     wx(j) = wx(j) * wsum
                     wx(j) = wx(j) - dwx(j)/ii0
 450              continue
                  wsum = wsum - dwsum/ii0
                  do 460 j = 1, p
                     wx(j) = wx(j)/wsum
 460              continue
                  ii1 = ii1 - 1
               end if
            end if
            z(ii4+i) = wsum
            do 500 j = 1, p
               g(j) = g(j) - z(ii + (j-1)*nobs + i)
               g(j) = g(j) + wx(j)
 500        continue
            f = f - mu(i) + log(wsum)
         end if
 600  continue
      z(6) = -f
      norm1 = zero
      if (kk .gt. 0) then
         do 700 j = 1, kk
            norm1 = norm1 + x(p+j)*z(ii5wu+j)
 700     continue
      end if
      norm2 = zero
      do 730 j = 1, p
         norm2 = norm2 + x(j)**2
 730  continue
      f = f + (norm1+half*lam2*norm2)*nobsdble
      do 750 j = 1, p
         g(j) = g(j) + lam2*x(j)*nobsdble
 750  continue
      if (kk .gt. 0) then
         do 780 j = 1, kk
            g(p+j) = z(ii5wu+j)*nobsdble
 780     continue
      end if
      return

      end

******************************************************************************
*     inform  says what happened; see Chapter 6.3 of the User's Guide.       *
*             A summary of possible values follows:                          *
*                                                                            *
*             inform   Meaning                                               *
*                                                                            *
*                0     Optimal solution found.                               *
*                1     The problem is infeasible.                            *
*                2     The problem is unbounded (or badly scaled).           *
*                3     Too many iterations.                                  *
*                4     Apparent stall.  The solution has not changed         *
*                      for a large number of iterations (e.g. 1000).         *
*                5     The Superbasics limit is too small.                   *
*                6     Subroutine funobj or funcon requested termination     *
*                      by returning mode < 0.                                *
*                7     Subroutine funobj seems to be giving incorrect        *
*                      gradients.                                            *
*                8     Subroutine funcon seems to be giving incorrect        *
*                      gradients.                                            *
*                9     The current point cannot be improved.                 *
*               10     Numerical error in trying to satisfy the linear       *
*                      constraints (or the linearized nonlinear              *
*                      constraints).  The basis is very ill-conditioned.     *
*               11     Cannot find a superbasic to replace a basic           *
*                      variable.                                             *
*               12     Basis factorization requested twice in a row.         *
*                      Should probably be treated as inform = 9.             *
*               13     Near-optimal solution found.                          *
*                      Should probably be treated as inform = 9.             *
*                                                                            *
*               20     Not enough storage for the basis factorization.       *
*               21     Error in basis package.                               *
*               22     The basis is singular after several attempts to       *
*                      factorize it (and add slacks where necessary).        *
*                                                                            *
*               30     An OLD BASIS file had dimensions that did not         *
*                      match the current problem.                            *
*               32     System error.  Wrong number of basic variables.       *
*                                                                            *
*               40     Fatal errors in the MPS file.                         *
*               41     Not enough storage to read the MPS file.              *
*               42     Not enough storage to solve the problem.              *
*                                                                            *
******************************************************************************
