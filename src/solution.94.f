      subroutine solutionwu ( k, n, nb, ne, 
     &                      hs, xn, zsmall, lenz, inform )

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
     &                   iprint, isumm, ispecs, i, ii, j
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
*     User workspace: 1 + 1 + 1  + 1  + 1 + 1 + 4 + (p+1)*nobs + nobs + nobs
*                   nobs lam lam2 dist  k  prec b   x, y data    w      off
*     ------------------------------------------------------------------
      call miopti( 'Workspace (user) ', lenz, 0, 0, inform )
      call miopti( 'LOG FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'PRINT LEVEL ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FILE ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'SUPERBASICS LIMIT ', n+30, 0, 0, inform )
      call miopti( 'PROBLEM NUMBER ', 1, 0, 0, inform )
      call mioptr( 'FUNCTION PRECISION ', zsmall(6), 0, 0, inform )
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
      do i = 1, fk
         ha(i) = 1
         a(i) = zero
      end do
      if (k .gt. 0) then
         do i = 1, k
            ii = fk+(i-1)*2+1
            ha(ii) = i*2
            ha(ii+1) = i*2+1
            a(ii) = -one
            a(ii+1) = one
         end do
         do i = 1, k
            ii = fk+(k+i-1)*2+1
            ha(ii) = i*2
            ha(ii+1) = i*2+1
            a(ii) = -one
            a(ii+1) = -one
         end do
      end if
      do i = 1, fk
         ka(i) = i
      end do
      if (k .gt. 0) then
         j = fk
         ii = fk-1
         do i = 1, 2*k
            j = j+1
            ii = ii+2
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

      end ! subroutine solution

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj( mode, n, x, f, g, nstate, nprob, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            mode, n, nstate, nprob, nwcore
      double precision   x(n), f, g(n), z(nwcore)

      integer            nobs, k, p, method
      double precision   lam, lam2, dstr
*     ------------------------------------------------------------------
*     User workspace: 1 + 1 + 1  + 1  + 1 + 1 + 4 + (p+1)*nobs + nobs + nobs
*                   nobs lam lam2 dist  k  prec b   x, y data    w      off
*     ------------------------------------------------------------------
*
      mode = mode
      nstate = nstate
      nprob = nprob
      nobs = int(z(1))
      lam = z(2)
      lam2 = z(3)
      k = int(z(5))
      p = n - k

      if (nprob .eq. 1) then
         dstr = z(4)
         method = 0
         call subfunobjwu( n, x, f, g, z, nwcore, nobs, 
     $                   lam, lam2, dstr, k, p )
      else
         dstr = 0.0d+0
         method = int(z(4))
         call subfunobjcox( n, x, f, g, z, nwcore, nobs, 
     $                      lam, lam2, method, k, p )
      endif 

      return

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine subfunobjwu( n, x, f, g, z, nwcore,
     $                      nobs, lam, lam2, dstr, k, p )  

      integer            n, nwcore, nobs, k, p
      double precision   x(n), g(n), f, z(nwcore),
     &                   lam, lam2, dstr

      integer            i, j, ii, ii1, ii2
      double precision   eta(nobs), mu(nobs), resid(nobs),
     &                   xi(nobs), y(nobs), wt(nobs),
     &                   loglik, ddot, norm1, norm2
      double precision   nobsdble
*     -------------------------------------------------------------

*     optimizing  LogLik+ nobs* wt_i*|b_i|+ nobs*0.5*b_i**2
      
*     -------------------------------------------------------------
      double precision zero,       half,       one,       two    
      parameter       (zero=0.0d+0,half=5.0d-1,one=1.0d+0,two=2.0d+0)
*     -------------------------------------------------------------

      ii = 10 
      ii1 = 10 + p*nobs
      ii2 = 10 + (p+2)*nobs

*     
      ii3wu = 10 + (p+3)*nobs
      nobsdble=dble(nobs)
*     
      
      do 200 i = 1, nobs
         eta(i) = z(ii2 + i)
         do 100 j = 1, p
            eta(i) = eta(i) + x(j) * z(ii + (j-1)*nobs + i)
 100     continue
         y(i) = z(ii1 + i)
 200  continue
      ii = 10 + (p+1)*nobs
      loglik = zero
      do 300 i = 1, nobs
         wt(i) = z(ii + i)
         if (dstr .eq. zero) then
            mu(i) = eta(i)
            loglik = loglik - wt(i)*((y(i)-mu(i))**2)/two
         else if (dstr .eq. one) then 
            mu(i) = one/(one+exp(-eta(i)))
            loglik = loglik + wt(i)*(y(i)*eta(i)-log(one+exp(eta(i))))
         else if (dstr .eq. two) then
            mu(i) = exp(eta(i))
            loglik = loglik + wt(i)*(y(i)*eta(i) - mu(i))
         end if
         resid(i) = wt(i)*(y(i) - mu(i))
 300  continue
      norm2 = zero
      if (p .gt. 1) then
         do 330 j = 2, p
            norm2 = norm2 + x(j)**2
 330     continue
      end if
      norm1 = zero
      if (k .gt. 0) then 
         do 360 j = 1, k
*     
*     
* 
            norm1 = norm1 + x(p+j)*z(ii3wu+j)
 360     continue
      end if
*      f = -loglik + lam*norm1 + half*lam2*norm2
*       changed by WU
      f = -loglik + (norm1 + half*lam2*norm2)*nobsdble
      do 500 j = 1, p
         ii = 10 + (j-1)*nobs
         do 400 i = 1, nobs
            xi(i) = -1*z(ii + i)
 400     continue
         g(j) = ddot (nobs, xi, 1, resid, 1)
 500  continue
      if (p .gt. 1) then
         do 530 j = 2, p
            g(j) = g(j) + lam2*x(j)*nobsdble
 530     continue
      end if
      if (k .gt. 0) then
         do 560 j = 1, k
* 
* 
            g(p+j) = z(ii3wu+j)*nobsdble
 560     continue
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
