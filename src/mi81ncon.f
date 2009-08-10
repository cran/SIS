************************************************************************
*
*     File  mi81ncon fortran.
*
*     m8aug1
*
*     12 Jul 2000: mi81ncon.f now separates m8aug1 from mi80ncon.f.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m8aug1( mode, ms, nncon, nnjac, njac, n, nb, inform,
     $                   ne, nka, a, ha, ka,
     $                   hs, kb, bl, bu, bbl, bbu, blslk, buslk,
     $                   gcon, gcon2, xn )

      implicit           double precision (a-h,o-z)
      integer            ha(ne), hs(nb)
      integer            ka(nka), kb(ms)
      double precision   a(ne), bl(nb), bu(nb), bbl(ms), bbu(ms),
     $                   blslk(nncon), buslk(nncon),
     $                   gcon (njac ), gcon2(njac ),
     $                   xn(nb)

*     ------------------------------------------------------------------
*     m8aug1  does the work for  m8augl.
*     11 Sep 1992: hs added as parameter for mode 5, to allow
*                  nonbasic slacks to be moved onto the relaxed bounds.
*     04 Mar 1994: Relax big bounds more (using relative and abs tols).
*     28 Nov 1996: mode = -1 and 0 added.
*                  The Jacobian part of a(*) is now random until the
*                  linear constraints are satisfied.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /m8al2 / radius,rhsmod,modpen,modrhs

      parameter        ( zero = 0.0d+0 )

      inform = 0
      bplus  = 0.9d+0*plinfy
      bminus = - bplus

      if (mode .eq. -2) then
*        ---------------------------------------------------------------
*        We want a(*) to be random until the linear constraints are
*        satisfied.  (Scale and Crash will use the random numbers.)
*        Generate random numbers in gcon and then load them into a(*).
*
*        Sometimes a(*) contains constant Jacobian elements.
*        Hence, save them in gcon and gcon2 before overwriting a(*).
*
*        m8rand generates numbers in the range (0,1).
*        To stop them being too small, we add 0.5 to the random a(*).
*        ---------------------------------------------------------------
         ia = 1547
         ib = 2671
         ic = 3770
         call m8rand( ia, ib, ic, njac, gcon, 1 )

         l      = 0
         do 40 j = 1, nnjac
            do 20 k = ka(j), ka(j+1) - 1
               ir       = ha(k)
               if (ir .gt. nncon) go to 40
               l        = l + 1
               ajac     = a(k)
               a(k)     = gcon(l) + 0.5d+0
               gcon(l)  = ajac
               gcon2(l) = ajac
   20       continue
   40    continue

      else if (mode .eq. -1) then
*        ---------------------------------------------------------------
*        Start of Major iteration 1 or 2.
*        Save the bounds on the nonlinear rows.
*        ---------------------------------------------------------------
         call dcopy ( nncon, bl(n+1), 1, blslk, 1 )
         call dcopy ( nncon, bu(n+1), 1, buslk, 1 )

      else if (mode .eq. 0) then
*        ---------------------------------------------------------------
*        Start of the Major iteration 1.
*        Initialize  modpen, rhsmod.
*        Make the nonlinear rows free (relax the bounds to infinity).
*        ---------------------------------------------------------------
         modpen = 0
         rhsmod = zero
         call dload ( nncon, (- plinfy), bl(n+1), 1 )
         call dload ( nncon,    plinfy , bu(n+1), 1 )

      else if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        Copy the Jacobian into  A  to form the next linearization.
*        ---------------------------------------------------------------
         l     = 0
         do 140 j = 1, nnjac
            do 120 k = ka(j), ka(j+1) - 1
               ir    = ha(k)
               if (ir .gt. nncon) go to 140
               l     = l + 1
               a(k)  = gcon(l)
  120       continue
  140    continue

      else if (mode .eq. 2) then
*        ---------------------------------------------------------------
*        Make sure the initial xn values are feasible.
*        Do the slacks too, in case the RHS has been perturbed.
*        ---------------------------------------------------------------
         do 220 j = 1, nb
            xn(j) = max( xn(j), bl(j) )
            xn(j) = min( xn(j), bu(j) )
  220    continue

      else if (mode .eq. 3) then
*        ---------------------------------------------------------------
*        This is a call from m5solv, when the linear constraints are
*        satisfied for the first time.
*        gcon2 contains the initial (unscaled) Jacobian,
*        possibly including constant Jacobian elements.
*        Copy them into a(*).   (gcon should already contain them.)
*        ---------------------------------------------------------------
         l     = 0
         do 340 j = 1, nnjac
            do 320 k = ka(j), ka(j+1) - 1
               ir      = ha(k)
               if (ir .gt. nncon) go to 340
               l       = l + 1
               a(k)    = gcon2(l)
  320       continue
  340    continue

      else if (mode .eq. 4) then
*        ---------------------------------------------------------------
*        Unbounded subproblem.
*        Increase the penalty parameter.
*        ---------------------------------------------------------------
         modpen = modpen + 1
         if (modpen .gt. 5  .or.  nlag .eq. 0) then
            inform = 1
         else
            t      = nncon
            if (penpar .le. zero) penpar = t / 100.0d+0
            penpar = 10.0d+0 * penpar
            if (iprint .gt. 0) write(iprint, 1400) penpar
            if (isumm  .gt. 0) write(isumm , 1400) penpar
         end if

      else if (mode .eq. 5) then
*        ---------------------------------------------------------------
*        Infeasible subproblem.
*        Relax the bounds on the linearized constraints.
*        Also, move nonbasic slacks onto the new bounds
*        so there won't be large numbers of them floating in between.
*
*        04 Mar 1994: Relax the bounds by both relative and absolute
*                     amounts.
*                     tolabs = sinf at first infeasible subproblem.
*                     tolrel = e-4, e-3, e-2, e-1, 1.
*                     Don't move slacks, but instead, terminate the
*                     subproblem at the first "feasible" point.
*        ---------------------------------------------------------------
         maxmod = 3
         if (modrhs .lt. maxmod) then
            modrhs = modrhs + 1
            if (modrhs .eq. 1) then
               rhsmod = sinf
            end if
            tolabs = rhsmod
            tolrel = 10.0d+0 ** (modrhs - maxmod)
            if (iprint .gt. 0) write(iprint, 1500) tolrel, tolabs
            if (isumm  .gt. 0) write(isumm , 1500) tolrel, tolabs

            call m5hs  ( 'External', nb, bl, bu, hs, xn )

            do 520 i = 1, nncon
               j     = n + i
               bnd   = blslk(i)
               if (bnd .gt. bminus) then
                   bl(j) = bnd - (abs( bnd ) + tolabs) * tolrel
*--                if (hs(j) .eq. 0) xn(j) = bl(j)
               end if
               bnd   = buslk(i)
               if (bnd .lt. bplus ) then
                   bu(j) = bnd + (abs( bnd ) + tolabs) * tolrel
*--                if (hs(j) .eq. 1) xn(j) = bu(j)
               end if
  520       continue

            call m5hs  ( 'Internal', nb, bl, bu, hs, xn )

            do 540 k  = 1, ms
               j      = kb(k)
               bbl(k) = bl(j)
               bbu(k) = bu(j)
  540       continue

         else
            inform = 1
         end if

      else if (mode .eq. 6) then
*        ---------------------------------------------------------------
*        Restore the bounds on the nonlinear constraints.
*        ---------------------------------------------------------------
         call dcopy ( nncon, blslk, 1, bl(n+1), 1 )
         call dcopy ( nncon, buslk, 1, bu(n+1), 1 )
      end if

      return

 1400 format(' Penalty parameter increased to', 1p, e12.2)
 1500 format(' XXX  Infeasible subproblem.  ',
     $   '  RHS perturbed by', 1p, e9.1, ' * (|RHS|  +', e9.1, ')')

*     end of m8aug1
      end

