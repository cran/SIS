************************************************************************
*
*     File  mi21amat fortran.
*
*     m2apr1   m2apr5
*
*     12 Jul 2000: mi21amat.f separates m2apr1, m2apr5 from mi20amat.f.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2apr1( mode, m, mbs, n, tolz,
     $                   ne, nka, a, ha, ka, kb, x, lenx, y, leny )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka), kb(mbs)
      double precision   a(ne)
      double precision   x(lenx), y(leny)

*     ------------------------------------------------------------------
*     m2apr1  computes various matrix-vector products involving
*     B  and  S,  the basic and superbasic columns of  A.
*
*     mode=1    y  =  y - B*x   or   y - (B S)*x  (depending on  lenx)
*     mode=2    y  =  y - S*x
*     mode=3    y  =  y - B(transpose)*x   or   y - (B S)(transpose)*x
*    (mode 3  is not used as yet)
*     mode=4    y  =  y - S(transpose)*x
*     ------------------------------------------------------------------

      kbase  = 0
      if (mode .le. 2) then
*        -------------
*        mode 1 and 2.
*        -------------
         if (mode .eq. 2) kbase = m

         do 200 k = 1, lenx
            xj    = x(k)
            if (abs(xj) .le. tolz) go to 200
            j     = kb(kbase + k)
            if (j .le. n) then
               do 150 i = ka(j), ka(j+1) - 1
                  ir    = ha(i)
                  y(ir) = y(ir) - a(i)*xj
  150          continue
            else
               ir    = j - n
               y(ir) = y(ir) - xj
            end if
  200    continue
      else
*        -------------
*        mode 3 and 4.
*        -------------
         if (mode .eq. 4) kbase = m

         do 400 k = 1, leny
            t     = y(k)
            j     = kb(kbase + k)
            if (j .le. n) then
               do 350 i = ka(j), ka(j+1) - 1
                  ir    = ha(i)
                  t     = t - x(ir)*a(i)
  350          continue
               y(k) = t
            else
               ir   = j - n
               y(k) = t - x(ir)
            end if
  400    continue
      end if

*     end of m2apr1
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m2apr5( mode, m, n, nb, nncon, nnjac, tolz,
     $                   ne, nka, a, ha, ka, xn, lenx, y, leny )

      implicit           double precision (a-h,o-z)
      integer            ha(ne)
      integer            ka(nka)
      double precision   a(ne)
      double precision   xn(lenx), y(leny)

*     ------------------------------------------------------------------
*     m2apr5  computes products involving parts of  A  and  xn.
*     ------------------------------------------------------------------

      if (mode .eq. 5) then
*        ----------------------------------
*        mode = 5.   Set  y = y - A*xn.
*        ----------------------------------
         do 500 j = 1, n
            xj    = xn(j)
            if (abs( xj ) .gt. tolz) then
               do 450 k = ka(j), ka(j+1) - 1
                  ir    = ha(k)
                  y(ir) = y(ir) - a(k)*xj
  450          continue
            end if
  500    continue

      else
*        ------------------------------------------------------------
*        If mode = 7, set  y  =  y - (Jacobian)*xn.
*        If mode = 8, set  y  =  y - (linear A)*xn.
*        Only the first nncon rows are involved, and no slacks.
*        ------------------------------------------------------------
         if (mode .eq. 7) then
            j1 = 1
            j2 = nnjac
         else
            j1 = nnjac + 1
            j2 = n
         end if

         do 800 j = j1, j2
            xj    = xn(j)
            if (abs( xj ) .gt. tolz) then
               do 750 k = ka(j), ka(j+1) - 1
                  ir    = ha(k)
                  if (ir .le. nncon) y(ir) = y(ir) - a(k)*xj
  750          continue
            end if
  800    continue
      end if

*     end of m2apr5
      end
