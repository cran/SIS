************************************************************************
*
*     File  mi35mps  fortran.
*
*     m3inpt
*
*     12 Jul 2000: mi36mps.f now contains most MPS routines.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3inpt( objadd, z, nwcore )

      implicit           double precision (a-h,o-z)
      double precision   z(nwcore)

*     ------------------------------------------------------------------
*     m3inpt  inputs constraint data in MPS format, and sets up
*     various quantities as follows:
*
*     objadd     (output) is minus the coefficient in row iobj of the
*                RHS section (zero by default).  MINOS adds it to the
*                objective function.
*               
*     m, n, ne   are the number of rows, columns and elements in A.
*
*     iobj       is the row number for the linear objective (if any).
*                It must come after any nonlinear rows.
*                iobj = 0 if there is no linear objective.
*
*     a, ha, ka  is the matrix a stored in z at locations la, lha, lka.
*     bl, bu     are the bounds  stored in z at locations lbl, lbu.
*     hs, xn     are states and values   stored in z at   lhs, lxn.
*
*     hs(j)      is set to  0, 1  to indicate a plausible initial state
*                (at lo or up bnd) for each variable j  (j = 1 to nb).
*                If crash is to be used, i.e., crash option gt 0 and
*                if no basis file will be supplied, the initial bounds
*                set may initialize hs(j) as follows to assist crash:
*
*     -1      if column or row j is likely to be in the optimal basis,
*      4      if column j is likely to be nonbasic at its lower bound,
*      5      if column j is likely to be nonbasic at its upper bound,
*      2      if column or row j should initially be superbasic,
*      0 or 1 otherwise.
*
*     xn(j)      is a corresponding set of initial values.
*                Safeguards are applied later by m4chek, so the
*                values of hs and xn are not desperately critical.
*
*     The arrays name(*), mobj(*), mrhs(*), mrng(*), mbnd(*) are loaded
*     with the appropriate names in 2a4 format.
*
*     m3inpt (and hence m3hash, m3mpsa, m3mpsb, m3mpsc and m3read)
*     may be replaced by routines that output the same information.
*
*     31 Oct 1991: Modified to be compatible with subroutine minoss.
*     10 Apr 1992: objadd added as an output parameter.
*     04 May 1992: m3mpsb now outputs pi.
*     16 May 1997: nint and hint(*) are new output parameters from
*                  m3mpsb.  For stand-alone MINOS, m2core allocates
*                  hint in z(*), but so far nint and hint are not used.
*     21 Oct 1997: nint and hint(*) are new I/O parameters for m3mpsc.
*                  Save them in iparm(30) and (29).
*     05 Feb 1998: Set njac here after MPS input if Jacobian = Dense.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m1word/ nwordr,nwordi,nwordh
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m2len / mrows ,mcols ,melms
      common    /m2mapa/ ne    ,nka   ,la    ,lha   ,lka
      common    /m2mapz/ maxw  ,maxz
      common    /m2parm/ dparm(30),iparm(30)
      common    /m3len / m     ,n     ,nb    ,nscl
      common    /m3loc / lascal,lbl   ,lbu   ,lbbl  ,lbbu  ,
     $                   lhrtyp,lhs   ,lkb
      common    /m3mps1/ lname1,lname2,lkeynm,nname
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m5len / maxr  ,maxs  ,mbs   ,nn    ,nn0   ,nr    ,nx
      common    /m5loc / lpi   ,lpi2  ,lw    ,lw2   ,
     $                   lx    ,lx2   ,ly    ,ly2   ,
     $                   lgsub ,lgsub2,lgrd  ,lgrd2 ,
     $                   lr    ,lrg   ,lrg2  ,lxn
      common    /m5log1/ idebug,ierr,lprint
      common    /m7len / fobj  ,fobj2 ,nnobj ,nnobj0
      common    /m8len / njac  ,nncon ,nncon0,nnjac
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor

      integer            ncard(6)
      logical            allnn
      character*4        key
      character*4        lenda
      data               lenda /'ENDA'/

*     ------------------------------------------------------------------
*     Here's a kludge to let users say they want all columns to be
*     nonlinear.  They should specify
*          Columns              n
*          Nonlinear variables  n
*     for any large enough n (ge the true n).
*     here we just test if the two n's are the same.
*     Later we reset nn once we know the true n.
*     ------------------------------------------------------------------
      allnn  = mcols .eq. nn

*     We may come back here to try again with more workspace.
*     key    retains the first 4 characters of the NAME, ROWS, COLUMNS
*            RHS, RANGES and BOUNDS cards.
*     ncard  counts the number of data records in each section.
*     m3getp finds a prime number for the length of the row hash table.

   10 ierr   = 0
      ncoll  = 0
      key    = '    '
      call iload1( 6, 0, ncard, 1 )
      call m3getp( mrows, lenh )

      call m2core( 2, mincor )
      if (maxz .lt. mincor) go to 600

*     ------------------------------------------------------------------
*     Input ROWS.
*     lrow   is the location of the first rowname in name1, name2.
*     lennm  is the initial length of name1, name2,
*            i.e. the maximum no. of names allowed for.
*     ------------------------------------------------------------------
      lrow   = mcols + 1
      lennm  = mcols + mrows
      call m3mpsa( mrows, mcols, melms, ncoll, m,
     $             lrow, lennm, lenh, nn, nncon, key, ncard,
     $             z(lhrtyp), z(lname1), z(lname2), z(lkeynm) )
      if (ierr .eq. 40) go to 400
      if (ierr .eq. 41) go to 500

*     ------------------------------------------------------------------
*     m  is now known.
*     Input COLUMNS, RHS, RANGES.
*     ------------------------------------------------------------------
      mrows  = m
      lhint  = lname2 + 1 + (nname/nwordi)
      call m3mpsb( mcols, melms, lrow, lennm, lenh, ncoll, objadd,
     $             m, n, nb, ne, nka, nint,
     $             nn, nncon, nnjac, nnobj, njac, key, ncard,
     $             z(lhint), z(lhrtyp), z(lname1), z(lname2), z(lkeynm),
     $             z(lka), z(lha), z(la), z(lbl), z(lbu),
     $             z(lkb), z(lpi) )
      if (ierr .eq. 40) go to 400
      if (ierr .eq. 41) go to 510

*     ------------------------------------------------------------------
*     n  and  ne  are now known.
*     Move the row names to be contiguous with the column names.
*     Input BOUNDS.
*     ------------------------------------------------------------------
      call m3imov( lennm, lrow, m, n, z(lname1) )
      call m3imov( lennm, lrow, m, n, z(lname2) )
      mcols  = n
      melms  = ne
      np1    = n + 1
      nb     = n + m
      if (maxs .gt. np1) maxs = np1
      if (maxr .gt. np1) maxr = np1
      if (nn   .ge. nb ) nn   = nb
      if (    allnn    ) then
          nn     = n
          if (nnobj .gt. 0) nnobj = n
          if (nnjac .gt. 0) nnjac = n
      end if

      call m3mpsc( m, n, nb, ne, nint, ns, lennm,
     $             key, ncard, z(lname1), z(lname2),
     $             z(lbl), z(lbu), z(lhint), z(lhs), z(lxn) )

      if (iprint .gt. 0) then
                            write(iprint, '(/)')
         if (lprint .gt. 0) write(iprint, 1300) lenh, ncoll
         if (na0    .gt. 0) write(iprint, 1320) na0
         if (nncon  .gt. 0) write(iprint, 1350) njac
         if (nncon  .gt. 0  .or.  ncard(5) .gt. 0)
     $                      write(iprint, 1400) ncard(5)
         if (nn     .gt. 0  .or.  ncard(6) .gt. 0)
     $                      write(iprint, 1420) ncard(6), ns
      end if

*     njac  counted the actual Jacobian entries in the MPS file,
*     but for Jacobian = Dense, we have to reset it.

      if (nden .eq. 1) njac = nnCon * nnJac

*     ------------------------------------------------------------------
*     Compress storage, now that we know the size of everything.
*     ------------------------------------------------------------------

*     Save current positions of  bl, bu, etc.

      kha    = lha
      kka    = lka
      kbl    = lbl
      kbu    = lbu
      kn1    = lname1
      kn2    = lname2
      khint  = lhint
      khs    = lhs
      kxn    = lxn
      kpi    = lpi

*     Redefine addresses in  z  in terms of the known dimensions.

      call m2core( 3, mincor )
      lhint  = lname2 + 1 + (nname/nwordi)
      if (maxz .lt. mincor) go to 800

*     Move bl, bu, etc. into their final positions.

      call hcopy ( ne, z(kha),   1, z(lha),   1 )
      call icopy ( nka,z(kka),   1, z(lka),   1 )
      call dcopy ( nb, z(kbl),   1, z(lbl),   1 )
      call dcopy ( nb, z(kbu),   1, z(lbu),   1 )
      if (nname .gt. 1) then
         call icopy ( nb, z(kn1), 1, z(lname1), 1 )
         call icopy ( nb, z(kn2), 1, z(lname2), 1 )
      end if
      call hcopy ( n , z(khint), 1, z(lhint), 1 )
      call hcopy ( nb, z(khs),   1, z(lhs),   1 )
      call dcopy ( nb, z(kxn),   1, z(lxn),   1 )
      call dcopy ( m , z(kpi),   1, z(lpi),   1 )

*     Save nint and khint in case anyone wants them later.

      iparm(30) = nint
      iparm(29) = lhint
      go to 900

*     ---------------------------
*     Fatal error in MPS file.
*     ---------------------------
  400 call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1100)
      if (isumm  .gt. 0) write(isumm , 1100)
      go to 700

*     ---------------------------
*     Too many rows.
*     ---------------------------
  500 mrows  = m
      go to 520

*     -----------------------------
*     Too many columns or elements.
*     -----------------------------
  510 mcols  = n
      melms  = ne

*     Try again.

  520 if (imps .ne. iread  .and.  imps .ne. ispecs) then
         rewind imps
         go to 10
      end if

*     ---------------------------------
*     Not enough core to read MPS file.
*     ---------------------------------
  600 call m1page(2)
      if (iprint .gt. 0) write(iprint, 1110)
      if (isumm  .gt. 0) write(isumm , 1110)
      ierr   = 41

*     ------------------------------------
*     Flush MPS file to the ENDATA card
*     if it is the same as the SPECS file.
*     ------------------------------------
  700 if (imps .eq. ispecs) then
         do 750 idummy = 1, 100000
            if (key .eq. lenda) go to 900
            read(imps, 1700, end=900) key
  750    continue
      end if
      go to 900

*     -------------------------------------
*     Not enough core to solve the problem.
*     -------------------------------------
  800 call m1page(2)
      if (iprint .gt. 0) write(iprint, 1120) mincor
      if (isumm  .gt. 0) write(isumm , 1120) mincor
      ierr   = 42

*     Exit.

  900 if (imps .ne. iread  .and.  imps .ne. ispecs) rewind imps
      return

 1100 format(' EXIT -- fatal errors in the MPS file')
 1110 format(' EXIT -- not enough storage to read the MPS file')
 1120 format(' EXIT -- not enough storage to start solving',
     $       ' the problem'
     $    // ' Workspace (total)  should be significantly',
     $       ' more than', i8)
 1300 format(/ ' Length of row-name hash table  ', i12
     $       / ' Collisions during table lookup ', i12)
 1320 format(  ' No. of rejected coefficients   ', i12)
 1350 format(  ' No. of Jacobian entries specified', i10)
 1400 format(  ' No. of LAGRANGE entries specified', i10)
 1420 format(  ' No. of INITIAL  bounds  specified', i10
     $       / ' No. of superbasics specified   ', i12)
 1700 format(a4)

*     end of m3inpt
      end
