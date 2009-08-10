************************************************************************
*
*     File  mi36mps  fortran.
*
*     m3getp   m3hash   m3imov   m3pint
*     m3mpsa   m3mpsb   m3mpsc   m3read
*
*     12 Jul 2000: mi36mps.f now separates m3inpt from other MPS routines.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3getp( maxm, lenh )

*     ------------------------------------------------------------------
*     m3getp finds a prime number lenh suitably larger than maxm.
*     It is used as the length of the hash table for the MPS row names.
*
*     20 Apr 1992: Implemented for use in m3inpt and matmps.
*     ------------------------------------------------------------------

      lenh   = maxm*2
      lenh   = max( lenh, 100 )
      lenh   = (lenh/2)*2 - 1
      k      = lenh/20 + 6

  100 k      = k + 1
      lenh   = lenh + 2
      do 120 i = 3, k, 2
         if (mod(lenh,i) .eq. 0) go to 100
  120 continue

*     end of m3getp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3hash( len  , nen  , ncoll,
     $                   key1 , key2 , mode , keytab,
     $                   name1, name2, ka   , found )

      integer            keytab(len), name1(nen), name2(nen)
      logical            found

*     ------------------------------------------------------------------
*     m3hash looks up and/or inserts integer keys in a table.
*     Reference:  R.P. Brent, CACM 16,2 (Feb 1973), pp. 105-109.
*     This version is simplified for the case where no entries are
*     deleted.
*     keytab is used as an index into a consecutive list of unique
*     identifiers name1 and name2.
*     Each pair   name1(i), name2(i)  is two 4-character identifiers,
*     treated as one 8-character identifier.
*     ------------------------------------------------------------------

      len2 = len - 2
      ic   = -1

*     Compute address of first probe (ir) and increment (iq).
*     ******************************************************************
*     NOTE -- the next statements are machine dependent.  The aim
*     is to produce a valid positive integer  key  out of the
*     two words  key1  and  key2 .  The latter contain four
*     characters left-justified (read under a4 format).  These
*     may turn on the sign bit (e.g. IBM 360/370) or give a
*     non-zero floating point exponent (e.g. Burroughs B6700).
*     ******************************************************************
*---+ DEC VAX
      k1  = abs(key1)
      k2  = abs(key2)

*---- On CDC systems (and possibly CRAY systems), change the previous 2
*-    lines to comments, and uncomment the next 3 lines.
*-    decode( 4, 10, key1 ) k1
*-    decode( 4, 10, key2 ) k2
*- 10 format( r4 )

      key = abs(k1 - k2)
      iq  = mod(key, len2) + 1
      ir  = mod(key, len)  + 1
      ka  = ir

*     Look in the table.
   20 kt  = keytab(ka)

*     Check for an empty space or a match.
      if (kt .eq. 0) go to 30
      if (key1 .eq. name1(kt)  .and.  key2 .eq. name2(kt)) go to 60
      ic    = ic + 1
      ncoll = ncoll + 1

*     Compute address of next probe.
      ka    = ka + iq
      if (ka .gt. len) ka = ka - len

*     See if whole table has been searched.
      if (ka .ne. ir ) go to 20

*     The key is not in the table.
   30 found = .false.

*     Return with KA = 0 unless an entry has to be made.
      if ((mode .eq. 2)  .and.  (ic .le. len2)) go to 70
      ka    = 0
      return

   60 found = .true.
      return

*     Look for the best way to make an entry.
   70 if (ic .le. 0) return
      ia  = ka
      is  = 0

*     Compute the maximum length to search along current chain.
   80 ix  = ic - is
      kt  = keytab(ir)

*     Compute increment JQ for current chain.
*     ******************************************************************
*     NOTE -- the next statements are machine dependent.  The same
*     transformation as discussed in the note above should be
*     applied to  name1(kt)  and  name2(kt)  to produce an integer  key.
*     ******************************************************************
*---+ DEC VAX
      k1  = abs(name1(kt))
      k2  = abs(name2(kt))

*---- On CDC systems (and possibly CRAY systems), change the previous 2
*-    lines to comments, and uncomment the next 2 lines.
*-    decode( 4, 10, name1(kt) ) k1
*-    decode( 4, 10, name2(kt) ) k2

      key = abs(k1 - k2)
      jq  = mod(key, len2) + 1
      jr  = ir

*     Look along the chain.
   90 jr  = jr + jq
      if (jr .gt. len) jr = jr - len

*     Check for a hole.
      if (keytab(jr) .eq. 0) go to 100
      ix  = ix - 1
      if (ix .gt. 0) go to 90
      go to 110

*     Save location of hole.
  100 ia  = jr
      ka  = ir
      ic  = ic - ix

*     Move down to the next chain.
  110 is  = is + 1
      ir  = ir + iq
      if (ir .gt. len) ir = ir - len

*     Go back if a better hole might still be found.
      if (ic .gt. is ) go to 80

*     If necessary move an old entry.
      if (ia .ne. ka ) keytab(ia) = keytab(ka)

*     end of m3hash
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3imov( lennm, lrow, m, n, name )

      integer            name(lennm)

*     ------------------------------------------------------------------
*     m3imov is needed by m3inpt and matmps to move the row names
*     to be contiguous with the column names inside the integer array
*     name(*), once the true number of columns is known.
*     On entry, the column names begin at name(1) and the row names
*               begin at name(lrow).
*     On exit,  the row names begin at name(n+1).
*
*     20 Apr 1992: m3imov implemented to help matmps.
*     ------------------------------------------------------------------

      if (lrow .gt. n+1) then
         call icopy ( m, name(lrow), 1, name(n+1), 1 )
      end if

*     end of m3imov
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3pint( n, nint, hint )

      integer            hint(n)

*     ------------------------------------------------------------------
*     m3pint  packs the "integer variable" array hint(*).
*
*     On entry,
*     nint    is the number of integer variables found by m3mpsb, m3mpsc.
*             It should be the number of positive entries in hint(*).
*     hint(j) = 0 if x(j) is continuous,
*             = 1 if x(j) is integer.
*
*     On exit,
*     hint(*) contains nint entries.  The k-th integer variable is
*             column number hint(k), k = 1, nint.
*
*     21 Oct 1997: First version of m3pint.
*     ------------------------------------------------------------------

      if (nint .gt. 0) then
         k      = 0
         do 100 j = 1, n
            if (hint(j) .gt. 0) then
               k       = k + 1
               hint(k) = j
            end if
  100    continue
      end if

*     end of m3pint
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3mpsa( mrows, mcols, melms, ncoll, m,
     $                   lrow, lennm, lenh, nn, nncon, key, ncard,
     $                   hrtype, name1, name2, keynam )

      implicit           double precision (a-h,o-z)
      character*4        key
      integer            hrtype(mrows)
      integer            ncard(6), name1(lennm), name2(lennm),
     $                   keynam(lenh)

*     ------------------------------------------------------------------
*     m3mpsa  inputs the name and rows sections of an MPS file.
*
*     Original version written by Keith Morris, Wellington, 1973.
*     1975: Use a hash table for the row names.
*     1979: Treat the (nncon by nnjac) principal submatrix
*           as a Jacobian for nonlinear constraints.
*     1980: Add Phantom columns to the end of  A.
*     1982: Store the rhs as bounds on the logicals,
*           instead of the last column of  A.  The constraints now have
*           the form    A*x  +  I*s  =  0,    bl .le. (x, s) .le. bu,
*           where A has m rows, n columns and ne nonzero elements.
*     1982  Treat * in column 1 correctly, and retry if the MPS file
*           is on disk and mrows, mcols or melms are too small.
*     Apr 1984: Added check for duplicate row entries in columns.
*     Mar 1985: Changes made to handle characters as in Fortran 77.
*     Oct 1985: m3mps split into m3mpsa, m3mpsb, m3mpsc.
*               Revisions made to handle characters more efficiently.
*     Oct 1991: More f77.  Row names now at end of name1, name2.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3mps5/ aelem(2), id(6), iblank
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint

      logical            found,gotnm
      character*4        lblank,
     $                   lname,lrows,lcolu,
     $                   lex,lgx,llx,lnx,lxe,lxg,lxl,lxn
      data               lblank            /'    '/
      data               lname,lrows,lcolu /'NAME','ROWS','COLU'       /
      data               lex,lgx,llx,lnx   /' E  ',' G  ',' L  ',' N  '/
      data               lxe,lxg,lxl,lxn   /'  E ','  G ','  L ','  N '/
*     ------------------------------------------------------------------

      call m1page( 1 )
      if (iprint .gt. 0) write(iprint, 1000)
      read (lblank, '(a4)') iblank
      inform = 0
      iobj   = 0
      line   = 0
      m      = 0
      gotnm  = mobj(1) .ne. iblank
      call iload1( 20  , 0, ier   , 1 )
      call iload1( lenh, 0, keynam, 1 )

*     Look for the NAME card.

   10 call m3read( 1, imps, line, 5, key, inform )
      if (key .ne. lname) then
         if (ier(1) .eq. 0) then
             ier(1) = 1
             if (iprint .gt. 0) write(iprint, 1100)
             if (isumm  .gt. 0) write(isumm , 1100)
         end if
         go to 10
      end if

      name(1) = id(3)
      name(2) = id(4)
      if (isumm  .gt. 0) write(isumm, 5000) name

*     Look for the ROWS card.

      call m3read( 1, imps, line, 5, key, inform )
      inform = 0
      if (key .ne. lrows) then
         ier(1) = ier(1) + 1
         if (iprint .gt. 0) write(iprint, 1120)
         if (isumm  .gt. 0) write(isumm , 1120)
         go to 35
      end if

*     ==================================================================
*     Read the row names and check if the relationals are valid.
*     ==================================================================
   30 call m3read( 1, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 110

   35 if      (key .eq. lgx  .or.  key .eq. lxg) then
         it  = -1
      else if (key .eq. lex  .or.  key .eq. lxe) then
         it  =  0
      else if (key .eq. llx  .or.  key .eq. lxl) then
         it  =  1
      else if (key .eq. lnx  .or.  key .eq. lxn) then
         it  =  2

*        Record objective name if we don't already have one.

         if (iobj .eq. 0) then
            if (.not. gotnm) then
               mobj(1) = id(1)
               mobj(2) = id(2)
               if (nn .gt. 0) then
                  if (iprint .gt. 0) write(iprint, 1170) mobj
                  if (isumm  .gt. 0) write(isumm , 1170) mobj
               end if
            end if

            if (id(1) .eq. mobj(1) .and. id(2) .eq. mobj(2)) then
               iobj     = m + 1
               ncard(1) = ncard(1) + 1
            end if
         end if
      else
         ier(3) = ier(3) + 1
         if (ier(3) .le. mer) then
            if (iprint .gt. 0) write(iprint, 1160) line,key,id(1),id(2)
            if (isumm  .gt. 0) write(isumm , 1160) line,key,id(1),id(2)
         end if
         go to 30
      end if

*     ..................................................................
*     Look up the row name  id(1), id(2)  in the hash table.
*     ..................................................................
      call m3hash( lenh, mrows, ncoll, id(1), id(2), 2,
     $             keynam, name1(lrow), name2(lrow), ia, found )

*     Error if the row name was already there.
*     Otherwise, enter the new name into the hash table.

      if (found) then
         ier(4) = ier(4) + 1
         if (ier(4) .le. mer) then
            if (iprint .gt. 0) write(iprint, 1200) id(1), id(2)
            if (isumm  .gt. 0) write(isumm , 1200) id(1), id(2)
         end if
      else
         m      = m + 1
         if (m .le. mrows) then
            jrow        = mcols + m
            keynam(ia)  = m
            name1(jrow) = id(1)
            name2(jrow) = id(2)
            hrtype(m)   = it
         end if
      end if
      go to 30

*     ==================================================================
*     Should be COLUMNS card.
*     ==================================================================
  110 if (key .ne. lcolu) then
         ier(1) = ier(1) + 1
         if (iprint .gt. 0) write(iprint, 1130)
         if (isumm  .gt. 0) write(isumm , 1130)
      end if

*     Error if no rows or too many rows.

      if (m .le. 0) then
         if (iprint .gt. 0) write(iprint, 1300)
         if (isumm  .gt. 0) write(isumm , 1300)
         ier(1) = ier(1) + 1
         ierr   = 40
         return
      else if (m .gt. mrows) then
         if (iprint .gt. 0) write(iprint, 3030) mrows, m
         if (isumm  .gt. 0) write(isumm , 3030) mrows, m
         ier(1) = ier(1) + 1
         ierr   = 41
         return
      end if

*     Warning if no objective row found.
*     Error if linear objective is ahead of nonlinear rows.

      if (iobj .eq. 0) then
         if (iprint .gt. 0) write(iprint, 1600)
         if (isumm  .gt. 0) write(isumm , 1600)
      else if (iobj .le. nncon) then
         if (iprint .gt. 0) write(iprint, 1180) mobj
         if (isumm  .gt. 0) write(isumm , 1180) mobj
         ierr   = 40
         return
      end if

      if (isumm  .gt. 0) write(isumm, 5100) m
      return

 1000 format(' MPS file' / ' --------')
 1100 format(' XXXX  Garbage before NAME card')
 1120 format(' XXXX  ROWS card not found')
 1130 format(' XXXX  COLUMNS card not found')
 1160 format(' XXXX  Illegal row type at line', i7, '... ', 3a4)
 1170 format(' ===>  Note:  row  ', 2a4,
     $   '  selected as linear part of objective.')
 1180 format(/ ' XXXX  The linear objective card      N ', 2a4
     $       / ' XXXX  is out of place.    Nonlinear constraints'
     $       / ' XXXX  must be listed first in the ROWS section.')
 1200 format(' XXXX  Duplicate row name --', 2a4, ' -- ignored')
 1300 format(' XXXX  No rows specified')
 1600 format(' ===>  Warning - no linear objective selected')
 2000 format(2a4)
 3030 format(' XXXX  Too many rows.  Limit was', i8,
     $   4x, '  Actual number is', i8)
 5000 format(' Name   ', 2a4)
 5100 format(' Rows   ',  i8)

*     end of m3mpsa
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3mpsb( mcols, melms, lrow, lennm, lenh, ncoll, objadd,
     $                   m, n, nb, ne, nka, nint,
     $                   nn, nncon, nnjac, nnobj, njac, key, ncard,
     $                   hint, hrtype, name1, name2, keynam,
     $                   ka, ha, a, bl, bu, kb, pi )

      implicit           double precision (a-h,o-z)
      character*4        key
      integer            hint(mcols), hrtype(m), ha(melms)
      integer            ncard(6) , name1(lennm), name2(lennm),
     $                   keynam(lenh), ka(nka)  , kb(m)
      double precision   a(melms) , bl(nb), bu(nb)
      double precision   pi(m)

*     ------------------------------------------------------------------
*     m3mpsb inputs the COLUMNS, RHS and RANGES sections of an MPS file.
*     10 Apr 1992: objadd added as an output parameter.  It is minus
*                  the coefficient in row iobj of the RHS.
*     04 May 1992: pi(m)  added as an output parameter to return the
*                  special RHS called LAGRANGE.
*     22 Apr 1993: Bug in loop 403 fixed.  The write and read are new.
*                  Needs character*8 phantm.
*     16 May 1997: nint, hint(*) now output if 'INTORG' and 'INTEND'
*                  markers are found in the COLUMNS section.
*     21 Oct 1997: nint, hint(*) now flag integer variables directly.
*                  This allows BV, LI, UI to be dealt with in BOUNDS.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3mps5/ aelem(2), id(6), iblank
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /cyclcm/ cnvtol,jnew,materr,maxcy,nephnt,nphant,nprint

      parameter        ( zero = 0.0d+0 )

      logical            dense, found, gotnm, intorg
      character*8        phantm
      character*4        lrhs , lrhsx, lrang,
     $                   llagr, lange,
     $                   lmar , lker ,
     $                   lint , lorg , lend
      data               lrhs , lrhsx, lrang / 'RHS ', 'RHS''', 'RANG'/
      data               llagr, lange        / 'LAGR', 'ANGE' /
      data               lmar , lker         /'''MAR', 'KER'''/
      data               lint , lorg , lend  /'''INT', 'ORG''', 'END'''/
*     ------------------------------------------------------------------

      read (llagr , '(a4)') ilagr
      read (lange , '(a4)') iange
      read (lmar  , '(a4)') imar
      read (lker  , '(a4)') iker
      read (lint  , '(a4)') iint
      read (lorg  , '(a4)') iorg
      read (lend  , '(a4)') iend
      objadd = zero
      dense  = nden .eq. 1
      bplus  = plinfy
      bminus = - bplus
      nmcol1 = 1234
      nmcol2 = 5678
      n      = 0
      na0    = 0
      ne     = 0
      ne1    = -1
      njac   = 0
      nint   = 0
      inform = 0
      intorg = .false.
      call iload1( m,    0, kb, 1 )
      call dload ( m, zero, pi, 1 )

*     ==================================================================
*     Read the next columns card.
*     ==================================================================
  210 call m3read( 2, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 310

  220 if (id(1) .ne. nmcol1  .or.  id(2) .ne. nmcol2) then

*        ---------------------------------------------------------------
*        Test for 'MARKER' in Field 3.
*        ---------------------------------------------------------------
         if (id(3) .eq. imar  .and.  id(4) .eq. iker) then
            if (id(5) .eq. iint  .and.  id(6) .eq. iorg) then

*              'INTORG' found.

*              if (intorg) then Error!
               intorg = .true.

            else if (id(5) .eq. iint  .and.  id(6) .eq. iend) then

*              'INTEND' found.

*              if (.not. intorg) then Error!
               intorg = .false.

            else
*              Error!
            end if
            go to 210
         end if

*        ---------------------------------------------------------------
*        Start a new column.
*        ---------------------------------------------------------------
         if (ne .le. ne1) go to 310
         n         = n + 1
         ne1       = ne
         nmcol1    = id(1)
         nmcol2    = id(2)
         if (n .le. mcols) then
            ka(n)    = ne + 1
            name1(n) = nmcol1
            name2(n) = nmcol2

            if (intorg) then
               nint       = nint + 1
*----          hint(nint) = n
               hint(n)    = 1
            else
               hint(n)    = 0
            end if

*           Make room for a dense Jacobian column.

            if (nncon .gt. 0) then
               ljac   = ne
               if (dense  .and.  n .le. nnjac) then
                  ne     = ne + nncon
                  if (ne .le. melms) then
                     ne     = ne - nncon
                     do 225 i  = 1, nncon
                        ne     = ne + 1
                        ha(ne) = i
                        a(ne)  = zero
  225                continue
                  end if
               end if
            end if
         end if
      end if

*     Process two row names and values.

      do 260 i  = 1, 2

*        Check for only one on the card.

         k      = i + i
         id1    = id(k+1)
         id2    = id(k+2)
         if (id1 .ne. iblank) go to 230
         if (id2 .eq. iblank) go to 260

*        Look up the row name.

  230    call m3hash( lenh, m, ncoll, id1, id2, 1,
     $                keynam, name1(lrow), name2(lrow), ia, found )

         if ( found ) then
            aij    = aelem(i)
            irow   = keynam(ia)

*           Test for a duplicate entry.

            if (kb(irow) .eq. n) then
               ier(8) = ier(8) + 1
               if (iprint .gt. 0  .and.  ier(8) .le. mer)
     $         write(iprint, 1420) nmcol1, nmcol2, id1, id2, aij, line
               go to 260
            end if

            kb(irow) = n
            if (irow .le. nncon  .and.  n .le. nnjac) then

*              Deal with Jacobian elements.

               njac   = njac + 1
               if ( dense ) then
                  a(ne1 + irow) = aij
                  go to 260
               end if

*              Sparse Jacobian -- make sure the new element is
*              squeezed in ahead of any linear-constraint elements.

               ljac   = ljac + 1
               if (ljac .le. ne) then
                  aij      = a(ljac)
                  irow     = ha(ljac)
                  a(ljac)  = aelem(i)
                  ha(ljac) = keynam(ia)
               end if

            else if (abs( aij ) .lt. aijtol) then

*              Ignore small aijs.

               na0    = na0 + 1
               go to 260
            end if

*           Pack the nonzero.

            ne     = ne + 1
            if (ne .le. melms) then
               ha(ne) = irow
               a(ne)  = aij
            end if
         else
            ier(5) = ier(5) + 1
            if (iprint .gt. 0  .and.  ier(5) .le. mer)
     $         write(iprint, 1400) id1, id2, line
         end if
  260 continue
      go to 210

*     Test for an empty column.

  310 if (ne .le. ne1) then

*        Column with no rows.   Warning unless variable is nonlinear.
*        Insert dummy column with zero in first row.

         if (n .gt. nn) then
            ier(6) = ier(6) + 1
            if (iprint .gt. 0  .and.  ier(6) .le. mer)
     $         write(iprint, 1500) nmcol1, nmcol2, line
         end if

         ne     = ne + 1
         if (ne .le. melms) then
            ha(ne) = 1
            a(ne)  = zero
         end if
         if (inform .eq. 0) go to 220
      end if

*     ==================================================================
*     See if we have hit the RHS.
*     ==================================================================
      if (key .ne. lrhs  .and.  key .ne. lrhsx) then

*        Nope sumpins rong.
*        Terminate the COLUMNS section anyway.

         ier(7) = ier(7) + 1
         if (iprint .gt. 0) write(iprint, 1140)
         if (isumm  .gt. 0) write(isumm , 1140)
      end if

*     Are there any columns at all?
*     Or too many columns or elements?
*     Include phantom columns and elements too.

      nephnt = max( nephnt, nphant )
      if (n .le. 0) then
         if (iprint .gt. 0) write(iprint, 1610)
         if (isumm  .gt. 0) write(isumm , 1610)
         ier(2) = ier(2) + 1
         ierr   = 40
         return
      else if (n + nphant .gt. mcols) then
         n      = n + nphant
         if (iprint .gt. 0) write(iprint, 3040) mcols, n
         if (isumm  .gt. 0) write(isumm , 3040) mcols, n
         ier(2) = ier(2) + 1
         ierr   = 41
         return
      else if (ne + nephnt .gt. melms) then
         ne     = ne + nephnt
         if (iprint .gt. 0) write(iprint, 3050) melms, ne
         if (isumm  .gt. 0) write(isumm , 3050) melms, ne
         ier(2) = ier(2) + 1
         ierr   = 41
         return
      end if

*     ------------------------------------------------------------------
*     Input the RHS.
*     ------------------------------------------------------------------
      if (isumm  .gt. 0) write(isumm, 5200) n, ne

*     Insert phantom columns for cycling algorithm.
*     Give them names  PHNT   1, PHNT   2, ...
*     (Trouble if more than 999999 phantom columns.)

      if (nphant .gt. 0) then
         do 402 i = 1, nephnt
            j     = ne + i
            ha(j) = 1
            a(j)  = zero
  402    continue

         do 403 k = 1, nphant
            n     = n    + 1
            ne    = ne   + 1
            write(phantm, '(a2, i6)') 'PH', k
            read (phantm, '( 2a4  )') name1(n), name2(n)
            ka(n) = ne
            if (k .eq. 1) ne = ne + nephnt - nphant
  403    continue
      end if

*     We finally know how big the problem is.

      ka(n+1) = ne + 1

*     Set bounds to default values.

      call dload ( n, bstruc(1), bl, 1 )
      call dload ( n, bstruc(2), bu, 1 )

      do 408 i  = 1, m
         k      = hrtype(i)
         jslack = n + i
         if (k .lt. 0) bl(jslack) = bminus
         if (k .le. 0) bu(jslack) = zero
         if (k .ge. 0) bl(jslack) = zero
         if (k .gt. 0) bu(jslack) = bplus
         if (k .eq. 2) bl(jslack) = bminus
  408 continue

*     Check for no RHS.

      if (key .ne. lrhs  .and.  key .ne. lrhsx) go to 600
      gotnm  = mrhs(1) .ne. iblank
      inform = 0

*     ==================================================================
*     Read next RHS card and see if it is the one we want.
*     ==================================================================
  410 call m3read( 2, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 600

*     A normal RHS is terminated if LAGRANGE is found.

      if (id(1) .eq. ilagr  .and.  id(2) .eq. iange) go to 490

      if (.not. gotnm) then
         gotnm   = .true.
         mrhs(1) = id(1)
         mrhs(2) = id(2)
      end if                                      

      if (id(1) .eq. mrhs(1)  .and.  id(2) .eq. mrhs(2)) then

*        Look at both halves of the record.

         do 440 i  = 1, 2
            k      = i + i
            id1    = id(k+1)
            id2    = id(k+2)
            if (id1 .eq. iblank  .and.  id2 .eq. iblank) go to 440
            call m3hash( lenh, m, ncoll, id1, id2, 1,
     $                   keynam, name1(lrow), name2(lrow), ia,found )

            if ( found ) then
               ncard(2) = ncard(2) + 1
               bnd      = aelem(i)
               irow     = keynam(ia)
               jslack   = n + irow
               k        = hrtype(irow)
               if (irow .eq. iobj) then
                  objadd = - bnd
               else if (k .ne. 2) then
                  if (k .le. 0) bu(jslack) = - bnd
                  if (k .ge. 0) bl(jslack) = - bnd
               end if
            else
               ier(5) = ier(5) + 1
               if (iprint .gt. 0  .and.  ier(5) .le. mer)
     $            write(iprint, 1400) id1, id2, line
            end if
  440    continue
      end if
      go to 410

*     LAGRANGE RHS found.

  490 if (ncard(2) .eq. 0) then
         mrhs(1) = iblank
         mrhs(2) = iblank
         if (iprint .gt. 0) write(iprint, 1720)
      end if
      go to 520

*     ==================================================================
*     Read next RHS card and see if it is a LAGRANGE one.
*     ==================================================================
  510 call m3read( 2, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 600

      if (id(1) .ne. ilagr  .or.   id(2) .ne. iange) go to 510

*     Find which row.
*     Look at both halves of the record.

  520 do 540 i  = 1, 2
         k      = i + i
         id1    = id(k+1)
         id2    = id(k+2)
         if (id1 .eq. iblank  .and.  id2 .eq. iblank) go to 540
         call m3hash( lenh, m, ncoll, id1, id2, 1,
     $                keynam, name1(lrow), name2(lrow), ia,found )
         
         if ( found ) then
            ncard(5) = ncard(5) + 1
            irow     = keynam(ia)
            pi(irow) = aelem(i)
         else
            ier(5) = ier(5) + 1
            if (iprint .gt. 0  .and.  ier(5) .le. mer)
     $         write(iprint, 1400) id1, id2, line
         end if
  540 continue
      go to 510

*     ------------------------------------------------------------------
*     RHS has been input.
*     ------------------------------------------------------------------
  600 if (ncard(2) .eq. 0) then
         if (iprint .gt. 0) write(iprint, 1620)
         if (isumm  .gt. 0) write(isumm , 1620)
      end if

      if (objadd .ne. zero) then
         if (iprint .gt. 0) write(iprint, 1630) objadd
         if (isumm  .gt. 0) write(isumm , 1630) objadd
      end if

*     ------------------------------------------------------------------
*     Input RANGES.
*     ------------------------------------------------------------------

*     Check for no RANGES.

      if (key .ne. lrang) go to 800
      gotnm  = mrng(1) .ne. iblank
      inform = 0
                                  
*     ==================================================================
*     Read card and see if it is the range we want.
*     ==================================================================
  610 call m3read( 2, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 800

      if (.not. gotnm) then
         gotnm   = .true.
         mrng(1) = id(1)
         mrng(2) = id(2)
      end if

      if (id(1) .eq. mrng(1)  .and.  id(2) .eq. mrng(2)) then

*        Look at both halves of the record.

         do 640 i  = 1, 2
            k      = i + i
            id1    = id(k+1)
            id2    = id(k+2)
            if (id1 .eq. iblank  .and.  id2 .eq. iblank) go to 640
            call m3hash( lenh, m, ncoll, id1, id2, 1,
     $                   keynam, name1(lrow), name2(lrow), ia,found )

            if ( found ) then
            ncard(3) = ncard(3)+1
            brng     = aelem(i)
            arng     = abs( brng )
            irow     = keynam(ia)
            jslack   = n + irow
            k        = hrtype(irow)
            if (k .ne. 2) then
               if (k .lt. 0) bl(jslack) = bu(jslack) - arng
               if (k .gt. 0) bu(jslack) = bl(jslack) + arng
               if (k .eq. 0) then
                  if (brng .gt. zero) bl(jslack) = bu(jslack) - arng
                  if (brng .lt. zero) bu(jslack) = bl(jslack) + arng
               end if
            end if
         else
            ier(5) = ier(5) + 1
            if (iprint .gt. 0  .and.  ier(5) .le. mer)
     $         write(iprint, 1400) id1, id2, line
            end if
  640    continue
      end if
      go to 610

*     RANGES have been input.

  800 return

 1140 format(' XXXX  RHS card not found')
 1400 format(' XXXX  Non-existent row    specified -- ', 2a4,
     $   ' -- entry ignored in line', i7)
 1420 format(' XXXX  Column  ', 2a4, '  has more than one entry',
     $       ' in row  ', 2a4
     $     / ' XXXX  Coefficient', 1p, e15.5, '  ignored in line', i10)
 1500 format(' XXXX  No valid row entries in column  ', 2a4,
     $       '   after line', i10)
 1610 format(' XXXX  No columns specified')
 1620 format(' ===>  Warning - the RHS is zero')
 1630 format(' ===>  Note:  constant', 1p, e15.7,
     $   '  is added to the objective.')
 1720 format(' ===>  Warning - first RHS is LAGRANGE.',
     $   '   Other RHS''s will be ignored.')
 2000 format(2a4)
 3040 format(' XXXX  Too many columns.   The limit was', i8,
     $   4x, '  Actual number is', i8)
 3050 format(' XXXX  Too many elements.  The limit was', i8,
     $   4x, '  Actual number is', i8)
 5200 format(' Columns', i8 / ' Elements', i7)

*     end of m3mpsb
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3mpsc( m, n, nb, ne, nint, ns, lennm,
     $                   key, ncard, name1, name2,
     $                   bl, bu, hint, hs, xn )

      implicit           double precision (a-h,o-z)
      character*4        key
      integer            hint(n), hs(nb)
      integer            ncard(6), name1(lennm), name2(lennm)
      double precision   bl(nb), bu(nb)
      double precision   xn(nb)

*     ------------------------------------------------------------------
*     m3mpsc inputs the BOUNDS section of an MPS file.
*
*     21 Oct 1997: nint and hint(*) are new I/O parameters.
*                  Now recognize bound types BV, LI, UI for integers.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m3mps5/ aelem(2), id(6), iblank
      common    /cyclcm/ cnvtol,jnew,materr,maxcy,nephnt,nphant,nprint

      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      logical            gotnm, ignore
      character*4        a4,
     $                   lboun, lenda,
     $                   lfr  , lfx, llo, lmi, lpl, lup,
     $                   lbv  , lli, lui,
     $                   linit, lial
      data               a4                   /'(a4)'/
      data               lboun, lenda         /'BOUN','ENDA'/
      data               lfr  , lfx, llo      /' FR ',' FX ',' LO '/
      data               lmi  , lpl, lup      /' MI ',' PL ',' UP '/
      data               lbv  , lli, lui      /' BV ',' LI ',' UI '/
      data               linit, lial          /'INIT', 'IAL '/
*     ------------------------------------------------------------------

      read (linit , '(a4)') iinit
      read (lial  , '(a4)') iial
      inform = 1
      bplus  = plinfy
      bminus = -bplus

*     Fix phantom variables at zero.

      j1     = n - nphant + 1
      call dload ( nphant, zero, bl(j1), 1 )
      call dload ( nphant, zero, bu(j1), 1 )

*     Check for no BOUNDS.

      if (key .ne. lboun) go to 700
      gotnm  = mbnd(1). ne. iblank
      inform = 0
      jmark  = 1

*     ==================================================================
*     Read and check BOUNDS cards.  Notice the double plural.
*     ==================================================================
  610 call m3read( 3, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 700

*     A normal bounds set is terminated if INITIAL is found.

      bnd     = aelem(1)
      if (id(1) .eq. iinit  .and.  id(2) .eq. iial ) go to 690

      if (.not. gotnm) then
         gotnm   = .true.
         mbnd(1) = id(1)
         mbnd(2) = id(2)
      end if

      if (id(1) .ne. mbnd(1)  .or.  id(2) .ne. mbnd(2)) go to 610

*     Find which column.  Input "nb" = n to search only column names.

      call m4name( m, n, n, lennm, name1, name2, id(3), id(4),
     $             line, ier(10), 0, 1, n, jmark, j )

      if (j .le. 0) then
         if (iprint .gt. 0  .and.  ier(10) .le. mer)
     $      write(iprint, 1400) id(3), id(4), line
      else

*        Select bound type for column j.

         ncard(4) = ncard(4) + 1
         if      (key .eq. lup) then
            bu(j)  = bnd
         else if (key .eq. llo) then
            bl(j)  = bnd
         else if (key .eq. lfx) then
            bu(j)  = bnd
            bl(j)  = bnd
         else if (key .eq. lfr) then
            bu(j)  = bplus
            bl(j)  = bminus
         else if (key .eq. lmi) then
            if (bu(j) .ge. bplus) bu(j)  = zero
            bl(j)  = bminus
         else if (key .eq. lpl) then
            bu(j)  = bplus
*---        bl(j)  = zero
*---        31-oct-1989: John Stone suggests delete this line
*---                     and add the "if" above.
         else if (key .eq. lbv) then
            if (hint(j) .eq. 0) then
               nint    = nint + 1
               hint(j) = 1
            end if
            bl(j)  = zero
            bu(j)  = one
         else if (key .eq. lli) then
            if (hint(j) .eq. 0) then
               nint    = nint + 1
               hint(j) = 1
            end if
            bl(j)  = bnd
         else if (key .eq. lui) then
            if (hint(j) .eq. 0) then
               nint    = nint + 1
               hint(j) = 1
            end if
            bu(j)  = bnd
         else
*           This lad didn't even make it to Form 1.

            ier(11) = ier(11) + 1
            if (iprint .gt. 0  .and.  ier(11) .le. mer)
     $         write(iprint, 1700) line, key, (id(i), i=1,4)
         end if
      end if
      go to 610

*     INITIAL bounds set found.

  690 if (ncard(4) .eq. 0) then
         mbnd(1) = iblank
         mbnd(2) = iblank
         if (iprint .gt. 0) write(iprint, 1720)
      end if

*     ------------------------------------------------------------------
*     End of normal bounds.
*     ------------------------------------------------------------------
  700 ns     = 0
      bplus  = 0.9*bplus
      bminus =   - bplus

*     Set variables to be nonbasic at zero (as long as that's feasible).

      do 706 j = 1, nb
         xn(j) = max( zero , bl(j) )
         xn(j) = min( xn(j), bu(j) )
         hs(j) = 0
         if (xn(j) .eq. bu(j)) hs(j) = 1
  706 continue

*     Ignore INITIAL bounds if a basis will be loaded.

      if (inform .ne. 0) go to 790
      ignore = ioldb .gt. 0  .or.  insrt .gt. 0  .or.  iload .gt. 0
      if (.not. ignore ) then
         jmark  = 1
         go to 720
      end if

*     ==================================================================
*     Read INITIAL bounds set.
*     ==================================================================
  710 call m3read( 3, imps, line, mlst, key, inform )
      if (inform .ne. 0) go to 790

      bnd    = aelem(1)
      if (ignore  .or.  id(1).ne.iinit  .or.  id(2).ne.iial) go to 710

*     Find which column.  Input "nb" = n to search only column names.

  720 call m4name( m, n, n, lennm, name1, name2, id(3), id(4),
     $             line, ier(12), 0, 1, n, jmark, j )

      if (j .le. 0) then
         if (iprint .gt. 0  .and.  ier(12) .le. mer)
     $      write(iprint, 1400) id(3), id(4), line
      else

*        Select bound type for column j.

         ncard(6) = ncard(6)+1
         if      (key .eq. lfr) then
            js  = -1
         else if (key .eq. lfx) then
            js  =  2
            ns  = ns + 1
         else if (key .eq. llo) then
            js  =  4
            bnd = bl(j)
         else if (key .eq. lup) then
            js  =  5
            bnd = bu(j)
         else if (key .eq. lmi) then
            js  =  4
         else if (key .eq. lpl) then
            js  =  5
         else
            ier(13) = ier(13) + 1
            if (iprint .gt. 0  .and.  ier(13) .le. mer)
     $         write(iprint, 1700) line, key, (id(i), i=1,4)
            go to 710
         end if
      end if

      if (abs( bnd ) .ge. bplus) bnd = zero
      xn(j)  = bnd
      hs(j)  = js
      go to 710

*     Should be ENDATA card.

  790 if (key .ne. lenda) then
         ier(14) = 1
         if (iprint .gt. 0) write(iprint, 1150)
         if (isumm  .gt. 0) write(isumm , 1150)
      end if

*     ------------------------------------------------------------------
*     Pass the buck - not got to Truman yet.
*     ------------------------------------------------------------------

*     Check that  bl .le. bu

      do 802 j = 1, n
         b1    = bl(j)
         b2    = bu(j)
         if (b1 .gt. b2) then
            ier(20) = ier(20) + 1
            if (iprint .gt. 0  .and.  ier(20) .le. mer)
     $         write(iprint, 1740) j, b1, b2
            bl(j) = b2
            bu(j) = b1
         end if
  802 continue

*     Count the errors.

      k      = 0
      do 804 i = 1, 20
         k     = k + ier(i)
  804 continue
      if (k .gt. 0) then
         if (iprint .gt. 0) write(iprint, 1900) k
         if (isumm  .gt. 0) write(isumm , 1900) k
      end if
      if (iprint .gt. 0) then
         write(iprint, 2100) mobj, minmax, ncard(1),
     $                       mrhs, ncard(2),
     $                       mrng, ncard(3),
     $                       mbnd, ncard(4)
      end if
      return

 1150 format(' XXXX  ENDATA card not found')
 1400 format(' XXXX  Non-existent column specified -- ', 2a4,
     $   ' -- entry ignored in line', i7)
 1700 format(' XXXX  Illegal bound type at line', i7, '... ',
     $   a1, a2, a1, 2a4, 2x, 2a4)
 1720 format(' ===>  Warning - first bounds set is  INITIAL .',
     $   '   Other bounds will be ignored.')
 1740 format(/' XXXX  Bounds back to front on column', i6,' :',
     $   1p, 2e15.5)
 1900 format(/' XXXX  Total no. of errors in MPS file', i6)
 2100 format(///
     $   ' Names selected' /
     $   ' --------------' /
     $   ' Objective', 6x, 2a4, ' (', a3, ')', i8 /
     $   ' RHS      ', 6x, 2a4, i14 /
     $   ' RANGES   ', 6x, 2a4, i14 /
     $   ' BOUNDS   ', 6x, 2a4, i14)

*     end of m3mpsc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m3read( mode, imps, line, mxlist, key, inform )

      implicit           double precision (a-h,o-z)
      character*4        key

*     ------------------------------------------------------------------
*     m3read  reads data from file imps and prints a listing on file
*     iprint.  The data is assumed to be in MPS format, with items of
*     interest in the following six fields...
*
*     Field:     1         2         3         4         5         6
*
*     Columns: 01-04     05-12     15-22     25-36     40-47     50-61
*
*     Format:    a4       2a4       2a4      e12.0      2a4      e12.0
*
*     Data:     key     id(1-2)   id(3-4)  aelem(1)   id(5-6)  aelem(2)
*
*
*     Comments may contain a * in column 1 and anything in columns 2-61.
*     They are listed and then ignored.
*
*
*     On entry,  mode    specifies which fields are to be processed.
*     On exit ,  inform  is set to 1 if column 1 is not blank.
*
*     04 Oct 1985: First version.
*     27 Sep 1991: More f77.  Each line read into a character buffer
*                  to simplify handling of comments and key.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m3mps5/ aelem(2), id(6), iblank

      character*61       buffer
      character*1        buff1
      character*1        lblank      , lstar
      data               lblank /' '/, lstar /'*'/

*     ------------------------------------------------------------------
*     Read a data card and look for keywords and comments.
*     ------------------------------------------------------------------
   10 read (imps, 1000) buffer
      buff1  = buffer(1:1)
      line   = line + 1

*     Print the buffer if column 1 is nonblank
*     or if a listing is wanted.

      if (buff1 .ne. lblank  .or.  line .le. mxlist) then

*        Find the last nonblank character.

         do 20 last = 61, 2, -1
            if (buffer(last:last) .ne. lblank) go to 30
   20    continue
         last   = 1

   30    if (iprint .gt. 0) write(iprint, 2000) line, buffer(1:last)
      end if

*     Ignore comments.

      if (buff1 .eq. lstar ) go to 10

*     If column 1 is nonblank, load key and exit.
*     The NAME card is unusual in having some data in field 3.
*     We have to load it into id(3)-id(4).

      if (buff1 .ne. lblank) then
         read(buffer, 1100) key, id(1), id(2), id(3), id(4)
         inform = 1
         return
      end if

*     ------------------------------------------------------------------
*     Process normal data cards.
*     ------------------------------------------------------------------
      if (mode .eq. 1) then

*        NAME or ROWS sections.

         read(buffer, 1100) key, id(1), id(2), id(3), id(4)
      else if (mode .eq. 2) then

*        COLUMNS, RHS or RANGES sections.

         read(buffer, 1100) key,   id(1), id(2),
     $                      id(3), id(4), aelem(1),
     $                      id(5), id(6), aelem(2)
      else

*        BOUNDS section.

         read(buffer, 1100) key, id(1), id(2), id(3), id(4), aelem(1)
      end if

      return

 1000 format(a61)
 1100 format(a4, 2a4, 2x, 2a4, 2x, bn, e12.0, 3x, 2a4, 2x, e12.0)
 2000 format(i7, 4x, a)

*     end of m3read
      end
