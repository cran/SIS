************************************************************************
*
*     File  mi02lib.f   f77 (compatible with f90)
*
*     mititle  Obtain MINOS version/date string if you want it.
*     mistart  Initialize MINOS.  CALL THIS BEFORE OTHER ROUTINES.
*     mispec   Read a SPECS file.
*     micore   Estimate length of work array z(*).
*     minoss   The primary solve routine.
*     minos    Alternative solve routine allowing constraint bounds
*              the same as SNOPT.  (But it calls minoss.)
*     minose   Embedded minos.  Same as minoss but some extra parameters
*              for use by environments like GAMS.
*
*     miopt    Set an option.
*     miopti   Set an integer option.
*     mioptr   Set a real option.
*
*     micjac   Count Jacobian elements.
*     micmps   Count MPS file rows, columns, elements.
*     mirmps   Read  an MPS file.
*     miwmps   Write an MPS file.
*     misolf   Set solution flags A D I N and LL UL SBS BS EQ FR.
*
*     12 Jul 2000: mi02lib.f  contains callable routines.
*     13 Aug 2000: mistart, mititle added to match SNOPT's
*                  snInit , sntitl.
*     18 Aug 2000: micjac , micmps  added.
*     20 Aug 2000: mirmps , miwmps  added.
*     23 Sep 2000: micore added.
*     02 Feb 2004: (At GAMS) minose implemented.
*     14 Feb 2004: misolf implemented.
*     10 Mar 2004: m5lpit and m7rgit prevent new nonbasics from being
*                  slightly inside their bound.  This fixes strange bug
*                  in "Max Dual infeas" following EXIT msg.
*     17 Jun 2004: misolf always flags infeasible jBinf1 as I.
*     22 Jul 2004: m5solv:
*                  Label 480 added after m5dgen sets nonbasics on bnd.
*                  Needed to ensure that m5frmc sets nphs correctly.
*     31 Jul 2004: m5solv:
*                  If infeasible and ns > 0, make sure tolrg is set
*                  using min of toldj(3) and 1e-6 (like in m5prc).
*     18 Aug 2005: minose: Added ncstate, vmax, vrel as output
*                  (responding to David Gay's suggestion about isfeas).
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mititle( title )

      implicit            none
      character*30        title

*     ==================================================================
*     mititle sets the title for MINOS.
*     ==================================================================
*     ----------123456789|123456789|123456789|--------------------------

      title  = 'M I N O S  5.51     (Aug 2005)'

      end ! subroutine mititle

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mistart( iprint, isumm, ispecs )

      implicit            none
      integer             iprint, isumm, ispecs

*     ==================================================================
*     mistart  must be called before any other MINOS library routine.
*     Actions:
*     1. Open default files (PRINT, SUMMARY and possibly SPECS).
*        It's ok to have ispecs = 0 and let mispec worry about it.
*     2. Initialize title.
*     3. Set options to default values.
*
*     13 Aug 2000: First version of mistart.
*     23 Sep 2000: Parameters iprint, isumm, ispecs have normal names.
*                  Now have   jprint, jsumm, jspecs in common blocks.
*     ==================================================================

      integer            iread,jprint,jsumm
      common    /m1file/ iread,jprint,jsumm

      integer            iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,jspecs,ireprt
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,jspecs,ireprt

      character*30       title

      jprint = iprint       ! Save specified file numbers.
      jsumm  = isumm
      jspecs = ispecs
      call mifile( 1 )      ! Open the PRINT, SUMMARY and SPECS files.
      call m1init( )        ! Set a few constants.
      call mititle( title ) ! Get title.

      call m1page( 1 )      ! Indicate new page, then print title.
      if (iPrint .gt. 0) write (iPrint, 1000) title
      if (iSumm  .gt. 0) write (iSumm , 1000) title

      call m3dflt( 1 )      ! Set the options to default values.
      return

 1000 format( 6x, '=============================='
     &      / 6x, a
     &      / 6x, '==============================')

      end ! subroutine mistart

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mispec( ispecs, inform )

*     ------------------------------------------------------------------
*     mispec  reads a SPECS file from unit ispecs.
*     If ispecs is not in the range [1,99], nothing happens.
*
*     01 Oct 1991: First version.
*     13 Aug 2000: mispec now reads a SPECS file and does nothing else.
*                  Other preliminary stuff moved to mistart.
*     23 Sep 2000: ispecs is now input parameter.
*                  jspecs is in common block.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,jspecs,ireprt

      external           m3key

      jspecs = ispecs
      ncalls = 1
      inform = 0

*     ------------------------------------------------------------------
*     Read the Specs file (if any).
*     ------------------------------------------------------------------
      if (ispecs .gt. 0) then
         call m3file( ncalls, ispecs, m3key,
     $                iprint, isumm, inform )
      end if

      end ! subroutine mispec

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine micore( m, n, ne, nscl, maxr, maxs,
     &                   nnobj, nncon, nnjac, nejac,
     &                   mincor )

      implicit           none
      integer            m, n, ne, nscl, maxr, maxs,
     &                   nnobj, nncon, nnjac, nejac,
     &                   mincor

*     ------------------------------------------------------------------
*     micore  estimates the amount of storage MINOS needs
*     for a problem with the specified dimensions.
*
*     On exit,
*     mincor  may be used to allocate an array real*8 z(1:mincor)
*             (or bigger) to be used in the call to minoss.
*
*     19 Sep 2000: First version of micore.  About time!
*                  (Requested long ago by Gerd Infanger.)
*     ------------------------------------------------------------------

      integer            mbs, nb, nn, nr, nx, lenlam, lenx
      integer            mina, minLU, necola

      ! We know the following things exactly (excluding LU arrays).

      mbs    = m + maxs
      nb     = n + m
      nn     = max( nnobj, nnjac, 1 )
      nr     = maxr*(maxr + 1)/2  +  (maxs - maxr)
      nx     = max( mbs, nn )

      if (nncon .eq. 0) then
         lenlam = 0
         lenx   = 0
      else
         lenlam = m
         lenx   = nb
      end if

      mincor =    6*mbs
     &         +    nn
     &         +  5*nx
     &         +    nscl
     &         +    nr
     &         +  2*maxs
     &         +  2*nnobj
     &         +  5*nncon
     &         +  2*nejac
     &         +    lenlam
     &         +  2*lenx

      ! Now estimate the LU storage.
      ! necola = estimate of nonzeros per column of A.
      ! We guess that the density of the basis factorization is
      ! 6 times as great.

      necola = max( (ne/n), 10 )
      mina   = 6 * min( m, n ) * necola

      minLU  =    3*mina    ! 1 real*8 + 2 integer*4 + 1 for elbow room
     &         +  4*mbs
     &         +  4*m

      ! Total storage estimate.

      mincor = mincor + minLU

      end ! subroutine micore

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine minoss( start, m, n, nb, ne, nname,
     $                   nncon, nnobj, nnjac, 
     $                   iobj, objadd, names,
     $                   a, ha, ka, bl, bu, name1, name2,
     $                   hs, xn, pi, rc, 
     $                   inform, mincor, ns, ninf, sinf, obj,
     $                   z, nwcore )

      implicit           double precision (a-h,o-z)
      character*(*)      start
      character*8        names(5)
      integer            ha(ne), hs(nb)
      integer            ka(n+1), name1(nname), name2(nname)
      double precision   a(ne), bl(nb), bu(nb)
      double precision   xn(nb), pi(m), rc(nb), z(nwcore)

*     ------------------------------------------------------------------
*     minoss (pronounced minos-s) is the subroutine version of MINOS.
*     It has all the data passed to it as parameters instead of reading
*     it from an MPS file.
*
*     ON ENTRY:
*
*     start   specifies how a starting basis (and certain other items)
*             are to be obtained.
*             start = 'Cold' means that Crash should be used to choose
*                      an initial basis, unless a basis file is given
*                      via Old basis, Insert or Load in the Specs file.
*             start = 'Basis file' means the same (but is more
*                      meaningful in the latter case).
*             start = 'Warm' means that a basis is already defined in hs
*                      (probably from an earlier call).
*             start = 'Hot' or 'Hot FHS' implies a hot start.
*                      hs defines a basis and an earlier call has
*                      defined certain other things that should also be
*                      kept.  The problem dimensions and the array z(*)
*                      must not have changed.
*                      F refers to the LU factors of the basis.
*                      H refers to the approximate reduced Hessian R.
*                      S refers to column and row scales.
*             start = 'Hot H' (for example) means that only the Hessian
*                      is defined.
*
*     m       is the number of general constraints.  For LP problems
*             this means the number of rows in the constraint matrix A.
*             m > 0 in principle, though sometimes m = 0 may be ok.
*             (Strictly speaking, Fortran declarations of the form
*                double precision   pi(m)
*             require m > 0.  In debug mode, compilers will probably
*             enforce m > 0, but optimized code may run ok with m = 0.)
*
*     n       is the number of variables, excluding slacks.
*             For LP problems, this is the number of columns in A.
*             n > 0.
*
*     nb      is n + m.
*
*     ne      is the number of nonzero entries in A (including the
*             Jacobian for any nonlinear constraints).
*             ne > 0 in principle, though again m = 0, ne = 0 may work
*             with some compilers.
*
*     nname   is the number of column and row names provided in the
*             arrays name1 and name2.  If nname = 1, there are NO names.
*             Generic names will be used in the printed solution.
*             Otherwise, nname = nb and all names must be provided.
*
*     nncon   is the number of nonlinear constraints.
*             nncon ge 0.
*
*     nnobj   is the number of nonlinear objective variables.
*             nnobj ge 0.
*
*     nnjac   is the number of nonlinear Jacobian variables.
*             If nncon = 0, nnjac = 0.
*             If nncon > 0, nnjac > 0.
*
*     iobj    says which row of A is a free row containing a linear
*             objective vector  c  (iobj = 0 if none).
*             iobj = 0  or  nncon < iobj le m.
*
*     objadd  is a constant that will be added to the objective.
*             Typically objadd = 0.0d+0.
*
*     names(5)is a set of 8-character names for the problem, the linear
*             objective, the rhs, the ranges and bounds.  (This is a
*             hangover from MPS files.  The names are used in the
*             printed solution and in some of the basis files.)
*
*     a(ne)   is the constraint matrix (Jacobian), stored column-wise.
*
*     ha(ne)  is the list of row indices for each nonzero in a(*).
*
*     ka(n+1) is a set of pointers to the beginning of each column of
*             the constraint matrix within a(*) and ha(*).
*             Must have ka(1) = 1 and ka(n+1) = ne+1.
*
*  NOTES:  1. If the problem has a nonlinear objective,
*             the first nnobj columns of a and ha belong to the
*             nonlinear objective variables.
*             Subroutine funobj deals with these variables.
*          
*          2. If the problem has nonlinear constraints,
*             the first nnjac columns of a and ha belong to the
*             nonlinear Jacobian variables, and
*             the first nncon rows of a and ha belong to the
*             nonlinear constraints.
*             Subroutine funcon deals with these variables and
*             constraints.
*          
*          3. If nnobj > 0 and nnjac > 0, the two sets of
*             nonlinear variables overlap.  The total number of
*             nonlinear variables is nn = max( nnobj, nnjac ).
*          
*          4. The Jacobian forms the top left corner of a and ha.
*             If a Jacobian column j (1 le j le nnjac) contains
*             any entries a(k), ha(k) associated with nonlinear
*             constraints (1 le ha(k) le nncon), those entries must
*             come before any other (linear) entries.
*          
*          5. The row indices ha(k) for a column may be in any order
*             (subject to Jacobian entries appearing first).
*             Subroutine funcon must define Jacobian entries in the
*             same order.
*          
*          6. If column j contains no entries, perhaps
*             ka(j) = ka(j+1) is acceptable.  (Must check this.
*             When MINOS reads an MPS with empty columns,
*             it inserts a dummy entry a(k) = 0.0d+0, ha(k) = 1.
*             This may not be necessary.)
*
*          7. To allocate storage, MINOS needs to know if the Jacobian
*             is dense or sparse.  The default is dense.  If this is
*             not appropriate, define
*                              Jacobian    Sparse
*             in the Specs file, or
*                call miopt ( 'Jacobian    Sparse', 0, 0, inform )
*             before calling minoss.
*                
*     bl(nb)  is the lower bounds on the variables and slacks (x, s).
*
*     bu(nb)  is the upper bounds on (x, s).
*
*     BEWARE: MINOS represents general constraints as Ax + s = 0.
*             Constraints of the form   l <= Ax <= u  
*             therefore mean            l <= -s <= u,
*             so that                  -u <=  s <= -l.
*             The last m components of bl and bu are -u and -l.
*
*     name1(nname), name2(nname) are two integer arrays.
*             If nname =  1, name1 and name2 are not used.  The printed
*             solution will use generic names for the columns and rows.
*             If nname = nb, name1(j) and name2(j) should contain the
*             name of the j-th variable in 2a4 format (j = 1, nb).
*             If j = n+i, the j-th variable is the i-th row.
*
*     hs(nb)  sometimes contains a set of initial states for each
*             variable x, or for each variable and slack (x, s).
*             See the following NOTES.
*
*     xn(nb)  sometimes contains a set of initial values for each
*             variable x, or for each variable and slack (x, s).
*             See the following NOTES.
*
*  NOTES:  1. If start = 'Cold' or 'Basis file' and a basis file
*             of some sort is to be input
*             (an OLD BASIS file, INSERT file or LOAD file),
*             hs and xn need not be set at all.
*
*          2. Otherwise, hs(j) and xn(j), j=1:n, must be defined for a
*             Cold start.  (The values for j=n+1:nb need not be set.)
*             If nothing special is known about the problem, or if
*             there is no wish to provide special information,
*             you may set hs(j) = 0, xn(j) = 0.0d+0 for all j=1:n.
*             All variables will be eligible for the initial basis.
*        
*             Less trivially, to say that variable j will probably
*             be equal to one of its bounds,
*             set hs(j) = 4 and xn(j) = bl(j)
*             or  hs(j) = 5 and xn(j) = bu(j) as appropriate.
*        
*          3. For Cold starts with no basis file, a Crash procedure
*             is used to select an initial basis.  The initial basis
*             matrix will be triangular (ignoring certain small
*             entries in each column).
*             The values hs(j) = 0, 1, 2, 3, 4, 5 have the following
*             meaning:
*                
*             hs(j)    State of variable j during Crash
*        
*             0, 1, 3  Eligible for the basis.  3 is given preference.
*             2, 4, 5  Ignored.
*        
*             After Crash, columns for which hs(j) = 2 are made superbasic.
*             Other columns not selected for the basis are made
*             nonbasic at the value xn(j) if bl(j) <= xn(j) <= bu(j),
*             or at the value bl(j) or bu(j) closest to xn(j).
*
*          4. For Warm or Hot starts, all of hs(1:nb) is assumed to be
*             set to the values 0, 1, 2 or 3 (probably from some
*             previous call) and all of xn(1:nb) must have values.
*        
*     pi(m)   contains an estimate of the vector of Lagrange multipliers
*             (shadow prices) for the NONLINEAR constraints.  The first
*             nncon components must be defined.  They will be used as
*             lambda in the subproblem objective function for the first
*             major iteration.  If nothing is known about lambda,
*             set pi(i) = 0.0d+0, i = 1 to nncon.
*
*     ns      need not be specified for Cold starts,
*             but should retain its value from a previous call
*             when a Warm or Hot start is used.
*
*
*     ON EXIT:
*
*     hs(nb)  is the final state vector:
*
*                hs(j)    State of variable j    Normal value of xn(j)
*
*                  0      nonbasic               bl(j)
*                  1      nonbasic               bu(j)
*                  2      superbasic             Between bl(j) and bu(j)
*                  3      basic                  ditto
*
*             Very occasionally there may be nonbasic variables for
*             which xn(j) lies strictly between its bounds.
*             If ninf = 0, basic and superbasic variables may be outside
*             their bounds by as much as the Feasibility tolerance.
*             Note that if Scale is specified, the Feasibility tolerance
*             applies to the variables of the SCALED problem. 
*             In this case, the variables of the original problem may be
*             as much as 0.1 outside their bounds, but this is unlikely
*             unless the problem is very badly scaled.
*
*     xn(nb)  is the final variables and slacks (x, s).
*
*     pi(m)   is the vector of Lagrange multipliers (shadow prices)
*             for the general constraints.
*
*     rc(nb)  is a vector of reduced costs: rc = g - (A I)'pi, where g
*             is the gradient of the objective function if xn is feasible
*             (or the gradient of the Phase-1 objective otherwise).
*             If ninf = 0, the last m entries are -pi (negative pi).
*
*     inform  says what happened; see Chapter 6.3 of the User's Guide.
*             A summary of possible values follows:
*
*             inform   Meaning
*
*                0     Optimal solution found.
*                1     The problem is infeasible.
*                2     The problem is unbounded (or badly scaled).
*                3     Too many iterations.
*                4     Apparent stall.  The solution has not changed
*                      for a large number of iterations (e.g. 1000).
*                5     The Superbasics limit is too small.
*                6     Subroutine funobj or funcon requested termination
*                      by returning mode < 0.
*                7     Subroutine funobj seems to be giving incorrect
*                      gradients.
*                8     Subroutine funcon seems to be giving incorrect
*                      gradients.
*                9     The current point cannot be improved.
*               10     Numerical error in trying to satisfy the linear
*                      constraints (or the linearized nonlinear
*                      constraints).  The basis is very ill-conditioned.
*               11     Cannot find a superbasic to replace a basic
*                      variable.
*               12     Basis factorization requested twice in a row.
*                      Should probably be treated as inform = 9.
*               13     Near-optimal solution found.
*                      Should probably be treated as inform = 9.
*
*               20     Not enough storage for the basis factorization.
*               21     Error in basis package.
*               22     The basis is singular after several attempts to
*                      factorize it (and add slacks where necessary).
*
*               30     An OLD BASIS file had dimensions that did not
*                      match the current problem.
*               32     System error.  Wrong number of basic variables.
*
*               40     Fatal errors in the MPS file.
*               41     Not enough storage to read the MPS file.
*               42     Not enough storage to solve the problem.
*
*     mincor  says how much storage is needed to solve the problem.
*             If inform = 42, the work array z(nwcore) was too small.
*             minoss may be called again with nwcore suitably larger
*             than mincor.  (The bigger the better, since it is
*             not certain how much storage the basis factors need.)
*
*     ns      is the final number of superbasics.
*
*     ninf    is the number of infeasibilities.
*
*     sinf    is the sum    of infeasibilities.
*
*     obj     is the value  of the objective function.
*             If ninf = 0, obj includes the nonlinear objective if any.
*             If ninf > 0, obj is just the linear objective if any.
*
*     30 Sep 1991: First version.
*     06 Dec 1991: A few more output parameters.
*     10 Apr 1992: Parameter  objadd added.  Parameters reordered.
*     20 Apr 1992: Parameters nname, name1, name2 added.
*     27 Apr 1992: Parameter  mincor added to allow reentry with more
*                  storage.
*     27 Jun 1992: Parameter  start  implemented.  Passed to misolv.
*     05 Feb 1998: Always set njac here.  Ignore Jacobian = dense.
*     12 May 1998: Use m3char to load names (ready for F90).
*     25 May 2000: Bug in minoss: Should not set nden = 2 (affects CUTE).
*     ------------------------------------------------------------------
*     NOTE:
*     In /m7len / and /m8len /, nnobjx, nnconx, nnjacx are normally
*                               nnobj , nncon , nnjac .
*     Here it is better to save those names for the minoss parameters.

      common    /m2len / mrows,mcols,melms
      common    /m2mapz/ maxw  ,maxz
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m7len / fobj  ,fobj2 ,nnobjx,nnobj0
      common    /m8len / njac  ,nnconx,nncon0,nnjacx
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /cycle2/ objtru,suminf,numinf

      external           micore     ! Any external for misolv to ignore.
*     ------------------------------------------------------------------

      ifuser = 0  ! So misolv won't call dummy external "micore".

*     Initialize timers.

      call m1time( 0,0 )

*     Load the Common variables with various problem dimensions.

      mrows  = m
      mcols  = n
      melms  = ne
      nnconx = nncon
      nnobjx = nnobj
      nnjacx = nnjac

*     Say how much z(*) we've got, in case mispec wasn't called,
*     or nwcore has been altered since mispec was called.
*     This means the Specs file can't set  Workspace (TOTAL) .

      maxz   = nwcore

*     The Specs file has been read (or the options have been
*     otherwise defined).  Check that the options have sensible values.

      call m3dflt( 2 )

*     ------------------------------------------------------------------
*     Determine storage requirements using the
*     following Common variables:
*        (m2len )   mrows, mcols, melms
*        (m3len )   nscl  (determined by lscale)
*        (m5len )   maxr, maxs, nn
*        (m7len )   nnobj
*        (m8len )   njac, nncon, nnjac
*     All have to be known exactly before calling m2core( ).
*     The only one in doubt is njac, the number of Jacobian elements.
*     Count them here.
*     ------------------------------------------------------------------
      call micjac( m, n, ne, nncon, nnjac, njac,
     &             ha, ka )

*     31 May 2000: Careful!

      if (njac .eq. nnCon*nnJac) then
         !!!!! nden = 1   might wreck CUTE, which assumes Sparse
      else
         nden = 2
      end if

      call m2core( 4, mincor )

      if (mincor .gt. nwcore) then
         inform = 42
         return
      end if

*     ------------------------------------------------------------------
*     Open files needed for this problem.
*     Print the options if iprint > 0, Print level > 0 and iparm(3) > 0.
*     ------------------------------------------------------------------
      call mifile( 2 )
      call m3dflt( 3 )

*     ------------------------------------------------------------------
*     Load names into the MINOS arrays.
*     It is laborious to use m3char, but F90 will require it.
*     ------------------------------------------------------------------
      call m3char( names(1)(1:4), name(1) )
      call m3char( names(1)(5:8), name(2) )
      call m3char( names(2)(1:4), mobj(1) )
      call m3char( names(2)(5:8), mobj(2) )
      call m3char( names(3)(1:4), mrhs(1) )
      call m3char( names(3)(5:8), mrhs(2) )
      call m3char( names(4)(1:4), mrng(1) )
      call m3char( names(4)(5:8), mrng(2) )
      call m3char( names(5)(1:4), mbnd(1) )
      call m3char( names(5)(5:8), mbnd(2) )
      
*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
      mimode = 2
      nka    = n + 1

      call misolv( mimode, start, m, n, nb, ne, nka, nname,
     $             iobj, objadd, 
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc, 
     $             inform, ns, z, nwcore,
     $             ifuser, micore, majitn, minitn )

      ninf   = numinf
      sinf   = suminf
      obj    = objtru

*     Print times for all clocks (if ltime > 0).

      call m1time( 0,2 )

      end ! subroutine minoss

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine minos ( start, m, n, nb, ne, nname,
     $                   nncon, nnobj, nnjac, 
     $                   iobj, objadd, names,
     $                   a, ha, ka, bl, bu, name1, name2,
     $                   hs, xn, pi, rc, 
     $                   inform, mincor, ns, ninf, sinf, obj,
     $                   z, nwcore )

      implicit           double precision (a-h,o-z)
      character*(*)      start
      character*8        names(5)
      integer            ha(ne), hs(nb)
      integer            ka(n+1), name1(nname), name2(nname)
      double precision   a(ne), bl(nb), bu(nb)
      double precision   xn(nb), pi(m), rc(nb), z(nwcore)

*     ------------------------------------------------------------------
*     minos  allows the constraint bounds bl and bu to be input in the
*     form
*                  bl <= (   x   ) <= bu.
*                        (  A*x  )
*     Hence, the constraints are
*                  Ax - s = 0,  bl <= (x,s) <= bu
*     rather than  Ax + s = 0,  bl <= (x,s) <= bu as in minoss.
*     pi and the last m components of bl, bu, hs, xn, rc are
*     "back to front", but the other parameters are the same as in
*     minoss.
*     ------------------------------------------------------------------

      call m2swap( 1, m, n, nb, bl, bu, hs, xn, pi, rc )
      call minoss( start, m, n, nb, ne, nname,
     $             nncon, nnobj, nnjac, 
     $             iobj, objadd, names,
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc, 
     $             inform, mincor, ns, ninf, sinf, obj,
     $             z, nwcore )
      call m2swap( 2, m, n, nb, bl, bu, hs, xn, pi, rc )

      end ! subroutine minos

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine minose( start, m, n, nb, ne, nname,
     $                   nncon, nnobj, nnjac, 
     $                   iobj, objadd, names,
     $                   a, ha, ka, bl, bu, name1, name2,
     $                   hs, xn, pi, rc, 
     $                   inform, mincor, ns, ninf, sinf, obj,
     $                   z, nwcore,
     $                   ifuser, m1user, majitn, minitn,
     &                   ncstate, vmax, vrel )

      implicit           double precision (a-h,o-z)
      character*(*)      start
      character*8        names(5)
      integer            ha(ne), hs(nb)
      integer            ka(n+1), name1(nname), name2(nname)
      double precision   a(ne), bl(nb), bu(nb)
      double precision   xn(nb), pi(m), rc(nb), z(nwcore)
      integer            ifuser, majitn, minitn, ncstate
      external           m1user
      double precision   vmax  , vrel

*     ------------------------------------------------------------------
*     minose (pronounced minos-e) is the same as minoss
*     but with extra parameters at the end.  The new parameters
*     are documented here.
*
*     ON ENTRY:
*
*     ifuser = 0 means m1user is a dummy routine, not to be called;
*            = 1 means m1user is to be called every major and minor itn.
*
*     m1user     is an external (dummy or otherwise as just explained).
*                See the calls in m5solv (mi50lp.f).
*
*
*     ON EXIT:
*
*     majitn     is the number of major iterations performed.
*
*     minitn     is the number of minor iterations performed.
*
*     ncstate    specifies the final state of the nonlinear constraints.
*     ncstate= 0 means the constraint functions were never evaluated.
*            = 1 means they are satisfied to the user-specified
*                tolerance (vrel <= rowtol) 
*            = 2 means they are NOT satisfied to the user-specified
*                tolerance (vrel >  rowtol).
*
*     vmax       is the largest infeasibility (violation) for the
*                nonlinear constraints.  (vmax = 0.0 if ncstate = 0)
*
*     vrel       is the largest relative infeasibility for the
*                nonlinear constraints.  vrel = vmax/(1 + xnorm).
*                This is the value that is compared to rowtol.
*
*     02 Feb 2004: (At GAMS) First version of minose.
*             1. GAMS has always had an m1user.  Now we can add a call
*                permanently for ifuser=1, so no edits are needed.
*             2. Erwin wanted jstate to show exactly what MINOS prints
*                (the flags computed for the SCALED solution).
*             3. majitn and minitn should have been output by minoss
*                from the beginning.  Now they save people from
*                digging them out of common /m8save/.
*
*     14 Feb 2004: Solution flags no longer returned in jstate(1:nb).
*                  misolf implemented instead.
*                  GAMS may call misolf once for every j = 1:nb
*                  following a normal call to minoss or minose.
*                  The flags will now match the UNSCALED solution!!
*                  (This was always a debatable issue anyway.
*                  It won't hurt to have flag matching what users see.)
*     18 Aug 2005: Added ncstate, vmax, vrel as output
*                  (responding to David Gay's suggestion about isfeas).

*     ------------------------------------------------------------------
*     NOTE:
*     In /m7len / and /m8len /, nnobjx, nnconx, nnjacx are normally
*                               nnobj , nncon , nnjac .
*     Here it is better to save those names for the minoss parameters.

      common    /m2len / mrows,mcols,melms
      common    /m2mapz/ maxw  ,maxz
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
      common    /m7len / fobj  ,fobj2 ,nnobjx,nnobj0
      common    /m8len / njac  ,nnconx,nncon0,nnjacx
      common    /m8al1 / penpar,rowtol,ncom,nden,nlag,nmajor,nminor
      common    /m8func/ nfcon(4),nfobj(4),nprob,nstat1,nstat2
      common    /m8save/ vimax ,virel ,maxvi ,majits,minits,nssave
      common    /cycle2/ objtru,suminf,numinf
*     ------------------------------------------------------------------

*     Initialize timers.

      call m1time( 0,0 )

*     Load the Common variables with various problem dimensions.

      mrows  = m
      mcols  = n
      melms  = ne
      nnconx = nncon
      nnobjx = nnobj
      nnjacx = nnjac

*     Say how much z(*) we've got, in case mispec wasn't called,
*     or nwcore has been altered since mispec was called.
*     This means the Specs file can't set  Workspace (TOTAL) .

      maxz   = nwcore

*     The Specs file has been read (or the options have been
*     otherwise defined).  Check that the options have sensible values.

      call m3dflt( 2 )

*     ------------------------------------------------------------------
*     Determine storage requirements using the
*     following Common variables:
*        (m2len )   mrows, mcols, melms
*        (m3len )   nscl  (determined by lscale)
*        (m5len )   maxr, maxs, nn
*        (m7len )   nnobj
*        (m8len )   njac, nncon, nnjac
*     All have to be known exactly before calling m2core( ).
*     The only one in doubt is njac, the number of Jacobian elements.
*     Count them here.
*     ------------------------------------------------------------------
      call micjac( m, n, ne, nncon, nnjac, njac,
     &             ha, ka )

*     31 May 2000: Careful!

      if (njac .eq. nnCon*nnJac) then
         !!!!! nden = 1   might wreck CUTE, which assumes Sparse
      else
         nden = 2
      end if

      call m2core( 4, mincor )

      if (mincor .gt. nwcore) then
         inform = 42
         return
      end if

*     ------------------------------------------------------------------
*     Open files needed for this problem.
*     Print the options if iprint > 0, Print level > 0 and iparm(3) > 0.
*     ------------------------------------------------------------------
      call mifile( 2 )
      call m3dflt( 3 )

*     ------------------------------------------------------------------
*     Load names into the MINOS arrays.
*     It is laborious to use m3char, but F90 will require it.
*     ------------------------------------------------------------------
      call m3char( names(1)(1:4), name(1) )
      call m3char( names(1)(5:8), name(2) )
      call m3char( names(2)(1:4), mobj(1) )
      call m3char( names(2)(5:8), mobj(2) )
      call m3char( names(3)(1:4), mrhs(1) )
      call m3char( names(3)(5:8), mrhs(2) )
      call m3char( names(4)(1:4), mrng(1) )
      call m3char( names(4)(5:8), mrng(2) )
      call m3char( names(5)(1:4), mbnd(1) )
      call m3char( names(5)(5:8), mbnd(2) )
      
*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
      mimode = 2
      nka    = n + 1
      call misolv( mimode, start, m, n, nb, ne, nka, nname,
     $             iobj, objadd, 
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc, 
     $             inform, ns, z, nwcore,
     $             ifuser, m1user, majitn, minitn )

      ninf   = numinf
      sinf   = suminf
      obj    = objtru

      if (nfcon(1) .le. 0) then
         ncstate = 0
         vmax    = 0.0d+0
         vrel    = 0.0d+0
      else
         if (vrel .le. rowtol) then
            ncstate = 1
         else
            ncstate = 2
         end if
         vmax    = vimax
         vrel    = virel
      end if

*     Print times for all clocks (if ltime > 0).

      call m1time( 0,2 )

      end ! subroutine minose

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine miopt ( buffer, iprint, isumm, inform )

      implicit           double precision (a-h,o-z)
      character*(*)      buffer

*     ------------------------------------------------------------------
*     miopt  decodes the option contained in  buffer.
*
*     The buffer is output to file iprint, minus trailing blanks.
*     Error messages are output to files iprint and isumm.
*     buffer is echoed to iprint but normally not to isumm.
*     It is echoed to isumm before any error msg.
*
*     On entry,
*     iprint is the Print   file.  No output occurs if iprint .le 0.
*     isumm  is the Summary file.  No output occurs if isumm  .le 0.
*     inform is the number of errors so far.
*
*     On exit,
*     inform is the number of errors so far.
*
*     27 Nov 1991: First version.
*     ------------------------------------------------------------------

      character*16       key

      call m3key ( buffer, key, iprint, isumm, inform )

      end ! subroutine miopt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine miopti( buffer, ivalue, iprint, isumm, inform )

      implicit           double precision (a-h,o-z)
      character*(*)      buffer
      integer            ivalue

*     ------------------------------------------------------------------
*     miopti decodes the option contained in  buffer // ivalue.
*     The parameters other than ivalue are as in miopt.
*
*     27 Nov 1991: First version.
*     17 Jan 1992: buff72 needed to comply with f77 standard.
*     ------------------------------------------------------------------

      character*16       key
      character*72       buff72

      write(key, '(i16)') ivalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      call m3key ( buff72, key, iprint, isumm, inform )

      end ! subroutine miopti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mioptr( buffer, rvalue, iprint, isumm, inform )

      implicit           double precision (a-h,o-z)
      character*(*)      buffer
      double precision   rvalue

*     ------------------------------------------------------------------
*     mioptr decodes the option contained in  buffer // rvalue.
*     The parameters other than rvalue are as in miopt.
*
*     27 Nov 1991: First version.
*     17 Jan 1992: buff72 needed to comply with f77 standard.
*     ------------------------------------------------------------------

      character*16       key
      character*72       buff72

      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      call m3key ( buff72, key, iprint, isumm, inform )

      end ! subroutine mioptr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine micjac( m, n, ne, nncon, nnjac, nejac,
     &                   ha, ka )

      implicit           none
      integer            m, n, ne, nncon, nnjac, nejac
      integer            ha(ne), ka(n+1)

!     ------------------------------------------------------------------
!     micjac (count Jacobian) counts how many elements are in the
!     nncon x nnjac top left-hand corner of the m x n sparse matrix
!     defined by ha, ka.
!
!     On entry:
!     All items except nejac are input data.
!
!     On exit:
!     nejac  returns the number of Jacobian elements.
!
!     18 Aug 2000: First version of micjac.
!                  Note that we use nejac here to match SNOPT,
!                  but MINOS uses   njac  internally.
!     ------------------------------------------------------------------

      integer            k, last

      nejac  = 0

      if (nncon .gt. 0) then
         last = ka(nnjac+1) - 1
         if (nncon .eq. m) then
            nejac  = last
         else
            do k = 1, last
               if (ha(k) .le. nncon) nejac = nejac + 1
            end do
         end if
      end if

      end ! subroutine micjac

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine micmps( imps, m, n, ne, nint, inform )

      implicit           none
      integer            imps, m, n, ne, nint, inform

!     ------------------------------------------------------------------
!     micmps  reads an MPS file from unit imps
!     and counts the important items in the file.
!
!     On entry:
!        imps   = unit number for MPS file
!
!     On exit:
!        m      = number of rows
!        n      = number of columns
!        ne     = number of entries in the columns section
!        nint   = number of integer bound types in FIRST BOUNDS section
!        inform = 0 if MPS file was processed correctly
!               = 1 if NAME    record wasn't found
!               = 2 if ROWS    record wasn't found
!               = 3 if COLUMNS record wasn't found
!               = 4 if ENDATA  record wasn't found
!               = 5 if unexpected END OF FILE occurred
!
!     29 Jul 2000: f77 version of micmps derived from f90 version.
!     ------------------------------------------------------------------

      character*1               c1            ! character in column 1
      character*2               c2            ! bound type
      character*47              buff          ! longer record
      character*8               name1, name2  ! column or bound names
      logical                   first         ! detect first bound name
      character*3               a             ! format
      integer                   i, infty      ! for infinite loops

      a      = '(a)'
      infty  = 1000000000

      m      = 0
      n      = 0
      ne     = 0
      nint   = 0
      inform = 0

      rewind imps

!     Find NAME record.

      do 100 i = 1, infty
         read(imps, a, end=900) c1
         if (c1 .eq. '*') go to 100
         if (c1 .ne. ' ') go to 110
  100 continue

  110 if (c1 .ne. 'N') then
         inform = 1
         return
      end if

!     Find ROWS record.

      do 200 i = 1, infty
         read(imps, a, end=900) c1
         if (c1 .eq. '*') go to 200
         if (c1 .ne. ' ') go to 210
  200 continue

  210 if (c1 .ne. 'R') then
         inform = 2
         return
      end if

!     Read ROWS section, counting rows.

      do 300 i = 1, infty
         read(imps, a, end=900) c1
         if (c1 .eq. '*') go to 300
         if (c1 .ne. ' ') go to 310
         m     = m + 1
  300 continue

  310 if (c1 .ne. 'C') then
         inform = 3
         return
      end if

!     Read COLUMNS sections, counting columns and elements.

      name1  = ' '
      do 400 i = 1, infty
         read(imps, a, end=900) buff
         c1     = buff(1:1)
         if (c1 .eq. '*') go to 400
         if (c1 .ne. ' ') go to 410
         name2  = buff(5:12)
         if (name1 .ne. name2) then
            n      = n + 1
            name1  = name2
         end if
         if (buff(15:22) .ne. '        ') ne    = ne + 1
         if (buff(40:47) .ne. '        ') ne    = ne + 1
  400 continue

!     Look for BOUNDS record.
!     Skip RHS or RANGE sections.
!     Stop if ENDATA.

  410 if (c1 .ne. 'B') then
         do 500 i = 1, infty
            read(imps, a, end=900) c1
            if (c1 .eq. '*') go to 500
            if (c1 .eq. 'B') go to 510
            if (c1 .eq. 'E') go to 510
  500    continue
      end if

!     Read BOUNDS section, counting integer bound types.
!     NB: For simplicity, process the FIRST bounds set only.

  510 if (c1 .eq. 'B') then
         first  = .true.
         do 600 i = 1, infty
            read(imps, a, end=900) buff(1:12)
            c1     = buff(1:1)
            if (c1 .eq. '*') go to 600
            if (c1 .ne. ' ') go to 610

            name2  = buff(5:12)
            if (first) then
                first  = .false.
                name1  = name2
            else
                if (name1 .ne. name2) go to 610
            end if

            c2     = buff(2:3)
            if (c2 .eq. 'BV') nint  = nint + 1
            if (c2 .eq. 'IV') nint  = nint + 1
  600    continue
      end if

!     We should have reached ENDATA.

  610 if (c1 .ne. 'E') then
         inform = 4
      end if

      go to 990
     
!     EOF encountered unexpectedly.

  900 inform = 5

!     Exit.

  990 rewind imps
      return

      end ! subroutine micmps

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mirmps( imps, maxm, maxn, maxnb, maxne,
     $                   nncon, nnjac, nnobj,
     $                   m, n, nb, ne, nint,
     $                   iobj, objadd, names,
     $                   a, ha, ka, bl, bu, name1, name2,
     $                   hint, hs, xn, pi,
     $                   inform, ns, z, nwcore )

      implicit           double precision (a-h,o-z)
      character*8        names(5)
      integer            ha(maxne) , hint(maxn)  , hs(maxnb)
      integer            ka(maxn+1), name1(maxnb), name2(maxnb)
      double precision   a(maxne)  , bl(maxnb)   , bu(maxnb)
      double precision   xn(maxnb) , pi(maxm)    , z(nwcore)

*     ------------------------------------------------------------------
*     mirmps  inputs constraint data for a linear or nonlinear program
*     in MPS format, consisting of NAME, ROWS, COLUMNS, RHS, RANGES and
*     BOUNDS sections in that order.  The RANGES and BOUNDS sections are
*     optional.
*
*     In the LP case, MPS format defines a set of constraints of the
*     form
*              l <= x <= u,      b1 <=  Ax  <= b2,
*     where l and u are specified by the BOUNDS section, and b1 and b2
*     are defined somewhat indirectly by the ROWS, RHS and RANGES
*     sections.  mirmps converts these constraints into the equivalent
*     form
*              Ax + s = 0,       bl <= ( x ) <= bu,
*                                      ( s )
*     where s is a set of slack variables.  This is the way MINOS deals
*     with the data.  The first n components of bl and bu are the same
*     as l and u.  The last m components are -b2 and -b1.
*
*     MPS format gives 8-character names to the rows and columns of A.
*     One of the rows of A may be regarded as a linear objective row.
*     This will be row iobj, where iobj = 0 means there is no such row.
*
*     The data defines a linear program if nncon = nnjac = nnobj = 0.
*     The nonlinear case is the same except for a few details.
*     1. If nncon = nnjac = 0 but nnobj > 0, the first nnobj columns
*        are associated with a nonlinear objective function.
*     2. If nncon > 0, then nnjac > 0 and nnobj may be zero or positive.
*        The first nncon rows and the first nnjac columns are associated
*        with a set of nonlinear constraints.
*     3. Let nn = max( nnjac, nnobj ).  The first nn columns correspond
*        to "nonlinear variables".
*     4. If an objective row is specified (iobj > 0), then it must be
*        such that iobj > nncon.
*     5. "Small" elements (below the Aij tolerance) are ignored only if
*        they lie outside the nncon by nnjac Jacobian, i.e. outside
*        the top-left corner of A.
*     6. No warning is given if some of the first nn columns are empty.
*
*
*     ON ENTRY
*     ========
*     imps   is the unit containing the MPS file.  On some systems, it
*            may be necessary to open file imps before calling mirmps.
*
*     maxm   is an overestimate of the number of rows in the ROWS
*            section of the MPS file.
*
*     maxn   is an overestimate of the number of columns in the COLUMNS
*            section of the MPS file.
*
*     maxnb  is maxm + maxn.
*
*     maxne  is an overestimate of the number of elements (matrix
*            coefficients) in the COLUMNS section.
*
*     nncon  is the no. of nonlinear constraints in the problem.
*            These must be the FIRST rows in the ROWS section.
*
*     nnjac  is the no. of nonlinear Jacobian variables in the problem.
*            These must be the FIRST columns in the COLUMNS section.
*
*     nnobj  is the no. of nonlinear objective variables in the problem.
*            These must be the FIRST columns in the COLUMNS section,
*            overlapping where necessary with the Jacobian variables.
*
*     names  is an array of five 8-character names.
*            names(1) need not be specified... it will be changed to
*            the name on the NAME card of the MPS file.
*            names(2) is the name of the objective row to be selected
*            from the ROWS section, or blank if mirmps should select
*            the first type N row encountered.
*            Similarly,
*            names(3), names(4) and names(5) are the names of the
*            RHS, RANGES and BOUNDS to be selected from the
*            RHS, RANGES and BOUNDS sections respectively, or blank
*            if mirmps should select the first ones encountered.
*
*     z(*)   is a workspace array of length nwcore.  It is needed to
*            hold the row-name hash table and a few other things.
*
*     nwcore is the length of z(*).  It should be at least 4*maxm.
*
*
*     ON EXIT
*     =======
*     m      is the number of rows in the ROWS section.
*
*     n      is the number of columns in the COLUMNS section.
*
*     nb     is n + m.
*
*     ne     is the number of matrix coefficients in COLUMNS section.
*
*     nint   is the number of integer variables detected.
*
*     iobj   is the row number of the specified objective row,
*            or zero if no such row was found.
*
*     objadd is a real constant extracted from row iobj of the RHS.
*            It is zero if the RHS contained no objective entry.
*            MINOS adds objadd to the objective function.
*
*     names(1)-names(5) contain the names of the
*               Problem, Objective row, RHS, RANGES and BOUNDS
*            respectively.
*
*     a(*)   contains the ne entries for each column of the matrix
*            specified in the COLUMNS section.
*
*     ha(*)  contains the corresponding row indices.
*
*     ka(j)  (j = 1 to n) points to the beginning of column j
*            in the parallel arrays a(*), ha(*).
*     ka(n+1) = ne+1.
*
*     bl(*)  contains nb lower bounds for the columns and slacks.
*            If there is no lower bound on x(j), then bl(j) = - 1.0d+20.
*
*     bu(*)  contains nb lower bounds for the columns and slacks.
*            If there is no upper bound on x(j), then bu(j) = + 1.0d+20.
*
*     name1(*), name2(*) contain nb column and row names in 2a4 format.
*            The j-th column name is stored in name1(j) and name2(j).
*            The i-th row    name is stored in name1(k) and name2(k),
*            where k = n + i.  
*
*     hint(j) = 0 if x(j) is continuous,
*             = 1 if x(j) is integer.
*
*     hs(*)  contains an initial state for each column and slack.
*
*     xn(*)  contains an initial value for each column and slack.
*
*            If there is no INITIAL bounds set,
*               xn(j) = 0 if that value lies between bl(j) and bu(j),
*                     = the bound closest to zero otherwise,
*               hs(j) = 0 if xn(j) < bu(j),
*                     = 1 if xn(j) = bu(j).
*
*            If there is an INITIAL bounds set, xn(j) and hs(j) are
*            set as follows.  Suppose the j-th variable has the name Xj,
*            and suppose any numerical value specified happens to be 3.
*                                                   xn(j)    hs(j)
*             FR INITIAL   Xj         3.0           3.0       -1
*             FX INITIAL   Xj         3.0           3.0        2
*             LO INITIAL   Xj                       bl(j)      4
*             UP INITIAL   Xj                       bu(j)      5
*             MI INITIAL   Xj         3.0           3.0        4
*             PL INITIAL   Xj         3.0           3.0        5
*
*     pi(*)  contains a vector defined by a special RHS called LAGRANGE.
*            If the MPS file contains no such RHS, pi(i) = 0.0, i=1:m.
*
*     inform =  0 if no fatal errors were encountered,
*            = 40 if the ROWS or COLUMNS sections were empty
*                 or iobj > 0 but iobj <= nncon,
*            = 41 if maxm, maxn or maxne were too small.
*
*     ns     is the no. of FX INITIAL entries in the INITIAL bounds set.
*              
*
*     19 Apr 1992: Original version, derived from m3inpt.
*     22 Jul 1993: Updated to match MINOS 5.4 Dec 1992.
*     17 Feb 1995: Fixed bug.  iobj is now output correctly.
*     16 May 1997: nint, hint(*) now output if 'INTORG' and 'INTEND'
*                  markers are found in the COLUMNS section.
*     21 Oct 1997: nint, hint(*) recognize bound types BV, LI, UI.
*                  m3pint may be called to pack hint(*).
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m1word/ nwordr,nwordi,nwordh
      common    /m2len / mrows,mcols,melms
      common    /m3mps3/ aijtol,bstruc(2),mlst,mer,
     $                   aijmin,aijmax,na0,line,ier(20)
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax
****  Note: iobjx is used in the next common block to avoid a clash
****        with the output parameter iobj (which returns the value
****        that m3mpsa gives to iobjx).
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobjx,jobj,kobj
      common    /m5log1/ idebug,ierr,lprint

      integer            ncard(6)
      character*4        key
      character*5        f1
      data               f1 /'(2a4)'/

*     Copy maximum dimensions into the global variables.

      mrows  = maxm
      mcols  = maxn
      melms  = maxne

*     The SPECS file has been read.
*     Make sure the MPS file is read from unit imps, and
*     make sure m3dflt knows if the problem is nonlinear.
*     The parameter 0 stops these two lines from being printed.

      nn     = max( nnjac, nnobj )
      call miopti( 'MPS file           ', imps  , 0, isumm, inform )
      call miopti( 'Nonlinear variables', nn    , 0, isumm, inform )

*     Set unspecified options to default values but don't print them
*     since minoss gets a chance later.

      iprnt  = iprint
      iprint = 0
      call m3dflt( 2 )
      iprint = iprnt

*     Load the specified obj, rhs, range and bound names into names(*).

      read (names(2), f1) mobj
      read (names(3), f1) mrhs
      read (names(4), f1) mrng
      read (names(5), f1) mbnd

*     key    retains the first 4 characters of the NAME, ROWS, COLUMNS
*            RHS, RANGES and BOUNDS cards.
*     ncard  counts the number of data records in each section.
*     m3getp finds a prime number for the length of the row hash table.
                                      
      ierr   = 0
      ncoll  = 0
      key    = '    '
      call iload ( 6, ncoll, ncard, 1 )
      call m3getp( maxm, lenh )

*     ------------------------------------------------------------------
*     Allocate workspace for the MPS input routines.
*     ------------------------------------------------------------------
      lhrtyp = 1
      lkb    = lhrtyp + 1 + (maxm/nwordh)
      lkeynm = lkb    + 1 + (maxm/nwordi)
      minmps = lkeynm + lenh
      if (nwcore .lt. minmps) go to 600

*     ------------------------------------------------------------------
*     Input ROWS.
*     lrow   is the location of the first rowname in name1, name2.
*            The column names will later go at the front.
*     lennm  is the initial length of name1, name2,
*            i.e. the maximum no. of names allowed for.
*     ------------------------------------------------------------------
      lrow   = maxn + 1
      lennm  = maxnb
      call m3mpsa( maxm, maxn, maxne, ncoll, m,
     $             lrow, lennm, lenh, nn, nncon, key, ncard,
     $             z(lhrtyp), name1, name2, z(lkeynm) )
      if (ierr .eq. 40) go to 400
      if (ierr .eq. 41) go to 600

*     ------------------------------------------------------------------
*     m  is now known.
*     Input COLUMNS, RHS, RANGES.
*     ------------------------------------------------------------------
      nb     = maxn + m 
      nka    = maxn + 1
      call m3mpsb( maxn, maxne, lrow, lennm, lenh, ncoll, objadd,
     $             m, n, nb, ne, nka, nint,
     $             nn, nncon, nnjac, nnobj, njac, key, ncard,
     $             hint, z(lhrtyp), name1, name2, z(lkeynm),
     $             ka, ha, a, bl, bu, z(lkb), pi )
      if (ierr .eq. 40) go to 400
      if (ierr .eq. 41) go to 600

*     ------------------------------------------------------------------
*     n  and  ne  are now known.
*     Move the row names to be contiguous with the column names.
*     Input BOUNDS.
*     ------------------------------------------------------------------
      call m3imov( lennm, lrow, m, n, name1 )
      call m3imov( lennm, lrow, m, n, name2 )
      nb     = n + m
      call m3mpsc( m, n, nb, ne, nint, ns, lennm,
     $             key, ncard, name1, name2,
     $             bl, bu, hint, hs, xn )

      if (iprint .gt. 0) then
         if (lprint .gt. 0) write(iprint, 1300) lenh, ncoll
         if (na0    .gt. 0) write(iprint, 1320) na0
         if (nncon  .gt. 0) write(iprint, 1350) njac
         if (nn     .gt. 0  .or.  ncard(6) .gt. 0)
     $                      write(iprint, 1400) ncard(6), ns
      end if

      go to 900

*     ---------------------------
*     Fatal error in MPS file.
*     ---------------------------
  400 call m1page( 2 )
      if (iprint .gt. 0) write(iprint, 1100)
      if (isumm  .gt. 0) write(isumm , 1100)
      go to 900

*     ------------------------------------
*     Not enough storage to read MPS file.
*     ------------------------------------
  600 call m1page(2)
      if (iprint .gt. 0) write(iprint, 1110)
      if (isumm  .gt. 0) write(isumm , 1110)

*     Exit.
*     Set iobj from m3mpsa.
*     Copy the various MPS names into the array names(*).

  900 inform = ierr
      iobj   = iobjx
      write(names, f1) name, mobj, mrhs, mrng, mbnd
      return

 1100 format(' EXIT -- fatal errors in the MPS file')
 1110 format(' EXIT -- not enough storage to read the MPS file')
 1300 format(/ ' Length of row-name hash table  ', i12
     $       / ' Collisions during table lookup ', i12)
 1320 format(  ' No. of rejected coefficients   ', i12)
 1350 format(  ' No. of Jacobian entries specified', i10)
 1400 format(  ' No. of INITIAL  bounds  specified', i10
     $       / ' No. of superbasics specified   ', i12)

      end ! subroutine mirmps

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine miwmps( mps, m, n, nb, ne, nname, names,
     $                   a, ha, ka, bl, bu, name1, name2 )

      implicit           double precision (a-h,o-z)
      character*8        names(5)
      integer            ha(ne)
      integer            ka(n+1), name1(nname), name2(nname)
      double precision   a(ne), bl(nb), bu(nb)

*     ------------------------------------------------------------------
*     miwmps  writes an MPS file to file number mps.
*     All parameters except mps are the same as for minoss in MINOS 5.4.
*     They are all input parameters.
*
*     25 Oct 1996: First version.  Not before time!
*     ------------------------------------------------------------------

      parameter        ( zero   = 0.0d+0,  one    =  1.0d+0  )
      parameter        ( bplus  = 1.0d+19, bminus = -1.0d+19 )

      logical            value
      character*1        rowtyp
      character*4        bndtyp
      character*29       form(6)

      data               form(1) /'(4x, 2a4, 2x, 2a4, f8.0     )'/
      data               form(2) /'(4x, 2a4, 2x, 2a4, 1p, e14.5)'/
      data               form(3) /'(4x,  a8, 2x, 2a4, f8.0     )'/
      data               form(4) /'(4x,  a8, 2x, 2a4, 1p, e14.5)'/
      data               form(5) /'(a4,  a8, 2x, 2a4, f8.0     )'/
      data               form(6) /'(a4,  a8, 2x, 2a4, 1p, e14.5)'/

      if (mps .le. 0) return

*     ------------------------------------------------------------------
*     ROWS section.
*     Note: b1 and b2 are bounds on ROWS, not slacks.
*           The objective row gets its name from name1(*), name2(*).
*           name(2) is ignored.
*     ------------------------------------------------------------------
      write(mps, '(a, 10x, a)') 'NAME', names(1)
      write(mps, '(a)'        ) 'ROWS'

      do i = 1, m
         j  =   n + i
         b1 = - bu(j)
         b2 = - bl(j)
         if (b1 .eq. b2) then
            rowtyp = 'E'
         else if (b1 .gt. bminus) then
            rowtyp = 'G'
         else if (b2 .lt. bplus ) then
            rowtyp = 'L'
         else
            rowtyp = 'N'
         end if

         call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
         write(mps, '(1x, a1, 2x, 2a4)') rowtyp, id1, id2
      end do

*     ------------------------------------------------------------------
*     COLUMNS section.
*     Note: Objective entries get their name from name1(*), name2(*).
*           name(2) is ignored.
*     ------------------------------------------------------------------
      write(mps, '(a)') 'COLUMNS'

      do j = 1, n
         call m4id  ( j, m, n, nb, nname, name1, name2, ic1, ic2 )

         do k = ka(j), ka(j+1) - 1
            i  = ha(k)
            call m4id  ( n+i, m, n, nb, nname, name1, name2, id1, id2 )

            ai    = a(k)
            kform = 2
            if (ai .eq. zero  .or.  abs(ai) .eq. one) kform = 1
            write(mps, form(kform)) ic1, ic2, id1, id2, ai
         end do
      end do

*     ------------------------------------------------------------------
*     RHS section.
*     Note: b1 and b2 are bounds on ROWS, not slacks.
*     ------------------------------------------------------------------
      write(mps, '(a)') 'RHS'

      do i = 1, m
         j   =   n + i
         b1  = - bu(j)
         b2  = - bl(j)
         bnd =   zero
         if (b1 .eq. b2) then
            bnd = b1
         else if (b1 .gt. bminus) then
            bnd = b1
         else if (b2 .lt. bplus ) then
            bnd = b2
         end if

         if (bnd .ne. zero) then
            call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
            kform = 4
            if (abs(bnd) .eq. one) kform = 3
            write(mps, form(kform)) names(3), id1, id2, bnd
         end if
      end do

*     ------------------------------------------------------------------
*     RANGES section.
*     ------------------------------------------------------------------
      write(mps, '(a)') 'RANGES'

      do i = 1, m
         j   =   n + i
         b1  = - bu(j)
         b2  = - bl(j)
         if (b1 .lt. b2  .and.  b1 .gt. bminus .and. b2 .lt. bplus) then
            rng   = b2 - b1
            call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
            kform = 4
            if (abs(rng) .eq. one) kform = 3
            write(mps, form(kform)) names(4), id1, id2, rng
         end if
      end do

*     ------------------------------------------------------------------
*     BOUNDS section.
*     ------------------------------------------------------------------
      write(mps, '(a)') 'BOUNDS'

      do j = 1, n
         b1     = bl(j)
         b2     = bu(j)
         bndtyp = '    '
         value  = .false.

*        Output lower bound, except for vanilla variables.

         if (b1 .eq. b2) then
            bndtyp = ' FX '
            value  = .true.
         else if (b1 .gt. bminus) then
            if (b1 .ne. zero) then
               bndtyp = ' LO '
               value  = .true.
            end if
         else if (b2 .lt. bplus) then
               bndtyp = ' MI '
         else
            bndtyp = ' FR '
         end if

         if (bndtyp .ne. '    ') then
            call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
            if (value) then
               kform = 6
               if (b1 .eq. zero  .or.  abs(b1) .eq. one) kform = 5
               write(mps, form(kform)) bndtyp, names(5), id1, id2, b1
            else
               write(mps, form(kform)) bndtyp, names(5), id1, id2
            end if
         end if

*        Output second bound if necessary.

         bndtyp = '    '
         value  = .false.

         if (b1 .eq. b2) then
*           do nothing
         else if (b1 .gt. bminus) then
            if (b2 .lt. bplus) then
               bndtyp = ' UP '
               value  = .true.
            end if
         else if (b2 .lt. bplus) then
            if (b2 .ne. zero) then
               bndtyp = ' UP '
               value  = .true.
            end if
         end if

         if (bndtyp .ne. '    ') then
            call m4id  ( j, m, n, nb, nname, name1, name2, id1, id2 )
            if (value) then
               kform = 6
               if (b2 .eq. zero  .or.  abs(b2) .eq. one) kform = 5
               write(mps, form(kform)) bndtyp, names(5), id1, id2, b2
            else
               write(mps, form(kform)) bndtyp, names(5), id1, id2
            end if
         end if
      end do

      write(mps, '(a)') 'ENDATA'

      end ! subroutine miwmps

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine misolf( m, n, nb, j, jkey, jstate,
     &                   hs, bl, bu, rc, xn )

      implicit           none
      integer            m, n, nb, j, jkey, jstate
      integer            hs(nb)
      double precision   bl(nb), bu(nb), rc(nb), xn(nb) 


*     ==================================================================
*     misolf sets the solution flags for the j-th variable:
*                      ' ' A  D  I  N    and   LL  UL SBS  BS  EQ  FR
*     by returning jkey=0  1  2  3  4,  jstate= 0   1   2   3   4   5
*
*     misolf is called by MINOS from m4soln.
*     misolf may also be called externally (e.g. by GAMS)
*     following a normal call of minoss.
*     At this stage the solution will be UNSCALED!!
*     Hence, MINOS (via m4soln) now outputs flags for the UNSCALED soln.
*
*     Input parameters m, n, nb, hs, bl, bu, rc, xn
*     are the same as for minoss.
*
*     j      (input ) is column j if j <= n;  otherwise row i = j - n.
*     jkey   (output) is one of 0 1 2 3 4.
*     jstate (output) is one of 0 1 2 3 4 5.
*
*     14 Feb 2004: First version of misolf.
*     18 Jun 2004: If the scaled problem was infeasible
*                  (with max inf at j = jBinf1), always flag that j
*                  as infeasible in the unscaled solution.
*     21 Jun 2004: Similarly, if the scaled problem wasn't optimal,
*                  (with max dual inf at j = jDinf1), always flag that j
*                  as nonoptimal in the unscaled solution.
*     ==================================================================

      ! Global variables

      double precision   dparm
      integer                      iparm
      common    /m2parm/ dparm(30),iparm(30)

      integer            name   ,mobj   ,mrhs   ,mrng   ,mbnd   ,minmax
      common    /m3mps4/ name(2),mobj(2),mrhs(2),mrng(2),mbnd(2),minmax

      double precision   sinf,wtobj
      integer                       minimz,ninf,iobj,jobj,kobj
      common    /m5lobj/ sinf,wtobj,minimz,ninf,iobj,jobj,kobj

      integer            idebug,ierr,lprint
      common    /m5log1/ idebug,ierr,lprint

      double precision   toldj   ,tolx,tolpiv,tolrow,rowerr,xnorm
      common    /m5tols/ toldj(3),tolx,tolpiv,tolrow,rowerr,xnorm

      double precision   xtol   ,ftol   ,gtol   ,pinorm,rgnorm,tolrg
      common    /m7tols/ xtol(2),ftol(2),gtol(2),pinorm,rgnorm,tolrg

      ! Local variables

      logical            feasbl, maximz
      integer            jBinf1, jDinf1, js
      double precision   b1, b2, d1, d2, dj, djtest, tolfea, tolopt, xj


      jBinf1 = iparm(27) ! Largest bound infeasibility (  scaled)
      jDinf1 = iparm(28) ! Largest dual  infeasibility (  scaled)

      tolfea = tolx
      tolopt = toldj(3) * pinorm
      feasbl = ninf   .eq. 0
      maximz = minimz .lt. 0

      js     = hs(j)
      b1     = bl(j)
      b2     = bu(j)
      xj     = xn(j)
      dj     = rc(j)

      d1     = b1 - xj
      d2     = xj - b2

      ! Change slacks into rows.
      ! Set djtest differently for row and cols.

      if (j .gt. n) then
         if (js .le. 1) js = 1 - js
         b1     = - b2
         b2     = - bl(j)
         djtest =   dj
      else
         djtest = - dj
      end if

      if (feasbl) then
         if (maximz) djtest = - djtest
         jBinf1 = 0
      end if

      if (ierr .eq. 0) then
         jDinf1 = 0
      end if

      ! Set keys and states.

      jkey   = 0   ! blank
      jstate = js  ! 0, 1, 2, 3

      if (js .le. 1) then            ! Nonbasic variables.
         if (b1 .eq. b2) jstate = 4
         if (- d1 .gt. tolfea  .and.     - d2 .gt. tolfea) jstate = 5
         if (jstate .eq. 1 ) djtest = - djtest
         if (jstate .ge. 4 ) djtest =   abs(djtest)
         if (                     abs(djtest) .le. tolopt) jkey = 1  ! A
         if (jstate .ne. 4     .and.  djtest  .gt. tolopt) jkey = 4  ! N

      else                           ! Basic and superbasic variables.
         if (abs(d1).le. tolfea  .or. abs(d2) .le. tolfea) jkey = 2  ! D
         if (jstate .eq. 2  .and. abs(djtest) .gt. tolopt) jkey = 4  ! N
         if (    d1 .gt. tolfea  .or.     d2  .gt. tolfea) jkey = 3  ! I
         if (     j .eq. jBinf1                          ) jkey = 3  ! I
      end if

      if (j .eq. jDinf1) jkey = 4  ! N

      end ! subroutine misolf
