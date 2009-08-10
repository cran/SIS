************************************************************************
*
*     File  mi10unix  fortran
*-->  Unix version of machine-dependent routines
*     (originally in file mi10mach  fortran).
*
*     mifile   m1open   m1clos   m1cpu
*
*     Machine dependencies are marked by the following line:
*-->  Machine dependency.
*
*     mi10pc.f, mi10unix.f, mi10vms.f are the same except for
*     m1open   m1cpu
*
*     14 Jul 2000: m1init, m1envt, m1page, m1time, m1timp moved
*                  into mi11sys.f.
*     30 Jan 2004: m1cpu: Added GAMS clock gfclck().
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mifile( mode )

*     ------------------------------------------------------------------
*     mifile  is a machine-dependent routine for opening various files.
*     It calls m1open (which is also machine-dependent).
*
*     MINOS uses sequential files only
*     and does not need to read and write to the same file.
*
*     ispecs, iprint, isumm  are defined in minos1 or mispec.
*     iread                  is defined here.
*
*     iread and iprint have the following use:
*        Input  files (MPS, Old Basis, Insert, Load)
*        are rewound after being read,
*        but not if they are the same as  iread.
*        Output files (Backup, New Basis, Punch, Dump, Solution, Report)
*        are rewound after being written,
*        but not if they are the same as  iprint.
*     
*     iread  = (conceptually) the Card Reader or the Keyboard that
*              can't be rewound.  MINOS does not use this file,
*              so there is no 'open'.
*     isumm  = the SUMMARY file.  Sometimes this is the screen.
*              If so, it may not need to be opened.
*     ispecs = the SPECS file, containing one or more problem specs.
*              This file is not rewound after use, because it may
*              contain another SPECS file.
*
*     Here are all the files used by MINOS.
*     The associated Index is passed to m1open
*     and must match the list of names in m1open, if that routine
*     uses method = 1.
*
*        Unit    Index    Description       Status
*        ispecs     1     Specs     file    In
*        iprint     2     Print     file    Out
*        isumm      3     Summary   file    Out
*        imps       4     MPS       file    In
*        ioldb      5     Old Basis file    In
*        insrt      6     Insert    file    In
*        iload      7     Load      file    In
*        iback      8     Backup    file    Out
*        inewb      9     New Basis file    Out
*        ipnch     10     Punch     file    Out
*        idump     11     Dump      file    Out
*        isoln     12     Solution  file    Out
*        ireprt    13     Report    file    Out
*        iread            Not opened, but used as described above
*        iprob            Not used
*        iscr             Scratch file -- no longer used
*
*     19 Jun 1989: Modified to call m1open.
*     11 Nov 1991: mode=1 now assumes that ispecs, iprint, isumm
*                  have been set by minos1 (for MINOS)
*                             or by mispec (for minoss).
*     03 Feb 1994: Pass integer index to m1open to help
*                  cater for various file-naming conventions.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
      common    /m2file/ iback,idump,iload,imps,inewb,insrt,
     $                   ioldb,ipnch,iprob,iscr,isoln,ispecs,ireprt

      integer            iprinx, isummx
      save               iprinx, isummx

*     ------------------------------------------------------------------
*-->  Machine dependency.
*     Set iread = some input unit number that should not be rewound.
*     Set iread = 0 if this is irrelevant.
*     ------------------------------------------------------------------
      iread  = 5

      if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        Mode 1: Open the Specs, Print and Summary files.
*        ispecs         remains the same throughout the run.
*        iprint, isumm  may be altered by the SPECS file.  They may
*                       need to be opened by both mode 1 and 2.
*        ---------------------------------------------------------------
         iprinx = iprint
         isummx = isumm
         call m1open( ispecs, 1, 'IN ' )
         call m1open( iprint, 2, 'OUT' )
         call m1open( isumm , 3, 'OUT' )

      else      
*        ---------------------------------------------------------------
*        Mode 2: Open files mentioned in the SPECS file just read.
*        Input files are opened first.  Only one basis file is needed.
*        ---------------------------------------------------------------
         if (imps .le. 0     ) imps = ispecs
         if (imps .ne. ispecs) call m1open( imps , 4, 'IN ' )

         if      (ioldb .gt. 0) then
            call m1open( ioldb, 5, 'IN ' )
         else if (insrt .gt. 0) then
            call m1open( insrt, 6, 'IN ' )
         else if (iload .gt. 0) then
            call m1open( iload, 7, 'IN ' )
         end if

         call m1open( iback ,  8, 'OUT' )
         call m1open( inewb ,  9, 'OUT' )
         call m1open( ipnch , 10, 'OUT' )
         call m1open( idump , 11, 'OUT' )
         call m1open( isoln , 12, 'OUT' )
         call m1open( ireprt, 13, 'OUT' )

*        Open new Print or Summary files if they were altered
*        by the Specs file.

         if (iprint .ne. iprinx) call m1open( iprint, 2, 'OUT' )
         if (isumm  .ne. isummx) call m1open( isumm , 3, 'OUT' )
      end if

*     Check that output files are different from Specs or MPS.

      if (iprint .gt. 0) then
         if (ispecs .gt. 0) then
            if (iback  .eq. ispecs) write(iprint, 1000) 'Backup'
            if (inewb  .eq. ispecs) write(iprint, 1000) 'New Basis'
            if (ipnch  .eq. ispecs) write(iprint, 1000) 'Punch'
            if (idump  .eq. ispecs) write(iprint, 1000) 'Dump'
            if (isoln  .eq. ispecs) write(iprint, 1000) 'Solution'
            if (ireprt .eq. ispecs) write(iprint, 1000) 'Report'
         end if

         if (imps .gt. 0) then
            if (iback  .eq. imps  ) write(iprint, 2000) 'Backup'
            if (inewb  .eq. imps  ) write(iprint, 2000) 'New Basis'
            if (ipnch  .eq. imps  ) write(iprint, 2000) 'Punch'
            if (idump  .eq. imps  ) write(iprint, 2000) 'Dump'
            if (isoln  .eq. imps  ) write(iprint, 2000) 'Solution'
            if (ireprt .eq. imps  ) write(iprint, 2000) 'Report'
         end if
      end if

      return

 1000 format(/ ' XXX  Warning:',
     $   ' the Specs file and ', a, ' file are on the same unit')
 2000 format(/ ' XXX  Warning:',
     $   ' the  MPS  file and ', a, ' file are on the same unit')

      end ! subroutine mifile

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1open( lun, index, state )

      integer            lun, index
      character*3        state

*     ------------------------------------------------------------------
*     m1open  is a machine-dependent routine.
*     In principal it opens a file with logical unit number lun
*     and positions it at the beginning.
*
*     Input files are treated that way.
*     An input file is opened with status='OLD'
*     (and F77 will terminate with an error if the file doesn't exist).
*
*     Output files are more machine-dependent.
*     With status='NEW', F77 would terminate if the file DID exist.
*     With status='UNKNOWN', existing files would be overwritten.
*     This is normal with Unix, but on systems that have file
*     version numbers (e.g. DEC OpenVMS), a new version is preferable.
*     It is then better not to open the file at all, but let a new
*     version be created by the first "write".
*
*     Nothing happens if
*     1. lun <= 0.  Saves us from testing lun before calling m1open.
*
*     2. Unit lun is already open.  Helps applications that call
*        minoss -- they can open files themselves if they want to.
*
*     3. lun = screen, where  screen  is a local machine-dependent
*        variable, typically 6.  It seems inadvisable to do an OPEN
*        on a file that is predefined to be an interactive screen.
*        With Unix on DECstations, open(6, file='minos.sum') sends
*        output to file minos.sum, not to the screen.
*     
*
*     lun     (input) is the unit number.
*
*     index   (input) points to one of the hardwired names below.
*             Used only if method = 1.
*
*     state   (input) is 'IN ' or 'OUT', indicating whether the file
*             is to be input or output.
*
*     15 Jul 1989: First version, follows some of the advice offered
*                  by David Gay, AT&T.
*     -- --- 1990: Added parameter "state".
*     03 Feb 1994: Added parameter "index" to help when method = 1.
*                  Local variable "method" must be set here to select
*                  various methods for naming files.
*                  Chris Jaensch, IFR Stuttgart, recommends not opening
*                  input files if they are already open.  This version
*                  ignores all open files (input or output), assuming
*                  that they are taken care of by the calling program.
*     13 May 1994: Local variable "screen" used to avoid opening screen.
*     ------------------------------------------------------------------

      logical            input, uopen
      integer            method, screen
      character*100      fname

*     ------------------------------------------------------------------
*-->  Machine dependency.
*     names(*) is needed if method = 1 below.
*     Make sure "character*n" sets n big enough below.
*     It doesn't matter if it is bigger than necessary, since
*     "open( lun, file=name )" allows name to have trailing blanks.
*     ------------------------------------------------------------------

      character*9        names(13)
      data   names( 1) /'MINOS.SPC'/
      data   names( 2) /'MINOS.PRN'/
      data   names( 3) /'MINOS.SUM'/
      data   names( 4) /'MINOS.MPS'/
      data   names( 5) /'MINOS.OLB'/
      data   names( 6) /'MINOS.INS'/
      data   names( 7) /'MINOS.LOD'/
      data   names( 8) /'MINOS.BAK'/
      data   names( 9) /'MINOS.NWB'/
      data   names(10) /'MINOS.PUN'/
      data   names(11) /'MINOS.DMP'/
      data   names(12) /'MINOS.SOL'/
      data   names(13) /'MINOS.RPT'/

*     ------------------------------------------------------------------
*-->  Machine-dependency.
*     Set "method" to suit your operating system.  See comments below.
*     Set "screen" to a unit number that never needs to be opened.
*     (Typically, screen = 6.  If unknown, set screen = 0.)
*     ------------------------------------------------------------------
      method = 2
      screen = 6

*     ------------------------------------------------------------------
*     Quit if lun<=0 or lun = iscreen or unit lun is already open.
*     ------------------------------------------------------------------
      if (lun .le. 0     ) go to 900
      if (lun .eq. screen) go to 900
      inquire( lun, opened=uopen )
      if (     uopen     ) go to 900

*     ------------------------------------------------------------------
*     Open file lun by a specified machine-dependent method.
*     ------------------------------------------------------------------
      input  = state .eq. 'IN '  .or.  state .eq. 'in '

      if (method .eq. 1) then
*        ---------------------------------------------------------------
*        Hardwired filename.
*        We use "index" to get it from names(*) above.
*        Typical machines: IBM PC under DOS.
*        ---------------------------------------------------------------
         if ( input ) then
            open( lun, file=names(index), status='OLD' )
         else
            open( lun, file=names(index), status='UNKNOWN' )
         end if

      else if (method .eq. 2) then
*        ---------------------------------------------------------------
*        External filename.
*        The operating system uses a default name
*        (e.g. fort.7 or for007)
*        or has already assigned a name to this unit
*        (e.g. via a command file or script).
*        Typical machines:  Unix with no command-line arguments,
*                           DEC OpenVMS,
*                           IBM Mainframes.
*        ---------------------------------------------------------------
         if ( input ) then
            open( lun, status='OLD' )
         else
*           Let the first "write" do it.
         end if

      else if (method .eq. 3) then
*        ---------------------------------------------------------------
*        Construct name from lun.
*        The default name used by the operating system
*        may not be quite what you wanted.
*        The code below constructs names like fort.7 and fort.15.
*        (This approach used by Chris Jaensch, IFR Stuttgart.)
*        Typical machines:  Arbitrary.
*        ---------------------------------------------------------------
         if (lun .le. 9) then
             write(fname, '(a,i1)') 'fort.', lun
         else
             write(fname, '(a,i2)') 'fort.', lun
         endif

         if ( input ) then
            open( lun, file=fname, status='OLD' )
         else
            open( lun, file=fname, status='UNKNOWN' )
         end if

      else if (method .eq. 4) then
*        ---------------------------------------------------------------
*        Lookup name.
*        Assume some routine "gfname" will provide a name at run-time.
*        (This approach is used by David Gay, AT&T.)
*        Typical machines:  Unix with command-line arguments.
*        Note that 'UNKNOWN' is equivalent to trying first with 'OLD',
*        and then with 'NEW' if the file doesn't exist.
*        ---------------------------------------------------------------
*-       call gfname( lun, fname )
         if ( input ) then
            open( lun, file=fname, form='formatted', status='OLD' )
         else
            open( lun, file=fname, form='formatted', status='UNKNOWN' )
         end if
      end if

*     ------------------------------------------------------------------
*     Rewind input files.
*     (Some systems position existing files at the end
*     rather than the beginning.)
*     err=900 covers files that have not yet been opened.
*     ------------------------------------------------------------------
      if ( input ) rewind( lun, err=900 )

  900 return
      
      end ! subroutine m1open

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1clos( lun )

*     ------------------------------------------------------------------
*     m1clos  closes the file with logical unit number lun.
*     This version is trivial and so far is not even used by MINOS.
*     Perhaps some implementations will need something fancier.
*     ------------------------------------------------------------------

      close( lun )

      end ! subroutine m1clos

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1cpu ( mode, time )

*-->  WinNT with DEC (Compaq) F90
*-->  USE DFPORT, ONLY: RTC  ! Must be before the other declarations

      integer            mode
      double precision   time

*     ------------------------------------------------------------------
*     m1cpu is a machine-dependent routine to return time = cpu time
*     in seconds, so that 2 consecutive calls will indicate the
*     time difference of operations between the 2 calls.
*     The parameter 'mode' indicates what function should be done
*     to the timer.  This allows necessary initialization for certain
*     machines.
*     mode =  1  indicates initialization,
*     mode =  0  indicates normal use,
*     mode = -1  indicates stop the timer.
*
*     1988:  Used Irv Lustig's approach...
*     On DEC VAX/VMS systems we need to call the correct library
*     routine to get the timer statistics.  These statistics are
*     found by using the times() function in the VAX C Runtime library.
*     To use this version of m1cpu, one must create an options file
*     called  vmsc.opt  with the line
*        SYS$LIBRARY:VAXCRTL/SHARE
*     in it.   Then link using the usual command and append ,vmsc/opt
*     to the end of the line.  The name vmsc can be anything.
*     
*     02 Apr 1993: Went back to VMS Fortran routines to avoid linking
*                  to the C library.  (On DEC AXP, the C runtime lib
*                  appears to be translated from the VAX executable,
*                  and therefore requires linking with /NONATIVE,
*                  which possibly adds a small overhead to all
*                  subroutine calls.
*     01 Mar 1998: Added timers for PC Lahey, IBM RS/6000 AIX, NeXT.
*     12 Jul 2000: Added timer for DEC (Compaq) F90 on WinNT.
*     30 Jan 2004: (At GAMS) Added GAMS clock gfclck().
*                  1. Uncomment two separate lines:
*                        external gfclck
*                           ...
*                        time   = gfclck()
*                  2. Following Unix (Sun, SGI, etc), comment out
*                        time   = etime ( tarray )
*     ------------------------------------------------------------------

*-->  Machine dependency.
*-->  DEC OpenVMS with Fortran runtime library.
*     external        lib$init_timer, lib$stat_timer, lib$free_timer
*     integer         itimad, istatu, idata
*     save            itimad

*-->  DEC VAX/VMS with C runtime library.
*     integer            itimad(4)

*-->  PC Lahey Fortran
*     integer            itimad(4)

*-->  WinNT with DEC (Compaq) F90
*-->  USE DFPORT, ONLY: RTC
*-->  REAL*8             dectime, decinit
*-->  SAVE               decinit

*-->  Unix (Sun, SGI, etc)
      real               tarray(2)

*-->  IBM RS/6000 method 1 (from Tim Barth, NASA Ames, Jan 1995)
*     REAL(4) etime_
*     TYPE TB_TYPE
*       SEQUENCE
*       REAL(4) USRTIME
*       REAL(4) SYSTIME
*     END TYPE
*     TYPE (TB_TYPE) ETIME_STRUCT

*-->  IBM RS/6000 AIX method 2
*     integer            mtime

*-->  NeXTstation M68040, Unix  (routine in ftime.c)
*     real               rtime

*-->  GAMS clock
*     external           gfclck
      double precision   gfclck

*     ------------------------------------------------------------------
*     Executable code.
*     ------------------------------------------------------------------

      if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        Initialize.
*        ---------------------------------------------------------------
         time   = 0.0d+0

*-->     DEC OpenVMS with Fortran library.
*        istatu = lib$init_timer( itimad )
*        if (.not. istatu) call lib$signal( %val(istatu) )

*-->     WinNT with DEC (Compaq) F90
*-->     decinit = RTC()

      else if (mode .eq. 0) then
*        ---------------------------------------------------------------
*        Normal call.
*        Return current timer value here.
*        ---------------------------------------------------------------

*-->     DEC OpenVMS with Fortran library.
*-->     istatu returns the number of  centiseconds.
*        istatu = lib$stat_timer( 2, idata, itimad )
*        if (.not. istatu) call lib$signal( %val(istatu) )
*        time   = dble( idata ) * 0.01d+0

*-->     DEC VAX/VMS with C library.
*-->     itimad(1) returns the number of  centiseconds.
*        call times ( itimad )
*        time   = dble( itimad(1) ) * 0.01d+0

*-->     PC Lahey Fortran.
*-->     itimad(1) returns the number of  centiseconds.
*        call timer ( itimad )
*        time   = dble( itimad(1) ) * 0.01d+0

*-->     WinNT with DEC (Compaq) F90
*-->     Time in secs since 00:00:00 GMT Jan 1st, 1970.
*-->     Bias must be subtracted to give adequate precision for real*4
*-->     dectime = RTC() - decinit
*-->     time    = real(dectime)

*-->     Unix (Sun, SGI, etc)
*-->     etime returns seconds.
*-->     time   = etime ( tarray )
         time   = etime ( tarray )

*-->     IBM RS/6000 method 1 (from Tim Barth, NASA Ames, Jan 1995)
*-->     Underscore must be added.
*        time   = etime_( tarray )

*-->     IBM RS/6000 AIX method 2
*-->     mclock returns hundredths of a second
*        time   = dble( mclock( ) ) * 0.01d+0

*-->     NeXTstation M68040, Unix  (routine in ftime.c)
*        call ftime ( rtime )
*        time   = dble( rtime )

*-->     GAMS clock
*        time   = gfclck()

*-->     On other machines, to forget about timing, just say
*        time   = -1.0d+0

      else if (mode .eq. -1) then
*        ---------------------------------------------------------------
*        Stop the clock.
*        ---------------------------------------------------------------
         time   = 0.0d+0

*-->     DEC OpenVMS with Fortran library.
*        istatu = lib$free_timer( itimad )
*        if (.not. istatu) call lib$signal( %val(istatu) )
      end if

      end ! subroutine m1cpu
