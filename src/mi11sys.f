************************************************************************
*
*     File  mi11sys   fortran
*
*     m1init   m1envt   m1page   m1time   m1timp
*
*     14 Jul 2000: mi11sys.f  collects system routines from mi10...
*                  that are the same for Unix, Linux, PC, OpenVMS, etc
*                  in the normal MINOS distribution.
*                  Some routines may be adapted for GAMS, AMPL, etc. 
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1init( )

      implicit           double precision (a-h,o-z)

*     ------------------------------------------------------------------
*     m1init defines certain machine-dependent constants.
*
*     eps    = floating-point precision (e.g., 2.0**(-47) for CDC)
*     nwordr = no. of reals      per word of  z(*)
*     nwordi = no. of integers   per word of  z(*)
*     nwordh = no. of integer*2s per word of  z(*)   IF ANY,
*            = no. of integer*4s per word of  z(*)   IF NO integer*2s
*        where  z(*)  is the main array of storage.
*
*     Original version:  integer*2  and  nwordh = 4  used throughout for
*                        certain arrays:
*                        ha, hs  in MINOS (and maybe a few others),
*                        indc, indr, ip, iq  in LUSOL.
*                        A quirk in LUSOL limits MINOS to 16383 rows
*                        when the limit should have been  32767.
*
*                        At present, nwordr is not used because there
*                        are no real*4 arrays.
*
*     22 May 1992:       integer*4  and  nwordh = 2  now used to allow
*                        essentially any number of rows.
*     13 Aug 2000: title moved into mititle.
*     ------------------------------------------------------------------

      common    /m1eps / eps,eps0,eps1,eps2,eps3,eps4,eps5,plinfy
      common    /m1word/ nwordr,nwordi,nwordh



*---+ IEEE standard: eps = 2**(-52) = 2.22e-16
      eps    = 2.0d+0**(-52)
      nwordr = 2
      nwordi = 2
      nwordh = 2

*     Set other machine-precision constants.

      eps0   = eps**(0.8 d+0)
      eps1   = eps**(0.67d+0)
      eps2   = eps**(0.5 d+0)
      eps3   = eps**(0.33d+0)
      eps4   = eps**(0.25d+0)
      eps5   = eps**(0.2 d+0)
      plinfy = 1.0d+20

*     Set the environment (for later use).

      call m1envt( 0 )

*     end of m1init
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1envt( mode )

*     ------------------------------------------------------------------
*     m1envt specifies the environment within which MINOS is being used.
*
*     When mode = 0, the various logicals should be initialized.
*     page1 says whether new pages are ever wanted on file iprint.
*     page2 says whether new pages are ever wanted on file isumm.
*
*     When mode is in the range 1 to 99, each environment does its
*     own thing.
*
*     When mode = 999, MINOS is asking if resource limits have been
*     reached.   To indicate YES, set ierr = 19.
*
*     The various environments are as follows:
*
*     ALONE:
*     This means MINOS is in stand-alone mode---the normal case.
*     Nothing special is done.
*
*     GAMS:
*     When mode = 1, 2, ..., 9 the characters =1, =2, ..., =9
*     are output to the print file.
*     When mode = 999, the resource limits are tested.
*
*     MINT:
*     Since branch-and-bound means a lot of CYCLES, we suppress
*     page ejects.  Otherwise, nothing special as yet.
*
*     AMPL:
*     Nothing special yet, but might want to test resource limits.
*
*     16 Sep 1987  Initial version.
*     24 Apr 1992  AMPL added.
*     ------------------------------------------------------------------

      logical            alone, AMPL, GAMS, MINT, page1, page2
      common    /m1env / alone, AMPL, GAMS, MINT, page1, page2
      common    /m1file/ iread,iprint,isumm
      common    /m5log1/ idebug,ierr,lprint

*     GAMS resource info:
*     ARESLM is limit on cpu time.
*     ATIME0 is starting time.
*     ATMGET is a function returning time so far.
*
*X    REAL               ARESLM, ATIME0, ADELT
*X    COMMON    /GAMS00/ ARESLM, ATIME0, ADELT
*X    REAL               ATMGET
*X    EXTERNAL           ATMGET

      if (mode .le. 0) then
*        ---------------------------------------------------------------
*        mode = 0.    Initialize.
*        page1 and page2 should be false for applications involving
*        many Cycles  (e.g. MINT).
*        ---------------------------------------------------------------
         alone  = .true.
         AMPL   = .false.
         GAMS   = .false.
         MINT   = .false.
         page1  = .true.
         page2  = .false.

      else if (mode .lt. 999) then
*        ---------------------------------------------------------------
*        mode = 1 or more.  Do what has to be done in each environment.
*        ---------------------------------------------------------------
         if (GAMS  .and.  iprint .gt. 0) then
            write(iprint, '(a1, i1)') '=', mode
         end if

      else if (mode .eq. 999) then
*        ---------------------------------------------------------------
*        mode = 999.  Test for excess time, etc.
*        ---------------------------------------------------------------
         if (GAMS) then
*X          IF (ATMGET() - ATIME0 .GE. ARESLM) IERR = 19
         else if (AMPL) then
*X
         end if
      end if

*     end of m1envt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1page( mode )

*     ------------------------------------------------------------------
*     m1page is an installation-dependent routine.  It is called at
*     points where some users might want output to files iprint or isumm
*     to begin on a new page.
*
*     page1 and page2 have already been set by m1envt.
*     If they are true, a page eject and a blank line are output.
*     Otherwise, just a blank line is output.
*
*     If mode = 1, just the page control is relevant.
*     If mode = 2, GAMS wants m1envt to print an =.
*     If mode = 0  and Summary level = 0, we don't want anything output
*                  to the Summary file.  At present, this is so m8setj
*                  will print just one line per major iteration, with
*                  no blank line in between.
*
*     16-Sep-1987:  First version.
*     20-Mar-1988:  mode 2 added.
*     12-Dec-1991:  mode 0 added.
*     ------------------------------------------------------------------

      logical            alone, AMPL, GAMS, MINT, page1, page2
      common    /m1env / alone, AMPL, GAMS, MINT, page1, page2
      common    /m1file/ iread,iprint,isumm
      logical            prnt0 ,prnt1 ,summ0 ,summ1 ,newhed
      common    /m5log4/ prnt0 ,prnt1 ,summ0 ,summ1 ,newhed

      if (iprint .gt. 0) then
         if ( page1 ) write(iprint, 1001)
                      write(iprint, 1002)
      end if

      if (mode   .eq. 2) call m1envt( 1 )

      if (isumm  .gt. 0) then
         if ( page2 ) write(isumm , 1001)
         if ( summ1  .or.  mode .ne. 0 )
     $                write(isumm , 1002)
      end if
      return

 1001 format('1')
 1002 format(' ')

*     end of m1page
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1time( clock, prtopt )

      implicit           double precision (a-h,o-z)
      integer            clock, prtopt

*     ------------------------------------------------------------------
*     m1time, m1timp and m1cpu are derived from timer, timout and nowcpu
*     written for DEC VAX/VMS systems by Irvin Lustig,
*     Department of Operations Research, Stanford University, 1987.
*
*     MINOS  calls m1time only.  m1time calls m1cpu  and  m1timp.
*     Only m1cpu is intrinsically machine dependent.
*
*     If a timer is available, call it in m1cpu  and arrange that
*     m1cpu  returns the current CPU time in seconds.
*
*     If a timer is not available, set time = -1.0d+0 in m1cpu.
*     Timing will be turned off and m1timp will not be called.
*     ------------------------------------------------------------------
*
*     m1time turns on or off a selected clock and optionally prints
*     statistics regarding all clocks or just the clock chosen.
*
*     The value of abs(clock) is which clock to use.
*     If clock = 0 and prtopt = 0, all clocks and statistics are reset.
*     If clock > 0, the clock is reset to start timing at the
*                   current time (determined by calling the
*                   machine-dependent subroutine m1cpu).
*     If clock < 0, the clock is turned off and the statistic is
*                   recorded for the amount of time since the clock
*                   was turned on.
*
*     prtopt is the print option.
*     If ltime < 0, nothing is printed.  Otherwise,
*     prtopt = 0 indicates print nothing,
*            = 1 indicates print last time for this clock,
*                only if clock < 0 (it has just been turned off),
*            = 2 indicates print total time for all clocks,
*            = 3 indicates print mean  time for all clocks.
*
*     The procedure for adding a new timer n is as follows:
*     1)  Change ntime to n in the parameter statement below (and in
*         all other routines referencing common block /m1tim /).
*     2)  Expand the array "label" to length n in subroutine m1timp.
*
*     04 Jun 1989: Irv's VMS/VAXC version of m1cpu installed,
*                  with changes to return time in seconds.
*     10 Jul 1992: More clocks added for use in AMPL (and elsewhere).
*     ------------------------------------------------------------------
*
*        Clock 1 is for input time.
*        Clock 2 is for solve time.
*        Clock 3 is for output time.
*        Clock 4 is for the nonlinear constraint functions.
*        Clock 5 is for the nonlinear objective.
*
*        numt(i)  is the number of times clock i has been turned on.
*        tlast(i) is the time at which clock i was last turned on.
*        tsum(i)  is the total time elapsed while clock i was on.
*        ltime    is the Timing level set in the Specs file.

      parameter        ( ntime = 5 )
      common    /m1tim / tlast(ntime), tsum(ntime), numt(ntime), ltime

      double precision   stat, time
      integer            iclock, ilo, ihi

      if (ltime .eq. 0) return
      iclock = iabs(clock)

      if (clock .eq. 0) then
         if (prtopt .eq. 0) then

*           clock = 0, prtopt = 0.  Reset everything.

            call m1cpu ( 1, time )
            call m1cpu ( 0, time )
            do 100 i = 1, ntime
               tlast(i) = time
               tsum(i)  = 0.0d+0
               numt(i)  = 0
  100       continue

*           If the m1cpu( 0, time ) gave time < 0.0, we assume that
*           the clock is a dummy.  Turn off future timing.

            if (time .lt. 0.0d+0) ltime = 0
         end if

      else
         call m1cpu ( 0, time )
         if (clock .gt. 0) then
            tlast(iclock) = time
         else
            stat         = time - tlast(iclock)
            tsum(iclock) = tsum(iclock) + stat
            numt(iclock) = numt(iclock) + 1
         end if
      end if

*     Now deal with print options.

      if (prtopt .eq. 0  .or.  ltime .lt. 0) then

*        Do nothing.

      else if (prtopt .eq. 1) then

*        Print statistic for last clock if just turned off.
      
         if (clock .lt. 0) then
            call m1timp( iclock, 'Last time', stat )
         end if

      else

*        prtopt >= 2.  Print all statistics if clock = 0,
*        or print statistic for individual clock.
      
         if (clock .eq. 0) then
            call m1cpu ( -1, time )
            ilo   = 1
            ihi   = ntime
         else
            ilo   = iclock
            ihi   = iclock
         end if
      
         do 400 i = ilo, ihi
            stat  = tsum(i)
            if (prtopt .eq. 2) then
               call m1timp( i, 'Time', stat )
            else if (prtopt .eq. 3) then
               istat = numt(i)
               if (istat .gt. 0) stat = stat / istat
               call m1timp( i, 'Mean time', stat )
            end if
  400    continue
      end if

*     end of m1time
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine m1timp( iclock, lstat, stat )

      integer            iclock
      character*(*)      lstat
      double precision   stat

*     ------------------------------------------------------------------
*     m1timp  prints CPU time for m1time on file iprint and/or isumm.
*     It is not intrinsically machine dependent.
*
*     iclock  selects the correct label.
*     lstat   is a string to print to tell which type of statistic.
*     stat    is the statistic to print out.
*
*     12 Jul 1992: Array of labels avoids multiple formats.
*     ------------------------------------------------------------------

      common    /m1file/ iread,iprint,isumm
*     common    /m7len / fobj  ,fobj2 ,nnobj ,nnobj0
*     common    /m8len / njac  ,nncon ,nncon0,nnjac

      character*24       label(5)
      data               label
     $                 / 'for MPS input',          
     $                   'for solving problem',    
     $                   'for solution output',    
     $                   'for constraint functions',
     $                   'for objective function' /

      if (iclock .eq. 1) then
         if (iprint .gt. 0) write(iprint, 1000)
         if (isumm  .gt. 0) write(isumm , 1000)
*      nncon and nnobj below are always 0 for some reason.
*      else
*         if (iclock .eq. 4  .and.  nncon .eq. 0) return
*         if (iclock .eq. 5  .and.  nnobj .eq. 0) return
      end if

      if (iprint .gt. 0) write(iprint, 1000) lstat, label(iclock), stat
      if (isumm  .gt. 0) write(isumm , 1000) lstat, label(iclock), stat
      return

 1000 format( 1x, a, 1x, a, t38, f13.2,' seconds')

*     end of m1timp
      end
