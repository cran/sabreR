!***********************************************************************
!*                                                                     *
!*                      MACHINE SPECIFIC ROUTINES                      *
!*                                                                     *
!***********************************************************************

      PROGRAM SABRE
!     Version: 23.10.07, 15:36
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'mpif.h'
!-----------------------------------------------------------------------
      include 'mpi_info.h'
!-----------------------------------------------------------------------
      include 'limitsf90.h'
      include 'chans.h'
      character start_date*8,finish_date*8,current_date*8,start_time*10, &
                finish_time*10,current_time*10,outbuf*80
      character(len=50), allocatable :: mnames(:), name(:)
      character(len=80) :: nextarg, arg_value, input_file
      double precision, allocatable :: x(:), beta(:)
      integer istat, iarg, nargs, iargc, length, ch, iostat, ivalue, &
              equals_position
      integer, allocatable :: ilev(:), ipos(:)
      integer clock1,clock2,hrs,mins,secs,timed
      logical stop, error
!
!    Assign default data space
!
      mxx = 100000000
      maxy = 500000
      maxvar = 200
      maxpar = 200
!
!    See if there are any program arguments
!
      terse = .false.
      iarg = 0
  1   continue
!
!       Look for argument prefix
!
         iarg = iarg + 1
         if( iarg > iargc() ) go to 2

!--------***********
         CALL GETARG( iarg, nextarg )
!--------***********

!       Convert the argument prefix to lower case
!       nextarg might actually be the prefix followed immediately by an
!       = sign and then the value with no separating spaces in which case only
!       convert the first part (up to the = sign)
!
         equals_position = index( nextarg, '=' )
         if( equals_position > 0 ) then
            CALL LCASE(nextarg(1:equals_position-1))
         else
            CALL LCASE(nextarg)
         end if

!
!       First look for version number argument.  If found, print and stop
!
         if( nextarg(1:7) == 'version' .or. nextarg(1:8) == '-version' .or. &
             nextarg(1:8) == '/version' ) then
!          Need to call begin to set the output unit number

!-----------**********
            CALL BEGIN
!-----------**********

            call wrtscr( 'Sabre version 5.0' )
            stop
         else if( nextarg(1:5) == 'terse' .or. nextarg(1:6) == '-terse' .or. &
             nextarg(1:6) == '/terse' ) then
            terse = .true.
         else if( nextarg(1:5) == 'input' .or. nextarg(1:6) == '-input' .or. &
                  nextarg(1:6) == '/input' ) then
!
!          Try to get the input file name which should follow
!
!-----------******************
            CALL GET_ARG_VALUE( iarg, nextarg, input_file )
!-----------******************
!
!          Attempt to open an input file with specified name
!
!          Need to call begin to set the input unit number
!
!-----------**********
            CALL BEGIN
!-----------**********

            open( unit=inch, file=input_file, status='old', iostat=iostat )
            if( iostat /= 0 ) then
               call wrtscr( 'Error in attempting to open input file:' )
               call wrtscr( input_file(1:length(input_file)) )
               stop
            end if     
         else if( nextarg(1:3) == 'dat' .or. nextarg(1:4) == '-dat' .or. &
             nextarg(1:4) == '/dat' ) then

!-----------******************
            CALL GET_ARG_VALUE( iarg, nextarg, arg_value )
!-----------******************

!-----------***********
            CALL CHAINT( arg_value, ivalue, length(arg_value), error )
!-----------***********

            if( error ) then
               call wrtscr( 'Invalid value given for DAT' )
               stop
            else
               if( ivalue <= 0 ) then
                  call wrtscr( 'DAT value must be greater than zero' )
                  stop
               else
                  mxx = ivalue
               end if
            end if
         else if( nextarg(1:3) == 'obs' .or. nextarg(1:4) == '-obs' .or. &
             nextarg(1:4) == '/obs' ) then

!-----------******************
            CALL GET_ARG_VALUE( iarg, nextarg, arg_value )
!-----------******************

!-----------***********
            CALL CHAINT( arg_value, ivalue, length(arg_value), error )
!-----------***********

            if( error ) then
               call wrtscr( 'Invalid value given for OBS' )
               stop
            else
               if( ivalue <= 0 ) then
                  call wrtscr( 'OBS value must be greater than zero' )
                  stop
               else
                  maxy = ivalue
               end if
            end if
         else if( nextarg(1:3) == 'var' .or. nextarg(1:4) == '-var' .or. &
             nextarg(1:4) == '/var' ) then

!-----------******************
            CALL GET_ARG_VALUE( iarg, nextarg, arg_value )
!-----------******************

!-----------***********
            CALL CHAINT( arg_value, ivalue, length(arg_value), error )
!-----------***********

            if( error ) then
               call wrtscr( 'Invalid value given for VAR' )
               stop
            else
               if( ivalue <= 0 ) then
                  call wrtscr( 'VAR value must be greater than zero' )
                  stop
               else
                  maxvar = ivalue
               end if
            end if
         else if( nextarg(1:3) == 'par' .or. nextarg(1:4) == '-par' .or. &
             nextarg(1:4) == '/par' ) then

!-----------******************
            CALL GET_ARG_VALUE( iarg, nextarg, arg_value )
!-----------******************

!-----------***********
            CALL CHAINT( arg_value, ivalue, length(arg_value), error )
!-----------***********

            if( error ) then
               call wrtscr( 'Invalid value given for PAR' )
               stop
            else
               if( ivalue <= 0 ) then
                  call wrtscr( 'PAR value must be greater than zero' )
                  stop
               else
                  maxpar= ivalue
               end if
            end if
         end if
!
      go to 1
  2   continue
!
!    Allocate array space
!
      allocate( mnames(MAXVAR), stat=istat )
      allocate( name(MAXVAR), stat=istat )
      allocate( x(MXX), stat=istat )
      if( istat /= 0 ) then
         print *, 'Failed to allocate space for X'
      end if
      allocate( beta(MAXPAR*(MAXPAR+5)), stat=istat )
      allocate( ilev(MAXVAR), stat=istat )
      allocate( ipos(MAXVAR), stat=istat )
!
!-----------------------------------------------------------------------
!---- Set MPI variables used in main Sabre code to values appropriate
!---- for the serial version
      num_processors = 1
      this_processor = 0

!-----*****************
      CALL SYSTEM_CLOCK(clock1)
!-----*****************

!---- Record start clock time
!-----******************
      CALL DATE_AND_TIME(date=start_date,time=start_time)
!-----******************

  3   continue

!-----***********
      CALL SABREM(x,beta,name,mnames,ilev,ipos,stop)
!-----***********

      if (.not. stop) then
          go to 3
      end if

!---- Display start and finish clock time
!-----******************
      CALL DATE_AND_TIME(date=finish_date,time=finish_time)
!-----******************

!-----*****************
      CALL SYSTEM_CLOCK(clock2)
!-----*****************

      write (outbuf,'(5(a,a2),a,a4)') 'Start:  ',start_time(1:2),':', &
      start_time(3:4),':',start_time(5:6),' on ',start_date(7:8),'/', &
      start_date(5:6),'/',start_date(1:4)
      if( .not. terse ) call wrtscr(outbuf)
      if (qflag) then
          call wrtflq(outbuf)
      end if

      write (outbuf,'(5(a,a2),a,a4)') 'Finish: ',finish_time(1:2),':', &
      finish_time(3:4),':',finish_time(5:6),' on ',finish_date(7:8),'/', &
      finish_date(5:6),'/',finish_date(1:4)
      if( .not. terse ) call wrtscr(outbuf)
      if (qflag) then
          call wrtflq(outbuf)
      end if

      timed = nint((clock2-clock1)/1000.0)
      hrs = timed/3600
      mins = (timed - hrs*3600)/60
      secs = timed - hrs*3600 - mins*60
      write (outbuf,'(a,i3,a,2(i2,a))') 'Time = ',hrs,'hr ',mins,'min ', &
      secs,'sec'
      if( .not. terse ) call wrtscr(outbuf)
      if (qflag) then
          call wrtflq(outbuf)
          call newlnq
      end if

      write (outbuf,'(i1,a)') num_processors,' processor used'
      if( .not. terse ) call wrtscr(outbuf)
      if (qflag) then
          call wrtflq(outbuf)
      end if

!---- close the log file
!-----***********
      CALL CLSFIL(ooutch)
!-----***********

!---- close the trace file
      if (tflag) then

!---------***********
          CALL CLSFIL(toutch)
!---------***********

      end if

!---- close the resid file
      if (rflag) then

!---------***********
          CALL CLSFIL(routch)
!---------***********

      end if

!---- close the time file
      if (qflag) then

!---------***********
          CALL CLSFIL(qoutch)
!---------***********

      end if

      stop
      end

!***********************************************************************

      SUBROUTINE BEGIN
!-----------------------------------------------------------------------
!     These are the sizes of various arrays, and other system dependent
!     structures for various things. Look at this subroutine carefully
!     to make sure that the settings are correct.
!-----------------------------------------------------------------------
!     Common block CHANS
!-----------------------------------------------------------------------
!     OUTCH  : output channel number
!     INCH   : terminal input channel number
!     CINCH  : current input channel
!     INPCH  : channel used by the input command
!     RINCH  : channel used by the read command
!     OOUTCH : file output channel number
!     tOUTCH : trace file output channel number
!     rOUTCH : resid file output channel number
!     qOUTCH : time file output channel number
!     tFLAG  : trace file output flag (value set in main program)
!     rFLAG  : resid file output flag (value set in main program)
!     qFLAG  : time file output flag (value set in main program)
!     terse  : If true causes reduced output to OUTCH
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     Common block ACCMAC
!-----------------------------------------------------------------------
!     Machine constants:
!     Z1     - largest double precision number possible on the machine
!     ZL1    - ln(Z1)
!     Z2     - smallest double precision number possible on the machine
!     ZL2    - ln(Z2)
!     IZ     - largest integer number possible on the machine
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
!     Initialisation of COMMON CHANS input and output channels. Note
!     that these should all be different numbers and they should also be
!     valid for use on your particular machine.
!-----------------------------------------------------------------------
      INCH  = 5
      OUTCH = 6
      RINCH = 15
      INPCH = 16
      OOUTCH = 21
      toutch = 22
      routch = 23
      qoutch = 24
      CINCH = INCH
!-----------------------------------------------------------------------
!     Initialisation of COMMON ACCMAC machine accuracy
!-----------------------------------------------------------------------
      Z1 = 1.0D+308
      Z2 = 1.0D-307
      ZL1 = DLOG(Z1)
      ZL2 = DLOG(Z2)
      IZ = 2147483647

      RETURN
      END

!***********************************************************************

      LOGICAL FUNCTION CALPHA(C)
!-----------------------------------------------------------------------
      CHARACTER C
!----------------------------------------------------------------------
!     This function is machine dependent. It should return .FALSE. if C
!     is not alphabetic and .TRUE. if C is alphabetic.
!     "alphabetic" = [a-z][A-Z][$]
!     It assumes an ASCII collating sequence, and will need to be
!     rewritten if moved to a machine using the EBCDIC character codes
!     or anything even more esoteric (BCD??)
!     BJF
!----------------------------------------------------------------------
      INTEGER ICHARC
!----------------------------------------------------------------------
      ICHARC = ICHAR(C)
      CALPHA = ((ICHARC .GE. 65 .AND. ICHARC .LE. 90) .OR. &
      (ICHARC .GE. 97 .AND. ICHARC .LE. 122) .OR. C .EQ. '$')

      RETURN
      END

!***********************************************************************

      SUBROUTINE LOWER(LINE)
!-----------------------------------------------------------------------
      CHARACTER LINE*82
!-----------------------------------------------------------------------
!     This routine
!     i) converts the input line to lower case if necessary
!     ii) removes all dollars from the command line
!     It assumes an ASCII collating sequence, and will need to be
!     rewritten if a different collating sequence for characters is used
!     by your compiler (eg EBCDIC). Alternatively, you can manually
!     switch into lower case when you use the program!
!     BJF 30/12/88
!-----------------------------------------------------------------------
      INTEGER I,J,ICH,DOLLAR
!-----------------------------------------------------------------------
      DATA DOLLAR /36/
!-----------------------------------------------------------------------
      J = 1

      DO 10 I = 1,82
          ICH = ICHAR(LINE(I:I))

          IF (ICH .EQ. DOLLAR) then
              GOTO 10
          end if

          IF (ICH .GE. 65 .AND. ICH .LE. 90) THEN
              LINE(J:J) = CHAR(ICH+32)
          ELSE
              LINE(J:J) = CHAR(ICH)
          END IF

          J = J+1
   10 CONTINUE

      DO 20 I = J,82
          LINE(I:I) = ' '
   20 continue


      RETURN
      END

!***********************************************************************

      SUBROUTINE LCASE(string)
!-----------------------------------------------------------------------
      CHARACTER string*(*)
!-----------------------------------------------------------------------
!     Version of LOWER which just converts letters to lower case, nothing else
!-----------------------------------------------------------------------
      INTEGER I,ICH
!-----------------------------------------------------------------------

      DO 10 I = 1, length(string)

          ICH = ICHAR(string(I:I))

          IF (ICH .GE. 65 .AND. ICH .LE. 90) THEN
              string(i:i) = CHAR(ICH+32)
          END IF

   10 CONTINUE

      RETURN
      END

!***********************************************************************

      SUBROUTINE SCRSET(smode)
!-----------------------------------------------------------------------
      integer smode
!-----------------------------------------------------------------------
!     This routine (re-)initialises the screen. It normally is necessary
!     only to clear the screen.
!     smode = 0 start of program
!     smode = 1 tidy up after end of program
!-----------------------------------------------------------------------
!     Brian Francis, University of Lancaster; 30.1.89
!-----------------------------------------------------------------------
      CHARACTER STR*7
!      integer rc
!-----------------------------------------------------------------------
!     VT100 terminal
!-----------------------------------------------------------------------
      STR = CHAR(27)//'[2J'//CHAR(27)//'[H'
      CALL PROMPT( STR )
!-----------------------------------------------------------------------
      if (smode .eq. 0) then
!          rc = dwfsetapptitle('SABRE')
!          rc = dwfsetaboutdlg('SABRE','SABRE for Windows 3.x, 95 & NT'//
!          char(13)//char(13)//'Version 3.1'//char(13)//char(13)//
!     &    'Copyright (c) 1999'//char(13)//
!     &    'Centre for Applied Statistics'//char(13)//
!     &    'Lancaster University'//char(13))
      else
!          rc = dwfshutdown()
      end if

      RETURN
      END

!***********************************************************************

      SUBROUTINE PROMPT(TEXT)
!-----------------------------------------------------------------------
      CHARACTER*(*) TEXT
!-----------------------------------------------------------------------
!     This routine should write a prompt to the screen, leaving the
!     cursor on the same line as the prompt and immediately to the right
!     of the prompt. This routine is highly machine dependent. Sometimes
!     it is possible to do this in FORTRAN - sometimes an assembler or
!     small routine written in another language is necessary.
!
!     No prompt if terse output is required
!-----------------------------------------------------------------------
!     Brian Francis, University of Lancaster; 30.1.89
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!
      if( .not. terse ) then
         WRITE ( OUTCH, '(a)', advance='no' ) TEXT
      end if
      RETURN
      END

!***********************************************************************

      SUBROUTINE WRTSCR(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the screen
!     moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (outch,'(a)') text
      call flush(outch)

      return
      end

!***********************************************************************

      SUBROUTINE WRTFLN(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the log file
!     moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (ooutch,'(a)') text
      call flush(ooutch)

      return
      end

!***********************************************************************

      SUBROUTINE WRTFLT(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the trace file
!     moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (toutch,'(a)') text
      call flush(toutch)

      return
      end

!***********************************************************************

      SUBROUTINE WRTFLR(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the resid file
!     moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (routch,'(a)') text
      call flush(routch)

      return
      end

!***********************************************************************

      SUBROUTINE WRTFLQ(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the time file
!     moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (qoutch,'(a)') text
      call flush(qoutch)

      return
      end

!***********************************************************************

      SUBROUTINE WRTFAQ(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the time file
!     without moving the cursor to the start of the next line.

!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (qoutch,'(a)',advance='no') text
      call flush(qoutch)

      return
      end

!***********************************************************************

      SUBROUTINE NEWLNS
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the screen
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (outch,'('' '')')
      call flush(outch)

      return
      end

!***********************************************************************

      SUBROUTINE NEWLNF
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the log file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      write (ooutch,'('' '')')
      call flush(ooutch)

      return
      end

!***********************************************************************

      SUBROUTINE NEWLNT
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     Write a new line to the trace file
!-----------------------------------------------------------------------
      write (toutch,'('' '')')
      call flush(toutch)

      return
      end

!***********************************************************************

      SUBROUTINE NEWLNR
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     Write a new line to the resid file
!-----------------------------------------------------------------------
      write (routch,'('' '')')
      call flush(routch)

      return
      end

!***********************************************************************

      SUBROUTINE NEWLNQ
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     Write a new line to the time file
!-----------------------------------------------------------------------
      write (qoutch,'('' '')')
      call flush(qoutch)

      return
      end

!***********************************************************************

      SUBROUTINE CLSFIL(unit)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer unit
!-----------------------------------------------------------------------
!     Close the file attached to unit
!-----------------------------------------------------------------------
      close (unit)

      return
      end

!***********************************************************************

      SUBROUTINE OPNFIL(unit,fname,status,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) fname,status
      integer unit
      logical error
!-----------------------------------------------------------------------
!     Open file fname on unit "unit".
!-----------------------------------------------------------------------
      error = .false.
      open (unit,file=fname,status=status,err=1)

!---- If no error then return here
      return

    1 continue

!---- Error in open statement
      error = .true.

      return
      end

!***********************************************************************

      SUBROUTINE READLN(unit,line,eof)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) line
      integer unit
      logical eof
!-----------------------------------------------------------------------
      eof = .false.
      read (unit,'(a)',end=1) line
      goto 2

    1 continue

!---- End of file reached
      eof = .true.
      line = ' '

    2 continue

      return
      end

!***********************************************************************

      SUBROUTINE READOB(x,maxobs,iobs,nvar,reinch,eof,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      double precision x(MXX)
      integer maxobs,iobs,nvar,reinch
      logical eof,error
!-----------------------------------------------------------------------
!     Read the iobs'th observation for nvar variables
!-----------------------------------------------------------------------
      double precision xbuff(MAXVAR)
      integer j
!-----------------------------------------------------------------------
      eof = .false.
      error = .false.
      read (reinch,*,end=1,err=2) (xbuff(j),j=1,nvar)
      goto 3

!---- End of file found
    1 continue
      eof = .true.
      goto 3

!---- Error in read
    2 continue
      error = .true.

    3 continue

!---- If data was read correctly then store in appropriate elements of x
      if (.not. eof .and. .not. error) then

          do 4 j = 1,nvar
              x(j*maxobs+iobs) = xbuff(j)
    4     continue

      end if

      return
      end
!
!***********************************************************************

      SUBROUTINE GET_ARG_VALUE( iarg, nextarg, arg_value )
!-----------------------------------------------------------------------
!
!    Look for argument value following nextarg
!    There may be an equals sign separating the two, with or without spaces
!
      implicit none
      integer :: iarg, iargc, length
      character(len=*) :: nextarg
      character(len=*) :: arg_value
      character(len=80) :: temp
!
!    Look at last character of argument name
!
      if( nextarg(length(nextarg):length(nextarg)) == '=' ) then
!       If it's an equals sign the next argument should be the value
         iarg = iarg + 1
         if( iarg > iargc() ) then
            call wrtscr( 'No value give for argument '// &
                         nextarg(1:length(nextarg)) )
            stop
         end if

!--------***********
         CALL GETARG( iarg, arg_value )
!--------***********

      else if( index( nextarg, '=' ) > 0 ) then
!       Argument value is actually part of the argument, separated
!       from the prefix by an equals sign
         arg_value = nextarg( index(nextarg,'=')+1:length(nextarg) )
      else
         iarg = iarg + 1
         if( iarg > iargc() ) then
            call wrtscr( 'No value give for argument '// &
                         nextarg(1:length(nextarg)) )
            stop
         end if

!--------***********
         CALL GETARG( iarg, arg_value )
!--------***********

         if( length(arg_value) == 1 .and. arg_value == '=' ) then
!          Equals sign by itself found so next argument should be value
            iarg = iarg + 1
            if( iarg > iargc() ) then
               call wrtscr( 'No value give for argument '// &
                            nextarg(1:length(nextarg)) )
               stop
            end if

!--------***********
         CALL GETARG( iarg, arg_value )
!--------***********

         else if( arg_value(1:1) == '=' ) then
!          Argument value is prefixed by equals sign
            temp = arg_value
            arg_value = temp(2:length(temp))
         else
!          No equals sign prefix so we have just read the argument value,
!          i.e. nothing more to do
         end if
      end if
      return
      end
!
!***********************************************************************
!     Dummy MPI routines for compatibility with parallel Sabre. They are
!     not called but have to be present to satisfy the linker
!***********************************************************************

      SUBROUTINE MPI_BCAST(buffer,count,data_type,root, &
                           communication_group,ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer buffer,count,data_type,root,communication_group,ierror
!-----------------------------------------------------------------------

      return
      end

!***********************************************************************

      SUBROUTINE MPI_ALLREDUCE(sendbuf,receivebuf,count,data_type, &
                               operation,communication_group,ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer sendbuf,receivebuf,count,data_type,operation, &
              communication_group,ierror
!-----------------------------------------------------------------------

      return
      end

!***********************************************************************

      SUBROUTINE MPI_SEND( buffer, count, data_type, destination, tag, &
                           communication_group, ierror )
!-----------------------------------------------------------------------
      implicit none
      integer buffer, count, data_type, destination, tag, &
              communication_group, ierror
      return
      end

!***********************************************************************

      SUBROUTINE MPI_RECV( buffer, count, data_type, source, tag, &
                           communication_group, status, ierror )
!-----------------------------------------------------------------------
      implicit none
      integer buffer, count, data_type, source, tag, &
              communication_group, status, ierror
      return
      end
