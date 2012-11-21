!
!***********************************************************************
!
      subroutine sabre_main( stop, mxx_size, maxy_size, maxvar_size, &
                             maxpar_size, versionflag, terseflag, &
                             inputflag, input_file )

!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character(len=10) :: start_time,finish_time,current_time
      character(len=8) ::  start_date,finish_date,current_date
      character(len=50), allocatable :: mnames(:), name(:)
      character(len=80) :: input_file
      double precision, allocatable :: x(:), beta(:)
      integer istat, length, iostat
      integer, allocatable :: ilev(:), ipos(:)
      integer clock1,clock2,hrs,mins,secs,timed,rate
      integer versionflag, terseflag, inputflag
      logical stop
      integer mxx_size, maxy_size, maxvar_size, maxpar_size
!
      terse = .false.
!
!    Copy the problem size arguments
!
      mxx = mxx_size
      maxy = maxy_size
      maxvar = maxvar_size
      maxpar = maxpar_size
!
!    If the version flag is set then just print version number, set stop flag
!    & return
!
      if( versionflag /= 0 ) then
!
!         Need to call begin to set the output unit number
!---------**********
          call begin
!---------**********
!
          call wrtscr('Sabre version 6.0')
          stop = .true.
          return
      end if
!
      if( terseflag /= 0 ) then
         terse = .true.
      end if
!
      if( inputflag /= 0 ) then
!
!         Attempt to open an input file with specified name
!
!         Need to call begin to set the input unit number
!---------**********
          call begin
!---------**********
!
          open(unit=inch,file=input_file(1:length(input_file)), &
               status='old',iostat=iostat)
!
          if (iostat /= 0) then
              call wrtscr('Error in attempting to open input file:')
              call wrtscr(input_file(1:length(input_file)))
              stop = .true.
              return
          end if
      end if
!
!     Allocate array space
!
      allocate(mnames(maxvar),stat=istat)
      allocate(name(maxvar),stat=istat)
      allocate(x(mxx),stat=istat)
!
      if (istat /= 0) then
          print *, 'Failed to allocate space for X'
      end if
!
      allocate(beta(maxpar*(maxpar+5)),stat=istat)
      allocate(ilev(maxvar),stat=istat)
      allocate(ipos(maxvar),stat=istat)
!
!
!-----------------------------------------------------------------------
!---- Find out how many processors there are and which one this is
!-----******************
      call mpi_comm_size(mpi_comm_world,num_processors,ierror)
!-----******************
!
      if (ierror /= 0) then
          write (outbuf,'(a,i9)') ' Error in mpi_comm_size, ierror = ', &
          ierror
          call wrtlin(outbuf)
          stop = .true.
          return
      end if
!
!-----******************
      call mpi_comm_rank(mpi_comm_world,this_processor,ierror)
!-----******************
!
      if (ierror /= 0) then
          write (outbuf,'(a,i9)') ' Error in mpi_comm_rank, ierror = ', &
          ierror
          call wrtlin(outbuf)
          stop = .true.
          return
      end if
!
!---- Record start clock time on master process
      if (this_processor == boss_processor) then
!
!---------*****************
!          call system_clock(clock1,rate)
!---------*****************
!
!---------******************
          call date_and_time(date=start_date,time=start_time)
!---------******************
!
      end if
!
    3 continue
!
!-----***********
      call sabrem(x,beta,name,mnames,ilev,ipos,stop)
!-----***********
!
      if (.not. stop) then
          go to 3
      end if
!
!---- Finished with MPI
!-----*****************
      call mpi_finalize(ierror)
!-----*****************
!
      if (ierror /= 0) then
          write (outbuf,'(a,i9)') ' Error in mpi_finalize, ierror = ', &
          ierror
          call wrtlin(outbuf)
          stop = .true.
          return
      end if
!
!---- Display start and finish clock time on master process
      if (this_processor == boss_processor) then
!
!---------******************
          call date_and_time(date=finish_date,time=finish_time)
!---------******************
!
!---------*****************
!          call system_clock(clock2)
!---------*****************
!
          write (outbuf,'(5(a,a2),a,a4)') 'Start:  ',start_time(1:2), &
          ':',start_time(3:4),':',start_time(5:6),' on ', &
          start_date(7:8),'/',start_date(5:6),'/',start_date(1:4)
!
          if (.not. terse) then
              call newlns
              call wrtscr(outbuf)
          end if
!
          if (qflag) then
              call wrtflq(outbuf)
          end if
!
          write (outbuf,'(5(a,a2),a,a4)') 'Finish: ',finish_time(1:2), &
          ':',finish_time(3:4),':',finish_time(5:6),' on ', &
          finish_date(7:8),'/',finish_date(5:6),'/',finish_date(1:4)
!
          if (.not. terse) then
              call wrtscr(outbuf)
          end if
!
          if (qflag) then
              call wrtflq(outbuf)
          end if
!
!          timed = nint((clock2-clock1)/real(rate))
!          hrs = timed/3600
!          mins = (timed - hrs*3600)/60
!          secs = timed - hrs*3600 - mins*60
!          write (outbuf,'(a,i3,a,2(i2,a))') 'Time = ',hrs,'hr ',mins, &
!          'min ',secs,'sec'
!
!          if (.not. terse) then
!              call wrtscr(outbuf)
!          end if
!
!          if (qflag) then
!              call wrtflq(outbuf)
!          end if
!
!         Removed the printing of the number of processors used since it's
!         misleading in the Openmp version (where in the code num_processors
!         is 1; all parallelism being done by directives instead of by code)
!
!          if (num_processors == 1) then
!              write (outbuf,'(i3,a)') num_processors,' processor used'
!          else
!              write (outbuf,'(i3,a)') num_processors,' processors used'
!          end if
!
          if (.not. terse) then
              call newlns
              call wrtscr(outbuf)
          end if
!
          if (qflag) then
              call newlnq
              call wrtflq(outbuf)
          end if
!
      end if

!----- close the input file
!-----***********
      call clsfil(inch)
!-----***********
!



!
!---- close the log file
!-----***********
      call clsfil(ooutch)
!-----***********
!
!---- close the trace file
      if (tflag) then
!
!---------***********
          call clsfil(toutch)
!---------***********
!
      end if
!
!---- close the resid file
      if (rflag) then
!
!---------***********
          call clsfil(routch)
!---------***********
!
      end if
!
!---- close the time file
      if (qflag) then
!
!---------***********
          call clsfil(qoutch)
!---------***********
!
      end if

!---- deallocate memory resources
!-----***********
      deallocate(mnames)
      deallocate(name)
      deallocate(x)
      deallocate(beta)
      deallocate(ilev)
      deallocate(ipos)


      return
!
      end subroutine sabre_main
!
!***********************************************************************
!
      subroutine begin
!-----------------------------------------------------------------------
      implicit none
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
      inch  = 17
      outch = 6
      rinch = 15
      inpch = 16
      ooutch = 21
      toutch = 22
      routch = 23
      qoutch = 24
      cinch = inch
!-----------------------------------------------------------------------
!     Initialisation of COMMON ACCMAC machine accuracy
!-----------------------------------------------------------------------
      z1 = 1.0d+308
      z2 = 1.0d-307
      zl1 = dlog(z1)
      zl2 = dlog(z2)
      iz = 2147483647
!
      return
!
      end subroutine begin
!
!***********************************************************************
!
      logical function calpha(c)
!-----------------------------------------------------------------------
      character c
!----------------------------------------------------------------------
!     This function is machine dependent. It should return .FALSE. if C
!     is not alphabetic and .TRUE. if C is alphabetic.
!     "alphabetic" = [a-z][A-Z][$]
!     It assumes an ASCII collating sequence, and will need to be
!     rewritten if moved to a machine using the EBCDIC character codes
!     or anything even more esoteric (BCD??)
!     BJF
!----------------------------------------------------------------------
      integer icharc
!----------------------------------------------------------------------
      icharc = ichar(c)
      calpha = ((icharc >= 65 .and. icharc <= 90) .or. &
      (icharc >= 97 .and. icharc <= 122) .or. c == '$')
!
      return
!
      end function calpha
!
!***********************************************************************
!
      subroutine lower(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character(len=*) line
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
      integer i,j,ich,dollar
      integer length
      data dollar /36/
!-----------------------------------------------------------------------
      j = 1
!
      do 10 i = 1,length(line)
          ich = ichar(line(i:i))
!
          if (ich == dollar) then
              cycle
          end if
!
          if (ich >= 65 .and. ich <= 90) then
              line(j:j) = char(ich+32)
          else
              line(j:j) = char(ich)
          end if
!
          j = j+1
   10 end do
!
      do 20 i = j,length(line)
          line(i:i) = ' '
   20 end do
!
      return
!
      end subroutine lower
!
!***********************************************************************
!
      subroutine lcase(string)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character string*(*)
!-----------------------------------------------------------------------
!     Version of LOWER which just converts letters to lower case, nothing else
!-----------------------------------------------------------------------
      integer i,ich
      integer length
!-----------------------------------------------------------------------
      do 10 i = 1, length(string)
          ich = ichar(string(i:i))
!
          if (ich >= 65 .and. ich <= 90) then
              string(i:i) = char(ich+32)
          end if
!
   10 end do
!
      return
!
      end subroutine lcase
!
!***********************************************************************
!
      subroutine scrset(smode)
!-----------------------------------------------------------------------
      implicit none
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
      character(len=7) :: str
!      integer rc
!-----------------------------------------------------------------------
!     VT100 terminal
!-----------------------------------------------------------------------
      str = char(27)//'[2J'//char(27)//'[H'

!-----***********
!      call prompt(str)
!-----***********

!-----------------------------------------------------------------------
      if (smode == 0) then
!          rc = dwfsetapptitle('SABRE')
!          rc = dwfsetaboutdlg('SABRE','SABRE for Windows 3.x, 95 & NT'//
!          char(13)//char(13)//'Version 3.1'//char(13)//char(13)//
!     &    'Copyright (c) 1999'//char(13)//
!     &    'Centre for Applied Statistics'//char(13)//
!     &    'Lancaster University'//char(13))
      else
!          rc = dwfshutdown()
      end if
!
      return
!
      end subroutine scrset
!
!***********************************************************************
!
      subroutine prompt(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
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
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor .and. .not. terse) then
          inquire(unit=outch,opened=is_open)
!
          if (is_open) then
              write (outch,'(a)',advance='no') text
          end if
!
      end if
!
      return
!
      end subroutine prompt
!
!***********************************************************************
!
      subroutine wrtscr(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the screen
!     moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=outch,opened=is_open)
!
          if (is_open) then
              write(outch,'(a)') text
              call flush(outch)
          end if
!
      end if
!
      return
!
      end subroutine wrtscr
!
!***********************************************************************
!
      subroutine wrtfln(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the log file
!     moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=ooutch,opened=is_open)
!
          if (is_open) then
              write (ooutch,'(a)') text
              call flush(ooutch)
          end if
!
      end if
!
      return
!
      end subroutine wrtfln
!
!***********************************************************************
!
      subroutine wrtflt(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the trace file
!     moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=toutch,opened=is_open)
!
          if (is_open) then
              write (toutch,'(a)') text
              call flush(toutch)
          end if
!
      end if
!
      return
!
      end subroutine wrtflt
!
!***********************************************************************
!
      subroutine wrtflr(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the resid file
!     moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=routch,opened=is_open)
!
          if (is_open) then
              write (routch,'(a)') text
              call flush(routch)
          end if
!
      end if
!
      return
!
      end subroutine wrtflr
!
!***********************************************************************
!
      subroutine wrtflq(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the time file
!     moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=qoutch,opened=is_open)
!
          if (is_open) then
              write (qoutch,'(a)') text
              call flush(qoutch)
          end if
!
      end if
!
      return
!
      end subroutine wrtflq
!
!***********************************************************************
!
      subroutine wrtfaq(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the time file
!     without moving the cursor to the start of the next line.
!
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
      logical is_open
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=qoutch,opened=is_open)
!
          if (is_open) then
              write (qoutch,'(a)',advance='no') text
              call flush(qoutch)
          end if
!
      end if
!
      return
!
      end subroutine wrtfaq
!
!***********************************************************************
!
      subroutine newlns
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the screen
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=outch,opened=is_open)
!
          if (is_open) then
              write (outch,'('' '')')
              call flush(outch)
          end if
!
      end if
!
      return
!
      end subroutine newlns
!
!***********************************************************************
!
      subroutine newlnf
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the log file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire( unit=ooutch, opened=is_open)
!
          if (is_open) then
              write (ooutch,'('' '')')
              call flush(ooutch)
          end if
!
      end if
!
      return
!
      end subroutine newlnf
!
!***********************************************************************
!
      subroutine newlnt
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the trace file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=toutch,opened=is_open)
!
          if (is_open) then
              write (toutch,'('' '')')
              call flush(toutch)
          end if
!
      end if
!
      return
!
      end subroutine newlnt
!
!***********************************************************************
!
      subroutine newlnr
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the resid file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=routch,opened=is_open)
!
          if (is_open) then
              write (routch,'('' '')')
              call flush(routch)
          end if
!
      end if
!
      return
!
      end subroutine newlnr
!
!***********************************************************************
!
      subroutine newlnq
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the time file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire( unit=qoutch, opened=is_open)
!
          if (is_open) then
              write (qoutch,'('' '')')
              call flush(qoutch)
          end if
!
      end if
!
      return
!
      end subroutine newlnq
!
!***********************************************************************
!
      subroutine clsfil(unit)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer unit
!-----------------------------------------------------------------------
!     Close the file attached to unit
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
      logical is_open
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          inquire(unit=unit,opened=is_open)
!
          if (is_open) then
              close(unit)
          end if
!
      end if
!
      return
!
      end subroutine clsfil
!
!***********************************************************************
!
      subroutine opnfil(unit,fname,status,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) fname,status
      integer unit
      logical error
!-----------------------------------------------------------------------
!     Open file fname on unit "unit".
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          error = .false.
          open(unit,file=fname,status=status,err=1)
!
!-------- If no error then return here
          return
!
    1     continue
!
!-------- Error in open statement
          error = .true.
          return
      end if
!
      return
!
      end subroutine opnfil
!
!***********************************************************************
!
      subroutine readln(unit,line,eof)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) line
      integer unit
      logical eof
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          eof = .false.
          read(unit,'(a)',end=1) line
          go to 2
!
    1     continue
!
!-------- End of file reached
          eof = .true.
          line = ' '
!
    2     continue
      end if
!
!---- Broadcast the data
!-----**************
      call mpi_bcast(line(1:len(line)),len(line),mpi_character, &
                     boss_processor, mpi_comm_world,ierror)
!-----**************
!
!-----**************
      call mpi_bcast(eof,1,mpi_logical,boss_processor,mpi_comm_world, &
                     ierror)
!-----**************
!
      return
!
      end subroutine readln
!
!***********************************************************************
!
      subroutine readob(x,maxobs,iobs,nvar,reinch,eof,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      double precision x(mxx)
      integer maxobs,iobs,nvar,reinch
      logical eof,error
!-----------------------------------------------------------------------
!     Read the iobs'th observation for nvar variables
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      double precision xbuff(maxvar)
      integer ierror,j
!-----------------------------------------------------------------------
!---- Only the BOSS does I/O
      if (this_processor == boss_processor) then
          eof = .false.
          error = .false.
          read(reinch,*,end=1,err=2) (xbuff(j),j=1,nvar)
          go to 3
!
!-------- End of file found
    1     continue
          eof = .true.
          go to 3
!
!-------- Error in read
    2     continue
          error = .true.
!
    3     continue
!
!-------- If data was read correctly then store in appropriate elements of x
          if (.not. eof .and. .not. error) then
!
              do 4 j = 1,nvar
                   x(j*maxobs+iobs) = xbuff(j)
  4           end do
!
          end if
!
      end if
!
      return
!
      end subroutine readob

