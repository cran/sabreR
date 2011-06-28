      subroutine fortfunc(sabfilename,lsabfilename)

      character(len=255) sabfilename
      integer lsabfilename
      character(len=80) :: input_file
      character(len=80) :: outbuf
      integer iarg, versionflag, terseflag, inputflag
      integer ierror
      integer mxx_size, maxy_size, maxvar_size, maxpar_size
      logical stop


!---- Initialise MPI
!-----*************
      call mpi_init(ierror)
!-----*************
!
      if (ierror /= 0) then
          write (outbuf,'(a,i9)') ' Error in mpi_init, ierror = ',ierror
          call wrtlin(outbuf)
          stop
      end if
!
!    Process any command line arguments
!
!
!     Assign default data space
!
      mxx_size = 100000000
      maxy_size = 500000
      maxvar_size = 200
      maxpar_size = 200
!
!    Set flags to argument not present
!
      versionflag = 0
      terseflag = 1
      inputflag = 1

      stop = .false.

      call sabre_main(stop, mxx_size, maxy_size, maxvar_size, maxpar_size, &
                       versionflag, terseflag, inputflag, sabfilename(1:lsabfilename))
      
      return
      end
