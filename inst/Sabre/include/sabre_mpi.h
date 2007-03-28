



! $RCSfile: sabre_mpi_serial.h,v $
! $Revision: 1.2 $
! $Date: 2008/02/15 14:37:19 $
!
!
!include 'mpif.h'
!
      integer :: mpi_comm_world, mpi_logical, mpi_integer, &
                 mpi_double_precision, mpi_lor, mpi_sum, mpi_character, &
                 mpi_any_tag
      integer, parameter :: mpi_status_size=3
!
      common /mpi_info/ num_processors,this_processor
      integer num_processors,this_processor
!
      integer boss_processor
      parameter(boss_processor=0)
