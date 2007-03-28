! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_init(ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!
      ierror = 0
!
      return
end subroutine mpi_init
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_finalize(ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!
      ierror = 0
!
      return
end subroutine mpi_finalize
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_comm_size(comm, size, ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer comm, size, ierror
!-----------------------------------------------------------------------
!
      size = 1
      ierror = 0
!
      return
end subroutine mpi_comm_size
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_comm_rank(comm, rank, ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer comm, rank, ierror
!-----------------------------------------------------------------------
!
      rank = 0
      ierror = 0
!
      return
end subroutine mpi_comm_rank
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_bcast(buffer, count, data_type, root, &
   & communication_group, ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer buffer, count, data_type, root, communication_group, &
         & ierror
!-----------------------------------------------------------------------
!
      ierror = 0
!
      return
end subroutine mpi_bcast
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_allreduce(sendbuf, receivebuf, count, data_type, &
   & operation, communication_group, ierror)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer sendbuf, receivebuf, count, data_type, operation, &
         & communication_group, ierror
!-----------------------------------------------------------------------
!
      ierror = 0
!
      return
end subroutine mpi_allreduce
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_send(buffer, count, data_type, destination, tag, &
   & communication_group, ierror)
!-----------------------------------------------------------------------
      implicit none
      integer buffer, count, data_type, destination, tag, &
         & communication_group, ierror
      return
!
      ierror = 0
!
end subroutine mpi_send
!
!***********************************************************************
!
! $RCSfile: sabre_mpi.f90,v $
! $Revision: 1.1 $
! $Date: 2008/02/12 11:36:14 $
!
subroutine mpi_recv(buffer, count, data_type, source, tag, &
   & communication_group, status, ierror)
!-----------------------------------------------------------------------
      implicit none
      integer buffer, count, data_type, source, tag, &
         & communication_group, status, ierror
!
      ierror = 0
!
      return
end subroutine mpi_recv

