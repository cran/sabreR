      module memory

      integer :: total_memory, max_memory = 0

      contains

      subroutine update_memory( n )

      total_memory = total_memory + n
      if( total_memory > max_memory ) max_memory = total_memory

      return
      end subroutine

      end module memory
