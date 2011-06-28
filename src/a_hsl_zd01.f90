! *******************************************************************
! COPYRIGHT (c) 2000 Council for the Central Laboratory
!                    of the Research Councils
! All rights reserved.
!
! None of the comments in this Copyright notice between the lines
! of asterisks shall be removed or altered in any way.
!
! This Package is intended for compilation without modification,
! so most of the embedded comments have been removed.
!
! ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
! Licence, see http://hsl.rl.ac.uk/acuk/cou.html
!
! Please note that for a UK ACADEMIC Licence:
!
! 1. The Packages may only be used for academic research or teaching
!    purposes by the Licensee, and must not be copied by the Licensee for
!    use by any other persons. Use of the Packages in any commercial
!    application shall be subject to prior written agreement between
!    Hyprotech UK Limited and the Licensee on suitable terms and
!    conditions, which will include financial conditions.
! 2. All information on the Package is provided to the Licensee on the
!    understanding that the details thereof are confidential.
! 3. All publications issued by the Licensee that include results obtained
!    with the help of one or more of the Packages shall acknowledge the
!    use of the Packages. The Licensee will notify the Numerical Analysis
!    Group at Rutherford Appleton Laboratory of any such publication.
! 4. The Packages may be modified by or on behalf of the Licensee
!    for such use in research applications but at no time shall such
!    Packages or modifications thereof become the property of the
!    Licensee. The Licensee shall make available free of charge to the
!    copyright holder for any purpose all information relating to
!    any modification.
! 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
!    direct or consequential loss or damage whatsoever arising out of
!    the use of Packages by the Licensee.
! *******************************************************************
!
! Original date 10 March 2000
! 8 Nov. 2000 Optional argument stat add to ZD01_put

! 12th July 2004 Version 1.0.0. Version numbering added.

module hsl_zd01_char
implicit none
contains
   subroutine zd01_put(array,string,stat)
     character, pointer :: array(:)
     character(*), intent(in), optional ::  string
     integer, intent(out), optional ::  stat
     integer :: i,l
     logical :: ok
     l = 0
     if (present(string)) l = len_trim(string)
     if (present(stat)) then
       allocate(array(l),stat=stat)
       ok = stat==0
     else
       allocate(array(l))
       ok = .true.
     end if
     if (ok) then
       do i = 1, l
         array(i) = string(i:i)
       end do
     end if
   end subroutine zd01_put
   function zd01_get(array)
     character :: array(:)
     character(size(array)) ::  zd01_get
     integer :: i
     do i = 1, size(array)
        zd01_get(i:i) = array(i)
     end do
   end function zd01_get
end module hsl_zd01_char
module hsl_zd01_double
      use hsl_zd01_char
      implicit none
      integer, private, parameter :: wp = kind( 1.0d+0 )
      type, public :: zd01_type
         integer :: m, n, ne
         character, pointer, dimension( : ) :: id, type
         integer, pointer, dimension( : ) :: row, col, ptr
         real ( kind = wp ), pointer, dimension( : ) :: val
      end type
   end module hsl_zd01_double

