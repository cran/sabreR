   module accurate_arithmetic
!  version: 23.10.07, 15:35

      implicit none

!    Set the range of values which can be exponentiated
!    (0.4343 is log10(exp(1d0)))
!    (2.30258509 is log(10))

!      double precision, parameter :: max_exp_exp = huge(0)/0.4343, &
!                                     min_exp_exp = -huge(0)/0.4343
      double precision, parameter :: max_exp_exp = huge(0)*2.30258509, &
                                     min_exp_exp = -huge(0)*2.30258509

      type accurate
         double precision :: man
         integer :: exp
      end type accurate

      interface operator (+)
         module procedure acc_acc_plus, int_acc_plus, acc_int_plus, &
                          double_acc_plus, acc_double_plus
      end interface

      interface operator (-)
         module procedure acc_acc_minus, int_acc_minus, acc_int_minus, &
                          double_acc_minus, acc_double_minus, unary_minus
      end interface

      interface operator (*)
         module procedure acc_acc_mult, int_acc_mult, acc_int_mult, &
                          double_acc_mult, acc_double_mult
      end interface

      interface operator (/)
         module procedure acc_acc_div, int_acc_div, acc_int_div, &
                          double_acc_div, acc_double_div
      end interface

      interface operator (**)
         module procedure acc_int_power
      end interface

      interface assignment(=)
         module procedure double_acc_assign, acc_double_assign, &
                          acc_real_assign, acc_int_assign
      end interface

      interface operator(.ge.)
         module procedure acc_double_ge
      end interface

      interface operator(.le.)
         module procedure acc_double_le
      end interface

      interface operator(.gt.)
         module procedure acc_double_gt, acc_int_gt
      end interface

      interface operator(.lt.)
         module procedure acc_double_lt
      end interface

      interface operator(.eq.)
         module procedure acc_double_eq, acc_real_eq, acc_int_eq
      end interface

      interface operator(.ne.)
         module procedure acc_double_ne, acc_real_ne, acc_int_ne
      end interface

      interface sqrt
         module procedure acc_sqrt
      end interface

      interface exp
         module procedure acc_exp
      end interface

      interface log
         module procedure acc_dlog
      end interface

      interface log10
         module procedure acc_dlog10
      end interface

      interface abs
         module procedure acc_dabs
      end interface

      interface dble
         module procedure acc_dble
      end interface

      interface fexp
         module procedure fexp, acc_fexp
      end interface

      interface erfc
         module procedure double_erfc, acc_erfc
      end interface


   contains


      function acc_acc_plus( left, right )

      implicit none
      type(accurate) :: acc_acc_plus
      type(accurate), intent(in) :: left, right

      call addnum( left%man, left%exp, right%man, right%exp, &
                   acc_acc_plus%man, acc_acc_plus%exp )

      end function acc_acc_plus


      function int_acc_plus( left, right )

      implicit none
      type(accurate) :: int_acc_plus
      integer, intent(in) :: left
      type(accurate), intent(in) :: right

      call addnum( real(left,kind(0d0)), 0, right%man, right%exp, &
                   int_acc_plus%man, int_acc_plus%exp )

      end function int_acc_plus


      function acc_int_plus( left, right )

      implicit none
      type(accurate) :: acc_int_plus
      type(accurate), intent(in) :: left
      integer, intent(in) :: right

      call addnum( left%man, left%exp, real(right,kind(0d0)), 0, &
                   acc_int_plus%man, acc_int_plus%exp )

      end function acc_int_plus


      function double_acc_plus( left, right )

      implicit none
      type(accurate) :: double_acc_plus
      double precision, intent(in) :: left
      type(accurate), intent(in) :: right

      call addnum( left, 0, right%man, right%exp, &
                   double_acc_plus%man, double_acc_plus%exp )

      end function double_acc_plus


      function acc_double_plus( left, right )

      implicit none
      type(accurate) :: acc_double_plus
      type(accurate), intent(in) :: left
      double precision, intent(in) :: right

      call addnum( left%man, left%exp, right, 0, &
                   acc_double_plus%man, acc_double_plus%exp )

      end function acc_double_plus


      function acc_acc_minus( left, right )

      implicit none
      type(accurate) :: acc_acc_minus
      type(accurate), intent(in) :: left, right

      call addnum( left%man, left%exp, -right%man, right%exp, &
                   acc_acc_minus%man, acc_acc_minus%exp )

      end function acc_acc_minus


      function int_acc_minus( left, right )

      implicit none
      type(accurate) :: int_acc_minus
      integer, intent(in) :: left
      type(accurate), intent(in) :: right

      call addnum( real(left,kind(0d0)), 0, -right%man, right%exp, &
                   int_acc_minus%man, int_acc_minus%exp )

      end function int_acc_minus


      function acc_int_minus( left, right )

      implicit none
      type(accurate) :: acc_int_minus
      type(accurate), intent(in) :: left
      integer, intent(in) :: right

      call addnum( left%man, left%exp, -real(right,kind(0d0)), 0, &
                   acc_int_minus%man, acc_int_minus%exp )

      end function acc_int_minus


      function double_acc_minus( left, right )

      implicit none
      type(accurate) :: double_acc_minus
      double precision, intent(in) :: left
      type(accurate), intent(in) :: right

      call addnum( left, 0, -right%man, right%exp, &
                   double_acc_minus%man, double_acc_minus%exp )

      end function double_acc_minus


      function acc_double_minus( left, right )

      implicit none
      type(accurate) :: acc_double_minus
      type(accurate), intent(in) :: left
      double precision, intent(in) :: right

      call addnum( left%man, left%exp, -right, 0, &
                   acc_double_minus%man, acc_double_minus%exp )

      end function acc_double_minus


      function unary_minus( right )

      implicit none
      type(accurate) :: unary_minus
      type(accurate), intent(in) :: right

      unary_minus%man = - right%man
      unary_minus%exp = right%exp

      end function unary_minus


      function acc_acc_mult( left, right )

      implicit none
      type(accurate) :: acc_acc_mult
      type(accurate), intent(in) :: left, right

      acc_acc_mult%man = left%man * right%man
      acc_acc_mult%exp = left%exp + right%exp

      call manexp( acc_acc_mult%man, acc_acc_mult%exp )

      end function acc_acc_mult


      function int_acc_mult( left, right )

      implicit none
      type(accurate) :: int_acc_mult
      integer, intent(in) :: left
      type(accurate), intent(in) :: right

      int_acc_mult%man = left * right%man
      int_acc_mult%exp = right%exp

      call manexp( int_acc_mult%man, int_acc_mult%exp )

      end function int_acc_mult


      function acc_int_mult( left, right )

      implicit none
      type(accurate) :: acc_int_mult
      type(accurate), intent(in) :: left
      integer, intent(in) :: right

      acc_int_mult%man = left%man * right
      acc_int_mult%exp = left%exp

      call manexp( acc_int_mult%man, acc_int_mult%exp )

      end function acc_int_mult


      function double_acc_mult( left, right )

      implicit none
      type(accurate) :: double_acc_mult
      double precision, intent(in) :: left
      type(accurate), intent(in) :: right

      double_acc_mult%man = left * right%man
      double_acc_mult%exp = right%exp

      call manexp( double_acc_mult%man, double_acc_mult%exp )

      end function double_acc_mult


      function acc_double_mult( left, right )

      implicit none
      type(accurate) :: acc_double_mult
      type(accurate), intent(in) :: left
      double precision, intent(in) :: right

      acc_double_mult%man = left%man * right
      acc_double_mult%exp = left%exp

      call manexp( acc_double_mult%man, acc_double_mult%exp )

      end function acc_double_mult


      function acc_acc_div( left, right )

      implicit none
      type(accurate) :: acc_acc_div
      type(accurate), intent(in) :: left, right

      acc_acc_div%man = left%man / right%man
      acc_acc_div%exp = left%exp - right%exp

      call manexp( acc_acc_div%man, acc_acc_div%exp )

      end function acc_acc_div


      function int_acc_div( left, right )

      implicit none
      type(accurate) :: int_acc_div
      integer, intent(in) :: left
      type(accurate), intent(in) :: right

      int_acc_div%man = left / right%man
      int_acc_div%exp = - right%exp

      call manexp( int_acc_div%man, int_acc_div%exp )

      end function int_acc_div


      function acc_int_div( left, right )

      implicit none
      type(accurate) :: acc_int_div
      type(accurate), intent(in) :: left
      integer, intent(in) :: right

      acc_int_div%man = left%man / right
      acc_int_div%exp = left%exp

      call manexp( acc_int_div%man, acc_int_div%exp )

      end function acc_int_div


      function double_acc_div( left, right )

      implicit none
      type(accurate) :: double_acc_div
      double precision, intent(in) :: left
      type(accurate), intent(in) :: right

      double_acc_div%man = left / right%man
      double_acc_div%exp = - right%exp

      call manexp( double_acc_div%man, double_acc_div%exp )

      end function double_acc_div


      function acc_double_div( left, right )

      implicit none
      type(accurate) :: acc_double_div
      type(accurate), intent(in) :: left
      double precision, intent(in) :: right

      acc_double_div%man = left%man / right
      acc_double_div%exp = left%exp

      call manexp( acc_double_div%man, acc_double_div%exp )

      end function acc_double_div


      function acc_int_power( left, right )

      implicit none
      type(accurate) :: acc_int_power
      type(accurate), intent(in) :: left
      integer, intent(in) :: right

      acc_int_power%man = left%man ** right
      acc_int_power%exp = left%exp * right

      call manexp( acc_int_power%man, acc_int_power%exp )

      end function acc_int_power


      function acc_double_ge( left, right )

      implicit none
      logical :: acc_double_ge
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_ge = left_double .ge. right

      end function acc_double_ge


      function acc_double_le( left, right )

      implicit none
      logical :: acc_double_le
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_le = left_double .le. right

      end function acc_double_le


      function acc_double_gt( left, right )

      implicit none
      logical :: acc_double_gt
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_gt = left_double .gt. right

      end function acc_double_gt


      function acc_int_gt( left, right )

      implicit none
      logical :: acc_int_gt
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      integer, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so make a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_int_gt = left_double .gt. right

      end function acc_int_gt


      function acc_double_lt( left, right )

      implicit none
      logical :: acc_double_lt
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_lt = left_double .lt. right

      end function acc_double_lt


      function acc_double_eq( left, right )

      implicit none
      logical :: acc_double_eq
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_eq = left_double .eq. right

      end function acc_double_eq


      function acc_real_eq( left, right )

      implicit none
      logical :: acc_real_eq
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      real, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_real_eq = left_double .eq. right

      end function acc_real_eq


      function acc_int_eq( left, right )

      implicit none
      logical :: acc_int_eq
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      integer, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_int_eq = left_double .eq. right

      end function acc_int_eq


      function acc_double_ne( left, right )

      implicit none
      logical :: acc_double_ne
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      double precision, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_double_ne = left_double .ne. right

      end function acc_double_ne


      function acc_real_ne( left, right )

      implicit none
      logical :: acc_real_ne
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      real, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_real_ne = left_double .ne. right

      end function acc_real_ne


      function acc_int_ne( left, right )

      implicit none
      logical :: acc_int_ne
      type(accurate), intent(in) :: left
      type(accurate) :: left_copy
      integer, intent(in) :: right
      double precision :: left_double

!    Make use of the fact that normal sets the double result to the
!    maximum double if the accurate argument is beyond that range

!    Normal may modify its input argument so take a copy and use that
      left_copy = left
      call normal( left_copy%man, left_copy%exp, left_double )

      acc_int_ne = left_double .ne. right

      end function acc_int_ne


      subroutine double_acc_assign( left, right )

      implicit none
      double precision, intent(out) :: left
      type(accurate), intent(in) :: right
      type(accurate) :: right_copy

!    Normal may modify its input argument so take a copy and use that
      right_copy = right
      call normal( right_copy%man, right_copy%exp, left )

      end subroutine double_acc_assign


      subroutine acc_double_assign( left, right )

      implicit none
      type(accurate), intent(out) :: left
      double precision, intent(in) :: right

      left%man = right
      left%exp = 0

      end subroutine acc_double_assign


      subroutine acc_real_assign( left, right )

      implicit none
      type(accurate), intent(out) :: left
      real, intent(in) :: right

      left%man = right
      left%exp = 0

      end subroutine acc_real_assign


      subroutine acc_int_assign( left, right )

      implicit none
      type(accurate), intent(out) :: left
      integer, intent(in) :: right

      left%man = right
      left%exp = 0

      end subroutine acc_int_assign


      function acc_sqrt( arg )

      implicit none
      type(accurate) :: acc_sqrt
      type(accurate), intent(in) :: arg

      if( mod( arg%exp, 2 ) .ne. 0 ) then
         acc_sqrt%man = arg%man * 1d1
         acc_sqrt%exp = arg%exp - 1
      else
         acc_sqrt%man = arg%man
         acc_sqrt%exp = arg%exp
      end if

      acc_sqrt%man = dsqrt( acc_sqrt%man )
      acc_sqrt%exp = acc_sqrt%exp / 2

      call manexp( acc_sqrt%man, acc_sqrt%exp )

      end function acc_sqrt


      function acc_exp( arg )

      implicit none
      type(accurate) :: arg_copy, acc_exp
      type(accurate), intent(in) :: arg
      double precision :: arg_double

      common /accmac/ z1,z2,zl1,zl2,iz
      double precision z1,z2,zl1,zl2
      integer iz

!    Check that argument is such that exp(arg) has an exponent which is
!    within the range of allowed integers

      if( arg .lt. max_exp_exp .and. arg .gt. min_exp_exp ) then

!       Convert to double & use emexp

!       Normal may modify its input argument so take a copy and use that
         arg_copy = arg
         call normal( arg_copy%man, arg_copy%exp, arg_double )

         call emexp( arg_double, acc_exp%man, acc_exp%exp )

      else if( arg .ge. max_exp_exp ) then

         call emexp( max_exp_exp, acc_exp%man, acc_exp%exp )

      else

         call emexp( min_exp_exp, acc_exp%man, acc_exp%exp )

      end if

      end function acc_exp


      function acc_dlog( arg )

      implicit none
      type(accurate) :: acc_dlog
      type(accurate), intent(in) :: arg

      acc_dlog%man = log( arg%man ) + arg%exp * log( 10d0 )
      acc_dlog%exp = 0

      call manexp( acc_dlog%man, acc_dlog%exp )

      end function acc_dlog


      function acc_dlog10( arg )

      implicit none
      type(accurate) :: acc_dlog10
      type(accurate), intent(in) :: arg

      acc_dlog10%man = dlog10( arg%man ) + arg%exp
      acc_dlog10%exp = 0

      call manexp( acc_dlog10%man, acc_dlog10%exp )

      end function acc_dlog10


      function acc_dabs( arg )

      implicit none
      type(accurate) :: acc_dabs
      type(accurate), intent(in) :: arg

      acc_dabs%man = dabs( arg%man )
      acc_dabs%exp = arg%exp

      end function acc_dabs


      function acc_dble( arg )

      implicit none
      double precision :: acc_dble
      type(accurate), intent(in) :: arg
      type(accurate) :: arg_copy
      double precision :: arg_double

!    Normal may modify its input argument so take a copy and use that
      arg_copy = arg
      call normal( arg_copy%man, arg_copy%exp, arg_double )

      acc_dble = arg_double

      end function acc_dble


      function acc_fexp( arg, iflow )

      implicit none
      type(accurate) :: arg_copy, acc_fexp
      type(accurate), intent(in) :: arg
      logical, intent(inout) :: iflow(2)
      double precision :: arg_double

!    Check that argument is such that exp(arg) has an exponent which is
!    less than the maximum integer allowed

      if( arg .lt. max_exp_exp .and. arg .gt. min_exp_exp ) then

!       Convert to double & use emexp

!       Normal may modify its input argument so take a copy and use that
          arg_copy = arg
         call normal( arg_copy%man, arg_copy%exp, arg_double )

         call emexp( arg_double, acc_fexp%man, acc_fexp%exp )

         iflow(1) = .false.
         iflow(2) = .false.
      else if( arg .le. min_exp_exp ) then
         acc_fexp%man = 0d0
         acc_fexp%exp = 0
!         iflow(1) = .true.
      else
         iflow(2) = .true.
      end if

      end function acc_fexp


      function acc_erfc( arg )

      implicit none
      type(accurate) :: arg_copy, acc_erfc
      type(accurate), intent(in) :: arg
      double precision :: arg_double

!    Convert to double & use normal erfc (normal adjusts to ensure that
!    arg_copy lies within double precision range

!    Normal may modify its input argument so take a copy and use that
      arg_copy = arg
      call normal( arg_copy%man, arg_copy%exp, arg_double )

      acc_erfc = double_erfc( arg_double )

      end function acc_erfc


!***********************************************************************

      function fexp(z,iflow)

      implicit none
      double precision fexp
!-----------------------------------------------------------------------
      double precision z
      logical iflow(2)
!-----------------------------------------------------------------------
!     Function : Takes exponentials, dealing with underflow & overflow.
!-----------------------------------------------------------------------
      common /accmac/ z1,z2,zl1,zl2,iz
      double precision z1,z2,zl1,zl2
      integer iz
!-----------------------------------------------------------------------
!---- prevent underflow
      if (z .lt. zl2) then
          fexp = 0
!          iflow(1) = .true.

!---- prevent overflow
      else if (z .gt. zl1) then
          iflow(2) = .true.

!---- calculate exp(z)
      else
          fexp = dexp(z)
      end if

      return

      end function fexp


!***********************************************************************

      function double_erfc(x)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision double_erfc
      double precision, intent(in) :: x
!-----------------------------------------------------------------------
! Complementary error function
! Taken from SUN's FDLIBM version 5.2 and translated from c to fortran.
!-----------------------------------------------------------------------
! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
!
! Developed at SunSoft, a Sun Microsystems, Inc. business.
! Permission to use, copy, modify, and distribute this
! software is freely granted, provided that this notice 
! is preserved.
!-----------------------------------------------------------------------
! Definition:
!------------
!                       x
!                2      |\
! erf(x)  =  ---------  | exp(-t*t)dt
!             sqrt(pi) \| 
!                       0
!
! erfc(x) =  1 - erf(x)
!
! Note that  erf(-x) = -erf(x)
!            erfc(-x) = 2 - erfc(x)
!
! Method:
!--------
!
! 1. For |x| in [0, 0.84375]
!    erf(x)  = x + x*R(x^2)
!    erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
!            = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
!    where R = P/Q where P is an odd poly of degree 8 and
!                        Q is an odd poly of degree 10.
!                                                 -57.90
!                        | R - (erf(x)-x)/x | <= 2
!
!
!    Remark. The formula is derived by noting
!            erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
!    and that
!            2/sqrt(pi) = 1.128379167095512573896158903121545171688
!    is close to one. The interval is chosen because the fix
!    point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
!    near 0.6174), and by some experiment, 0.84375 is chosen to
!    guarantee the error is less than one ulp for erf.
!
! 2. For |x| in [0.84375,1.25], let s = |x| - 1, and c = 0.84506291151
!    rounded to single (24 bits)
!    erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
!    erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
!              1+(c+P1(s)/Q1(s))    if x < 0
!              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
!
!    Remark: here we use the taylor series expansion at x=1.
!            erf(1+s) = erf(1) + s*Poly(s)
!                     = 0.845.. + P1(s)/Q1(s)
!            That is, we use rational approximation to approximate
!            erf(1+s) - (c = (single)0.84506291151)
!    Note that |P1/Q1|< 0.078 for x in [0.84375,1.25] where 
!    P1(s) = degree 6 poly in s
!    Q1(s) = degree 6 poly in s
!
! 3. For x in [1.25,1/0.35(~2.857143)], 
!    erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
!    erf(x)  = 1 - erfc(x)
!    where 
!    R1(z) = degree 7 poly in z, (z=1/x^2)
!    S1(z) = degree 8 poly in z
!
! 4. For x in [1/0.35,28]
!    erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
!            = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
!            = 2.0 - tiny(if x <= -6)
!    erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
!    erf(x)  = sign(x)*(1.0 - tiny)
!    where
!    R2(z) = degree 6 poly in z, (z=1/x^2)
!    S2(z) = degree 7 poly in z
!
!    Note1: To compute exp(-x*x-0.5625+R/S), let s be a single
!           precision number and s := x; then
!           -x*x = -s*s + (s-x)*(s+x)
!            exp(-x*x-0.5626+R/S) = exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S)
!    Note2: Here 4 and 5 make use of the asymptotic series
!                      exp(-x*x)
!           erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
!                     x*sqrt(pi)
!           We use rational approximation to approximate
!           g(s) = f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
!           Here is the error bound for R1/S1 and R2/S2
!           |R1/S1 - f(x)|  < 2**(-62.57)
!           |R2/S2 - f(x)|  < 2**(-61.52)
!
! 5. For inf > x >= 28
!    erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
!    erfc(x) = tiny*tiny (raise underflow) if x > 0
!            = 2 - tiny if x<0
!
! 7. Special case:
!    erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
!    erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2, 
!    erfc/erf(NaN) is NaN
!-----------------------------------------------------------------------
      double precision             :: ax,p,q,r,s,y,z
      double precision, parameter  :: zero =  0.0d0
      double precision, parameter  :: half =  0.5d0
      double precision, parameter  :: one  =  1.0d0
      double precision, parameter  :: two  =  2.0d0
      double precision, parameter  :: erx  =  8.45062911510467529297d-01
!-----------------------------------------------------------------------
!     Coefficients for approximation to erf on [0,0.84375]
!-----------------------------------------------------------------------
      double precision, parameter  :: efx  =  1.28379167095512586316d-01
      double precision, parameter  :: efx8 =  1.02703333676410069053d+00
      double precision, parameter  :: pp0  =  1.28379167095512558561d-01
      double precision, parameter  :: pp1  = -3.25042107247001499370d-01
      double precision, parameter  :: pp2  = -2.84817495755985104766d-02
      double precision, parameter  :: pp3  = -5.77027029648944159157d-03
      double precision, parameter  :: pp4  = -2.37630166566501626084d-05
      double precision, parameter  :: qq1  =  3.97917223959155352819d-01
      double precision, parameter  :: qq2  =  6.50222499887672944485d-02
      double precision, parameter  :: qq3  =  5.08130628187576562776d-03
      double precision, parameter  :: qq4  =  1.32494738004321644526d-04
      double precision, parameter  :: qq5  = -3.96022827877536812320d-06
!-----------------------------------------------------------------------
!     Coefficients for approximation to erf in [0.84375,1.25] 
!-----------------------------------------------------------------------
      double precision, parameter  :: pa0  = -2.36211856075265944077d-03
      double precision, parameter  :: pa1  =  4.14856118683748331666d-01
      double precision, parameter  :: pa2  = -3.72207876035701323847d-01
      double precision, parameter  :: pa3  =  3.18346619901161753674d-01
      double precision, parameter  :: pa4  = -1.10894694282396677476d-01
      double precision, parameter  :: pa5  =  3.54783043256182359371d-02
      double precision, parameter  :: pa6  = -2.16637559486879084300d-03
      double precision, parameter  :: qa1  =  1.06420880400844228286d-01
      double precision, parameter  :: qa2  =  5.40397917702171048937d-01
      double precision, parameter  :: qa3  =  7.18286544141962662868d-02
      double precision, parameter  :: qa4  =  1.26171219808761642112d-01
      double precision, parameter  :: qa5  =  1.36370839120290507362d-02
      double precision, parameter  :: qa6  =  1.19844998467991074170d-02
!-----------------------------------------------------------------------
!     Coefficients for approximation to erfc in [1.25,1/0.35]
!-----------------------------------------------------------------------
      double precision, parameter  :: ra0  = -9.86494403484714822705d-03
      double precision, parameter  :: ra1  = -6.93858572707181764372d-01
      double precision, parameter  :: ra2  = -1.05586262253232909814d+01
      double precision, parameter  :: ra3  = -6.23753324503260060396d+01
      double precision, parameter  :: ra4  = -1.62396669462573470355d+02
      double precision, parameter  :: ra5  = -1.84605092906711035994d+02
      double precision, parameter  :: ra6  = -8.12874355063065934246d+01
      double precision, parameter  :: ra7  = -9.81432934416914548592d+00
      double precision, parameter  :: sa1  =  1.96512716674392571292d+01
      double precision, parameter  :: sa2  =  1.37657754143519042600d+02
      double precision, parameter  :: sa3  =  4.34565877475229228821d+02
      double precision, parameter  :: sa4  =  6.45387271733267880336d+02
      double precision, parameter  :: sa5  =  4.29008140027567833386d+02
      double precision, parameter  :: sa6  =  1.08635005541779435134d+02
      double precision, parameter  :: sa7  =  6.57024977031928170135d+00
      double precision, parameter  :: sa8  = -6.04244152148580987438d-02
!-----------------------------------------------------------------------
!     Coefficients for approximation to erfc in [1/.35,28]
!-----------------------------------------------------------------------
      double precision, parameter  :: rb0  = -9.86494292470009928597d-03
      double precision, parameter  :: rb1  = -7.99283237680523006574d-01
      double precision, parameter  :: rb2  = -1.77579549177547519889d+01
      double precision, parameter  :: rb3  = -1.60636384855821916062d+02
      double precision, parameter  :: rb4  = -6.37566443368389627722d+02
      double precision, parameter  :: rb5  = -1.02509513161107724954d+03
      double precision, parameter  :: rb6  = -4.83519191608651397019d+02
      double precision, parameter  :: sb1  =  3.03380607434824582924d+01
      double precision, parameter  :: sb2  =  3.25792512996573918826d+02
      double precision, parameter  :: sb3  =  1.53672958608443695994d+03
      double precision, parameter  :: sb4  =  3.19985821950859553908d+03
      double precision, parameter  :: sb5  =  2.55305040643316442583d+03
      double precision, parameter  :: sb6  =  4.74528541206955367215d+02
      double precision, parameter  :: sb7  = -2.24409524465858183362d+01
!-----------------------------------------------------------------------
      ax = dabs(x)

      if (ax < 0.84375d0) then

          if (ax < epsilon(x)) then
              double_erfc = one-x
          else
              z = x**2
              r = pp0 + z*(pp1 + z*(pp2 + z*(pp3 + z*pp4)))
              s = one + z*(qq1 + z*(qq2 + z*(qq3 + z*(qq4 + z*qq5))))
              y = r/s

              if (x < 0.25d0) then
                  double_erfc =  one - (x + x*y)
              else
                  r = x*y
                  r = r+x-half
                  double_erfc = half-r 
              end if

          end if

      else if (ax < 1.25d0) then
          s = ax-one
          p = pa0 + s*(pa1 + s*(pa2 + s*(pa3 + s*(pa4 + &
          s*(pa5 + s*pa6)))))
          q = one + s*(qa1 + s*(qa2 + s*(qa3 + s*(qa4 + &
          s*(qa5 + s*qa6)))))

          if (x > zero) then
              z  = one-erx
              double_erfc = z - p/q
          else
              z = erx + p/q
              double_erfc = one+z
          end if

      else if (ax < 28.0d0) then
          s = one/(ax**2)

          if (ax < 2.857143d0) then
              p = ra0 + s*(ra1 + s*(ra2 + s*(ra3 + s*(ra4 + &
              s*(ra5 + s*(ra6 + s*ra7))))))
              q = one + s*(sa1 + s*(sa2 + s*(sa3 + s*(sa4 + &
              s*(sa5 + s*(sa6 + s*(sa7 + s*sa8)))))))
          else

              if (x < -6.0d0) then
                  double_erfc = two
                  return
              end if

              p = rb0 + s*(rb1 + s*(rb2 + s*(rb3 + s*(rb4 + &
              s*(rb5 + s*rb6)))))
              q = one + s*(sb1 + s*(sb2 + s*(sb3 + s*(sb4 + &
              s*(sb5 + s*(sb6 + s*sb7))))))
          end if

          z = ax
          r = exp(-z**2 - 0.5625d0)*exp((z-ax)*(z+ax) + p/q)

          if (x > zero) then
              double_erfc = r/x
          else
              double_erfc = two + r/x
!             double_erfc = two - r/x
          end if

      else

          if (x > zero) then
              double_erfc = zero
          else
              double_erfc = two
          end if

      end if

      end function double_erfc


   end module accurate_arithmetic

