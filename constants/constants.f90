module constants_mod

use prec_mod

implicit none

real(dp), parameter :: zero = 0.0_dp
real(dp), parameter :: one = 1.0_dp
real(dp), parameter :: two = 2.0_dp
real(dp), parameter :: ten = 10.0_dp
real(dp), parameter :: pi = acos(-one)
real(dp), parameter :: two_pi = two*pi
real(dp), parameter :: pi_per_two = pi/two

integer, parameter :: length = 10000
integer, parameter :: is = selected_int_kind(2)

end module constants_mod
