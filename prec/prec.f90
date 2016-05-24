module prec_mod

implicit none

! precision
#ifndef QP
integer, parameter :: dp = 8
#else
integer, parameter :: dp = 16
#endif

end module prec_mod
