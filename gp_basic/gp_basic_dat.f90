module gp_basic_dat_mod

use prec_mod
use constants_mod
use fparser_dat_mod

implicit none

! GP basic and matrix related variables
integer, parameter :: NOT_FACTORISED = 0
integer, parameter :: CH             = 1
integer, parameter :: QR             = 2
integer, parameter :: BK             = 3

integer, parameter :: USER_kernel = 0
integer, parameter :: SE_kernel = 1

type GP_Matrix
    logical :: initialised = .false.
    logical :: equilibrated = .false.
    integer :: factorised = NOT_FACTORISED
    integer :: n
    integer :: m
    real(dp), dimension(:,:), allocatable :: matrix
    real(dp), dimension(:,:), allocatable :: factor
    real(dp), dimension(:), allocatable :: s
    real(dp), dimension(:), allocatable :: tau
    integer, dimension(:), allocatable :: ipiv
end type GP_Matrix

type gp_basic
    logical :: initialised = .false.
    logical :: sparsified = .false.
    character(len=length) :: mode = "full"  ! possibilities: partial, compact, full
    character(len=length) :: sparsemethod = "dtc" ! sparsification approximation: dic, dtc, fitc
    integer :: kernel_type
    character(len=length) :: kernel_string = ""
    type(tComp), pointer :: fcomp(:) => null()
    integer :: n_dof
    integer :: m_f
    integer :: m_g
    integer :: m_teach
    integer :: n_tot
    real(dp) :: jitter
    real(dp) :: logL
    real(dp) :: delta_sq
    real(dp), allocatable :: theta_sq(:)
    real(dp), allocatable :: periodicity(:)
    character(len=length), allocatable :: var_names(:)
    real(dp), allocatable :: f_r(:,:)
    real(dp), allocatable :: g_r(:,:)
    real(dp), allocatable :: k(:)
    real(dp), allocatable :: k_grad(:,:)
    real(dp), allocatable :: k_hess(:,:,:)
    real(dp), allocatable :: mat_inv_k(:)
    real(dp), allocatable :: mat_inv_k_grad(:,:)
    real(dp), allocatable :: Cmat_inv_v(:)
    real(dp) :: y_siginvsq_y
    real(dp), allocatable :: Kmn_siginvsq_y(:)
    real(dp), allocatable :: Kmn_siginvsq_Knm(:,:)
    real(dp) :: y_siginvsq_y_error
    real(dp), allocatable :: Kmn_siginvsq_y_error(:)
    real(dp), allocatable :: Kmn_siginvsq_Knm_error(:,:)
    type(GP_Matrix) :: Kmm
    type(GP_Matrix) :: noise_Kmm
    type(GP_Matrix) :: Cmat
end type gp_basic

end module gp_basic_dat_mod
