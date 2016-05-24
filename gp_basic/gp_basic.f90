!% A simple Gaussian process module for n-D functions learned from samples of
!% their values or gradients, with sparsification

module gp_basic_mod

use prec_mod
#ifndef EXTERNAL_BLAS_LAPACK
use blas_lapack_mod 
#endif
use gp_basic_dat_mod
use fparser_mod

implicit none
private

public :: gp_basic_teach, finalise
public :: f_predict, f_predict_grad, f_predict_hess, f_predict_int
public :: f_predict_var, f_predict_grad_var, f_predict_int_var
public :: f_predict_var_grad
public :: gp_basic_change, gp_basic_merge, gp_basic_complete
public :: gp_basic_write, gp_basic_read

!!% initialise (and teach) a gp_basic
!interface initialise
!   module procedure gp_basic_initialise_nd
!end interface initialise

!% finalise and deallocate a gp_basic
interface finalise
   module procedure gp_basic_finalise
end interface finalise

!% predict a function value from a gp
interface f_predict
   module procedure f_predict_r
end interface f_predict

!% predict a gradient of a function from a gp
interface f_predict_grad
   module procedure f_predict_grad_r
end interface f_predict_grad

!% predict a hessian of a function from a gp
interface f_predict_hess
   module procedure f_predict_hess_r
end interface f_predict_hess

!% predict the integral of a function from a gp
interface f_predict_int
   module procedure f_predict_int_r
end interface f_predict_int

!% predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r
end interface f_predict_var

!% predict the variance of a gradient of a function from a gp
interface f_predict_grad_var
   module procedure f_predict_grad_var_r
end interface f_predict_grad_var

!% predict the variance of an integral of a function from a gp
interface f_predict_int_var
   module procedure f_predict_int_var_r
end interface f_predict_int_var

!% predict the gradient of a function variance from a gp
interface f_predict_var_grad
   module procedure f_predict_var_grad_r
end interface f_predict_var_grad

interface Matrix_Solve
   module procedure GP_Matrix_Solve_Vector, GP_Matrix_Solve_Matrix
end interface Matrix_Solve

interface assignment(=)
   module procedure GP_Matrix_Assign
end interface

interface assignment(=)
   module procedure gp_basic_assign
end interface

interface operator(+)
   module procedure gp_basic_add
end interface

contains

subroutine gp_basic_teach(self, kernel_type, &
                          f_r, f_v, f_s, g_r, g_v, g_s, &
                          f_sparse_r, g_sparse_r, &
                          kernel_string, var_names, &
                          delta, theta, periodicity, &
                          sparsemethod, jitter, factorization, partial, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self !% object to store GP
   integer, intent(in) :: kernel_type
   real(dp), optional, intent(in) :: f_r(:,:), f_v(:), f_s(:) !% arrays of function positions, values, sigma 
   real(dp), optional, intent(in) :: g_r(:,:), g_v(:,:), g_s(:,:) !% arrays of function gradient positions, values, sigma
   real(dp), optional, intent(in) :: f_sparse_r(:,:), g_sparse_r(:,:) !% pseudo positions to use for sparsifcation for values and gradients
   character(len=length), optional, intent(in) :: kernel_string !% LAM
   character(len=length), optional, intent(in) :: var_names(:) !% LAM
   real(dp), optional, intent(in) :: delta !% delta
   real(dp), optional, intent(in) :: theta(:) !% theta 
   real(dp), optional, intent(in) :: periodicity(:) !% periodicity for GP
   character(len=length), optional, intent(in) :: sparsemethod
   real(dp), optional, intent(in) :: jitter !% jitter for sparse calculation
   integer, optional, intent(in) :: factorization
   logical, optional, intent(in) :: partial
   logical, optional, intent(in) :: verbose

   real(dp), allocatable :: Cmat(:,:), Kmn(:,:), y(:), siginvsq_y(:), siginvsq_Knm(:,:), Kmm(:,:)
   real(dp), allocatable :: Kmn_siginvsq_y(:), Kmn_siginvsq_Knm(:,:)
   real(dp), allocatable :: Kmm_inv_kn(:), Lnn(:)
   integer :: i, n_f, n_g, n_tot, i_glob, ii
   integer :: factor
   logical :: do_partial, do_verbose
   real(dp) :: logdet
   real(dp) :: maxdiag, sjitter

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (self%initialised) then
       if (trim(self%mode) .eq. 'compact') then
           return
       end if
   end if

   if (present(factorization)) then
      if(factorization /= CH .and. factorization /= QR .and. factorization /= BK) then
        print *,"unknown factorization, factorization= ", factorization
        stop
      end if 
      factor = factorization
   else
      factor = CH
   end if

   if (present(partial)) then
      do_partial = partial
   else
      do_partial = .false.
   end if

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   ! check input consistency
   i = count( (/ present(f_r), present(f_v), present(f_s) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of f_r, f_v, f_s"
      stop
   endif
   i = count( (/ present(g_r), present(g_v), present(g_s) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of g_r, g_v, g_s"
      stop
   endif

   ! compare f_r to f_[gn]
   if (present(f_r)) then
      if (size(f_r,2) /= size(f_v) .or. size(f_r,2) /= size(f_s)) then
         print *,"size(f_r,2)=",size(f_r,2)," size(f_v)=",size(f_v)," size(f_s)=",size(f_s)," not all equal"
         stop
      endif
   endif
   ! compare g_r to g_[gn]
   if (present(g_r)) then
      if (any(shape(g_r) /= shape(g_v)) .or. any(shape(g_r) /= shape(g_s))) then
         print *,"shape(g_r)=",shape(g_r)," shape(g_v)=",shape(g_v)," shape(g_s)=",shape(g_s)," not all equal"
         stop
      endif
   endif
   ! compare f_r to g_r
   if (present(f_r) .and. present(g_r)) then
      if (size(f_r,1) /= size(g_r,1)) then
         print *,"size(f_r,1)=",size(f_r,1)," size(g_r,1)=",size(g_r,1)," not equal"
         stop
      endif
   endif

   if (present(f_r)) then
      n_f = size(f_r,2)
   else
      n_f = 0
   end if

   if (present(g_r)) then
      n_g = size(g_r,2)
   else
      n_g = 0
   end if

   if ( (n_f == 0) .and. (n_g == 0)) then
      print *,"no teaching data provided"
      stop
   endif

   if ( (.not. self%initialised) .or. (.not. self%sparsified) ) then
      call finalise(self)

      if (present(f_r)) then
         self%n_dof = size(f_r,1)
      end if
      if (present(g_r)) then
         self%n_dof = size(g_r,1)
      end if

      if (present(f_sparse_r)) then
         if (size(f_sparse_r,1) /= self%n_dof) then
            print *,"size(f_sparse_r,1)=",size(f_sparse_r,1)," n_dof=",self%n_dof," not equal"
            stop
         end if
         self%sparsified = .true.
         self%m_f = size(f_sparse_r,2)
         allocate(self%f_r(self%n_dof,self%m_f))
         self%f_r(:,:) = f_sparse_r(:,:)
      endif
      if (present(g_sparse_r)) then
         if (size(g_sparse_r,1) /= self%n_dof) then
            print *,"size(g_sparse_r,1)=",size(g_sparse_r,1)," n_dof=",self%n_dof," not equal"
            stop
         end if
         self%sparsified = .true.
         self%m_g = size(g_sparse_r,2)
         allocate(self%g_r(self%n_dof,self%m_g))
         self%g_r(:,:) = g_sparse_r(:,:)
      endif

      if ( (self%m_f == 0) .and. (self%m_g == 0) ) then
         if (present(f_r)) then
            self%m_f = n_f
            allocate(self%f_r(self%n_dof,self%m_f))
            self%f_r(:,:) = f_r(:,:)
         end if
         if (present(g_r)) then
            self%m_g = n_g
            allocate(self%g_r(self%n_dof,self%m_g))
            self%g_r(:,:) = g_r(:,:)
         end if
      end if

      if (present(jitter)) self%jitter = jitter

      if (present(sparsemethod)) self%sparsemethod = trim(sparsemethod)

      self%kernel_type = kernel_type

      select case(self%kernel_type)
         case(USER_kernel)

            if (.not. present(kernel_string)) then
               print *,"kernel_string is not present"
               stop
            end if

            if (.not. present(var_names)) then
               print *,"var_names is not present"
               stop
            end if

            if (size(var_names) /= 2*self%n_dof) then
               print *,"size(var_names)=",size(var_names)," , 2*n_dof=",2*self%n_dof
               stop
            end if

            allocate(self%var_names(2*self%n_dof))
            self%kernel_string = kernel_string
            self%var_names = var_names

            call initf(self%fcomp, 1)
            call parsef(self%fcomp, 1, self%kernel_string, self%var_names)

         case(SE_kernel)
      
           if (.not. present(delta)) then
              print *,"delta is not present"
              stop
           endif

           if (.not. present(theta)) then 
              print *,"theta is not present"
              stop
           end if

           if (.not. present(periodicity)) then
              print *,"periodicity is not present"
              stop
           end if

           if (delta <= 0.0_dp) then
              print *,"invalid delta=",delta
              stop
           endif

           if (size(theta) /= self%n_dof) then
              print *,"size(theta)=",size(theta)," , n_dof=",self%n_dof
              stop
           endif

           if (size(periodicity) /= self%n_dof) then
              print *,"size(periodicity)=",size(periodicity)," /= n_dof=",self%n_dof
              stop
           endif

           do i=1,self%n_dof
              if (theta(i) <= 0.0_dp) then
                 print *,"invalid theta(", i, ")=",theta(i)
                 stop
              end if
           end do

           allocate(self%theta_sq(self%n_dof))
           allocate(self%periodicity(self%n_dof))
           self%theta_sq = theta**2
           self%periodicity = periodicity
           self%delta_sq = delta**2

      case default
         print *,"kernel type does not exist"
         stop
      end select

      if (self%m_f + self%m_g == 0) then
         print *,"no sparsified teaching data provided"
         stop
      endif

      self%m_teach = self%m_f + self%m_g*self%n_dof
   end if

   n_tot = n_f + n_g*self%n_dof

   ! everyone needs Kmm - calculate it only if required and not already available
   if ( ((.not. self%sparsified) .or. (.not. do_partial) .or. (self%sparsemethod .eq. 'fitc')) &
         .and. (.not. self%Kmm%initialised) ) then
      allocate(Kmm(self%m_teach, self%m_teach))
      if(do_verbose) print *,"Calculation of Kmm matrix has started..."
      call kernel_mat(self, self%m_f, self%f_r, self%m_g, self%g_r, &
                      self%m_f, self%f_r, self%m_g, self%g_r, Kmm)
      if(do_verbose) print *,"Calculation of Kmm matrix has finished..."
   end if

   if (self%sparsified) then
      allocate(Kmn(self%m_teach, n_tot))
      allocate(siginvsq_Knm(n_tot, self%m_teach))
      allocate(y(n_tot))
      allocate(siginvsq_y(n_tot))

      ! calculate Kmn
      if(do_verbose) print *,"Calculation of Kmn matrix has started..."
      call kernel_mat(self, self%m_f, self%f_r, self%m_g, self%g_r, & 
                      n_f, f_r, n_g, g_r, Kmn)
      if(do_verbose) print *,"Calculation of Kmn matrix has finished..."

      if ( ((.not. do_partial) .or. (trim(self%sparsemethod) .eq. 'fitc')) .and. (.not. self%Kmm%initialised) ) then
         ! we'll need Kmm^{-1} for variance
         call gp_matrix_initialise(self%Kmm, Kmm, factor)
      end if

      if (trim(self%sparsemethod) .eq. 'fitc') then
          if(do_verbose) print *,"Calculation of Lnn matrix has started..."
          allocate(Lnn(n_tot))
          allocate(Kmm_inv_kn(self%m_teach))

          call gp_matrix_initialise(self%noise_Kmm, self%Kmm%matrix, factor)
          maxdiag = self%noise_Kmm%matrix(1,1)
          do i=2,self%m_teach
             if(self%noise_Kmm%matrix(i,i) .gt. maxdiag) maxdiag = self%noise_Kmm%matrix(i,i)
          end do
          sjitter = self%jitter*10.0_dp**(int(log10(maxdiag))+1)
!$omp parallel do
          do i=1,self%m_teach
             self%noise_Kmm%matrix(i,i) = self%noise_Kmm%matrix(i,i) + sjitter
          end do
!$omp end parallel do

          call GP_Matrix_QR_Factorise(self%noise_Kmm)

!$omp parallel do private(Kmm_inv_kn)
          do i=1, n_f
             call Matrix_Solve(self%noise_Kmm, Kmn(:,i), Kmm_inv_kn)
             Lnn(i) = f_kernel_sca(self, f_r(:,i), f_r(:,i)) - ddot(self%m_teach, Kmm_inv_kn(:), 1, Kmn(:,i), 1)
          end do
!$omp end parallel do
!$omp parallel do private(Kmm_inv_kn)
          do i=1, n_g
             Lnn(n_f+(i-1)*self%n_dof+1:i*self%n_dof) = g_kernel_sca(self, g_r(:,i), g_r(:,i))
             do ii=1, self%n_dof
                call Matrix_Solve(self%noise_Kmm, Kmn(:,n_f+(i-1)*self%n_dof+ii), Kmm_inv_kn)
                Lnn(n_f+(i-1)*self%n_dof+ii) = Lnn(n_f+(i-1)*self%n_dof+ii) &
                                               - ddot(self%m_teach, Kmm_inv_kn(:), 1, Kmn(:,n_f+(i-1)*self%n_dof+ii), 1)
             end do
          end do
!$omp end parallel do
          deallocate(Kmm_inv_kn)
          call gp_matrix_finalise(self%noise_Kmm)
          if(do_verbose) print *,"Calculation of Lnn matrix has finished..."
      end if

      if(do_verbose) print *,"Calculation of siginvsq_y has started..."
      if (trim(self%sparsemethod) .eq. 'fitc') then
!$omp parallel do
          do i=1, n_f
             siginvsq_y(i) = f_v(i)/(Lnn(i)+f_s(i)**2)
          end do
!$omp end parallel do
!$omp parallel do
          do i=1, n_g
             siginvsq_y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = &
                                                       g_v(:,i)/(Lnn(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof)+g_s(:,i)**2)
          end do
!$omp end parallel do
      else
!$omp parallel do
          do i=1, n_f
             y(i) = f_v(i)
             siginvsq_y(i) = f_v(i)/f_s(i)**2
          end do
!$omp end parallel do
!$omp parallel do
          do i=1, n_g
             y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = g_v(:,i)
             siginvsq_y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = g_v(:,i)/g_s(:,i)**2
          end do
!$omp end parallel do
      end if
      if(do_verbose) print *,"Calculation of siginvsq_y has finished..."

      if(do_verbose) print *,"Calculation of Kmn_siginvsq_y has started..."
      if (allocated(self%Kmn_siginvsq_y)) then
          allocate(Kmn_siginvsq_y(self%m_teach))
!#ifndef QP
          call dgemv('N', size(Kmn,1), size(Kmn,2), 1.0_dp, Kmn, size(Kmn,1), &
                     siginvsq_y, 1, 0.0_dp, Kmn_siginvsq_y, 1)
!$omp parallel do
          do i=1, self%m_teach
             self%Kmn_siginvsq_y(i) = self%Kmn_siginvsq_y(i) + Kmn_siginvsq_y(i)
          end do
!$omp end parallel do
          deallocate(Kmn_siginvsq_y)
      else
          allocate(self%Kmn_siginvsq_y(self%m_teach))
          call dgemv('N', size(Kmn,1), size(Kmn,2), 1.0_dp, Kmn, size(Kmn,1), &
                     siginvsq_y, 1, 0.0_dp, self%Kmn_siginvsq_y, 1)
      end if
      if(do_verbose) print *,"Calculation of Kmn_siginvsq_y has finished..."

      if(do_verbose) print *,"Calculation of siginvsq_Knm has started..."
!$omp parallel workshare
      siginvsq_Knm = transpose(Kmn)
!$omp end parallel workshare
      if (trim(self%sparsemethod) .eq. 'fitc') then
!$omp parallel do
          do i=1, n_f
             siginvsq_Knm(i,:) = siginvsq_Knm(i,:) / (Lnn(i)+f_s(i)**2)
          end do
!$omp end parallel do
!$omp parallel do
          do i=1, n_g
             do ii=1, self%n_dof
                siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) = &
                                        siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) / (Lnn(n_f+(i-1)*self%n_dof+ii)+g_s(ii,i)**2)
             end do
          end do
!$omp end parallel do
      else
!$omp parallel do
          do i=1, n_f
             siginvsq_Knm(i,:) = siginvsq_Knm(i,:) / f_s(i)**2
          end do
!$omp end parallel do
!$omp parallel do
          do i=1, n_g
             do ii=1, self%n_dof
                siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) = siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) / g_s(ii,i)**2
             end do
          end do
!$omp end parallel do
      end if
      if(do_verbose) print *,"Calculation of siginvsq_Knm has finished..."

      if(do_verbose) print *,"Calculation of Kmn_siginvsq_Knm has started..."
      if (allocated(self%Kmn_siginvsq_Knm)) then
          allocate(Kmn_siginvsq_Knm(self%m_teach, self%m_teach))
          call dgemm('N', 'N', size(Kmn,1), size(siginvsq_Knm,2), size(Kmn,2), &
                     1.0_dp, Kmn, size(Kmn,1), siginvsq_Knm, size(siginvsq_Knm,1), &
                     0.0_dp, Kmn_siginvsq_Knm, size(Kmn_siginvsq_Knm,1))
!$omp parallel workshare
          self%Kmn_siginvsq_Knm = self%Kmn_siginvsq_Knm + Kmn_siginvsq_Knm
!$omp end parallel workshare
          deallocate(Kmn_siginvsq_Knm)
      else
          allocate(self%Kmn_siginvsq_Knm(self%m_teach, self%m_teach))
          call dgemm('N', 'N', size(Kmn,1), size(siginvsq_Knm,2), size(Kmn,2), &
                     1.0_dp, Kmn, size(Kmn,1), siginvsq_Knm, size(siginvsq_Knm,1), &
                     0.0_dp, self%Kmn_siginvsq_Knm, size(self%Kmn_siginvsq_Knm,1))
      end if
      if(do_verbose) print *,"Calculation of Kmn_siginvsq_Knm has finished..."

      self%y_siginvsq_y = self%y_siginvsq_y + ddot(n_tot, y, 1, siginvsq_y, 1)
      self%n_tot = self%n_tot + n_tot

      if (do_partial) then
         self%mode = 'partial'
         if (self%Cmat%initialised) call gp_matrix_finalise(self%Cmat)
      else
         allocate(Cmat(self%m_teach, self%m_teach))
!$omp parallel workshare
         Cmat = self%Kmm%matrix + self%Kmn_siginvsq_Knm
!$omp end parallel workshare
         ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})
         maxdiag = Cmat(1,1)
         do i=2,self%m_teach
            if(Cmat(i,i) .gt. maxdiag) maxdiag = Cmat(i,i)
         end do
         sjitter = self%jitter*10.0_dp**(int(log10(maxdiag))+1)
!$omp parallel do
         do i=1,self%m_teach
            Cmat(i,i) = Cmat(i,i) + sjitter
         end do
!$omp end parallel do

         call gp_matrix_initialise(self%Cmat, Cmat, factor)

         if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
         if(allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
         allocate(self%Cmat_inv_v(self%m_teach))
         call Matrix_Solve(self%Cmat, self%Kmn_siginvsq_y, self%Cmat_inv_v)
         if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."

         if(do_verbose) print *,"Calculation of logL has started..."
         ! logL = -0.5*log(det(C)) - 0.5*y^{T} C^{-1} y - N/2*log(2*pi), where C = \Sigma^{2} + K_{nm} K_{mm}^{-1} K_{mn}
         ! -0.5*log(det(C)) = -0.5*log(R_{11}*R_{22}*...*R_{mm})
         logdet = 0.0_dp
         do i=1, self%Cmat%m
            logdet = logdet + log(abs(self%Cmat%factor(i,i)))
         end do
         self%logL = -0.5_dp * logdet

         ! -0.5*y^{T} C^{-1} y = -0.5*(y^{T} \Sigma^{-2} y - y^{T} \Sigma^{-2} K_{nm} Cmat^{-1}  K_{mn} \Sigma^{-2} y)
         ! where Cmat = K_{mm} + K_{mn} \Sigma^{-2} K_{nm}
         self%logL = self%logL - 0.5_dp * (self%y_siginvsq_y - ddot(self%m_teach, self%Kmn_siginvsq_y, 1, self%Cmat_inv_v, 1))

         ! - N/2*log(2*pi)
         self%logL = self%logL - 0.5_dp * real(self%n_tot) * log(TWO_PI)
         if(do_verbose) print *,"Calculation of logL has finished..."

         deallocate(Cmat)
      end if

      deallocate(Kmn)
      deallocate(siginvsq_Knm)
      deallocate(y)
      deallocate(siginvsq_y)
      if (trim(self%sparsemethod) .eq. 'fitc') deallocate(Lnn)
   else 
      ! not sparsified
      allocate(y(n_tot))

      ! we'll need self%noise_Kmm, which is Kmm shifted by noise, for variance
      if(do_verbose) print *,"Calculation of Cmat has started..."
!$omp parallel do
      do i=1, self%m_f
         Kmm(i,i) = Kmm(i,i) + f_s(i)**2
      end do
!$omp end parallel do
!$omp parallel do
      do i=1, self%m_g
         do ii=1, self%n_dof
            i_glob = self%m_f + (i-1)*self%n_dof+ii
            Kmm(i_glob, i_glob)= Kmm(i_glob, i_glob) + g_s(ii,i)**2
         end do
      end do
!$omp end parallel do
      maxdiag = Kmm(1,1)
      do i=2,self%m_teach
         if(Kmm(i,i) .gt. maxdiag) maxdiag = Kmm(i,i)
      end do
      sjitter = self%jitter*10.0_dp**(int(log10(maxdiag))+1)
!$omp parallel do
      do i=1,self%m_teach
         Kmm(i,i) = Kmm(i,i) + sjitter
      end do
!$omp end parallel do

      call gp_matrix_initialise(self%noise_Kmm, Kmm, factor)
      if(do_verbose) print *,"Calculation of Cmat has finished..."

!$omp parallel do
      do i=1,n_f
         y(i) = f_v(i)
      end do
!$omp end parallel do
!$omp parallel do
      do i=1, n_g
         i_glob = n_f + (i-1)*self%n_dof
         y(i_glob+1:i_glob+self%n_dof) = g_v(1:self%n_dof,i)
      end do
!$omp end parallel do

      if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
      allocate(self%Cmat_inv_v(self%m_teach))
      call Matrix_Solve(self%noise_Kmm, y, self%Cmat_inv_v)
      if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."

      self%n_tot = n_tot
      if(do_verbose) print *,"Calculation of logL has started..."
      ! logL = -0.5*log(det(C)) - 0.5*y^{T} C^{-1} y - N/2*log(2*pi), where C = noise_Kmm
      ! -0.5*log(det(C)) = -0.5*log(R_{11}*R_{22}*...*R_{mm})
      logdet = 0.0_dp
      do i=1, self%noise_Kmm%m
         logdet = logdet + log(abs(self%noise_Kmm%factor(i,i)))
      end do
      self%logL = -0.5_dp * logdet

      ! -0.5*y^{T} C^{-1} y
      self%logL = self%logL - 0.5_dp * ddot(n_tot, y, 1, self%Cmat_inv_v, 1)

      ! - N/2*log(2*pi)
      self%logL = self%logL - 0.5_dp * real(self%n_tot) * log(TWO_PI)
      if(do_verbose) print *,"Calculation of logL has finished..."

      deallocate(y)
   endif

   if(allocated(Kmm)) deallocate(Kmm)
   if(.not. allocated(self%k)) allocate(self%k(self%m_teach))
   if(.not. allocated(self%k_grad)) allocate(self%k_grad(self%n_dof, self%m_teach))
   if(.not. allocated(self%k_hess)) allocate(self%k_hess(self%n_dof, self%n_dof, self%m_teach))
   if(.not. allocated(self%mat_inv_k)) allocate(self%mat_inv_k(self%m_teach))
   if(.not. allocated(self%mat_inv_k_grad)) allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))

   self%initialised = .true.
end subroutine gp_basic_teach

subroutine gp_basic_finalise(self)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP

   self%initialised = .false.
   self%sparsified = .false.
   self%kernel_type = -1
   self%mode = 'full'
   self%sparsemethod = 'dtc'

   self%n_dof = 0
   self%m_f = 0
   self%m_g = 0
   self%m_teach = 0
   self%n_tot = 0
   self%jitter = 0.0_dp
   self%delta_sq = 0.0_dp
   self%y_siginvsq_y = 0.0_dp
   self%logL = 0.0_dp

   self%kernel_string = ''

   if (associated(self%fcomp)) deallocate(self%fcomp)
   if (allocated(self%var_names)) deallocate(self%var_names)
   if (allocated(self%theta_sq)) deallocate(self%theta_sq)
   if (allocated(self%periodicity)) deallocate(self%periodicity)
   if (allocated(self%f_r)) deallocate(self%f_r)
   if (allocated(self%g_r)) deallocate(self%g_r)
   if (allocated(self%k)) deallocate(self%k)
   if (allocated(self%k_grad)) deallocate(self%k_grad)
   if (allocated(self%k_hess)) deallocate(self%k_hess)
   if (allocated(self%mat_inv_k)) deallocate(self%mat_inv_k)
   if (allocated(self%mat_inv_k_grad)) deallocate(self%mat_inv_k_grad)
   if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
   if (allocated(self%Kmn_siginvsq_y)) deallocate(self%Kmn_siginvsq_y)
   if (allocated(self%Kmn_siginvsq_Knm)) deallocate(self%Kmn_siginvsq_Knm)
   call gp_matrix_finalise(self%Kmm)
   call gp_matrix_finalise(self%Cmat)
   call gp_matrix_finalise(self%noise_Kmm)

end subroutine gp_basic_finalise

subroutine gp_basic_assign(left_self, right_self)

   implicit none

   type(gp_basic), intent(inout) :: left_self
   type(gp_basic), intent(in) :: right_self

   call gp_basic_finalise(left_self)

   left_self%initialised = right_self%initialised
   if (.not. left_self%initialised) return

   left_self%sparsified = right_self%sparsified
   left_self%mode = right_self%mode
   left_self%sparsemethod = right_self%sparsemethod

   left_self%kernel_type = right_self%kernel_type
       
   left_self%n_dof = right_self%n_dof
   left_self%m_f = right_self%m_f
   left_self%m_g = right_self%m_g
   left_self%m_teach = right_self%m_teach
   left_self%n_tot = right_self%n_tot
   left_self%jitter = right_self%jitter
   left_self%delta_sq = right_self%delta_sq
   left_self%y_siginvsq_y = right_self%y_siginvsq_y
   left_self%logL = right_self%logL

   left_self%kernel_string = right_self%kernel_string

   if (allocated(right_self%var_names)) then
      allocate(left_self%var_names(size(right_self%var_names)))
      left_self%var_names = right_self%var_names
   end if

   if (left_self%kernel_type == USER_kernel) then
      call initf(left_self%fcomp, 1)
      call parsef(left_self%fcomp, 1, left_self%kernel_string, left_self%var_names)
   end if

   if (allocated(right_self%theta_sq)) then
      allocate(left_self%theta_sq(left_self%n_dof))
      left_self%theta_sq = right_self%theta_sq
   end if
   if (allocated(right_self%periodicity)) then
      allocate(left_self%periodicity(left_self%n_dof))
      left_self%periodicity = right_self%periodicity
   end if
   if (allocated(right_self%f_r)) then
      allocate(left_self%f_r(left_self%n_dof,left_self%m_f))
      left_self%f_r = right_self%f_r
   end if
   if (allocated(right_self%g_r)) then
      allocate(left_self%g_r(left_self%n_dof,left_self%m_g))
      left_self%g_r = right_self%g_r
   end if
   if (allocated(right_self%k)) then
      allocate(left_self%k(left_self%m_teach))
      left_self%k = right_self%k
   end if
   if (allocated(right_self%k_grad)) then
      allocate(left_self%k_grad(left_self%n_dof,left_self%m_teach))
      left_self%k_grad = right_self%k_grad
   end if
   if (allocated(right_self%k_hess)) then
      allocate(left_self%k_hess(left_self%n_dof,left_self%n_dof,left_self%m_teach))
      left_self%k_hess = right_self%k_hess
   end if
   if (allocated(right_self%mat_inv_k)) then
      allocate(left_self%mat_inv_k(left_self%m_teach))
      left_self%mat_inv_k = right_self%mat_inv_k
   end if
   if (allocated(right_self%mat_inv_k_grad)) then
      allocate(left_self%mat_inv_k_grad(left_self%m_teach,left_self%n_dof))
      left_self%mat_inv_k_grad = right_self%mat_inv_k_grad
   end if
   if (allocated(right_self%Cmat_inv_v)) then
      allocate(left_self%Cmat_inv_v(left_self%m_teach))
      left_self%Cmat_inv_v = right_self%Cmat_inv_v
   end if
   if (allocated(right_self%Kmn_siginvsq_y)) then
      allocate(left_self%Kmn_siginvsq_y(left_self%m_teach))
      left_self%Kmn_siginvsq_y = right_self%Kmn_siginvsq_y
   end if
   if (allocated(right_self%Kmn_siginvsq_Knm)) then
      allocate(left_self%Kmn_siginvsq_Knm(left_self%m_teach,left_self%m_teach))
      left_self%Kmn_siginvsq_Knm = right_self%Kmn_siginvsq_Knm
   end if

   left_self%Kmm = right_self%Kmm
   left_self%noise_Kmm = right_self%noise_Kmm
   left_self%Cmat = right_self%Cmat

end subroutine gp_basic_assign

function gp_basic_add(self1, self2) result(self3)

   implicit none

   type(gp_basic), intent(in) :: self1
   type(gp_basic), intent(in) :: self2
   type(gp_basic) :: self3

   integer :: i, j

   if ( .not. self1%initialised ) then
      print *,"self1 is not initialized!"
      stop
   end if

   if ( .not. self2%initialised ) then
      print *,"self2 is not initialized!"
      stop
   end if

   ! check compatibility of self1 and self2
   if (self1%kernel_type .ne. self2%kernel_type) then
      print *,"self1%kernel_type=", self1%kernel_type, " is not equal to self2%kernel_type=", self2%kernel_type
      stop
   end if

   if (self1%n_dof .ne. self2%n_dof) then
      print *,"self1%n_dof=", self1%n_dof, " is not equal to self2%n_dof=", self2%n_dof
      stop
   end if

   if (self1%m_f .ne. self2%m_f) then
      print *,"self1%m_f=", self1%m_f, " is not equal to self2%m_f=", self2%m_f
      stop
   end if

   if (self1%m_g .ne. self2%m_g) then
      print *,"self1%m_g=", self1%m_g, " is not equal to self2%m_g=", self2%m_g
      stop
   end if

   if (self1%m_teach .ne. self2%m_teach) then
      print *,"self1%m_teach=", self1%m_teach, " is not equal to self2%m_teach=", self2%m_teach
      stop
   end if

   if (self1%jitter .ne. self2%jitter) then
      print *,"self1%jitter=", self1%jitter, " is not equal to self2%jitter=", self2%jitter
      stop
   end if

   if (self1%delta_sq .ne. self2%delta_sq) then
      print *,"self1%delta_sq=", self1%delta_sq, " is not equal to self2%delta_sq=", self2%delta_sq
      stop
   end if

   if (self1%kernel_string .ne. self2%kernel_string) then
      print *,"self1%kernel_string=", self1%kernel_string, " is not equal to self2%kernel_string=", self2%kernel_string
      stop
   end if

   if (.not. self1%sparsified) then
      print *,"self1 is not sparsified!"
      stop
   end if

   if (.not. self2%sparsified) then
      print *,"self2 is not sparsified!"
      stop
   end if

   select case(self1%kernel_type)
      case(USER_kernel)

         do i=1,2*self1%n_dof
            if (self1%var_names(i) .ne. self2%var_names(i)) then
               print *,"self1%var_names(",i,")=",self1%var_names(i), &
                       " is not equal to self2%var_names(",i,")=",self2%var_names(i)
            end if
         end do

      case(SE_kernel) 
         do i=1,self1%n_dof
            if (self1%theta_sq(i) .ne. self2%theta_sq(i)) then
               print *,"self1%theta_sq(",i,")=",self1%theta_sq(i), &
                       " is not equal to self2%theta_sq(",i,")=",self2%theta_sq(i)
               stop
            end if
         end do

         do i=1,self1%n_dof
            if (self1%periodicity(i) .ne. self2%periodicity(i)) then
               print *,"self1%periodicity(",i,")=",self1%periodicity(i), &
                       " is not equal to self2%periodicity(",i,")=",self2%periodicity(i)
               stop
            end if
         end do

   end select

   if (self1%m_f .ne. 0) then
      do j=1,self1%m_f
         do i=1,self1%n_dof
            if (self1%f_r(i,j) .ne. self2%f_r(i,j)) then
               print *,"self1%f_r(",i,j,")=",self1%f_r(i,j)," is not equal to self2%f_r(",i,j,")=",self2%f_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if (self1%m_g .ne. 0) then
      do j=1,self1%m_g
         do i=1,self1%n_dof
            if (self1%g_r(i,j) .ne. self2%g_r(i,j)) then
               print *,"self1%g_r(",i,j,")=",self1%g_r(i,j)," is not equal to self2%g_r(",i,j,")=",self2%g_r(i,j)
               stop
            end if
         end do
      end do
   end if

   call finalise(self3)

   self3%sparsified = self1%sparsified
   self3%mode = 'partial'
   self3%sparsemethod = self1%sparsemethod

   self3%kernel_type = self1%kernel_type

   self3%n_dof = self1%n_dof
   self3%m_f = self1%m_f
   self3%m_g = self1%m_g
   self3%m_teach = self1%m_teach
   self3%n_tot = self1%n_tot + self2%n_tot
   self3%jitter = self1%jitter
   self3%delta_sq = self1%delta_sq
   self3%kernel_string = self1%kernel_string

   self3%y_siginvsq_y = self1%y_siginvsq_y + self2%y_siginvsq_y

   select case(self3%kernel_type)
      case(USER_kernel)
         allocate(self3%var_names(size(self1%var_names)))
         self3%var_names = self1%var_names
         call initf(self3%fcomp, 1)
         call parsef(self3%fcomp, 1, self3%kernel_string, self3%var_names)
      case(SE_kernel)
         allocate(self3%theta_sq(self3%n_dof))
         self3%theta_sq = self1%theta_sq
         allocate(self3%periodicity(self3%n_dof))
         self3%periodicity = self1%periodicity
   end select
   if (self3%m_f .ne. 0) then
      allocate(self3%f_r(self3%n_dof,self3%m_f))
      self3%f_r = self1%f_r
   end if
   if (self3%m_g .ne. 0) then
      allocate(self3%g_r(self3%n_dof,self3%m_g))
      self3%g_r = self1%g_r
   end if
   allocate(self3%Kmn_siginvsq_y(self3%m_teach))
   self3%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
   allocate(self3%Kmn_siginvsq_Knm(self3%m_teach,self3%m_teach))
   self3%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm

end function gp_basic_add

subroutine gp_basic_write(self,fileid,filename,mode)

   implicit none

   type(gp_basic), intent(in) :: self !% object for GP
   integer, intent(in) :: fileid !% fileid
   character(len=*), intent(in) :: filename !% filename to write GP out
   character(len=*), intent(in), optional :: mode !% what kind of information should be written out

   character(len=length) :: do_mode

   integer :: i, j

   if (.not. self%initialised) then
      return
   endif

   if(present(mode)) then
      do_mode = trim(mode)
      if(trim(self%mode) .eq. 'compact') then
         do_mode = 'compact'
      else if(trim(self%mode) .eq. 'partial') then
         do_mode = 'partial'
      end if
   else
      do_mode = self%mode
   end if

   open(fileid,file=filename)

   ! obligatory variables
   write(fileid,*) '# self%mode'
   write(fileid,*) trim(do_mode)
   write(fileid,*) '# self%n_dof'
   write(fileid,*) self%n_dof
   write(fileid,*) '# self%m_f'
   write(fileid,*) self%m_f
   write(fileid,*) '# self%m_g'
   write(fileid,*) self%m_g
   write(fileid,*) '# self%m_teach'
   write(fileid,*) self%m_teach
   write(fileid,*) '# self%sparsified'
   write(fileid,*) self%sparsified
   if(self%sparsified) then
      write(fileid,*) '# self%sparsemethod'
      write(fileid,*) trim(self%sparsemethod)
   end if
   write(fileid,*) '# self%kernel_type'
   write(fileid,*) self%kernel_type
   select case(self%kernel_type)
      case(USER_kernel)
         write(fileid,*) '# self%kernel_string'
         write(fileid,*) trim(self%kernel_string)
         write(fileid,*) '# self%var_names(2*self%n_dof)'
         do i=1, 2*self%n_dof
            write(fileid,*) trim(self%var_names(i))
         end do
      case(SE_kernel)
         write(fileid,*) '# self%delta_sq'
         write(fileid,100) self%delta_sq
         write(fileid,*) '# self%theta_sq(self%n_dof)'
         do i=1, self%n_dof
            write(fileid,100,advance='no') self%theta_sq(i)
         end do
         write(fileid,*)
         write(fileid,*) '# self%periodicity(self%n_dof)'
         do i=1, self%n_dof
            write(fileid,100,advance='no') self%periodicity(i)
         end do
         write(fileid,*)
   end select
   if(self%m_f .ne. 0) then
      write(fileid,*) '# self%f_r(self%n_dof,self%m_f)'
      do j=1, self%m_f
         do i=1, self%n_dof
            write(fileid,100,advance='no') self%f_r(i,j)
         end do
         write(fileid,*)
      end do
   end if
   if(self%m_g .ne. 0) then
      write(fileid,*) '# self%g_r(self%n_dof,self%m_g)'
      do j=1, self%m_g
         do i=1, self%n_dof
            write(fileid,100,advance='no') self%g_r(i,j)
         end do
         write(fileid,*)
      end do
   end if
   ! compact mode
   if(trim(do_mode) .eq. 'compact') then
      write(fileid,*) '# self%Cmat_inv_v(self%m_teach)'
      do i=1, self%m_teach
         write(fileid,100,advance='no') self%Cmat_inv_v(i)
      end do
      write(fileid,*)
      write(fileid,*) '# self%logL'
      write(fileid,100) self%logL
   ! partial mode
   else
      if(self%sparsified) then
         write(fileid,*) '# self%n_tot'
         write(fileid,*) self%n_tot
         write(fileid,*) '# self%jitter'
         write(fileid,100) self%jitter
         write(fileid,*) '# self%y_siginvsq_y'
         write(fileid,100) self%y_siginvsq_y
         write(fileid,*) '# self%Kmn_siginvsq_y(self%m_teach)'
         do i=1, self%m_teach
            write(fileid,100,advance='no') self%Kmn_siginvsq_y(i)
         end do
         write(fileid,*)
         write(fileid,*) '# self%Kmn_siginvsq_Knm(self%m_teach,self%m_teach)'
         do j=1, self%m_teach
            do i=1, self%m_teach 
               write(fileid,100,advance='no') self%Kmn_siginvsq_Knm(i,j)
            end do
            write(fileid,*)
         end do
      else
      end if
   end if
   ! full mode
   if(trim(do_mode) .eq. 'full') then
      if(self%sparsified) then
         write(fileid,*) '# self%Cmat_inv_v(self%m_teach)'
         do i=1, self%m_teach
            write(fileid,100,advance='no') self%Cmat_inv_v(i)
         end do
         write(fileid,*)
         write(fileid,*) '# self%Kmm%n'
         write(fileid,*) self%Kmm%n
         write(fileid,*) '# self%Kmm%m'
         write(fileid,*) self%Kmm%m
         write(fileid,*) '# self%Kmm%matrix(self%Kmm%n,self%Kmm%m)'
         do j=1, self%Kmm%m
            do i=1, self%Kmm%n
               write(fileid,100,advance='no') self%Kmm%matrix(i,j)
            end do
            write(fileid,*)
         end do
         write(fileid,*) '# self%Cmat%n'
         write(fileid,*) self%Cmat%n
         write(fileid,*) '# self%Cmat%m'
         write(fileid,*) self%Cmat%m
         write(fileid,*) '# self%Cmat%matrix(self%Cmat%n,self%Cmat%m)'
         do j=1, self%Cmat%m
            do i=1, self%Cmat%n
               write(fileid,100,advance='no') self%Cmat%matrix(i,j)
            end do
            write(fileid,*)
         end do
         write(fileid,*) '# self%logL'
         write(fileid,100) self%logL
      else
         write(fileid,*) '# self%Cmat_inv_v(self%m_teach)'
         do i=1, self%m_teach
            write(fileid,100,advance='no') self%Cmat_inv_v(i)
         end do
         write(fileid,*)
         write(fileid,*) '# self%noise_Kmm%n'
         write(fileid,*) self%noise_Kmm%n
         write(fileid,*) '# self%noise_Kmm%m'
         write(fileid,*) self%noise_Kmm%m
         write(fileid,*) '# self%noise_Kmm%matrix(self%noise_Kmm%n,self%noise_Kmm%m)'
         do j=1, self%noise_Kmm%m
            do i=1, self%noise_Kmm%n
               write(fileid,100,advance='no') self%noise_Kmm%matrix(i,j)
            end do
         end do
         write(fileid,*)
         write(fileid,*) '# self%logL'
         write(fileid,100) self%logL
      end if
   end if

   close(fileid)

#ifndef QP
100 format(1X,E24.16E3)
#else
100 format(1X,E43.34E4)
#endif 

end subroutine gp_basic_write

subroutine gp_basic_read(self, fileid, filename, factorization)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   integer, intent(in) :: fileid !% fileid
   character(len=*), intent(in) :: filename !% filename to read GP from
   integer, optional, intent(in) :: factorization

   integer :: factor
   integer :: n, m
   real(dp), allocatable :: matrix(:,:)

   integer :: i, j

   if (present(factorization)) then
      if(factorization /= CH .and. factorization /= QR .and. factorization /= BK) then
        print *,"unknown factorization, factorization= ", factorization
        stop
      end if
      factor = factorization
   else
      factor = CH
   end if

   open(fileid,file=filename)

   call finalise(self)

   self%initialised=.true.

   ! obligatory variables
   read(fileid,*)
   read(fileid,*) self%mode
   read(fileid,*)
   read(fileid,*) self%n_dof
   read(fileid,*)
   read(fileid,*) self%m_f 
   read(fileid,*)
   read(fileid,*) self%m_g
   read(fileid,*)
   read(fileid,*) self%m_teach
   read(fileid,*)
   read(fileid,*) self%sparsified
   if(self%sparsified) then
      read(fileid,*)
      read(fileid,*) self%sparsemethod
   end if
   read(fileid,*)
   read(fileid,*) self%kernel_type
   select case(self%kernel_type)
      case(USER_kernel)
         read(fileid,*)
         read(fileid,*) self%kernel_string
         read(fileid,*)
         allocate(self%var_names(2*self%n_dof))
         do i=1, 2*self%n_dof
            read(fileid,*) self%var_names(i)
         end do
         call initf(self%fcomp, 1)
         call parsef(self%fcomp, 1, self%kernel_string, self%var_names)
      case(SE_kernel) 
         read(fileid,*)
         read(fileid,100) self%delta_sq
         read(fileid,*)
         allocate(self%theta_sq(self%n_dof))
         do i=1, self%n_dof
            read(fileid,100,advance='no') self%theta_sq(i)
         end do
         read(fileid,*)
         read(fileid,*)
         allocate(self%periodicity(self%n_dof))
         do i=1, self%n_dof
            read(fileid,100,advance='no') self%periodicity(i)
         end do
         read(fileid,*)
   end select
   if(self%m_f .ne. 0) then
      allocate(self%f_r(self%n_dof,self%m_f))
      read(fileid,*)
      do j=1, self%m_f
         do i=1, self%n_dof
            read(fileid,100,advance='no') self%f_r(i,j)
         end do
         read(fileid,*)
      end do
   endif
   if(self%m_g .ne. 0) then
      allocate(self%g_r(self%n_dof,self%m_g))
      read(fileid,*)
      do j=1, self%m_g
         do i=1, self%n_dof
            read(fileid,100,advance='no') self%g_r(i,j)
         end do
         read(fileid,*)
      end do
   endif
   ! compact mode
   if(trim(self%mode) .eq. 'compact') then
      allocate(self%Cmat_inv_v(self%m_teach))
      read(fileid,*)
      do i=1, self%m_teach
         read(fileid,100,advance='no') self%Cmat_inv_v(i)
      end do
      read(fileid,*)
      allocate(self%k(self%m_teach))
      allocate(self%k_grad(self%n_dof, self%m_teach))
      allocate(self%k_hess(self%n_dof, self%n_dof, self%m_teach))
      allocate(self%mat_inv_k(self%m_teach))
      allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))
      read(fileid,*)
      read(fileid,100) self%logL
   ! partial mode
   else
      if(self%sparsified) then
         read(fileid,*)
         read(fileid,*) self%n_tot
         read(fileid,*)
         read(fileid,100) self%jitter
         read(fileid,*)
         read(fileid,100) self%y_siginvsq_y
         allocate(self%Kmn_siginvsq_y(self%m_teach))
         read(fileid,*)
         do i=1, self%m_teach
            read(fileid,100,advance='no') self%Kmn_siginvsq_y(i)
         end do
         read(fileid,*)
         allocate(self%Kmn_siginvsq_Knm(self%m_teach,self%m_teach))
         read(fileid,*)
         do j=1, self%m_teach
            do i=1, self%m_teach
               read(fileid,100,advance='no') self%Kmn_siginvsq_Knm(i,j)
            end do
            read(fileid,*)
         end do
      else
      end if
   end if
   ! full mode
   if(trim(self%mode) .eq. 'full') then
      if(self%sparsified) then
         allocate(self%Cmat_inv_v(self%m_teach))
         read(fileid,*)
         do i=1, self%m_teach
            read(fileid,100,advance='no') self%Cmat_inv_v(i)
         end do
         read(fileid,*)
         allocate(self%k(self%m_teach))
         allocate(self%k_grad(self%n_dof, self%m_teach))
         allocate(self%k_hess(self%n_dof, self%n_dof, self%m_teach))
         allocate(self%mat_inv_k(self%m_teach))
         allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))
         read(fileid,*)
         read(fileid,*) n
         read(fileid,*)
         read(fileid,*) m
         allocate(matrix(n,m))
         read(fileid,*)
         do j=1, m
            do i=1, n
               read(fileid,100,advance='no') matrix(i,j)
            end do
            read(fileid,*)
         end do
         call gp_matrix_initialise(self%Kmm, matrix, factor)
         deallocate(matrix)
         read(fileid,*)
         read(fileid,*) n
         read(fileid,*)
         read(fileid,*) m
         allocate(matrix(n,m))
         read(fileid,*)
         do j=1, m
            do i=1, n
               read(fileid,100,advance='no') matrix(i,j)
            end do
            read(fileid,*)
         end do         
         call gp_matrix_initialise(self%Cmat, matrix, factor)
         deallocate(matrix)
         read(fileid,*)
         read(fileid,100) self%logL
      else
         allocate(self%Cmat_inv_v(self%m_teach))
         read(fileid,*)
         do i=1, self%m_teach
            read(fileid,100,advance='no') self%Cmat_inv_v(i)
         end do
         read(fileid,*) 
         allocate(self%k(self%m_teach))
         allocate(self%k_grad(self%n_dof, self%m_teach))
         allocate(self%k_hess(self%n_dof, self%n_dof, self%m_teach))
         allocate(self%mat_inv_k(self%m_teach))
         allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))
         read(fileid,*)
         read(fileid,*) n
         read(fileid,*)
         read(fileid,*) m
         allocate(matrix(n,m))
         read(fileid,*)
         do j=1, m
            do i=1, n
               read(fileid,100,advance='no') matrix(i,j)
            end do
            read(fileid,*)
         end do
         call gp_matrix_initialise(self%noise_Kmm, matrix, factor)
         deallocate(matrix)
         read(fileid,*)
         read(fileid,100) self%logL
      end if
   end if

   close(fileid)

#ifndef QP
100 format(1X,E24.16E3)
#else
100 format(1X,E43.34E4)
#endif

end subroutine gp_basic_read

! this subroutine changes the position of sparse points
subroutine gp_basic_change(self, f_set, g_set, factorization, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), optional, intent(in) :: f_set(:,:)
   real(dp), optional, intent(in) :: g_set(:,:)
   integer, optional, intent(in) :: factorization
   logical, optional, intent(in) :: verbose

   integer :: factor
   integer :: i
   type(gp_basic) :: newself
   logical :: do_verbose
   real(dp), allocatable :: Kll(:,:), Klm(:,:), Qmm(:,:)
   real(dp), allocatable :: Kmm_inv_Kml(:,:), Qmm_Kmm_inv_Kml(:,:)

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if ( .not. self%initialised ) then
      return
   end if

   if ( .not. self%sparsified ) then
      print *,"self is not sparsified!"
      stop
   end if

   if (present(f_set)) then
      if (size(f_set,1) .ne. self%n_dof) then
         print *,"size(f_set,1)=",size(f_set,1)," is not equal to self%n_dof=",self%n_dof
      end if
   end if

   if (present(g_set)) then
      if (size(g_set,1) .ne. self%n_dof) then
         print *,"size(g_set,1)=",size(g_set,1)," is not equal to self%n_dof=",self%n_dof
      end if
   end if

   if (present(factorization)) then
      if(factorization /= CH .and. factorization /= QR .and. factorization /= BK) then
        print *,"unknown factorization, factorization= ", factorization
        stop
      end if
      factor = factorization
   else
      factor = CH
   end if

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   call finalise(newself)

   newself%initialised = .true.
   newself%sparsified = .true.
   newself%mode = 'partial'
   newself%sparsemethod = self%sparsemethod

   newself%kernel_type = self%kernel_type

   newself%n_dof = self%n_dof
   if(present(f_set)) newself%m_f = size(f_set,2)
   if(present(g_set)) newself%m_g = size(g_set,2)
   if (newself%m_f + newself%m_g == 0) then
      print *,"no new set of sparse points provided"
      stop
   end if
   newself%m_teach = newself%m_f + newself%m_g*newself%n_dof
   newself%n_tot = self%n_tot

   newself%jitter = self%jitter

   select case(newself%kernel_type)
      case(USER_kernel)
         newself%kernel_string = self%kernel_string
         allocate(newself%var_names(size(self%var_names)))
         newself%var_names = self%var_names
         call initf(newself%fcomp, 1)
         call parsef(newself%fcomp, 1, newself%kernel_string, newself%var_names)
      case(SE_kernel)
         newself%delta_sq = self%delta_sq
         allocate(newself%theta_sq(newself%n_dof))
         newself%theta_sq = self%theta_sq
         allocate(newself%periodicity(newself%n_dof))
         newself%periodicity = self%periodicity
   end select

   if (newself%m_f .ne. 0) then
      allocate(newself%f_r(newself%n_dof,newself%m_f))
      newself%f_r = f_set
   end if
   if (newself%m_g .ne. 0) then
      allocate(newself%g_r(newself%n_dof,newself%m_g))
      newself%g_r = g_set
   end if

   newself%y_siginvsq_y = self%y_siginvsq_y

   ! complete self if necessary
   if(do_verbose) print *,"Completion of original teching has started..."
   call gp_basic_complete(self, factor, do_verbose)
   if(do_verbose) print *,"Completion of original teching has finished..."

   ! calculate Kll
   allocate(Kll(newself%m_teach, newself%m_teach))
   if(do_verbose) print *,"Calculation of Kll matrix has started..."
   call kernel_mat(self, newself%m_f, newself%f_r, newself%m_g, newself%g_r, &
                   newself%m_f, newself%f_r, newself%m_g, newself%g_r, Kll)
   if(do_verbose) print *,"Calculation of Kll matrix has finished..."
   ! "jitter" (ABP e-mail 14 Aug)
!$omp parallel do
   do i=1, size(Kll,1)
      Kll(i,i) = Kll(i,i) + newself%jitter
   end do
!$omp end parallel do
   call gp_matrix_initialise(newself%Kmm, Kll, factor)

   ! calculate Klm
   allocate(Klm(newself%m_teach, self%m_teach))
   if(do_verbose) print *,"Calculation of Klm matrix has started..."
   call kernel_mat(newself, newself%m_f, newself%f_r, newself%m_g, newself%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, Klm)
   if(do_verbose) print *,"Calculation of Klm matrix has finished..."

   ! calculate Kmm_inv_Kml
   allocate(Kmm_inv_Kml(self%m_teach, newself%m_teach))
   if(do_verbose) print *,"Calculation of Kmm_inv_Kml matrix has started..."
   call Matrix_Solve(self%Kmm, transpose(Klm), Kmm_inv_Kml)
   if(do_verbose) print *,"Calculation of Kmm_inv_Kml matrix has finished..."

   ! calculation of new Kmn_siginvsq_y
   allocate(newself%Kmn_siginvsq_y(newself%m_teach))
   newself%Kmn_siginvsq_y = 0.0_dp
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_y has started..."
   call dgemv('T', size(Kmm_inv_Kml,1), size(Kmm_inv_Kml,2), 1.0_dp, Kmm_inv_Kml, size(Kmm_inv_Kml,1), &
              self%Kmn_siginvsq_y, 1, 0.0_dp, newself%Kmn_siginvsq_y, 1)
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_y has finished..."

   ! calculate Qmm
   allocate(Qmm(self%m_teach,self%m_teach))
!  Qmm = self%Kmm%matrix + self%Kmn_siginvsq_Knm !!! LAM why?
   Qmm = self%Kmn_siginvsq_Knm

   ! calculate Qmm_Kmm_inv_Kml
   allocate(Qmm_Kmm_inv_Kml(self%m_teach, newself%m_teach))
   Qmm_Kmm_inv_Kml = 0.0_dp
   if(do_verbose) print *,"Calculation of Qmm_Kmm_inv_Kml matrix has started..."
   call dgemm('N', 'N', size(Qmm,1), size(Kmm_inv_Kml,2), size(Qmm,2), 1.0_dp, &
              Qmm, size(Qmm,1), Kmm_inv_Kml, size(Qmm,2), 0.0_dp, &
              Qmm_Kmm_inv_Kml, size(Qmm,1))
   if(do_verbose) print *,"Calculation of Qmm_Kmm_inv_Kml matrix has finished..."

   ! calculation of new Kmn_siginvsq_Knm
   allocate(newself%Kmn_siginvsq_Knm(newself%m_teach, newself%m_teach))
   newself%Kmn_siginvsq_Knm = 0.0_dp
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_Knm has started..."
   call dgemm('T', 'N', size(Kmm_inv_Kml,2), size(Qmm_Kmm_inv_Kml,2), size(Kmm_inv_Kml,1), 1.0_dp, &
              Kmm_inv_Kml, size(Kmm_inv_Kml,1), Qmm_Kmm_inv_Kml, size(Kmm_inv_Kml,1), 0.0_dp, &
              newself%Kmn_siginvsq_Knm, size(Kmm_inv_Kml,2))
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_Knm has finished..."

   deallocate(Kll)
   deallocate(Klm)
   deallocate(Kmm_inv_Kml)
   deallocate(Qmm)
   deallocate(Qmm_Kmm_inv_Kml)

   self = newself

   call finalise(newself)

end subroutine gp_basic_change

! this subroutine merges two teachings that have the same teaching points
! if self3 is specified then the result is written into self3, otherwise into self2 
subroutine gp_basic_merge(self1, self2, self3)

   implicit none

   type(gp_basic), intent(inout) :: self1, self2
   type(gp_basic), optional, intent(inout) :: self3

   integer :: i, j

   if ( .not. self1%initialised ) then
      print *,"self1 is not initialized!"
      stop
   end if

   if ( .not. self2%initialised ) then
      print *,"self2 is not initialized!"
      stop
   end if

   ! check compatibility of self1 and self2
   if (self1%kernel_type .ne. self2%kernel_type) then
      print *,"self1%kernel_type=", self1%kernel_type, " is not equal to self2%kernel_type=", self2%kernel_type
      stop
   end if

   if (self1%n_dof .ne. self2%n_dof) then
      print *,"self1%n_dof=", self1%n_dof, " is not equal to self2%n_dof=", self2%n_dof
      stop
   end if

   if (self1%m_f .ne. self2%m_f) then
      print *,"self1%m_f=", self1%m_f, " is not equal to self2%m_f=", self2%m_f
      stop
   end if

   if (self1%m_g .ne. self2%m_g) then
      print *,"self1%m_g=", self1%m_g, " is not equal to self2%m_g=", self2%m_g
      stop
   end if

   if (self1%m_teach .ne. self2%m_teach) then
      print *,"self1%m_teach=", self1%m_teach, " is not equal to self2%m_teach=", self2%m_teach
      stop
   end if

   if (self1%delta_sq .ne. self2%delta_sq) then
      print *,"self1%delta_sq=", self1%delta_sq, " is not equal to self2%delta_sq=", self2%delta_sq
      stop
   end if

   if (.not. self1%sparsified) then
      print *,"self1 is not sparsified!"
      stop
   end if

   if (.not. self2%sparsified) then
      print *,"self2 is not sparsified!"
      stop
   end if

   if (trim(self1%sparsemethod) .ne. trim(self2%sparsemethod)) then
      print *,"self1%sparsemethod=", trim(self1%sparsemethod), " is not equal to self2%sparsemethod=", trim(self1%sparsemethod)
      stop
   end if

   select case(self1%kernel_type)
      case(USER_kernel)
         if (self1%kernel_string .ne. self2%kernel_string) then
               print *,"self1%kernel_string=",self1%kernel_string, &
                       " is not equal to self2%kernel_string=",self2%kernel_string
               stop
         end if

         do i=1,2*self1%n_dof
            if (self1%var_names(i) .ne. self2%var_names(i)) then
               print *,"self1%var_names(",i,")=",self1%var_names(i), &
                       " is not equal to self2%var_names(",i,")=",self2%var_names(i)
               stop  
            end if
         end do

      case(SE_kernel)
         do i=1,self1%n_dof
            if (self1%theta_sq(i) .ne. self2%theta_sq(i)) then
               print *,"self1%theta_sq(",i,")=",self1%theta_sq(i), &
                       " is not equal to self2%theta_sq(",i,")=",self2%theta_sq(i)
               stop  
            end if
         end do

         do i=1,self1%n_dof
            if (self1%periodicity(i) .ne. self2%periodicity(i)) then
               print *,"self1%periodicity(",i,")=",self1%periodicity(i), &
                       " is not equal to self2%periodicity(",i,")=",self2%periodicity(i)
               stop
            end if
         end do

   end select

   if (self1%m_f .ne. 0) then
      do j=1,self1%m_f
         do i=1,self1%n_dof
            if (self1%f_r(i,j) .ne. self2%f_r(i,j)) then
               print *,"self1%f_r(",i,j,")=",self1%f_r(i,j)," is not equal to self2%f_r(",i,j,")=",self2%f_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if (self1%m_g .ne. 0) then
      do j=1,self1%m_g
         do i=1,self1%n_dof
            if (self1%g_r(i,j) .ne. self2%g_r(i,j)) then
               print *,"self1%g_r(",i,j,")=",self1%g_r(i,j)," is not equal to self2%g_r(",i,j,")=",self2%g_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if ( present(self3) ) then
      call finalise(self3)

      self3%n_dof = self1%n_dof
      self3%m_f = self1%m_f
      self3%m_g = self1%m_g 
      self3%m_teach = self1%m_teach
      self3%n_tot = self1%n_tot + self2%n_tot
      self3%sparsified = self1%sparsified
      self3%kernel_type = self1%kernel_type
      select case(self3%kernel_type)
         case(USER_kernel)
            self3%kernel_string = self1%kernel_string
            allocate(self3%var_names(size(self1%var_names)))
            self3%var_names = self1%var_names
            call initf(self3%fcomp, 1)
            call parsef(self3%fcomp, 1, self3%kernel_string, self3%var_names)
         case(SE_kernel)
            self3%delta_sq = self1%delta_sq
            allocate(self3%theta_sq(self3%n_dof))
            self3%theta_sq = self1%theta_sq
            allocate(self3%periodicity(self3%n_dof))
            self3%periodicity = self1%periodicity
      end select
      if (self3%m_f .ne. 0) then
         allocate(self3%f_r(self3%n_dof,self3%m_f))
         self3%f_r = self1%f_r
      end if
      if (self3%m_g .ne. 0) then
         allocate(self3%g_r(self3%n_dof,self3%m_g))
         self3%g_r = self1%g_r
      end if

      self3%mode = 'partial'
      self3%sparsemethod = self1%sparsemethod
      self3%y_siginvsq_y = self1%y_siginvsq_y + self2%y_siginvsq_y
      allocate(self3%Kmn_siginvsq_y(self3%m_teach))
      self3%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
      allocate(self3%Kmn_siginvsq_Knm(self3%m_teach,self3%m_teach))
      self3%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm

      self3%initialised = .true.
   else
      self2%n_tot = self1%n_tot + self2%n_tot

      self2%mode = 'partial'
      self2%sparsemethod = self1%sparsemethod
      self2%y_siginvsq_y = self1%y_siginvsq_y + self2%y_siginvsq_y
      self2%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
      self2%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm
   end if
    
end subroutine gp_basic_merge 

! this subroutine completes a teaching (i.e. calculates Kmm and Cmat) if gp_basic is not complete
subroutine gp_basic_complete(self, factorization, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self
   integer, optional, intent(in) :: factorization
   logical, optional, intent(in) :: verbose

   integer :: factor
   logical :: do_verbose 
   integer :: i
   real(dp), allocatable :: Kmm(:,:), Cmat(:,:)
   real(dp) :: logdet
   real(dp) :: maxdiag, sjitter

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      print *,"self is not initialised!"
      stop
   end if

   if (trim(self%mode) .ne. 'partial') return

   if (present(factorization)) then                               
      if(factorization /= CH .and. factorization /= QR .and. factorization /= BK) then
        print *,"unknown factorization, factorization= ", factorization
        stop
      end if
      factor = factorization
   else
      factor = CH
   end if

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   allocate(Kmm(self%m_teach, self%m_teach))
   if(do_verbose) print *,"Calculation of Kmm has started..."
   call kernel_mat(self, self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, Kmm)
   if(do_verbose) print *,"Calculation of Kmm has finished..."

   call gp_matrix_initialise(self%Kmm, Kmm, factor)
   
   allocate(Cmat(self%m_teach, self%m_teach))

!$omp parallel workshare
   Cmat = self%Kmm%matrix + self%Kmn_siginvsq_Knm
!$omp end parallel workshare
   ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})
   maxdiag = Cmat(1,1)
   do i=2,self%m_teach
      if(Cmat(i,i) .gt. maxdiag) maxdiag = Cmat(i,i)
   end do
   sjitter = self%jitter*10.0_dp**(int(log10(maxdiag))+1)
!$omp parallel do
   do i=1,self%m_teach
      Cmat(i,i) = Cmat(i,i) + sjitter
   end do
!$omp end parallel do

   call gp_matrix_initialise(self%Cmat, Cmat, factor)

   if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
   if(allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
   allocate(self%Cmat_inv_v(self%m_teach))
   call Matrix_Solve(self%Cmat, self%Kmn_siginvsq_y, self%Cmat_inv_v)
   if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."

   if(do_verbose) print *,"Calculation of logL has started..."
   ! logL = -0.5*log(det(C)) - 0.5*y^{T} C^{-1} y - N/2*log(2*pi), where C = \Sigma^{2} + K_{nm} K_{mm}^{-1} K_{mn}
   ! -0.5*log(det(C)) = -0.5*log(R_{11}*R_{22}*...*R_{mm})
   logdet = 0.0_dp
   do i=1, self%Cmat%m
      logdet = logdet + log(abs(self%Cmat%factor(i,i)))
   end do
   self%logL = -0.5_dp * logdet

   ! -0.5*y^{T} C^{-1} y = -0.5*(y^{T} \Sigma^{-2} y - y^{T} \Sigma^{-2} K_{nm} Cmat^{-1}  K_{mn} \Sigma^{-2} y)
   ! where Cmat = K_{mm} + K_{mn} \Sigma^{-2} K_{nm}
   self%logL = self%logL - 0.5_dp * (self%y_siginvsq_y - ddot(self%m_teach, self%Kmn_siginvsq_y, 1, self%Cmat_inv_v, 1))

   ! - N/2*log(2*pi)
   self%logL = self%logL - 0.5_dp * real(self%n_tot) * log(TWO_PI)
   if(do_verbose) print *,"Calculation of logL has finished..."

   self%mode = 'full'

   deallocate(Kmm)
   deallocate(Cmat)

   if(.not. allocated(self%k)) allocate(self%k(self%m_teach))
   if(.not. allocated(self%k_grad)) allocate(self%k_grad(self%n_dof, self%m_teach))
   if(.not. allocated(self%k_hess)) allocate(self%k_hess(self%n_dof, self%n_dof, self%m_teach))
   if(.not. allocated(self%mat_inv_k)) allocate(self%mat_inv_k(self%m_teach))
   if(.not. allocated(self%mat_inv_k_grad)) allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))

end subroutine gp_basic_complete

function f_predict_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to predict value
   real(dp) :: f_predict_r

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   end if

   call f_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k)
   f_predict_r = ddot(size(self%k), self%k, 1, self%Cmat_inv_v, 1)

end function f_predict_r

function f_predict_grad_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to predict gradient
   real(dp) :: f_predict_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_grad_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   end if

   call g_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k_grad)
   call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 1.0_dp, self%k_grad, size(self%k_grad,1), &
      self%Cmat_inv_v, 1, 0.0_dp, f_predict_grad_r, 1)

end function f_predict_grad_r

function f_predict_hess_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to predict hessian
   real(dp) :: f_predict_hess_r(size(r), size(r))
   integer :: i, j

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_hess_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   end if

   call h_kernel_mat(self, r, self%m_f, self%f_r, self%k_hess)
   do i=1, size(r)
      f_predict_hess_r(i,i) = ddot(self%m_f, self%k_hess(i,i,:), 1, self%Cmat_inv_v, 1)
      do j=i+1, size(r)
         f_predict_hess_r(i,j) = ddot(self%m_f, self%k_hess(i,j,:), 1, self%Cmat_inv_v, 1)
         f_predict_hess_r(j,i) = f_predict_hess_r(i,j)
      end do
   end do

end function f_predict_hess_r

function f_predict_int_r(self, r, interval)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp), intent(in) :: interval(:,:)
   real(dp) :: f_predict_int_r

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_int_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   end if

   call e_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, interval, self%k)
   f_predict_int_r = ddot(size(self%k), self%k, 1, self%Cmat_inv_v, 1)

end function f_predict_int_r

function f_predict_var_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: f_predict_var_r

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   else if (trim(self%mode) .eq. 'compact') then
      f_predict_var_r = 0.0_dp
      return
   end if

   call f_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k)

   if (self%sparsified) then
      if (trim(self%sparsemethod) .eq. 'dic') then
          f_predict_var_r = 0.0_dp
      else
!         f_predict_var_r = self%delta_sq
          f_predict_var_r = f_kernel_sca(self, r, r)
          call Matrix_Solve(self%Kmm, self%k, self%mat_inv_k)
          f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
      end if
      call Matrix_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r + ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
   else
!     f_predict_var_r = self%delta_sq
      f_predict_var_r = f_kernel_sca(self, r, r)
      call Matrix_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
   endif

end function f_predict_var_r

function f_predict_grad_var_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: f_predict_grad_var_r(size(r))

   integer :: ii

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_grad_var_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   else if (trim(self%mode) .eq. 'compact') then
      f_predict_grad_var_r = 0.0_dp
      return
   end if

   ! From Lockwood and Anitescu preprint ANL/MCS-P1808-1110 "Gradient-Enhanced Universal Kriging for Uncertainty Propagation"
   ! Eq. 2.22, not including last term which is relevant only for underlying polynomial basis (which we don't have)
   call g_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k_grad)

   if (self%sparsified) then
      if (trim(self%sparsemethod) .eq. 'dic') then
          f_predict_grad_var_r = 0.0_dp
      else
          f_predict_grad_var_r = g_kernel_sca(self, r, r)
          ! -k' K_{mm}^{-1} k'
          call Matrix_Solve(self%Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
          do ii=1, self%n_dof
             f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), &
                                                                        size(self%k_grad,1), self%mat_inv_k_grad(1,ii), 1)
          end do
      end if
      ! +k' C^{-1} k'
      call Matrix_Solve(self%Cmat, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
         f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) + ddot(size(self%k_grad,2), self%k_grad(ii,1), &
                                                                    size(self%k_grad,1), self%mat_inv_k_grad(1,ii), 1)
      end do
   else
      f_predict_grad_var_r = g_kernel_sca(self, r, r)
      call g_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k_grad)
      ! -k' C^{-1} k'
      call Matrix_Solve(self%noise_Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
         f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), &
                                                                    size(self%k_grad,1), self%mat_inv_k_grad(:,ii), 1)
      end do
   endif

end function f_predict_grad_var_r

function f_predict_int_var_r(self, r, interval)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp), intent(in) :: interval(:,:)
   real(dp) :: f_predict_int_var_r

#ifdef EXTERNAL_BLAS_LAPACK
   real(dp), external :: ddot
#endif

   if (.not. self%initialised) then
      f_predict_int_var_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   else if (trim(self%mode) .eq. 'compact') then
      f_predict_int_var_r = 0.0_dp
      return
   end if

   call e_kernel_vec(self, r , self%m_f, self%f_r, self%m_g, self%g_r, interval, self%k)

   if (self%sparsified) then
      if (trim(self%sparsemethod) .eq. 'dic') then
          f_predict_int_var_r = 0.0_dp
      else
          f_predict_int_var_r = e_kernel_sca(self, r, r, interval)
          call Matrix_Solve(self%Kmm, self%k, self%mat_inv_k)
          f_predict_int_var_r = f_predict_int_var_r - ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
      end if
      call Matrix_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_int_var_r = f_predict_int_var_r + ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
   else
      f_predict_int_var_r = e_kernel_sca(self, r, r, interval)
      call Matrix_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_int_var_r = f_predict_int_var_r - ddot(size(self%k), self%k, 1, self%mat_inv_k, 1)
   endif

end function f_predict_int_var_r

function f_predict_var_grad_r(self, r)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: f_predict_var_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_var_grad_r = 0.0_dp
      return
   endif

   if (trim(self%mode) .eq. 'partial') then
      call gp_basic_complete(self)
   else if (trim(self%mode) .eq. 'compact') then
      f_predict_var_grad_r = 0.0_dp
      return 
   end if

   call f_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k)
   call g_kernel_vec(self, r, self%m_f, self%f_r, self%m_g, self%g_r, self%k_grad)

   f_predict_var_grad_r = 0.0_dp
   if (self%sparsified) then
      if (trim(self%sparsemethod) .eq. 'dic') then
          f_predict_var_grad_r = 0.0_dp
      else
          call Matrix_Solve(self%Kmm, self%k, self%mat_inv_k)
          call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad, size(self%k_grad,1), &
                     self%mat_inv_k, 1, 1.0_dp, f_predict_var_grad_r, 1)
      end if
      call Matrix_Solve(self%Cmat, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 2.0_dp, self%k_grad, size(self%k_grad,1), &
                 self%mat_inv_k, 1, 1.0_dp, f_predict_var_grad_r, 1)
   else
      call Matrix_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad, size(self%k_grad,1), &
                 self%mat_inv_k, 1, 1.0_dp, f_predict_var_grad_r, 1)
   endif

end function f_predict_var_grad_r

subroutine kernel_mat(self, l_n_f, l_f_r, l_n_g, l_g_r, r_n_f, r_f_r, r_n_g, r_g_r, mat)

   implicit none

   type(gp_basic), intent(in) :: self
   integer, intent(in) :: l_n_f, l_n_g
   real(dp), optional, intent(in) :: l_f_r(:,:), l_g_r(:,:)
   integer, intent(in) :: r_n_f, r_n_g
   real(dp), optional, intent(in) :: r_f_r(:,:), r_g_r(:,:)
   real(dp), intent(out) :: mat(:,:)

   integer :: i, i_glob, n_dof

   if (l_n_f > 0) n_dof=size(l_f_r,1)
   if (l_n_g > 0) n_dof=size(l_g_r,1)

   select case(self%kernel_type)
      case(USER_kernel)
!$omp parallel do
        do i=1, l_n_f
           ! f on f
           if (r_n_f > 0) call USER_kernel_r_rr(self, l_f_r(1:n_dof,i), r_f_r(1:n_dof,:), f_f=mat(i,1:r_n_f))
           ! f on g
           if (r_n_g > 0) call USER_kernel_r_rr(self, l_f_r(1:n_dof,i), r_g_r(1:n_dof,:), f_g=mat(i,r_n_f+1:r_n_f+r_n_g*n_dof))
        end do
!$omp end parallel do
!$omp parallel do private(i_glob)
        do i=1, l_n_g
           i_glob = l_n_f + (i-1)*n_dof + 1
           ! g on f
           if (r_n_f > 0) call USER_kernel_r_rr(self, l_g_r(1:n_dof,i), r_f_r(1:n_dof,:), g_f=mat(i_glob:i_glob+n_dof-1,1:r_n_f))
           ! g on g
           if (r_n_g > 0) call USER_kernel_r_rr(self, l_g_r(1:n_dof,i), r_g_r(1:n_dof,:), &
                                                                        g_g=mat(i_glob:i_glob+n_dof-1,r_n_f+1:r_n_f+r_n_g*n_dof))
        end do
!$omp end parallel do
      case(SE_kernel)
!$omp parallel do
        do i=1, l_n_f
           ! f on f
           if (r_n_f > 0) call SE_kernel_r_rr(self, l_f_r(1:n_dof,i), r_f_r(1:n_dof,:), f_f=mat(i,1:r_n_f))
           ! f on g
           if (r_n_g > 0) call SE_kernel_r_rr(self, l_f_r(1:n_dof,i), r_g_r(1:n_dof,:), f_g=mat(i,r_n_f+1:r_n_f+r_n_g*n_dof))
        end do
!$omp end parallel do
!$omp parallel do private(i_glob)
        do i=1, l_n_g
           i_glob = l_n_f + (i-1)*n_dof + 1
           ! g on f
           if (r_n_f > 0) call SE_kernel_r_rr(self, l_g_r(1:n_dof,i), r_f_r(1:n_dof,:), g_f=mat(i_glob:i_glob+n_dof-1,1:r_n_f))
           ! g on g
           if (r_n_g > 0) call SE_kernel_r_rr(self, l_g_r(1:n_dof,i), r_g_r(1:n_dof,:), &
                                                                      g_g=mat(i_glob:i_glob+n_dof-1,r_n_f+1:r_n_f+r_n_g*n_dof))
        end do
!$omp end parallel do 
      end select

end subroutine kernel_mat

subroutine h_kernel_mat(self, r, n_f, f_r, mat)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   real(dp), intent(out) :: mat(:,:,:)

   integer :: n_dof

   n_dof=size(r)

   select case(self%kernel_type)
      case(USER_kernel)
        if (n_f > 0) call USER_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), h_f=mat(1:n_dof,1:n_dof,1:n_f))
      case(SE_kernel)
        if (n_f > 0) call SE_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), h_f=mat(1:n_dof,1:n_dof,1:n_f))
   end select

end subroutine

subroutine e_kernel_vec(self, r, n_f, f_r, n_g, g_r, interval, vec)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: interval(:,:)
   real(dp), intent(out) :: vec(:)

   integer :: n_dof

   n_dof=size(r)

   select case(self%kernel_type)
      case(USER_kernel)
        print *,"integral computation for USER kernel is not implemented!"
        stop
      case(SE_kernel)
        if (n_f > 0) call SE_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), interval, e_f=vec(1:n_f))
        if (n_g > 0) call SE_kernel_r_rr(self, r(1:n_dof), g_r(1:n_dof,1:n_g), interval, e_g=vec(n_f+1:n_f+n_g*n_dof))
   end select

end subroutine e_kernel_vec

subroutine f_kernel_vec(self, r, n_f, f_r, n_g, g_r, vec)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(out) :: vec(:)

   integer :: n_dof

   n_dof=size(r)

   select case(self%kernel_type)
      case(USER_kernel)
        if (n_f > 0) call USER_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), f_f=vec(1:n_f))
        if (n_g > 0) call USER_kernel_r_rr(self, r(1:n_dof), g_r(1:n_dof,1:n_g), f_g=vec(n_f+1:n_f+n_g*n_dof))
      case(SE_kernel)
        if (n_f > 0) call SE_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), f_f=vec(1:n_f))
        if (n_g > 0) call SE_kernel_r_rr(self, r(1:n_dof), g_r(1:n_dof,1:n_g), f_g=vec(n_f+1:n_f+n_g*n_dof))
   end select

end subroutine f_kernel_vec

subroutine g_kernel_vec(self, r, n_f, f_r, n_g, g_r, vec)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(out) :: vec(:,:)

   integer n_dof

   n_dof = size(r)

   select case(self%kernel_type)
      case(USER_kernel)
        if (n_f > 0) call USER_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), g_f=vec(1:n_dof,1:n_f))
        if (n_g > 0) call USER_kernel_r_rr(self, r(1:n_dof), g_r(1:n_dof,1:n_g), g_g=vec(1:n_dof,n_f+1:n_f+n_g*n_dof))
      case(SE_kernel)
        if (n_f > 0) call SE_kernel_r_rr(self, r(1:n_dof), f_r(1:n_dof,1:n_f), g_f=vec(1:n_dof,1:n_f))
        if (n_g > 0) call SE_kernel_r_rr(self, r(1:n_dof), g_r(1:n_dof,1:n_g), g_g=vec(1:n_dof,n_f+1:n_f+n_g*n_dof))
   end select

end subroutine g_kernel_vec

function e_kernel_sca(self, r, f_r, interval)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   real(dp), intent(in) :: f_r(:)
   real(dp), intent(in) :: interval(:,:)
   real(dp) :: e_kernel_sca

   real(dp) :: sca(1)
   integer n_dof

   n_dof = size(r)

   select case(self%kernel_type)
      case(USER_kernel)
      case(SE_kernel)
        call SE_kernel_r_rr(self, r(1:n_dof), reshape(f_r(1:n_dof), (/n_dof,1/)), interval, e_e=sca)
   end select

   e_kernel_sca = sca(1)

end function e_kernel_sca

function f_kernel_sca(self, r, f_r)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   real(dp), intent(in) :: f_r(:)
   real(dp) :: f_kernel_sca

   real(dp) :: sca(1)
   integer n_dof

   n_dof = size(r)

   select case(self%kernel_type)
      case(USER_kernel)
         call USER_kernel_r_rr(self, r(1:n_dof), reshape(f_r(1:n_dof), (/n_dof,1/)), f_f=sca)
      case(SE_kernel)
         call SE_kernel_r_rr(self, r(1:n_dof), reshape(f_r(1:n_dof), (/n_dof,1/)), f_f=sca)
   end select

   f_kernel_sca = sca(1)

end function f_kernel_sca

function g_kernel_sca(self, r, g_r)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: r(:)
   real(dp), intent(in) :: g_r(:)
   real(dp) :: g_kernel_sca(size(r))

   real(dp) :: sca(size(r))
   integer n_dof

   n_dof = size(r)

   select case(self%kernel_type)
      case(USER_kernel)
         call USER_kernel_r_rr(self, r(1:n_dof), reshape(g_r(1:n_dof), (/n_dof,1/)), g_g_diag=sca)
      case(SE_kernel)
         call SE_kernel_r_rr(self, r(1:n_dof), reshape(g_r(1:n_dof), (/n_dof,1/)), g_g_diag=sca)
   end select

   g_kernel_sca = sca

end function g_kernel_sca

subroutine USER_kernel_r_rr(self, x1, x2, f_f, f_g, g_f, g_g, g_g_diag, h_f)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: x1(:), x2(:,:)
   real(dp), optional, intent(out) :: f_f(:), f_g(:), g_f(:,:)
   real(dp), optional, intent(out) :: g_g(:,:), g_g_diag(:)
   real(dp), optional, intent(out) :: h_f(:,:,:)

   integer :: i, j, k, l, nv, n_dof

   real(dp), allocatable :: val(:)

   n_dof = size(x1)
   nv = size(x2,2)

   allocate(val(2*n_dof))

   if (present(f_f)) then
      val(1:n_dof) = x1(1:n_dof)
      do k=1, nv
         val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
         f_f(k) = evalf(self%fcomp, 1, val)
      end do
   end if

   if (present(f_g)) then
      val(1:n_dof) = x1(1:n_dof)
      l = 0
      do k=1, nv
         val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
         do i=1, n_dof
            l = l + 1
            f_g(l) = evald(self%fcomp, 1, n_dof+i, val)
         end do
      end do
   end if

   if (present(g_f)) then
      val(1:n_dof) = x1(1:n_dof)
      do k=1, nv
         val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
         do i=1, n_dof
            g_f(i,k) = evald(self%fcomp, 1, i, val)
         end do
      end do
   end if

   if (present(g_g)) then
      val(1:n_dof) = x1(1:n_dof)
      do i=1, n_dof
         l = 0
         do k=1, nv
            val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
            do j=1, n_dof
               l = l + 1
               g_g(i,l) = evaldd(self%fcomp, 1, i, n_dof+j, val)
            end do
         end do
      end do
   end if

   if (present(g_g_diag)) then
      val(1:n_dof) = x1(1:n_dof)
      l = 0
      do k=1, nv
         val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
         do i=1, n_dof
            l = l + 1
            g_g_diag(l) = evaldd(self%fcomp, 1, i, n_dof+i, val)
         end do
      end do
   end if

   if (present(h_f)) then
      val(1:n_dof) = x1(1:n_dof)
      do k=1, nv
         val(n_dof+1:2*n_dof) = x2(1:n_dof,k)
         do j=1, n_dof
            do i=1, n_dof
               h_f(i,j,k) = evaldd(self%fcomp, 1, i, j, val)
            end do
         end do
      end do
   end if

   deallocate(val)

end subroutine USER_kernel_r_rr

subroutine SE_kernel_r_rr(self, x1, x2, interval, e_e, e_f, f_e, e_g, g_e, f_f, f_g, g_f, g_g, g_g_diag, h_f)

   implicit none

   type(gp_basic), intent(in) :: self
   real(dp), intent(in) :: x1(:), x2(:,:)
   real(dp), optional, intent(in) :: interval(:,:)
   real(dp), optional, intent(out) :: e_e(:), e_f(:), f_e(:), e_g(:), g_e(:,:)
   real(dp), optional, intent(out) :: f_f(:), f_g(:), g_f(:,:)
   real(dp), optional, intent(out) :: g_g(:,:), g_g_diag(:)
   real(dp), optional, intent(out) :: h_f(:,:,:)

   real(dp), allocatable :: iexp_arg(:), iexp_argj(:), iexp_argi(:)
   real(dp), allocatable :: iexp_arg_j(:,:), iexp_arg_i(:,:)
   real(dp), allocatable :: idexp_arg_j(:,:), idexp_arg_i(:,:)
   real(dp), allocatable :: exp_arg(:)
   real(dp), allocatable :: dexp_arg_i(:,:), ddexp_arg_ij(:,:), dexp_arg_j(:,:)
   real(dp), allocatable :: ddexp_arg_ii(:,:)
   integer :: i, j, nv, n_dof

   real(dp), allocatable :: pi_per_periodicity(:), two_pi_per_periodicity(:)
   real(dp), allocatable :: diff(:)

   if( present(interval) .neqv. (present(e_e) .or. present(e_f) .or. present(f_e) .or. present(e_g) .or. present(g_e)) ) then
      print *,"interval and e_e/e_f/f_e/g_e/e_g must be present together"
      stop
   end if

   n_dof = size(x1)
   nv = size(x2,2)

   allocate(pi_per_periodicity(n_dof))
   allocate(two_pi_per_periodicity(n_dof))

   do i=1, n_dof
      if (self%periodicity(i) > 0.0_dp) then
         pi_per_periodicity(i) = PI/self%periodicity(i)
         two_pi_per_periodicity(i) = TWO_PI/self%periodicity(i)
      else
         pi_per_periodicity(i) = 0.0_dp
         two_pi_per_periodicity(i) = 0.0_dp
      end if
   end do

   allocate(diff(nv))

   if(present(e_e)) allocate(iexp_arg(nv))
   if(present(e_f) .or. present(e_g)) allocate(iexp_argj(nv))
   if(present(e_g)) allocate(iexp_arg_j(n_dof,nv))
   if(present(e_g)) allocate(idexp_arg_j(n_dof,nv))
   if(present(f_e) .or. present(g_e)) allocate(iexp_argi(nv))
   if(present(g_e)) allocate(iexp_arg_i(n_dof,nv))
   if(present(g_e)) allocate(idexp_arg_i(n_dof,nv))
   if(present(f_f) .or. present(g_f) .or. present(f_g) .or. &
      present(g_g) .or. present(g_g_diag) .or. present(h_f)) allocate(exp_arg(nv))
   if(present(g_f) .or. present(g_g) .or. present(g_g_diag) .or. present(h_f)) allocate(dexp_arg_i(n_dof,nv))
   if(present(g_g) .or. present(g_g_diag)) allocate(ddexp_arg_ij(n_dof,nv))
   if(present(f_g) .or. present(g_g)) allocate(dexp_arg_j(n_dof,nv))
   if(present(h_f)) allocate(ddexp_arg_ii(n_dof,nv))

   if (present(e_e)) then
      iexp_arg = 0.0_dp
      do i=1, n_dof
         if (interval(i,1) /= interval(i,2)) then
            if (self%periodicity(i) > 0.0_dp) then
               diff(:) = log(self%periodicity(i)**2*bessi0(1.0_dp/self%theta_sq(i))) - 1.0_dp/self%theta_sq(i)
            else if (self%periodicity(i) < 0.0_dp) then ! LAM
               diff(:) = log( sqrt(two_pi*self%theta_sq(i)) * &
                             ( sqrt(self%theta_sq(i)/pi_per_two)* &
                               (exp(-0.5_dp*(interval(i,2)-interval(i,1))**2/self%theta_sq(i))-1.0_dp) &
                               + 2.0_dp*(interval(i,2)-interval(i,1))* &
                                 erf((interval(i,2)-interval(i,1))/sqrt(2.0_dp*self%theta_sq(i))) ) )
            else
               diff(:) = log( sqrt(two_pi*self%theta_sq(i)) * &
                             ( sqrt(self%theta_sq(i)/pi_per_two)* &
                               (exp(-0.5_dp*(interval(i,2)-interval(i,1))**2/self%theta_sq(i))-1.0_dp) &
                               + 2.0_dp*(interval(i,2)-interval(i,1))* &
                                 erf((interval(i,2)-interval(i,1))/sqrt(2.0_dp*self%theta_sq(i))) ) )
            end if
            iexp_arg(:) = iexp_arg(:) + diff(:)
         else
            diff(:) = x2(i,:)-x1(i)
            if (self%periodicity(i) > 0.0_dp) then
               iexp_arg(:) = iexp_arg(:) - 2.0_dp*sin(pi_per_periodicity(i)*diff(:))**2/self%theta_sq(i)
            else if (self%periodicity(i) < 0.0_dp) then
               diff(:) = diff(:) + nint(-diff(:)/self%periodicity(i))*self%periodicity(i)
               iexp_arg(:) = iexp_arg(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
            else
               iexp_arg(:) = iexp_arg(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
            end if
         end if
      end do
   end if

   if (present(e_f) .or. present(e_g)) then
      iexp_argj = 0.0_dp
      do j=1, n_dof
         if (interval(j,1) /= interval(j,2)) then
            if (self%periodicity(j) > 0.0_dp) then
               diff(:) = log(self%periodicity(j)*bessi0(1.0_dp/self%theta_sq(j))) - 1.0_dp/self%theta_sq(j)
            else if (self%periodicity(j) < 0.0_dp) then ! LAM
               diff(:) = log( sqrt(pi_per_two*self%theta_sq(j))* &
                         ( erf((x2(j,:)-interval(j,1))/sqrt(two*self%theta_sq(j))) &
                          -erf((x2(j,:)-interval(j,2))/sqrt(two*self%theta_sq(j))) ) )
            else
               diff(:) = log( sqrt(pi_per_two*self%theta_sq(j))* &
                         ( erf((x2(j,:)-interval(j,1))/sqrt(two*self%theta_sq(j))) &
                          -erf((x2(j,:)-interval(j,2))/sqrt(two*self%theta_sq(j))) ) )
            end if
            iexp_argj(:) = iexp_argj(:) + diff(:)
            if (present(e_g)) iexp_arg_j(j,:) = diff(:)
         else
            diff(:) = x2(j,:)-x1(j)
            if (self%periodicity(j) > 0.0_dp) then
               iexp_argj(:) = iexp_argj(:) - 2.0_dp*sin(pi_per_periodicity(j)*diff(:))**2/self%theta_sq(j)
            else if (self%periodicity(j) < 0.0_dp) then
               diff(:) = diff(:) + nint(-diff(:)/self%periodicity(j))*self%periodicity(j)
               iexp_argj(:) = iexp_argj(:) - 0.5_dp*diff(:)**2/self%theta_sq(j)
            else
               iexp_argj(:) = iexp_argj(:) - 0.5_dp*diff(:)**2/self%theta_sq(j)
            end if
         end if
      end do
   end if

   if (present(e_g)) then
      do j=1, n_dof
         if (interval(j,1) /= interval(j,2)) then
            if (self%periodicity(j) > 0.0_dp) then
               idexp_arg_j(j,:) = exp(-2.0_dp*sin(pi_per_periodicity(j)*(x2(j,:)+self%periodicity(j)/2.0_dp))**2 &
                                                                                                       /self%theta_sq(j)) &
                                - exp(-2.0_dp*sin(pi_per_periodicity(j)*(x2(j,:)-self%periodicity(j)/2.0_dp))**2 &
                                                                                                       /self%theta_sq(j))
            else if (self%periodicity(j) < 0.0_dp) then ! LAM
               idexp_arg_j(j,:) = exp(-0.5_dp*((x2(j,:)-interval(j,1))**2/self%theta_sq(j))) &
                                - exp(-0.5_dp*((x2(j,:)-interval(j,2))**2/self%theta_sq(j)))
            else
               idexp_arg_j(j,:) = exp(-0.5_dp*((x2(j,:)-interval(j,1))**2/self%theta_sq(j))) &
                                - exp(-0.5_dp*((x2(j,:)-interval(j,2))**2/self%theta_sq(j)))
            end if
         else
            diff(:) = x2(j,:)-x1(j)
            if (self%periodicity(j) > 0.0_dp) then
               idexp_arg_j(j,:) = -2.0_dp*sin(pi_per_periodicity(j)*diff(:)) &
                                         *cos(pi_per_periodicity(j)*diff(:))*two_pi_per_periodicity(j) / self%theta_sq(j)
            else if (self%periodicity(j) < 0.0_dp) then
               diff(:) = diff(:) + nint(-diff(:)/self%periodicity(j))*self%periodicity(j)
               idexp_arg_j(j,:) = -diff(:)/self%theta_sq(j)
            else
               idexp_arg_j(j,:) = -diff(:)/self%theta_sq(j)
            end if
         end if
      end do
   end if

   if (present(f_e) .or. present(g_e)) then
      iexp_argi = 0.0_dp
      do i=1, n_dof
         if (interval(i,1) /= interval(i,2)) then
             if (self%periodicity(i) > 0.0_dp) then
                diff(:) = log(self%periodicity(i)*bessi0(1.0_dp/self%theta_sq(i))) - 1.0_dp/self%theta_sq(i)
             else if (self%periodicity(i) < 0.0_dp) then ! LAM
                diff(:) = log( sqrt(pi_per_two*self%theta_sq(i))* &
                          ( erf((x1(i)-interval(i,1))/sqrt(two*self%theta_sq(i))) &
                           -erf((x1(i)-interval(i,2))/sqrt(two*self%theta_sq(i))) ) )
             else
                diff(:) = log( sqrt(pi_per_two*self%theta_sq(i))* &
                          ( erf((x1(i)-interval(i,1))/sqrt(two*self%theta_sq(i))) &
                           -erf((x1(i)-interval(i,2))/sqrt(two*self%theta_sq(i))) ) )
             end if
             iexp_argi(:) = iexp_argi(:) + diff(:)
             if (present(g_e)) iexp_arg_i(i,:) = diff(:)
         else
             diff(:) = x2(i,:)-x1(i)
             if (self%periodicity(i) > 0.0_dp) then
                iexp_argi(:) = iexp_argi(:) - 2.0_dp*sin(pi_per_periodicity(i)*diff(:))**2/self%theta_sq(i)
             else if (self%periodicity(i) < 0.0_dp) then
                diff(:) = diff(:) + nint(-diff(:)/self%periodicity(i))*self%periodicity(i)
                iexp_argi(:) = iexp_argi(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
             else
                iexp_argi(:) = iexp_argi(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
             end if
         end if
      end do
   end if

   if (present(g_e)) then
      do i=1, n_dof
         if (interval(i,1) /= interval(i,2)) then
            if (self%periodicity(i) > 0.0_dp) then
               idexp_arg_i(i,:) = exp(-2.0_dp*sin(pi_per_periodicity(i)*(x1(i)+self%periodicity(i)/2.0_dp))**2 &
                                                                                                   /self%theta_sq(i)) &
                                - exp(-2.0_dp*sin(pi_per_periodicity(i)*(x1(i)-self%periodicity(i)/2.0_dp))**2 &
                                                                                                   /self%theta_sq(i))
            else if (self%periodicity(i) < 0.0_dp) then ! LAM
               idexp_arg_i(i,:) = exp(-0.5_dp*((x1(i)-interval(i,1))**2/self%theta_sq(i))) &
                                - exp(-0.5_dp*((x1(i)-interval(i,2))**2/self%theta_sq(i)))
            else
               idexp_arg_i(i,:) = exp(-0.5_dp*((x1(i)-interval(i,1))**2/self%theta_sq(i))) &
                                - exp(-0.5_dp*((x1(i)-interval(i,2))**2/self%theta_sq(i)))
            end if
         else
            diff(:) = x2(i,:)-x1(i)
            if (self%periodicity(i) > 0.0_dp) then
               idexp_arg_i(i,:) = -2.0_dp*sin(pi_per_periodicity(i)*diff(:)) &
                                         *cos(pi_per_periodicity(i)*diff(:))*two_pi_per_periodicity(i) / self%theta_sq(i)
            else if (self%periodicity(i) < 0.0_dp) then
               diff(:) = diff(:) + nint(-diff(:)/self%periodicity(i))*self%periodicity(i)
               idexp_arg_i(i,:) = -diff(:)/self%theta_sq(i)
            else
               idexp_arg_i(i,:) = -diff(:)/self%theta_sq(i)
            end if
         end if
      end do 
   end if

   if (present(f_f) .or. present(g_f) .or. present(f_g) .or. present(g_g) .or. present(g_g_diag) .or. present(h_f)) then
      exp_arg = 0.0_dp
      do i=1, n_dof
         diff(:) = x2(i,:)-x1(i)
         if (self%periodicity(i) > 0.0_dp) then
            exp_arg(:) = exp_arg(:) - 2.0_dp*sin(pi_per_periodicity(i)*diff(:))**2/self%theta_sq(i)
         else if (self%periodicity(i) < 0.0_dp) then
            diff(:) = diff(:) + nint(-diff(:)/self%periodicity(i))*self%periodicity(i)
            exp_arg(:) = exp_arg(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
         else
            exp_arg(:) = exp_arg(:) - 0.5_dp*diff(:)**2/self%theta_sq(i)
         endif
      end do
   end if

   if (present(g_f) .or. present(g_g) .or. present(g_g_diag) .or. present(h_f)) then
      do i=1, n_dof
         diff(:) = x2(i,:)-x1(i)
         if (self%periodicity(i) > 0.0_dp) then
            dexp_arg_i(i,:) = 2.0_dp*sin(pi_per_periodicity(i)*diff(:)) &
                            * cos(pi_per_periodicity(i)*diff(:))*two_pi_per_periodicity(i) / self%theta_sq(i)
            if (present(g_g) .or. present(g_g_diag)) then
                ddexp_arg_ij(i,:) = (-sin(pi_per_periodicity(i)*diff(:))**2 &
                            + cos(pi_per_periodicity(i)*diff(:))**2)*two_pi_per_periodicity(i)**2/self%theta_sq(i)
            end if
            if (present(h_f)) then
                ddexp_arg_ii(i,:) = (sin(pi_per_periodicity(i)*diff(:))**2 &
                            - cos(pi_per_periodicity(i)*diff(:))**2)*two_pi_per_periodicity(i)**2/self%theta_sq(i)
            end if
         else if (self%periodicity(i) < 0.0_dp) then
            diff(:) = diff(:) + nint(-diff(:)/self%periodicity(i))*self%periodicity(i)
            dexp_arg_i(i,:) = diff(:)/self%theta_sq(i)
            if (present(g_g) .or. present(g_g_diag)) then
                ddexp_arg_ij(i,:) = 1.0_dp/self%theta_sq(i)
            end if
            if (present(h_f)) then
                ddexp_arg_ij(i,:) = -1.0_dp/self%theta_sq(i)
            end if
         else
            dexp_arg_i(i,:) = diff(:)/self%theta_sq(i)
            if (present(g_g) .or. present(g_g_diag)) then
                ddexp_arg_ij(i,:) = 1.0_dp/self%theta_sq(i)
            end if
            if (present(h_f)) then
                ddexp_arg_ij(i,:) = -1.0_dp/self%theta_sq(i)
            end if
         endif
      end do
   endif
   if (present(f_g) .or. present(g_g)) then
      do j=1, n_dof
         diff(:) = x2(j,:)-x1(j)
         if (self%periodicity(j) > 0.0_dp) then
            dexp_arg_j(j,:) = -2.0_dp*sin(pi_per_periodicity(j)*diff(:)) & 
                            * cos(pi_per_periodicity(j)*diff(:))*two_pi_per_periodicity(j) / self%theta_sq(j)
         else if (self%periodicity(j) < 0.0_dp) then
            diff(:) = diff(:) + nint(-diff(:)/self%periodicity(j))*self%periodicity(j)
            dexp_arg_j(j,:) = -diff(:)/self%theta_sq(j)
         else
            dexp_arg_j(j,:) = -diff(:)/self%theta_sq(j)
         endif
      end do
   endif

   if (present(e_e)) then
      e_e(:) = self%delta_sq * exp(iexp_arg(:))
   end if

   if (present(e_f)) then
      e_f(:) = self%delta_sq * exp(iexp_argj(:))
   end if

   if (present(f_e)) then
      f_e(:) = self%delta_sq * exp(iexp_argi(:))
   end if

   if (present(e_g)) then
      do j=1, n_dof
         if (interval(j,1) /= interval(j,2)) then
            e_g(j:nv*n_dof:n_dof) = self%delta_sq * exp(iexp_argj(:)-iexp_arg_j(j,:))*(idexp_arg_j(j,:))
         else
            e_g(j:nv*n_dof:n_dof) = self%delta_sq * exp(iexp_argj(:))*(idexp_arg_j(j,:))
         end if
      end do
   end if

   if (present(g_e)) then
      do i=1, n_dof
         if (interval(i,1) /= interval(i,2)) then
            g_e(i,1:nv) = self%delta_sq * exp(iexp_argi(:)-iexp_arg_i(i,:))*(idexp_arg_i(i,:))
         else
            g_e(i,1:nv) = self%delta_sq * exp(iexp_argi(:))*(idexp_arg_i(i,:))
         end if
      end do
   end if

   if (present(f_f)) then
      f_f(:) = self%delta_sq * exp(exp_arg(:))
   end if

   if (present(f_g)) then
      do j=1, n_dof
         f_g(j:nv*n_dof:n_dof) = self%delta_sq * exp(exp_arg(:))*(dexp_arg_j(j,:))
      end do
   endif

   if (present(g_f)) then
      do i=1, n_dof
         g_f(i,1:nv) = self%delta_sq * exp(exp_arg(:))*(dexp_arg_i(i,:))
      end do
   endif

   if (present(g_g)) then
      do i=1, n_dof
         do j=1, n_dof
            if (i /= j) then
               g_g(i,j:nv*n_dof:n_dof) = self%delta_sq * exp(exp_arg(:))*(dexp_arg_i(i,:))*(dexp_arg_j(j,:))
            else
               g_g(i,j:nv*n_dof:n_dof) = self%delta_sq * exp(exp_arg(:))*(ddexp_arg_ij(i,:)-dexp_arg_i(i,:)**2)
            endif
         end do
      end do
   end if

   if (present(g_g_diag)) then
      do i=1, n_dof
         g_g_diag(i:nv*n_dof:n_dof) = self%delta_sq * exp(exp_arg(:))*(ddexp_arg_ij(i,:)-dexp_arg_i(i,:)**2)
      end do
   end if

   if (present(h_f)) then
      do i=1, n_dof
         do j=1, n_dof
            if (i /= j) then
               h_f(i,j,1:nv) = self%delta_sq * exp(exp_arg(:))*(dexp_arg_i(i,:))*(dexp_arg_i(j,:))
            else
               h_f(i,j,1:nv) = self%delta_sq * exp(exp_arg(:))*(ddexp_arg_ii(i,:)+dexp_arg_i(i,:)**2)
            end if
         end do
      end do
   end if

   deallocate(pi_per_periodicity)
   deallocate(two_pi_per_periodicity)
   deallocate(diff)

   if(allocated(iexp_arg)) deallocate(iexp_arg)
   if(allocated(iexp_argj)) deallocate(iexp_argj)
   if(allocated(iexp_arg_j)) deallocate(iexp_arg_j)
   if(allocated(idexp_arg_j)) deallocate(idexp_arg_j) 
   if(allocated(iexp_argi)) deallocate(iexp_argi)
   if(allocated(iexp_arg_i)) deallocate(iexp_arg_i)
   if(allocated(idexp_arg_i)) deallocate(idexp_arg_i)
   if(allocated(exp_arg)) deallocate(exp_arg)
   if(allocated(dexp_arg_i)) deallocate(dexp_arg_i)
   if(allocated(ddexp_arg_ij)) deallocate(ddexp_arg_ij)
   if(allocated(dexp_arg_j)) deallocate(dexp_arg_j)
   if(allocated(ddexp_arg_ii)) deallocate(ddexp_arg_ii)

end subroutine SE_kernel_r_rr

subroutine GP_Matrix_Initialise(this,matrix,factorization)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  integer, intent(in) :: factorization

  if(this%initialised) call gp_matrix_finalise(this)

  this%n = size(matrix,1)
  this%m = size(matrix,2)
  allocate(this%matrix(this%n,this%m), this%factor(this%n,this%m), this%s(this%n), this%tau(this%m), this%ipiv(this%n))

  this%matrix = matrix
  this%initialised = .true.
  this%factorised = -factorization

end subroutine GP_Matrix_Initialise

subroutine GP_Matrix_Update(this,matrix,factorization)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  integer, optional, intent(in) :: factorization

  integer :: factor

  if (present(factorization)) then
      if(factorization /= CH .and. factorization /= QR .and. factorization /= BK) then
        print *,"unknown factorization, factorization= ", factorization
        stop
      end if
      factor = factorization
  else
      factor = CH
  end if

  if(this%initialised) then
     if( all(shape(matrix) == (/this%n,this%m/)) ) then
        this%matrix = matrix
     else
        call gp_matrix_initialise(this, matrix, factor)
     endif
  else
     call gp_matrix_initialise(this, matrix, factor)
  endif

  select case(-this%factorised)
     case(CH)
        call GP_Matrix_CH_Factorise(this)
     case(QR)
        call GP_Matrix_QR_Factorise(this)
     case(BK)
        call GP_Matrix_BK_Factorise(this)
  endselect

end subroutine GP_Matrix_Update

subroutine GP_Matrix_Finalise(this)

  implicit none

  type(GP_Matrix), intent(inout) :: this

  if(.not. this%initialised) return

  this%initialised = .false.
  this%equilibrated = .false.
  this%factorised = NOT_FACTORISED

  this%n = 0
  this%m = 0
  if(allocated(this%matrix) ) deallocate(this%matrix)
  if(allocated(this%factor) ) deallocate(this%factor)
  if(allocated(this%s) ) deallocate(this%s)
  if(allocated(this%tau) ) deallocate(this%tau)
  if(allocated(this%ipiv) ) deallocate(this%ipiv)

end subroutine GP_Matrix_Finalise

subroutine GP_Matrix_Assign(left_this, right_this)

  implicit none

  type(GP_Matrix), intent(inout) :: left_this
  type(GP_Matrix), intent(in) :: right_this

  call gp_matrix_finalise(left_this)

  left_this%initialised = right_this%initialised
  if( .not. left_this%initialised) return

  left_this%n = right_this%n
  left_this%m = right_this%m

  allocate(left_this%matrix(left_this%n,left_this%m))
  allocate(left_this%factor(left_this%n,left_this%m))
  allocate(left_this%s(left_this%n))
  allocate(left_this%tau(left_this%m))
  allocate(left_this%ipiv(left_this%n))

  left_this%matrix = right_this%matrix
  left_this%factorised = right_this%factorised

  if(left_this%factorised > NOT_FACTORISED) then
     left_this%factor = right_this%factor
     left_this%s = right_this%s
     left_this%tau = right_this%tau
     left_this%ipiv = right_this%ipiv
  end if

end subroutine GP_Matrix_Assign

subroutine GP_Matrix_Factorise(this)

  implicit none

  type(GP_Matrix), intent(inout) :: this

  select case(-this%factorised)
     case(CH)
       call GP_Matrix_CH_Factorise(this)
     case(QR)
       call GP_Matrix_QR_Factorise(this)
     case(BK)
       call GP_Matrix_BK_Factorise(this)
  end select

end subroutine GP_Matrix_Factorise

subroutine GP_Matrix_Solve_Matrix(this,matrix,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  if(this%factorised < NOT_FACTORISED) then
     call GP_Matrix_Factorise(this)
  endif

  select case(this%factorised)
     case(CH)
       call GP_Matrix_CH_Solve_Matrix(this,matrix,result)
     case(QR)
       call GP_Matrix_QR_Solve_Matrix(this,matrix,result)
     case(BK)
       call GP_Matrix_BK_Solve_Matrix(this,matrix,result)
  end select

end subroutine GP_Matrix_Solve_Matrix

subroutine GP_Matrix_Solve_Vector(this,vector,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:), intent(in) :: vector
  real(dp), dimension(:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: n, m

  n = size(vector)
  m = size(result)

  allocate(my_result(m,1))

  call GP_Matrix_Solve_Matrix(this,reshape(vector,(/n,1/)),my_result)
  result = my_result(:,1)

  deallocate(my_result)

end subroutine GP_Matrix_Solve_Vector

! Bunch-Kaufman factorization: A = UDU^{T}
subroutine GP_Matrix_BK_Factorise(this,u,d)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(out), optional :: u, d

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: info

  this%factor = this%matrix

  allocate(work(1))
  lwork = -1
  ! enquiry about optimal lwork for factorization
  call dsytrf('U', this%m, this%factor, this%n, this%ipiv, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! do the BK factorization
  call dsytrf('U', this%m, this%factor, this%n, this%ipiv, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_BK_Factorise: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  this%factorised = BK

! if( present(u) .and. present(d) ) then
!     call GP_Matrix_GetBK(this,u=u,d=d)
! else if( present(u) ) then
!     call GP_Matrix_GetBK(this,u=u)
! else if( present(d) ) then
!     call GP_Matrix_GetBK(this,d=d)
! end if

end subroutine GP_Matrix_BK_Factorise

! computes A^{-1}B matrix product using BK factorization of A: A^{-1}B=(UDU^T)^{-1}B
! this contains the factorization of A, matrix is B
subroutine GP_Matrix_BK_Solve_Matrix(this,matrix,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: info, n, o

  integer :: lwork
  real(dp), allocatable :: work(:)
  real(dp) :: rcond
  real(dp), allocatable :: ferr(:), berr(:)
  integer, allocatable :: iwork(:)

  n = size(matrix,1)
  o = size(matrix,2)
  if (size(result,1) /= this%m .or. size(result,2) /= o) then
     print *, "GP_Matrix_BK_Solve_matrix: shape(result) ",shape(result),"does not match",this%m,o
  endif

  if( n /= this%n ) then
     print *,'GP_Matrix_BK_Solve_Matrix: dimensions of U and matrix do not match.'
     stop
  endif

  allocate(my_result(n,o))
  lwork = -1
  allocate(work(1))
  allocate(ferr(o), berr(o), iwork(this%n))
  ! enquiry about the optimal lwork for calculating (UDU^T)^{-1}B
  call dsysvx('F', 'U', this%n, o, this%matrix, this%n, this%factor, this%n, this%ipiv, matrix, this%n, my_result, this%n, &
              rcond, ferr, berr, work, lwork, iwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! calculate (UDU^T)^{-1}B
  call dsysvx('F', 'U', this%n, o, this%matrix, this%n, this%factor, this%n, this%ipiv, matrix, this%n, my_result, this%n, &
              rcond, ferr, berr, work, lwork, iwork, info)
  deallocate(work)
  deallocate(ferr, berr, iwork)

  if( info /= 0 ) then
     print *,'GP_Matrix_BK_Solve_Matrix: dsysvx ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  result = my_result(1:this%m,:)
  deallocate(my_result)

end subroutine GP_Matrix_BK_Solve_Matrix

! QR factorization: A=QR
subroutine GP_Matrix_QR_Factorise(this,q,r)

  implicit none

  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: info

  this%factor = this%matrix

  allocate(work(1))
  lwork = -1
  ! enquiry about optimal lwork for factorization
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! do the QR factorization
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Factorise: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  this%factorised = QR

  if( present(q) .and. present(r) ) then
      call GP_Matrix_GetQR(this,q=q,r=r)
  else if( present(q) ) then
      call GP_Matrix_GetQR(this,q=q)
  else if( present(r) ) then
      call GP_Matrix_GetQR(this,r=r)
  end if

end subroutine GP_Matrix_QR_Factorise

subroutine GP_Matrix_GetQR(this,q,r)

  implicit none

  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: j, info

  if( this%factorised /= QR ) then
     print *,'GP_Matrix_GetQR: not QR-factorised, call GP_Matrix_QR_Factorise first.'
     stop
  endif

  if(present(q)) then
     if (size(q,1) /= this%n .or. size(q,2) /= this%m) then
        print *, "GT_Matrix_GetQR: shape(q) ",shape(q),"does not match",this%n,this%m
     endif
     q = this%factor

     allocate(work(1))
     lwork = -1
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     deallocate(work)
  endif

  if(present(r)) then
     if (size(r,1) /= this%n .or. size(r,2) /= this%m) then
        print *, "GP_Matrix_GetQR: shape(r) ",shape(r),"does not match",this%n,this%m
     endif
     r = this%factor(1:this%m,1:this%m)
     do j = 1, this%m - 1
        r(j+1:this%m,j) = 0.0_dp
     enddo
  endif

  if( info /= 0 ) then
     print *,'GP_Matrix_GetQR: ',(info),'-th parameter had an illegal value.'
     stop
  endif

end subroutine GP_Matrix_GetQR

! computes A^{-1}B matrix product using QR factorization of A: A^{-1}B=R^{-1}Q^{T}B
! this contains the factorization of A, matrix is B
subroutine GP_Matrix_QR_Solve_Matrix(this,matrix,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: info, n, o

  integer :: lwork
  real(dp), allocatable :: work(:)

  n = size(matrix,1)
  o = size(matrix,2)
  if (size(result,1) /= this%m .or. size(result,2) /= o) then
     print *, "GP_Matrix_QR_Solve_matrix: shape(result) ",shape(result),"does not match",this%m,o
  endif

  if( n /= this%n ) then
     print *,'GP_Matrix_QR_Solve_Matrix: dimensions of Q and matrix do not match.'
     stop
  endif

  allocate(my_result(n,o))
  my_result = matrix
  lwork = -1
  allocate(work(1))
  ! enquiry about the optimal lwork for calculating Q^{T}B
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! calculate Q^{T}B
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Solve_Matrix: dormqr ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  ! calculate R^{-1}Q^{T}B using back substitution
  call dtrtrs('U', 'N', 'N', this%n, o, this%factor, this%n, my_result, this%n, info)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Solve_Matrix: dtrtrs ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  result = my_result(1:this%m,:)
  deallocate(my_result)

end subroutine GP_Matrix_QR_Solve_Matrix

! Cholesky factorization: A=LL^{T}
subroutine GP_Matrix_CH_Factorise(this,l)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(out), optional :: l

  integer :: info

  this%factor = this%matrix

  call dpotrf('L', this%n, this%factor, this%n, info)

  if( info /= 0 ) then
     print *,'GP_Matrix_CH_Factorise: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  this%factorised = CH

! if( present(l) ) then
!     call GP_Matrix_GetU(this,l=l)
! end if

end subroutine GP_Matrix_CH_Factorise

! computes A^{-1}B matrix product using Cholesky factorization of A: A^{-1}B=(LL^{T})^{-1}B
! this contains the factorization of A, matrix is B
subroutine GP_Matrix_CH_Solve_Matrix(this,matrix,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: info, n, o

  integer :: lwork
  real(dp), allocatable :: work(:)

  n = size(matrix,1)
  o = size(matrix,2)
  if (size(result,1) /= this%m .or. size(result,2) /= o) then
     print *, "GP_Matrix_QR_Solve_matrix: shape(result) ",shape(result),"does not match",this%m,o
  endif

  if( n /= this%n ) then
     print *,'GP_Matrix_CH_Solve_Matrix: dimensions of L and matrix do not match.'
     stop
  endif

  allocate(my_result(n,o))
  my_result = matrix  

  call dpotrs('L', this%n, o, this%factor, this%n, my_result, this%n, info)

  if( info /= 0 ) then
     print *,'GP_Matrix_CH_Solve_Matrix: dsysvx ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  result = my_result(1:this%m,:)
  deallocate(my_result)

end subroutine GP_Matrix_CH_Solve_Matrix

real(dp) function bessi0(x)

  implicit none

  real(dp), intent(in) :: x
  real(dp) :: ax
  real(dp) :: y,p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
  DATA p1,p2,p3,p4,p5,p6,p7/1.0_dp,3.5156229_dp,3.0899424_dp,1.2067492_dp, &
                            0.2659732_dp,0.360768e-1_dp,0.45813e-2_dp/
  DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228_dp,0.1328592e-1_dp, &
                                  0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,-0.2057706e-1_dp, &
                                  0.2635537e-1_dp,-0.1647633e-1_dp,0.392377e-2_dp/
  if (abs(x) .lt. 3.75_dp) then
      y=(x/3.75_dp)**2
         bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
  else
      ax=abs(x)
      y=3.75_dp/ax
      bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
  endif
  return
end function bessi0

end module gp_basic_mod
