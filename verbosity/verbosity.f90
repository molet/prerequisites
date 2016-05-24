module verbosity_mod

use prec_mod

private

integer :: start_time, stop_time, clock_rate, clock_max

public :: verbosity_start_task, verbosity_finish_task

contains

!--------------------------------------------------------------------------------------------------

subroutine verbosity_start_task(task)

  implicit none

  character(len=*), intent(in) :: task

  write(6,*) trim(task), ' has started...'

  call system_clock ( start_time, clock_rate, clock_max )

end subroutine verbosity_start_task

!--------------------------------------------------------------------------------------------------

subroutine verbosity_finish_task(task)

  implicit none

  character(len=*), intent(in) :: task

  call system_clock ( stop_time, clock_rate, clock_max )
 
  write(6,*) 'elapsed time: ', real(stop_time-start_time,dp)/real(clock_rate,dp)
  write(6,*) trim(task), ' has finished...'

end subroutine verbosity_finish_task

!--------------------------------------------------------------------------------------------------

end module verbosity_mod
