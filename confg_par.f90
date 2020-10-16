module parameters_confg
use constants_and_parameters

implicit none


contains
!================================================================================================================
subroutine parameters_module( type_potential , nl , n_of_traj , pos_in , pos_f , dt_i , statein , vel , kxin )
!================================================================================================================
real*8 , intent(out) :: dt_i , vel , pos_in , pos_f
integer , intent(out) :: type_potential , nl , n_of_traj , statein
real*8 , intent(in)  :: kxin

train_IA = .false.
statein         =   1
! Select a type of potential
! Select 1 if u are interested in Simple Avoided Crossing
! Select 2 if u are interested in Dual Avoided Crossing
! Select 3 if u are interested in Extended Coupling with Reflection
type_potential  =   1
nl              =   2
n_of_traj       =   2000
pos_in          =   -20.0d0
pos_f           =   20.0d0
!kxin            =   4d0
vel             =   kxin/massp
dt_i            =   20.d0 * ( 10.d0 / kxin )
max_steps       =   8000
!4.512d0*(vel**half) - 0.2028d0!-1.065d0*(vel**half) + 0.1715d0    !0.11d0*( nine/kxin )**(half)
nsteps_el       = 1

!0.11d0 au para k=9a

end subroutine

end module
