module constants_and_parameters
!	Declarar Unidades

real*8 :: N

real*8 , parameter  :: half = 0.5d0, zero = 0.0d0, one = 1.0d0, two = 2.0d0, five = 5.0d0, four = 4.0d0, three = 3.0d0, eight = 8.0d0 , six = 6.0d0 , nine = 9.0d0

real*8 , parameter :: ev_au             = 27.211
real*8 , parameter :: j_to_ev           = 6.242e18
real*8 , parameter :: ang_au            = 5.2917721092e-9
real*8 , parameter :: hbar             = 6.626070040e-34*j_to_ev
real*8 , parameter :: t_au              = 2.42e-17
real*8 , parameter :: kvec              = 21947463e0
real*8 , parameter :: pi                = 3.141592e0
real*8 , parameter :: a0_to_comp        = 1.158826e-3
real*8 , parameter :: a0                = 0.528e-10 
real*8 , parameter :: eps               = 1.0d-8
real*8 , parameter :: massp             = 2000.d0!*9.10938356e-31 !mass in a.u For comparasion, the mass of the hydrogen atom is 1836 a.u!

integer , parameter :: pos_index = 2

integer :: grid_size , passo , npt , ni , loop , nsteps_el , nrk0 , mkf , k_step
integer :: n_all_roots , R_index , n_of_traj , type_potential , max_steps , size_hnn , n_collapse , size_inn

real*8  , allocatable   :: x(:) , coef_coupl_af(:,:) , eigenenergy(:) , enp(:) , potential(:,:) , nh(:) , vs1(:), vs2(:),vs3(:),vs4(:)
 
real*8 , allocatable :: coef_coupl(:,:) , rho_coupl(:,:) , dcoef_coupl(:,:) , force_HF(:) , vec_stat(:,:) , mat_prob(:,:) 
real*8 , allocatable :: vec_P_tr1(:,:),vec_P_tr2(:,:) , vec_P_refl1(:,:), vec_P_refl2(:,:), vec_P_conf(:,:)

real*8 , allocatable :: matriz_estocastica(:,:) , prob_evento(:) , dcoast(:,:) , x_toll(:)

complex*16 , parameter :: zi = (0.0d0,1.d0)

real*8 :: dt ,dx, kx, tfinal , nr12  , nt12 , nnt1 , nnr1 , energyconserv ,time , de_i , dx_i , de_max , de_min , dde_max , dde_min,count_inter_conf , konst , konst2 , transition_n , noftrans12 , noftrans21 , nofNtrans1 , nofNtrans2 , t_and_t_wot ,t_and_nt_wot ,r_and_t_wot , r_norm , find_n1 , find_n2 , find_confined , beta_bn_par , gama_bn_par 
real*8 :: pos, vel, acel , theta_af , var_phi_before , inflex_denerg , dt_el , dx_el, vs_p1, vs_p2, vs_p3, vs_p4 , dx_pos , kx_ini ,p_norm , x_inter , shannon_11, shannon_12r, shannon_21i, rho_N_1,rho_N_2

real*8 , allocatable :: Ntr2(:) , Ntr1(:) , Nrefl1(:), Nrefl2(:), Nconf(:) 
real*8 :: Tr_R_rho , Tr_P_rho , Tr_F_rho ,time_collapse

real*8    , allocatable         :: rho_Real(:,:)                  ! Real part of rho Re(rho) 
real*8    , allocatable         :: rho_Imag(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: weight0_tmp(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: input_tmp(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: input_neuron0(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: input_neuron2(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: a_parameter(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: output_desired(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: output_desired_wot(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: sigmoid_arg(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: bias(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: weight1(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: m_vec(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: input_neuron1(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: nneural_vec(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: impulse_tully(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: nneural_vec_tmp(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: f_loss(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: input_neuron_epoch(:,:,:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: weight0(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: dw1(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: dw0(:,:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: db1(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: db0(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: k_rand(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: traj_rand(:)                  ! Imaginary part of rho Im(rho)
real*8    , allocatable         :: dR_fg(:,:)
real*8    , allocatable         :: dP_fg(:,:)
real*8    , allocatable         :: dF_fg(:,:)



character(len=255) :: pwd

real*8 :: posi , posf , t_rate_el , dbeta , dgama

logical :: check_time_initial , trans_prob , train_IA , batch_norm , finish_tully , end_step

complex*16, allocatable         :: c_phi_dot(:) , rho_m(:,:) , nonadiabatic_coupling(:,:) , rho_mold(:,:) , j_m(:,:) ,rho_collaps(:,:)
! rho dot d(rho)/dt

integer :: count_step , count_train

integer :: nl
! System size
integer :: neqn
! Number of equations
integer :: dimtot
! Total Dimension


!**********************************************************************
! ode parameters
real*8 , parameter :: relerr = 1.d-12                   ! ode parameters
real*8 , parameter :: abserr = 1.d-12                   ! ode parameters

real*8     , allocatable :: work(:)                    ! ode parameters
integer::  iwork(5)                                    ! ode parameters
integer::  iflag                                       ! ode parameters
!**********************************************************************



integer, dimension(8) :: values

end module constants_and_parameters

