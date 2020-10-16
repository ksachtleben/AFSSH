program dwell

use constants_and_parameters
use griding                  , only : grade , random_normal , phase_gauss
use nuclear_din              , only : dinamic_nuclei , draft_psi1
use el_din                   
use energy_surface
use mediadat
use parameters_confg


use mkl_dfti

implicit none


real*8 :: start, finish 
integer :: io1 , i , j , io4 , a , iounit , nii , nif , m , gap , n_of_pictures , imag_index
real*8  :: pos_in , pos_f , kxin , velocity , E_i , dk , par_eki , coast , anya

real*8 , allocatable :: coast_k(:) , rn(:) , temp1(:) ,temp2(:) , temp3(:) ,temp4(:) , k_vec(:) , montante(:)

integer , dimension(:) , allocatable :: seed
call system( 'rm trained_weight/*.txt' )

call getcwd(pwd)

call cpu_time(start)

finish_tully = .false.

kxin = sqrt(4000.d0*2.72d0)
mkf  = 90
par_eki = (float(mkf) - one)*5.d0/float(mkf) - 0.99d0

dgama = 0.d0
dbeta = 0.d0

size_hnn = 100
size_inn = 8

n_of_pictures = 300 ! max number of pictures is mkf*n_of_traj

call parameters_module( type_potential , nl , n_of_traj , pos_in , pos_f , dt , ni , velocity , kxin )
if(.not. allocated( bias ))     allocate( bias( size_hnn + 3 ) )
if(.not. allocated( weight1 ))  allocate( weight1(3, size_hnn ) )
if(.not. allocated( weight0 ))  allocate( weight0( size_hnn ,size_inn) )
if(.not. allocated( rn ) )      allocate( rn( size(bias) + size(weight1) + size(weight0) + 2) )
if(.not. allocated( traj_rand ))allocate( traj_rand( n_of_pictures ) )
if(.not. allocated( k_rand ))   allocate( k_rand( n_of_pictures ) )
if(.not. allocated( temp1 ))    allocate( temp1( n_of_pictures ) , source=0.d0 )
if(.not. allocated( temp2 ))    allocate( temp2( n_of_pictures ) , source=0.d0 )
if(.not. allocated( temp3 ))    allocate( temp3( n_of_pictures ) , source=0.d0 )
if(.not. allocated( temp4 ))    allocate( temp4( n_of_pictures ) , source=0.d0 )
call grade
if(.not. allocated( montante )) allocate( montante( grid_size ) )
if(.not. allocated( k_vec ))    allocate( k_vec(mkf) ) 

call random_number( traj_rand )
call random_number( k_rand )
call phase_gauss( rn )

coast = 1.d0
count_train = 1
do !i = 1 , 500
            if(.not. allocated( vec_P_tr1 ))            allocate( vec_P_tr1(mkf,grid_size) , source= 0.d0) 
            if(.not. allocated( vec_P_tr2 ))            allocate( vec_P_tr2(mkf,grid_size) , source= 0.d0) 
            if(.not. allocated( vec_P_refl1 ))          allocate( vec_P_refl1(mkf,grid_size) , source= 0.d0) 
            if(.not. allocated( vec_P_refl2 ))          allocate( vec_P_refl2(mkf,grid_size) , source= 0.d0) 
            if(.not. allocated( vec_P_conf ))           allocate( vec_P_conf(mkf,grid_size) , source= 0.d0) 

        do k_step = mkf , 1 , -1 

        kxin = 50.0d0 - 40.0d0 * float( k_step - 1 ) / float( mkf - 1 )
        !kxin = sqrt( two*massp* ( 0.625d0 - 0.6d0 * float( k_step ) / float( mkf )  )) 
        kxin = 35.0d0
        if( finish_tully == .false. ) print*, 'Axis X' , kxin
        kx_ini  = kxin
        
        call parameters_module( type_potential , nl , n_of_traj , pos_in , pos_f , dt , ni , velocity , kxin )
            vel  = velocity
            posi = pos_in
            posf = pos_f
            passo = 0
            time = 0.d0
            transition_n = 0.d0
            check_time_initial = .false.
        
            Ntr1        = 0.d0
            Ntr2        = 0.d0
            Nrefl1      = 0.d0
            Nrefl2      = 0.d0
            Nconf       = 0.d0
        
            p_norm = 1.d0
            r_norm = 1.d0
        
            noftrans12 = 0d0
            noftrans21 = 0.d0
            nofNtrans1 = 0.d0
            nofNtrans2 = 0.d0
            find_n2       = 0.d0
            find_n1       = 0.d0
            find_confined = 0.d0
        
            batch_norm = .false.
        
            if(.not. allocated( vec_stat ))      allocate( vec_stat(2,2) , source= 0.d0) 

            energyconserv = 0.d0
            if( finish_tully == .true. ) n_of_traj = 1

        do loop = 1 , n_of_traj
            call parameters_module( type_potential , nl , n_of_traj , pos_in , pos_f , dt , ni , velocity , kxin )
            vel = velocity
            pos = pos_in
            posf = pos_f
            nii = ni
            if( allocated( impulse_tully ) )deallocate( impulse_tully )
            if(.not. allocated(rho_m))      allocate( rho_m(nl,nl) , source=(0.d0,0.d0) )
            if(.not. allocated(j_m))        allocate( j_m(nl,nl) , source=(0.d0,0.d0) )
            if(.not. allocated(rho_mold))   allocate( rho_mold(nl,nl) , source=(0.d0,0.d0) )
            if(.not. allocated(dR_fg))      allocate( dR_fg(nl,nl) , source=0.d0 )
            if(.not. allocated(dP_fg))      allocate( dP_fg(nl,nl) , source=0.d0 )
            if( allocated( rho_collaps ) ) rho_collaps = (0.d0,0.d0)
            rho_m            = (0.d0,0.d0)
            j_m              = (0.d0,0.d0)
            rho_mold         = (0.d0,0.d0)
            dR_fg            = 0.d0
            dP_fg            = 0.d0
            rho_m( ni , ni ) = (1.d0,0.d0)
            rho_mold( ni , ni ) = (1.d0,0.d0)
            time = 0.d0
            shannon_11  = 0.d0
            shannon_12r = 0.d0
            shannon_21i = 0.d0
            rho_N_1 = 0.d0
            rho_N_2 = 0.d0
            n_collapse = 0
            
            end_step = .false.
            call dinamic_nuclei( ni )
stop
            if( allocated( rho_collaps ) ) deallocate(rho_collaps)

        
            if(.not. allocated( coast_k ))              allocate( coast_k(mkf) , source= 0.d0) 
        
!
        !deallocate( coef_coupl , rho_coupl, coef_coupl_af , eigenenergy , enp , potential , rho_m , impulse_tully )
        enddo
                write( 200 + k_step , 200 ) "started:" , posi , "tpoints:" , grid_size
                do i = 1 , grid_size
                montante(i)               = Ntr1(i) + Ntr2(i) + Nrefl2(i) + Nrefl1(i) + Nconf(i)
                vec_P_tr1(k_step,i)       = Ntr1(i) / montante(i)
                vec_P_tr2(k_step,i)       = Ntr2(i) / montante(i)
                vec_P_refl1(k_step,i)     = Nrefl1(i) / montante(i)
                vec_P_refl2(k_step,i)     = Nrefl2(i) / montante(i)
                vec_P_conf(k_step,i)      = Nconf(i) / montante(i)
       !         print 16, x(i) , vec_p_tr1(k_step,i),vec_p_tr2(k_step,i),vec_p_refl1(k_step,i),vec_p_refl2(k_step,i),montante(i)  , float(i)
                write( 200 + k_step , 16 ) kxin , x(i) , vec_p_tr1(j,i),vec_p_tr2(j,i),vec_p_refl1(j,i),vec_p_refl2(j,i),vec_p_conf(j,i)
                enddo
            !    print*, x( grid_size )
        enddo
        
        call draft_psi1( x , k_vec , vec_p_tr1 , 1 )
        call draft_psi1( x , k_vec , vec_p_tr2 , 2 )
        call draft_psi1( x , k_vec , vec_p_refl1 , 3 )
        call draft_psi1( x , k_vec , vec_p_refl2 , 4 )
        call draft_psi1( x , k_vec , vec_p_conf , 5 )

        if( finish_tully == .true. ) coast_k( k_step ) = sum((input_neuron2 - output_desired)**2)/two


        if( count_train == n_of_pictures ) count_train = 0
        count_train = count_train + 1
        if( finish_tully == .true. ) coast = sum( coast_k )/mkf
!        print*, coast , ' <=== custo ', count_train
        if( coast < 5e-3 ) exit
        if( train_ia == .true. ) write(45 , 17) coast
        if( train_ia == .false. ) exit
        finish_tully = .true.
enddo
18  format(a16,a16,a16,a16,a16,a16)
13  format(3es16.4E3)
14  format(4es16.4E3)
15  format(2es16.4E3)
16  format(7f20.12)
17  format(1es16.4E3)
110 format(a16,a16,a16,a16,a16,a16,i3)
100 format(3es16.4E3,a5,3es16.4E3)
200 format(a8,1f20.12,a8,1f20.12)

call cpu_time(finish)
PRINT '("Time =",f8.3," hours." )', (finish - start)/3600.0
PRINT '("Time =",f8.3," seconds." )', finish - start
!PRINT '("Time =",f8.3," minutes." )', (finish - start)/60.0

end program
