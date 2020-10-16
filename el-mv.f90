module el_din
use f95_precision
use blas95
use lapack95

use constants_and_parameters
use griding                    , only : grade , sumtrap , phase_gauss
use energy_surface
use ordinary_equation_solution

implicit none

public 


contains
!===============================================================================
 Subroutine dinamic_el( dtime )
!===============================================================================
real*8  , intent(inout) :: dtime


integer  :: i , j , n , m , k , nvec , io , l , istep , i_prob 
real*8   :: norma , a ,b
complex*16 :: norm

real*8 :: trans_hop                                            ! Grid increment
real*8 :: tout , dte  , tau_kl  , exp_var  , dt_call                                 ! Evolution time
real*8 :: t  , t_i
real*8     , allocatable :: vec(:) , rho_it(:) , par_hopping(:,:,:) ,alg_eign(:)
real*8 , allocatable :: id_matrix(:,:) , diag_m_En(:,:)
complex*16 , allocatable :: U(:,:) , matrix_alloc(:,:), diag_m_En2(:,:)
logical :: transition_pr


nl = nl
dimtot = nl*nl
neqn = two*dimtot

if(.not. allocated(work))       allocate( work( 100 + 21*neqn ) )
if(.not. allocated(vec))        allocate( vec(neqn) , source = (0.d0) )                 ! building rho vector  
!if(.not. allocated(rho_m))      allocate( rho_m(nl,nl) , source=(0.d0,0.d0) )
!if(.not. allocated(rho_mold))   allocate( rho_mold(nl,nl) , source=(0.d0,0.d0) )
if(.not. allocated(U))          allocate( U(nl,nl)  )
if(.not. allocated(diag_m_En))  allocate( diag_m_En(nl,nl) , source=0.d0)
if(.not. allocated(diag_m_En2))  allocate( diag_m_En2(nl,nl) , source=(0.d0,0.d0))
if(.not. allocated(alg_eign))   allocate( alg_eign(nl) )
if(.not. allocated(matrix_alloc))  allocate( matrix_alloc(nl,nl) , source=(0.d0,0.d0))

!rho_m(ni,ni) = (1.d0,0.d0)


t_rate_el       = dtime/float(nsteps_el)

diag_m_En(1,1) = enp(1)
diag_m_En(2,2) = enp(2)

!!if( dtime == dt ) vel = half*( vel_d + vel_a )

U = diag_m_En - zi*nonadiabatic_coupling*vel 
call heevd( U , alg_eign , 'V' , 'U' )
!estava comentado ja \/
!print*, U(1,1) , U(1,2) !eigenvector of eigenvalue 1
!print*, U(2,1) , U(2,2) !eigenvector of eigenvalue 2

!estava comentado ja \/
!print*, matmul( U , transpose(conjg(U)) )

!write( 26 , 13 ) pos , alg_eign(1)
!write( 27 , 13 ) pos , alg_eign(2)

forall( i = 1 : nl ) diag_m_En2(i,i) = exp( -zi*alg_eign(i)*t_rate_el )
matrix_alloc = matmul( U , matmul( diag_m_En2 , transpose(conjg(U)) ) )

iflag = 1
    
    rho_m = matmul( matmul( matrix_alloc , rho_m ) , transpose(conjg( matrix_alloc ) ) )
    j_m = j_m - rho_mold / vel 
    rho_mold = matmul( matmul( matrix_alloc , rho_mold ) , transpose(conjg( matrix_alloc ) ) )
    j_m = j_m + rho_mold / vel 


    vs_p1 = rho_m(1,1)!/float(n_of_traj)
    vs_p2 = rho_m(2,2)!/float(n_of_traj)
    vs_p3 = aimag(rho_m(1,2))!/float(n_of_traj)
    vs_p4 = real(rho_m(1,2))!/float(n_of_traj)


!if( loop == n_of_traj ) then
!    write( 36 , 13 ) pos , vs_p1
!    write( 37 , 13 ) pos , vs_p2
!    write( 38 , 13 ) pos , vs_p3
!    write( 39 , 13 ) pos , vs_p4
!end if

13 format(4es16.4E3)
deallocate( diag_m_en )
end subroutine dinamic_el
!==============================================================================
 subroutine hopping( rho , rho2 , cur_rho , save_ , toll_i )
!==============================================================================
complex*16 , intent(in)     :: rho(:,:) 
complex*16 , intent(in)     :: rho2(:,:) 
complex*16 , intent(in)     :: cur_rho(:,:) 
logical    , intent(in)     :: save_
integer    , intent(in)     :: toll_i

integer :: i,j,k,m , nr , nii , nif , n , nall_step , loop_sh , nrk , nrki , io4 , iounit, a_int , b_int
real*8  :: trace , tau_kl , dvel , a_r , b_r , c_r , reescale_comp , m_direc , rand , end_loop_training , a_par , b_par , weight_weight , friction_rate ,a 
real*8  :: c_par,d_par,abs_vel , t_and_t , t_and_nt , r_and_t , norm1112 , norm2221 , coast ,learn_rate , variance , mean , sizem , dvariance , direction

real*8  , allocatable :: nad_coupl(:,:) ,par_hop(:) , dir_nadcoupl(:) , aux(:,:) , x_grid(:) , Matriz_input(:,:) ,mat_temp(:,:) ,error_d2j(:) , error_d1j(:) , rand_n2(:)
real*8  , allocatable :: var_lay1(:) , var_lay2(:) , gama_BN(:) , beta_BN(:) , BN_lay1(:), error_w0(:,:) , input_neuron1_norm(:) , rand_n(:)
integer , allocatable :: nh(:) 

integer , dimension(:) , allocatable :: seed

character :: fmt 
character(len=2) :: x1

nall_step = 120005!nsteps_el*npt

if(.not. allocated(nad_coupl) )     allocate( nad_coupl( nl,nl ) )
if(.not. allocated(aux) )           allocate( aux( nl,nl ) )
if(.not. allocated( nh )    )       allocate( nh( nl ) )
if(.not. allocated( rand_n )    )   allocate( rand_n( nl ) )
if(.not. allocated( rand_n2)    )   allocate( rand_n2(nl ) )
if(.not. allocated( dir_nadcoupl )) allocate( dir_nadcoupl( nl ) )
if(.not. allocated( par_hop ))      allocate( par_hop( nl ) )
if(.not. allocated( x_grid ))       allocate( x_grid( 2001 ) )
if(.not. allocated( mat_temp ) )    allocate( mat_temp(2,2) , source=0.d0)
    

nad_coupl = nonadiabatic_coupling


   ! rn = exp(-half*rn**2)/sqrt(two*pi)


call random_number( rand_n )

m= 1
do i = 1 , nl
        dir_nadcoupl(:) = two*dt*(-vel*real(rho(ni,:))*nad_coupl(:,ni))/rho(ni,ni)
        forall( i = 1 : nl ) par_hop(i) = sum(dir_nadcoupl(1:i))
        nh(m) = i 

        if( par_hop(nh(m))*heaviside( two*(eigenenergy(ni) - eigenenergy(nh(m)))/(massp*vel**2) + one ) > rand_n(nh(m)) ) then 

        direction = sign( one , nad_coupl( ni , nh(m) ) )

        a_r = ( massp * direction ) * ( massp * direction )
        b_r = two * ( massp * vel ) * ( massp * direction )
        c_r = -two * massp * ( eigenenergy( ni ) - eigenenergy( nh(m) ) )

        dvel = roots_square( a_r , b_r , c_r )

        vel = vel + dvel * direction 

        !if( b_r**2 - four * a_r * c_r >= 0.d0 ) then
!            nr      = ni
!            ni      = nh(m)
!            nh(m)   = nr
        !endif

        !rho_N_2 = rho_N_2 + one
                forall( i =1:nl ) dP_fg(i,i) = dP_fg(i,i) - dP_fg(ni,ni)
                forall( i =1:nl ) dR_fg(i,i) = dR_fg(i,i) - dR_fg(ni,ni)

        dP_fg(ni,ni) = 0.d0
        dR_fg(ni,ni) = 0.d0
            
        trans_prob = .true.
!        if( ni == 2 ) transition_n = transition_n + one
!        exit
!        stop "transition can be occur"
        endif
        if( trans_prob == .true. ) exit
        m = m + 1
enddo

write( 55+type_potential , * ) real(two*dt*(vel*real(rho(1,2))*nad_coupl(1,2))/rho(1,1) )

        if(.not. allocated( gama_BN ))  allocate( gama_BN( size_hnn ) )
        if(.not. allocated( beta_BN ))  allocate( beta_BN( size_hnn ) )
        if(.not. allocated( db1 ))     allocate( db1( 3 ) , source=0.d0 )
        if(.not. allocated( db0 ))     allocate( db0( size_hnn ) , source=0.d0 )
        if(.not. allocated( dw1 ))  allocate( dw1(3,size_hnn) , source=0.d0 )
        if(.not. allocated( dw0 ))  allocate( dw0(size_hnn,size_inn) , source=0.d0)

        if(.not. allocated( m_vec ))            allocate( m_vec( size_hnn ) )
        if(.not. allocated( input_neuron1 ))    allocate( input_neuron1(size_hnn) , source = 0.d0 )
        if(.not. allocated( sigmoid_arg ))      allocate( sigmoid_arg( 3 ) )
        if(.not. allocated( input_neuron2 ))    allocate( input_neuron2( 3 ) )
        if(.not. allocated( error_d2j ))        allocate( error_d2j( 3 ) )
        if(.not. allocated( error_d1j ))        allocate( error_d1j( size_hnn ) )
        if(.not. allocated( error_w0 ))         allocate( error_w0( size_hnn , size_hnn ) )
        if(.not. allocated( f_loss ))           allocate( f_loss( 3 ) )
        if(.not. allocated( BN_lay1 ))          allocate( BN_lay1( size_hnn ) )
        if(.not. allocated( var_lay2 ))         allocate( var_lay2( 3 ) )
        if(.not. allocated( prob_evento ))      allocate( prob_evento( 3 ) , source = 0.d0 )
        if(.not. allocated( impulse_tully ))    allocate( impulse_tully( 4 ) , source = 0.d0 )
        if(.not. allocated( output_desired ))   allocate( output_desired( 3 ) , source= 0.d0 )
        if(.not. allocated( dcoast ))           allocate( dcoast( 3, mkf ) , source= 0.d0 )
        if(.not. allocated( input_neuron0 ))    allocate( input_neuron0(size_inn) , source = 0.d0 )
 
        if(.not. allocated( input_neuron1_norm ))    allocate( input_neuron1_norm(size_hnn) , source = 0.d0 )
        
        dir_nadcoupl(:) = two*dt*(-vel*real(rho(1,:))*nad_coupl(:,1))/rho(1,1)
        par_hop(2) = sum(dir_nadcoupl)
        dir_nadcoupl(:) = two*dt*(-vel*real(rho(2,:))*nad_coupl(2,:))/rho(2,2)
        par_hop(1) = sum(dir_nadcoupl)

        impulse_tully(1) = impulse_tully(1) + (par_hop(2))*heaviside( vel**2 - two*(eigenenergy(1) - eigenenergy(2)/massp ))
 !       impulse_tully(1) = impulse_tully(1) + par_hop(2)
        impulse_tully(2) = impulse_tully(2) + (par_hop(1))*heaviside( vel**2 - two*(eigenenergy(2) - eigenenergy(1))/massp )
        impulse_tully(3) = impulse_tully(3) + ( one - par_hop(2) )*heaviside( vel**2 - two*(eigenenergy(1) - eigenenergy(2))/massp )
        impulse_tully(4) = impulse_tully(4) + ( one - par_hop(1) )*heaviside( vel**2 - two*(eigenenergy(2) - eigenenergy(1))/massp )

        shannon_11   = shannon_11  - rho2(1,1)*log( rho2(1,1) + eps ) - rho2(2,2)*log(rho2(2,2) + eps )

        if( ni == 1 ) rho_N_1 = rho_N_1 + one
        if( ni == 2 ) rho_N_2 = rho_N_2 + one


!shannon_21i  = shannon_21i + a 
if( save_ == .true. ) then

        input_neuron0(1)   = real(rho( 1 , 1 ))

        input_neuron0(2)   = real(rho( 2 , 2 ))

        input_neuron0(3)   = (kx_ini - massp*vel) / sqrt( two * massp )
        
        input_neuron0(4)   = (kx_ini + massp*vel) / sqrt( two * massp )

        input_neuron0(5)   = eigenenergy( 2 ) - eigenenergy( 1 )
                
        input_neuron0(6)  = (dR_fg(1,1) - dR_fg(2,2))*sign( one , ( dR_fg(1,1) - dR_fg(2,2) ) / ( dP_fg(1,1) - dP_fg(2,2) ) )
        
        input_neuron0(7)  = half*(Force_HF(1) - Force_HF(2))

        input_neuron0(8)  = nad_coupl(1,2)

endif

call random_number( rand_n2 )
a = ( dR_fg(1,1) - dR_fg(2,2) ) / ( dP_fg(1,1) - dP_fg(2,2) )
a = half*dt*( Force_hf(1) - Force_hf(2) )*( dR_fg(1,1) - dR_fg(2,2) )*sign( one , a )

if( rand_n2(1) < a ) then
        dR_fg = 0.d0
        dP_fg = 0.d0
        rho_m = 0.d0
        rho_m(ni,ni) = one
endif


13 format(3es16.4E3)
17 format(1es16.4E3)
18 format(6es16.4E3)
15 format(i5.3,i5.3,3es16.4E3)
end subroutine hopping

!====================================================================
 function roots_square( a , b , c )
!====================================================================
real*8 , intent(in) :: a , b , c

real*8 , dimension(2) :: min_value_abs 
real*8 :: roots_square

    min_value_abs(1) = -b/(two*a) - sqrt( b**2 - four*a*c )/(two*a)
    min_value_abs(2) = -b/(two*a) + sqrt( b**2 - four*a*c )/(two*a)
    
    roots_square = minval( abs(min_value_abs) )  


end function

!====================================================================
character(len=20) function str(k)
!====================================================================
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!====================================================================

!====================================================================
real*8 function heaviside(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x

    heaviside = sign(0.5d0,x) + 0.5d0
end function heaviside
!====================================================================

!====================================================================
function relu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8, dimension(size(x)) :: relu    

    relu = (sign(x(:),x(:)) + abs(x(:)))/two
end function relu
!====================================================================

!====================================================================
function Smoothrelu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8, dimension(size(x)) :: smoothrelu    

    real*8 :: y

    y = exp(-maxval(x))

    smoothrelu = log( exp( x(:) + log(y) ) + y ) - log(y)
end function smoothrelu
!====================================================================

!====================================================================
real*8 function absolute(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x

    absolute = (sign(0.5d0,x) + 0.5d0)*x
end function absolute
!====================================================================

!====================================================================
function drelu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8 , dimension(size(X)) :: drelu

    drelu = sign(0.5d0,x(:)) + 0.5d0
end function drelu
!====================================================================
!====================================================================
function Dsmoothrelu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8, dimension(size(x)) :: dsmoothrelu    
    
    real*8 :: y

    y = exp(-maxval(x))

    dsmoothrelu = exp(x + log(y))/( exp(x + log(y)) + y )
end function dsmoothrelu
!====================================================================
!====================================================================
function elu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8, dimension(size(x)) :: elu    

    integer :: n
    real*8  :: alph

    alph = one

    do n = 1 , size(x)
        if( x(n) >= 0.d0 )      elu     = x
        if( x(n) < 0.d0 )       elu(n)  = alph*(exp(x(n)) - one)
    enddo

end function elu
!====================================================================
!====================================================================
function delu(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8, dimension(size(x)) :: delu    

    integer :: n
    real*8  :: alph

    alph = one

    do n = 1 , size(x)
        if( x(n) >= 0.d0 )     delu     = one
        if( x(n) < 0.d0 )      delu(n)  = alph*(exp(x(n)) - one) + alph
    enddo

end function delu
!====================================================================
!====================================================================
function softmax(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8 , dimension(size(X)) :: softmax
    
    real*8 :: y

    y = exp(-maxval(x))

    softmax = exp(x + log(y))/sum(exp(x(:)+log(y)))

end function softmax
!====================================================================
function dsoftmax(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8 , dimension(size(X)) :: dsoftmax

    real*8 :: y

    y = exp(-maxval(x))

    dsoftmax = (one - exp(x + log(y))/sum(exp(x(:)+log(y))))*(exp(x + log(y))/sum(exp(x(:)+log(y))))
end function dsoftmax
!====================================================================
function sigmoid(x)
!====================================================================
!   "Convert an integer to string."
    real*8, intent(in) :: x(:)
    real*8 , dimension(size(X)) :: sigmoid

    real*8 :: y

    y = exp(-maxval(x))

    sigmoid = one/( one + exp(-x) )
end function sigmoid
!====================================================================
function sign_matriz(x)
!====================================================================
    real*8, intent(in) :: x(:,:)
    real*8 , allocatable :: sign_matriz(:,:)
    real*8 , allocatable :: sign_Is(:,:)

    if(.not. allocated( sign_is ))      allocate( sign_is , mold=x )
    if(.not. allocated( sign_matriz ))  allocate( sign_matriz , mold=x )
    sign_is = one
    sign_matriz = sign( sign_is , x )

    end function
!====================================================================
function sign_vec(x)
!====================================================================
    real*8, intent(in) :: x(:)
    real*8 , allocatable :: sign_vec(:)
    real*8 , allocatable :: sign_Is(:)

    if(.not. allocated( sign_is ))      allocate( sign_is , mold=x )
    if(.not. allocated( sign_vec ))  allocate( sign_vec , mold=x )
    sign_is = one
    sign_vec = sign( sign_is , x )

    end function

end module


