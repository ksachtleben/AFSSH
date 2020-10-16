module energy_surface

use constants_and_parameters
use griding                     , only : grade , sumtrap
use fitting

public

contains
!==============================================================================
 Subroutine potential_grid( eigenenergy , potential , type_potential )
!==============================================================================
integer , intent(in)  :: type_potential
real*8 , intent(out)  :: eigenenergy(:) , potential(:,:) 
!------------------------------------------------------------------------------
integer , parameter :: degree = 5

integer :: i , j , k, l , n , io4 , io5 , iounit

real*8 :: Ca
real*8 :: Cb
real*8 :: Cc
real*8 :: Cd
real*8 :: E0

real*8 , allocatable    :: norm_vec(:) , phi_ao(:,:) , soma(:) , base_loc(:,:) , base_des(:,:), aloc_matriz(:,:) , dpotential(:,:) 
real*8 , allocatable    :: var_x(:) , var_y(:)
real*8 :: theta , var_phi , dtheta , dvar_phi , dE , energia , va , vb

real*8 , dimension(degree+1) :: a , pot

if(.not. allocated( dpotential ))               allocate( dpotential(nl,nl) )
if(.not. allocated( aloc_matriz ) )             allocate( aloc_matriz( nl , nl ) )
if(.not. allocated( base_loc ) )                allocate( base_loc( nl , nl ) )
if(.not. allocated( base_des ) )                allocate( base_des( nl , nl ) )
if(.not. allocated( norm_vec ) )                allocate( norm_vec( grid_size ) )
if(.not. allocated( coef_coupl ))               allocate( coef_coupl( nl , nl ) )
if(.not. allocated( coef_coupl_af ))            allocate( coef_coupl_af( nl , nl  ) )
if(.not. allocated( rho_coupl ))                allocate( rho_coupl( nl , nl ) )
if(.not. allocated( nonadiabatic_coupling) )    allocate( nonadiabatic_coupling( nl,nl ) )
if(.not. allocated( enp ))                      allocate( enp(nl) )
if(.not. allocated( var_y ))                    allocate( var_y(10) )
if(.not. allocated( var_x ))                    allocate( var_x(10) )
if(.not. allocated( force_HF ))                 allocate( force_HF(nl) )
if(.not. allocated( dF_fg ))                 allocate( dF_fg(nl,nl) )


select case( type_potential )
    case(1)
        Ca = 0.01d0
        Cb = 1.6d0
        Cc = 0.005d0
        Cd = 1.0d0

        if( pos >= 0.d0 ) potential(1,1)    = Ca*( 1.d0 - exp( -Cb*pos ))
        if( pos < 0.d0 )  potential(1,1)    = -Ca*( 1.d0 - exp(  Cb*pos ))
        potential(2,2)    = -potential(1,1)
        potential(1,2)    = Cc*exp( -Cd*pos**2 )
        potential(2,1)    = potential(1,2)

        if( pos >= 0.d0 ) dpotential(1,1)   = Ca*(exp( -Cb*pos ))*Cb
        if( pos < 0.d0 )  dpotential(1,1)   = Ca*(exp(  Cb*pos ))*Cb
        dpotential(2,2)   = -dpotential(1,1)
        dpotential(1,2)   = -Cc*exp( -Cd*(pos)**2 )*two*Cd*pos
        dpotential(2,1)   = dpotential(1,2)
                                                
    case(2)
        Ca = 0.10d0
        Cb = 0.28d0
        Cc = 0.015d0
        Cd = 0.06d0
        E0 = 0.05d0

        potential(1,1)    = 0.d0
        potential(2,2)    = -Ca*exp( -Cb*(pos)**2 ) + E0
        potential(1,2)    = Cc*exp( -Cd*(pos)**2 )
        potential(2,1)    = potential(1,2)
    
        dpotential(1,1)   = 0.d0
        dpotential(2,2)   = Ca*exp( -Cb*(pos)**2 )*Cb*two*pos
        dpotential(1,2)   = -Cc*exp( -Cd*(pos)**2 )*Cd*pos*two
        dpotential(2,1)   = dpotential(1,2)

    case(3)
        Ca = 6.e-4
        Cb = 0.10d0
        Cc = 0.9d0

        if( pos >= 0.d0 ) potential(1,2) = Cb*( 2.d0 - exp( -Cc*pos ))
        if( pos < 0.d0 )  potential(1,2) = Cb*( exp( Cc*pos ))
        potential(1,1) = -Ca
        potential(2,2) = Ca
        potential(2,1) = potential(1,2)

        if( pos >= 0.d0 ) dpotential(1,2) = Cb*( exp( -Cc*pos ))*Cc
        if( pos < 0.d0 )  dpotential(1,2) = Cb*( exp( Cc*pos ))*Cc
        dpotential(1,1) = 0.d0
        dpotential(2,2) = 0.d0
        dpotential(2,1) = dpotential(1,2)

    case(4)
        Ca = 6.0d-4
        Cb = 1.0d-1
        Cc = 9.0d-1
        Cd = 10.0d0
        potential(1,1) = -Ca
        potential(2,2) = Ca

        if( pos <= -Cd )                        potential(1,2) = Cb * dexp( Cc * ( pos - Cd ) ) + Cb * ( two - dexp( Cc * ( pos + Cd ) ) )
        if( -cd < pos .and. pos < cd )          potential(1,2) = Cb * dexp( Cc * ( pos - Cd ) ) + Cb * dexp( -Cc * ( pos + Cd ) )
        if( pos >= cd )                        potential(1,2) = Cb * ( two - dexp( -Cc * ( pos - Cd ) ) ) + Cb * dexp( -Cc * ( pos + Cd ) )  
        potential(2,1) = potential(1,2)

        dpotential(1,1) = 0.d0
        dpotential(2,2) = 0.d0

        if(  pos <= - Cd )                      dpotential(1,2) = Cb * Cc * ( dexp( Cc * ( pos - Cd ) ) - dexp( Cc * ( pos + Cd ) ) )
        if( -Cd < pos .and. pos <  Cd )         dpotential(1,2) = Cb * Cc * dexp( -Cc * ( pos + Cd ) ) * ( dexp( two * Cc * pos ) - one )
        if(  pos >=  Cd )                       dpotential(1,2) = Cb * Cc * dexp( -Cc * ( pos + Cd ) ) * ( dexp( two * Cc * Cd ) - one )

        dpotential(2,1) = dpotential(1,2)

    case(5)
        Ca = 6.0d-4
        Cb = 1.0d-1
        Cc = 9.0d-1
        Cd = 4.0d0
        potential(1,1) = -Ca
        potential(2,2) = Ca

        if( pos <= -Cd )                        potential(1,2) = -Cb * dexp( Cc * ( pos - Cd ) ) + Cb * dexp( Cc * ( pos + Cd ) ) 
        if( -cd < pos .and. pos < cd )          potential(1,2) = -Cb * dexp( Cc * ( pos - Cd ) ) - Cb * dexp( -Cc * ( pos + Cd ) ) + two*Cb
        if( pos >= cd )                         potential(1,2) = Cb * dexp( -Cc * ( pos - Cd ) ) - Cb * dexp( -Cc * ( pos + Cd ) )  
        potential(2,1) = potential(1,2)
       
        dpotential(1,1) = 0.d0
        dpotential(2,2) = 0.d0

        if(  pos <= - Cd )                      dpotential(1,2) = Cb * Cc * two * dexp( Cc * pos ) * dsinh( Cc * Cd )
        if( -Cd < pos .and. pos <  Cd )         dpotential(1,2) = -Cb * Cc * dexp( -Cc * ( pos + Cd ) ) * ( dexp( two * Cc * pos ) - one )
        if(  pos >=  Cd )                       dpotential(1,2) = -Cb * Cc * dexp( -Cc * ( pos + Cd ) ) * ( dexp( two * Cc * Cd ) - one )

        dpotential(2,1) = dpotential(1,2)

end select
    
    eigenenergy(1)= (potential(1,1) + potential(2,2))/two -sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)/two
    eigenenergy(2)= (potential(1,1) + potential(2,2))/two +sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)/two


    coef_coupl(1,1) = sqrt( half*( one + ( potential(2,2) - potential(1,1) )/sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)))
    coef_coupl(2,1) = -sqrt( half*( one - ( potential(2,2) - potential(1,1) )/sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)))
    coef_coupl(1,2) = sqrt( half*( one - ( potential(2,2) - potential(1,1) )/sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)))
    coef_coupl(2,2) = sqrt( half*( one + ( potential(2,2) - potential(1,1) )/sqrt( (potential(1,1) - potential(2,2) )**2+ four*(potential(1,2))**2)))

aloc_matriz = matmul(transpose(coef_coupl),matmul(dpotential,coef_coupl))

do i = 1 , nl
    do j = 1 , nl
    dE = eigenenergy(j) - eigenenergy(i)
    nonadiabatic_coupling(i,j) = aloc_matriz(i,j)/dE
    enddo
    nonadiabatic_coupling(i,i) = 0.d0
    enp(i) = eigenenergy(i)
enddo

forall( i = 1 : nl) force_HF(i) = -aloc_matriz(i,i)


dF_fg(1,1) = force_HF(1) 
dF_fg(2,2) = force_HF(2)
dF_fg(1,2) = 0.d0
dF_fg(2,1) = 0.d0

!if( finish_tully == .false. ) then
!if( check_time_initial == .true. ) then
!write( 24 , 13 ) pos , nonadiabatic_coupling(1,2)
!write( 26 , 13 ) pos , force_HF(1)
!write( 27 , 13 ) pos , force_HF(2)
!write( 25 , 13 ) pos , nonadiabatic_coupling(2,1)
!write( 90  , 13 ) pos , eigenenergy(1)
!write( 91  , 13 ) pos , eigenenergy(2)
!write( 92  , 13 ) x(pos_index) , potential(1,2,pos_index)/(vel*abs(dpotential(2,2) - dpotential(1,1)))
!endif

!write( 94  , 13 ) time , cos(theta/two)**2*potential(1,1,pos_index) + sin(theta/two)**2*potential(2,2,pos_index) + two*cos(theta/two)*sin(theta/two)*var_phi**2*potential(1,2,pos_index)
write( 105 , 13 ) pos , eigenenergy(1)
write( 106 , 13 ) pos , eigenenergy(2)
write( 101 , 13 ) pos , potential(1,1)
write( 102 , 13 ) pos , potential(1,2)
write( 103 , 13 ) pos , potential(2,1)
write( 104 , 13 ) pos , potential(2,2)
!write( 111 , 13 ) x(pos_index) , dpotential(1,1)
!write( 112 , 13 ) x(pos_index) , dpotential(1,2)
!write( 113 , 13 ) x(pos_index) , dpotential(2,1)
!write( 114 , 13 ) x(pos_index) , dpotential(2,2)

!endif


13 format(3es16.4E3)
end subroutine potential_grid

end module
