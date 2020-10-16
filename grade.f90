module griding

use constants_and_parameters

implicit none

!   module variables 

    public :: grade , sumtrap

    contains

!========================
 subroutine Grade
!========================
! variaveis locais ...
integer                 :: i , n , k
real*8  :: xl , xr


! domain 1 ...
!dx = 20.0d0
xl = posi
xr = posf

grid_size = 4001
!grid_size = aint(abs(xr - xl)/dx)  + 1
dx = aint( abs(xr - xl) )  / ( grid_size - 1 )

if(.not. allocated( x )) allocate( x(grid_size) )

do i = 1 , grid_size
        x(i) = xl + float(i-1)*dx
enddo

end subroutine Grade

!============================
 function sumtrap(i1,i2,y,yp)
!============================
integer , intent(in) :: i1 , i2
real*8  , intent(in) :: y(:)
real*8  , intent(in) :: yp(:)

real*8  :: sumtrap

!------------------------------------------------------------------------------
! CALCULA A INTEGRAL DA FUNCAO Y(I) PELO METODO DO TRAPEZIO COM PASSO VARIAVEL
!------------------------------------------------------------------------------

sumtrap  = sum( (y(i1+1:i2)-y(i1:i2-1)) * (yp(i1+1:i2)+yp(i1:i2-1)) ) / two

end function sumtrap

!==============================
subroutine phase_gauss( vec )
!==============================
real*8 , intent(inout) :: vec(:)

integer , dimension(:) , allocatable :: seed

integer :: size_vec , i , n
real*8 , parameter :: gama = 1.7253d0
real*8 , parameter :: variance = 1.0d0   ! sigma²
real*8 , parameter :: mean = 0.d0       ! mu

real*8 :: gauss_f , phase_rn , u , v 

size_vec = size(vec)

call date_and_time( values=values )
call random_seed(size = n)
allocate(seed(1:n))
seed(:) = values(8)
call random_seed(put=seed)

do i = 1 , size_vec
        do 

        call random_number(u)
        !phase_rn = (2.d0*u - 1.d0)
        phase_rn = u
        
        call random_number(v)
        gauss_f = v
        
        if( gauss_f <= exp( -half*(phase_rn - mean)**2/variance )/( sqrt(two*pi*variance) ) ) exit

        enddo
        vec(i) = phase_rn
enddo 

end subroutine

!==============================
subroutine random_normal( vec )
!==============================
! Adapted from the following Fortran 77 code
! ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
! THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
! The function random_normaL() returns a normally distributed pseudo-random number with zero mean and unit variance.
! The algorithm uses the ratio of uniforms method of A.J. Kinderman and J.F. Monahan augmented with quadratic bounding curves.
implicit none
real*8  , intent(inout) :: vec(:)

! local variables ...
real*8  :: q , u , v , x , y
integer :: i , n , m

! local seed 
integer , dimension(:) , allocatable :: seed

! local parameter ...
real*8  , parameter :: s =   0.449871d0
real*8  , parameter :: t = - 0.386595d0
real*8  , parameter :: a =   0.196000d0
real*8  , parameter :: b =   0.254720d0
real*8  , parameter :: r1 =  0.275970d0
real*8  , parameter :: r2 =  0.278460d0

m = size( vec )

call date_and_time( values=values )
call random_seed(size = n)
if(.not. allocated(seed) ) allocate(seed(1:n))
seed(:) = values(8) !put 8 means take the seconds as seed
call random_seed(put=seed)

! Generate P = (u,v) uniform in rectangle enclosing acceptance region
do i = 1 , m

    do

        CALL random_number( u )
        CALL random_number( v )

        v = 1.7156 * ( v - half )

!       evaluate the quadratic form
        x = u - s
        y = Abs( v ) - t
        q = x * x + y * ( a * y - b * x )

!       accept P if inside inner ellipse
        if( q < r1 ) exit

!       reject P if outside outer ellipse
        if( q > r2 ) cycle

!       reject P if outside acceptance region
        if( v * v < - 4.0d0 * dlog( u ) * u * u ) exit

!       reject values greater than one
!        if( abs(v/u) < 1.d0 ) exit
    end do

!   return ratio of P's coordinates as the normal deviate
    vec( i ) = v / u

end do
end subroutine random_normal

end module griding
