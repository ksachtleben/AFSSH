module griding

use constants_and_parameters

implicit none

!   module variables 

    public :: grade , sumtrap

    contains

!=====================
 subroutine Grade(pos)
!=====================
real*8 , intent(in) :: pos
! variaveis locais ...
integer                 :: i , n , k
real*8  :: xl , xr


! domain 1 ...
xl = pos - 3.d0!*a0
xr = pos + 3.d0!*a0
dx = 0.01d0

grid_size = int(abs(xr - xl)/dx)  + 1

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

end module griding
