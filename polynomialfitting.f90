module polynomialfitting  

use fitting
    implicit none
 
  contains
  !===========================================
  Subroutine fitting_pol( pot )
  !===========================================
  real*8 , intent(inout) :: pot(:,:)

  ! let us test it
  integer, parameter      :: degree = 3
  integer                 :: i , io1 , j
  real*8                  :: k , va , vb, vc, vd 
  real*8 , allocatable :: x(:) 
  real*8 , allocatable :: y(:)
  real*8 , allocatable :: funcx(:)
  real*8 , allocatable :: guncx(:)

  real*8 , dimension(degree+1) :: a
  
  integer  , parameter :: nlines = 10

  if( .not. allocated(y) ) allocate( y(nlines) , source=0.d0)
  if( .not. allocated(x) ) allocate( x(nlines) , source=0.d0)
  if( .not. allocated(funcx) ) allocate( funcx(degree+1) )
  if( .not. allocated(guncx) ) allocate( guncx(degree+1) )

  open( 31 , file='prob_trans_up.dat' )
  i = 1
  do
    read( 31 , * , iostat=io1 ) va , vb !, vc !, vd
    if( io1 < 0 ) exit
    y(i) = vb
    !x(i) = log(va*va/4000.d0)
    x(i) = va
  write( 14 , * ) x(i) , y(i)
    i = i + 1
    print*, va,vb
  enddo
  a = polyfit(x, y, degree)
 
  write ( * , *  ) a

    
  do j = 1 , 1050
   k = -4.d0 + 5.d0*float(j)/1050.d0
   !k = 1d0 + 1.d0*float(j)/1050.d0
    do i = 1 , degree+1
      funcx(i) =  a(i)*k**(i-1)  
      guncx(i) = -float(i-1)*a(i)*k**(i-2)
    enddo
   !write( 15 , * )  k , a(6)*k**5 + a(5)*k**4 + a(4)*k**3 + a(3)*k**2 + a(2)*k**1 + a(1)
   write( 16 , * ) k , sum(funcx(:))
   write( 17 , * ) k , sum(guncx(:))
   enddo

close(31)
100 format(f5.2,f5.1)
13 format(f17.16)
end subroutine

end module
