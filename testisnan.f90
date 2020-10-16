program test_nan
  implicit none
  real , allocatable :: x(:)
  integer               :: i
  complex*16            :: zi = ( 0.d0 , 1.d0 )
        allocate( x( 10 ) )
        do i =1 , 10
        x(i)  = (zi)**i
        x(i) = sqrt(x(i))
        if (isnan(x(i))) stop '"x" is a NaN'
        print*, i
        enddo
end program test_nan
