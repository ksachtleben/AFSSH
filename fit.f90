module fitting
contains

  function polyfit(vx, vy, d)
    implicit none
    integer, intent(in)                   :: d
    real*8, dimension(d+1)              :: polyfit
    real*8  , dimension(:), intent(in)    :: vy , vx

    real*8  , dimension(:,:), allocatable :: X
    real*8 , dimension(:,:), allocatable :: XT
    real*8  , dimension(:,:), allocatable :: XTX

    integer :: i, j

    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real*8  , dimension(:), allocatable :: work

    n = d+1
    lda = n
    lwork = n

    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))

    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem1"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem2"
       return
    end if

    polyfit = matmul( matmul(XTX, XT), vy)
    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  end function

end module

