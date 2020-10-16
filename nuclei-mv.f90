module nuclear_din

use f95_precision
use blas95
use lapack95

use constants_and_parameters
use griding                     , only : grade , sumtrap
use energy_surface
use el_din
use parameters_confg

implicit none 

public !:: dinamic_nuclei

contains
!=======================================================================================================
 Subroutine dinamic_nuclei( ni )
!=======================================================================================================
integer , intent(inout) :: ni
!-------------------------------------------------

! local variables
integer :: i , n  , j , m , x_in , k , nii , nif , io1, io2 , nf_test , iounit

real*8     , allocatable :: coefi(:) , rho2_t(:) , temp(:) , time_step(:) , temp_de(:) , dtime_step(:)

complex*16 , allocatable ::  c_phi(:)

real*8 :: norm , t_zero ,a,b,c,d , p_lz , inflex_denerg , pos_in , dx_test , velocity , pos_f , kxin , last_velocity , dv

logical :: interaction_e_N = .false.


!real*8  :: matrix1(2,2) , matrix2(2,2) , matrix3(2,2)
!matrix1(1,1) = 1.d0
!matrix1(1,2) = 2.d0
!matrix1(2,1) = -1.d0
!matrix1(2,2) = 3.d0

!matrix2(1,1) = 2.d0
!matrix2(1,2) = -1.d0
!matrix2(2,1) = -1.d0
!matrix2(2,2) = 1.d0
!call gemm( matrix1 , matrix2 , matrix3 )

!print*, matrix3(1,1)
!print*, matrix3(1,2)
!print*, matrix3(2,1)
!print*, matrix3(2,2)
!stop

call grade( posi )
if(.not. allocated( coefi ))            allocate( coefi(grid_size) )
if(.not. allocated( rho2_t ))           allocate( rho2_t(npt) )
if(.not. allocated( c_phi       ))      allocate( c_phi(nl) )
if(.not. allocated( eigenenergy ))      allocate( eigenenergy(nl,grid_size) )
if(.not. allocated( potential ))        allocate( potential(nl,nl,grid_size) )
!First step is calculate the acelaration
select case (check_time_initial)
!========================
! only classical dynamic
!========================
    case( .true. )
            pos_in = posi
            if(.not. allocated( temp ))             allocate( temp( 5000000 ) , source=0.d0 )
            call potential_grid( eigenenergy , potential , type_potential , inflex_denerg)
                acel = -( eigenenergy(ni,pos_index+1 ) - eigenenergy(ni,pos_index-1) )/(massp*two*dx)
                j = 0
                do while( abs(posi) <= posf )
                call grade( posi ) 
                call potential_grid( eigenenergy , potential , type_potential , inflex_denerg ) 
                j = j + 1
                temp(j) = abs(eigenenergy(2,pos_index) - eigenenergy(1,pos_index))
                   time = time + dt
                   posi = posi + vel*dt + acel*dt**2/two
                call grade( posi ) 
                call potential_grid( eigenenergy , potential , type_potential , inflex_denerg ) 
                   if( abs(posi) > posf ) tfinal = time
                   vel  = vel + acel*dt/two
                   acel = -( eigenenergy(ni,pos_index+1 ) - eigenenergy(ni,pos_index-1) )/(massp*two*dx)
            
                   vel  = vel + acel*dt/two
                if( j > 10000 ) then
                    print*, "ok"
                    exit
                endif
               enddo
            posi = pos_in
            npt = j
            if(.not. allocated( time_step ) )       allocate( time_step( npt ) , source=temp )
            deallocate( temp )
            time = 0.d0
            de_min = minval( time_step )
            de_max = maxval( time_step )
            de_i   = time_step(1)
            deallocate( time_step )
!===========================
! classical/quantum dynamic
!===========================
    case( .false. )
        call grade( pos ) 
        call potential_grid( eigenenergy , potential , type_potential , inflex_denerg ) 

            acel = -( eigenenergy(ni,pos_index+1 ) - eigenenergy(ni,pos_index-1) )/(massp*two*dx)
            
        !do passo = 1 , npt
        dx_pos = 0.d0
        do while( pos <= posf .and. pos >= posi )
        count_step = 1 
               
!        dt = abs(eigenenergy(2,pos_index) - eigenenergy(1,pos_index))/(two*pi*de_max*de_min)       
!         dt   = abs( 0.05d0/vel )
            dt = 1.d0
            nii = ni
            
            if( pos > -5.d0 .and. pos < 5.d0 ) then
                interaction_e_N = .true.
                call dinamic_el
            else
                interaction_e_N = .false.
            endif

            nif = ni
            if( nii /= nif ) then
                ni = nif
            else
                ni = nii
            endif

            energyconserv = energyconserv + (eigenenergy(ni,2) + massp*vel**2/(two))/float(npt*n_of_traj)
            pos  = pos + vel*dt + acel*dt**2/two 

        call grade( pos ) 
        call potential_grid( eigenenergy , potential , type_potential , inflex_denerg ) 

            vel  = vel + acel*dt/two
            acel = -( eigenenergy(ni,pos_index+1 ) - eigenenergy(ni,pos_index-1) )/(massp*two*dx)
            vel  = vel + acel*dt/two
            time = time + dt

                write( 12, 13 ) time , pos
                write( 13, 13 ) pos , vel
                write( 14, 13 ) pos , energyconserv
            
            if( time > max_steps ) then
                select case( interaction_e_N )
                case( .true.  )
                    count_inter_conf = count_inter_conf + 1
                    pos = sign( posi , pos )
                case( .false. )
                    max_steps = two*max_steps
                endselect
            endif
        enddo
endselect

13 format(4es16.4E3)
end subroutine

!============================================
 Subroutine draft_psi1( psi , step )
!============================================
!complex*16 , intent(in) :: psi(:)
real*8     , intent(in) :: psi(:)
integer    , intent(in) :: step
integer                 :: i 
!-------------------------------------------------

!========================================================================
! ESCREVE AS FUNÇÕES QUE SERÃO DESENHADAS NO GNUPLOT.
!========================================================================
open(20 , file="plot1.gnu" )
 write(20,fmt='(a)')        'clear'
 write(20,fmt='(a)')        'reset'
 write(20,fmt='(a)')        'unset key'
 write(20,fmt='(a)')        "set title'Dinâmica do Pacote'"
 write(20,fmt='(a)')        "set terminal pngcairo size 1000,600 enhanced font 'Verdana,20'"
 write(20,fmt='(a,f6.1,a,f6.1,a)') "set xrange[",x(1),":",x(grid_size),"]"
! write(20,fmt='(a,f6.1,a,f6.1,a)') "set xrange[-6.0e-7:6.0e-7]"
 write(20,fmt='(a)')        "set yrange[0:0.4]"
 write(20,fmt='(a)')        "set y2range[0:0.4]"
 write(20,fmt='(a)')        "set xlabel'Comprimento'"
 write(20,fmt='(a)')        "set y2label'Energy'"
 write(20,fmt='(a)')        "set ylabel '20*{/Symbol y} + E_{package}' textcolor rgb 'blue'"
 write(20,fmt='(a)')        "set xtics 100"
!\/ ja estava comentado
! write(20,fmt='(a)')        'set ytics 300 textcolor rgb "blue"'
!/\ ja estava comentado
 write(20,fmt='(a)')        "set y2tics 250"
 write(20,fmt='(a)')        "set pointsize 0.005"
 write(20,fmt='(a)')        'set border linewidth 1.0'
 write(20,fmt='(a)')        "set samples 300000"
 write(20,fmt='(a,i0,a)')        "system('mkdir -p 1state')"


!\/ ja estava comentado
!    write(*,*)  '   Frame = ' , tlinha
!/\ ja estava comentado
    write(20,fmt='(a,i0,a,i0,a)') "outfile = sprintf('1state/geral%05.0f.png',",step,")"
    write(20,fmt='(a)')      "set output outfile"
    write(20,fmt='(a,i3,a,i3,a)') "plot '-' u 1:2 w l lw 4 lc 3 , '-' u 1:4 w l lw 4 lc 0"

!=================================================
!=================================================
do i = 1 , grid_size
    write( 20 , 13 ) x(i) , psi(i)!*conjg(psi(i))
enddo

 write(20,fmt='(a)') "unset output outfile"
 write(20,fmt='(a)') 'e'
 write(20,fmt='(a)') ' '

13 format(3es16.4E3)
!==========================================
end subroutine

!============================================
 Subroutine draft_psi2( psi , step )
!============================================
!complex*16 , intent(in) :: psi(:)
real*8     , intent(in) :: psi(:)
integer    , intent(in) :: step
integer                 :: i 

!========================================================================
! ESCREVE AS FUNÇÕES QUE SERÃO DESENHADAS NO GNUPLOT.
!========================================================================
open(21 , file="plot2.gnu" )
 write(21,fmt='(a)')        'clear'
 write(21,fmt='(a)')        'reset'
 write(21,fmt='(a)')        'unset key'
 write(21,fmt='(a)')        "set title'Dinâmica do Pacote'"
 write(21,fmt='(a)')        "set terminal pngcairo size 1000,600 enhanced font 'Verdana,20'"
 write(21,fmt='(a,f6.1,a,f6.1,a)') "set xrange[",x(1),":",x(grid_size),"]"
 !write(21,fmt='(a,f6.1,a,f6.1,a)') "set xrange[-6.0e-7:6.0e-7]"
 write(21,fmt='(a)')        "set yrange[0:0.4]"
 write(21,fmt='(a)')        "set y2range[0:0.4]"
 write(21,fmt='(a)')        "set xlabel'Comprimento'"
 write(21,fmt='(a)')        "set y2label'Energy'"
 write(21,fmt='(a)')        "set ylabel '20*{/Symbol y} + E_{package}' textcolor rgb 'blue'"
 write(21,fmt='(a)')        "set xtics 100"
!\/ ja estava comentado
! write(20,fmt='(a)')        'set ytics 300 textcolor rgb "blue"'
!/\ ja estava comentado
 write(21,fmt='(a)')        "set y2tics 250"
 write(21,fmt='(a)')        "set pointsize 0.005"
 write(21,fmt='(a)')        'set border linewidth 1.0'
 write(21,fmt='(a)')        "set samples 300000"
 write(21,fmt='(a,i0,a)')        "system('mkdir -p 2state')"


!\/ ja estava comentado
!    write(*,*)  '   Frame = ' , tlinha
!/\ ja estava comentado
    write(21,fmt='(a,i0,a,i0,a)') "outfile = sprintf('2state/geral%05.0f.png',",step,")"
    write(21,fmt='(a)')      "set output outfile"
    write(21,fmt='(a,i3,a,i3,a)') "plot '-' u 1:2 w l lw 4 lc 3 , '-' u 1:4 w l lw 4 lc 0"

!=================================================
!=================================================
do i = 1 , grid_size
    write( 21 , 13 ) x(i) , psi(i)!*conjg(psi(i))
enddo

 write(21,fmt='(a)') "unset output outfile"
 write(21,fmt='(a)') 'e'
 write(21,fmt='(a)') ' '

13 format(3es16.4E3)
!==========================================
end subroutine

end module
