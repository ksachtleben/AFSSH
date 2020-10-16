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
 Subroutine dinamic_nuclei( ni )!, interaction_e_N )
!=======================================================================================================
integer , intent(inout) :: ni 
!logical , intent(out)   :: interaction_e_N 
!-------------------------------------------------

! local variables
integer :: i , n  , j , m , x_in , k , nii , nif , io1, io2 , nf_test , iounit , regiao , regiao_ant

real*8     , allocatable :: coefi(:) , rho2_t(:) , temp(:) , time_step(:) , temp_de(:) , dtime_step(:) , alg_eign(:) , E_M(:,:) , local(:)

complex*16 , allocatable ::  c_phi(:)

real*8 :: norm , t_zero ,a,b,c,d , p_lz , inflex_denerg , pos_in , dx_test , velocity , pos_f , kxin , velo , dv , dtime
logical :: transmission , reflection , int_e_N , save_


!interaction_e_N = .false.

call grade
if(.not. allocated( coefi ))            allocate( coefi(grid_size) )
if(.not. allocated( rho2_t ))           allocate( rho2_t(npt) )
if(.not. allocated( c_phi       ))      allocate( c_phi(nl) )
if(.not. allocated( eigenenergy ))      allocate( eigenenergy(nl) )
if(.not. allocated( potential ))        allocate( potential(nl,nl) )
if(.not. allocated( E_m ))        allocate( E_M(nl,nl) )
if(.not. allocated(alg_eign))   allocate( alg_eign(nl) )
if(.not. allocated(Ntr2))       allocate( NTr2(grid_size) )
if(.not. allocated(Ntr1))       allocate( NTr1(grid_size) )
if(.not. allocated(Nrefl2))     allocate( Nrefl2(grid_size) )
if(.not. allocated(Nrefl1))     allocate( Nrefl1(grid_size) )
if(.not. allocated(Nconf))      allocate( Nconf(grid_size) )
if(.not. allocated(local))      allocate( local(grid_size) , source = 0.d0 )
if(.not. allocated(x_toll))     allocate( x_toll(grid_size) , source = 0.d0 )
!===========================
! classical/quantum dynamic
!===========================

        !First step is calculate the acelaration

        call potential_grid( eigenenergy , potential , type_potential )
 
        acel    =   force_HF(ni)/massp

        dv      =   half*acel*dt

        vel     =   vel + dv

        dtime = half * dt
        call dinamic_el( dtime )
        count_step = 1
!        do while( pos <= posf .and. pos >= posi )
        transmission    = .false.
        reflection      = .false.
        int_e_n         = .false.
        regiao          = 1
        regiao_ant      = 0

        do while( end_step == .false. )

            if( count_step > max_steps ) end_step = .true.

                do i = regiao_ant + 1 , regiao

                        if ( ni == 2 )  Ntr2( i ) = Ntr2( i ) + one
                        if ( ni == 1 )  Ntr1( i ) = Ntr1( i ) + one

                enddo


           if( vel < 0.d0 ) then
                
                do i = regiao + 1 , grid_size

                        if ( ni == 2 )  Nrefl2( i ) = Nrefl2( i ) + one
                        if ( ni == 1 )  Nrefl1( i ) = Nrefl1( i ) + one

                end do

                end_step = .true.

           end if

           if( regiao == grid_size ) exit
           pos  = pos + vel*dt 

           regiao_ant = regiao
           do while( pos > x( regiao ) ) 
                regiao = regiao + 1 
                if( regiao == grid_size + 1 ) exit
           enddo
           regiao = regiao - 1 

           if( regiao_ant - regiao - 1 == 0 .and. end_step == .false. ) print*, "careful"

            call potential_grid( eigenenergy , potential , type_potential )
            
            acel    = force_HF(ni)/massp
            
            vel     = vel + acel*dt*half
            
            E_M(1,1) = enp(1)
            E_M(2,2) = enp(2)
            E_M(2,1) = 0.0d0
            E_M(1,2) = 0.0d0

            dR_fg   = dR_fg + dt * dP_fg / massp
            dP_fg   = dP_fg + dt * half * ( matmul( dF_fg , rho_m ) + matmul( rho_m , dF_fg ) )

            count_step = count_step + 1 

                    trans_prob = .false. 
                    call dinamic_el( dt )
                    
                    nii = ni
                    call hopping( rho_m , rho_mold , j_m , save_ , regiao )
                    nif = ni
        
                    select case( trans_prob )
                        case( .true. )
                            dtime   = float(i)*dt/float(nsteps_el)
                            ni      = nif
                        case( .false. )
                            dtime   = dt
                            ni      = nii
                    endselect

            vel     = vel + acel*dt*half
            time    = time + dt

        enddo

13 format(4es16.4E3)
end subroutine

!==============================================
 Subroutine draft_psi1( x_i , k_i , psi , step)
!==============================================
!complex*16 , intent(in) :: psi(:)
integer    , intent(in) :: step
real*8     , intent(in) :: psi(:,:) , x_i(:) , k_i(:)
!-------------------------------------------------
integer :: m , n , i, j
character( len=255 ) :: potent

m = size(k_i)
n = size(x_i) - 1
potent = str( step )
!========================================================================
! ESCREVE AS FUNÇÕES QUE SERÃO DESENHADAS NO GNUPLOT.
!========================================================================
open(20 , file="plot"//trim(potent)//".dat" )
do j = 1 , m ! k é m
    do i = 1 , n ! x é n
        write( 20 , 13 ) x_i(i+1) , k_i(j) , psi(j,i)
    end do
end do
close(20)
13  format(3f20.12)
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

13  format(3f20.12)
!==========================================
end subroutine

end module
