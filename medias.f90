module mediadat
use constants_and_parameters


contains
!=======================================
subroutine medias_sub_dat(npt,n_of_traj)
!=======================================
integer, intent(in) :: npt , n_of_traj

!local variables
integer :: j , io1 , io2 
real*8  :: va , vb , vd
complex*16 :: ve , vf

if( .not. allocated(vs1)) allocate(vs1(100) ,source=0.d0)
if( .not. allocated(vs2)) allocate(vs2(100) ,source=0.d0)
if( .not. allocated(vs3)) allocate(vs3(100) ,source=0.d0)
if( .not. allocated(vs4)) allocate(vs4(100) ,source=0.d0)


!    open( 31 , file='pop1.dat' , action='read' )
    open( 32 , file='sumpop1.dat' , action='write' )
!    open( 35 , file='pop2.dat' , action='read' )
    open( 36 , file='sumpop2.dat' , action='write'  )
!    open( 12 , file='coeh12.dat' , action='read' )
    open( 13 , file='sumcoeh12.dat' , action='write'  )
!    open( 21 , file='coeh21.dat' , action='read' )
    open( 22 , file='sumcoeh21.dat' , action='write' )
    j = 1
    do
        read( 31 , 100 , iostat=io1 ) va , vb
        if( io1 < 0 ) exit
        vs1(j) = vs1(j) + vb/float(n_of_traj)
        write( 32 , 13 ) va , vs1(j)
        j = j + 1
    enddo
    j = 1
    do
        read( 35 , 100 , iostat=io2 ) va , vd
        if( io2 < 0 ) exit
        vs2(j) = vs2(j) + vd/float(n_of_traj)
        write( 36 , 13 ) va , vs2(j)
        j = j + 1
    enddo
    j = 1
    do
        read( 12 , 100 , iostat=io2 ) va , ve
        if( io2 < 0 ) exit
        vs3(j) = vs3(j) + ve/float(n_of_traj)
        write( 13 , 13 ) va , vs3(j)
        j = j + 1
    enddo
    j = 1
    do
        read( 21 , 100 , iostat=io2 ) va , vf
        if( io2 < 0 ) exit
        vs4(j) = vs4(j) + vf/float(n_of_traj)
        write( 22 , 13 ) va , vs4(j)
        j = j + 1
    enddo
!    print*, vb + vd , vb , vd !"aqui"
close(31)
close(32)
close(35)
close(36)
close(12)
close(13)
close(21)

13  format(3es16.4E3)
100 format(3es16.4E3,a5,3es16.4E3)


end subroutine

!=====================================================================================
 subroutine Trans_or_Reflex( pos , coef_down ,coef_up , state_in )
!=====================================================================================
real*8      , intent(in) :: pos , coef_down , coef_up
integer                  :: state_in

real*8  , allocatable    :: checkup(:) , checkdown(:) , temp(:)
integer                  :: upper_state, down_state

character (len=6)        :: statei , statef

Logical , parameter      :: T_ = .true. , F_ = .false.
logical                  :: state_end

if(.not. allocated(temp) ) allocate( temp(ntraj) )
!if(.not. allocated(ch_state_count) ) allocate(
!if(.not. allocated(re_state_count) ) allocate(
!if(.not. allocated(tr_state_count) ) allocate(

upper_state = nint( coef_up )
down_state  = nint( coef_down )

select case( state_in )
    case( 1 )
        statei = "down"
        select case( upper_state )
            case( 1 )
                statef = "upper"

                statef = "down"

                print*, pos !pos > 0 transmission, pos < 0 reflexion
                print*, "reflexion or transmission"!can be reflected or transmiison check you pos
        endselect
    case( 2 )
        statei = "upper"
        select case( down_state )
            case( 1 )
                statef = "down"
                print*, "change state" !change of state
            case( 0 )
                statef = "upper"

                print*, pos !pos > 0 transmission, pos < 0 reflexion
                print*, "reflexion or transmission"!can be reflected or transmiison check you pos
        endselect
endselect

13 format(3es16.4E3)
end subroutine Trans_or_Reflex


end module
