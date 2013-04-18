
program wavepacket

  implicit none


!! HARDCODED DIMENSIONS, AXES, and STEPS
  integer, parameter :: d = 1          ! Dimensionality of the space
  integer, parameter :: LENGTH = 10    ! Total length of axes considered
  integer :: R0(d) = 0                 ! D-dimensional location of origin

  real(8) :: dt = .0001, dr = .01      ! time step and space step
  integer :: TIMEINIT = 0              ! Initial time of numerical method
  integer :: TIMEFINAL = 1000          ! Length of time to run numerical method

  real(8),parameter :: MU=1, SIGMA=1   ! norm coefficient, variance
  real(8),parameter :: k(d) = 1        ! wavevector with dimension D



!! Other Declarations
  integer :: i,j                       ! loop variables
  integer :: SIZE                      ! Size of position lattice
  real(8) :: r(d)                      ! D-dimensional position vector (a particular location in the position lattice
  real(8), allocatable :: init(:,:)    ! initial condition lattice (initial form of wave at t=0)
  complex :: complex = (0,1)           ! complex unit i



!! Initialization

  SIZE = LENGTH/dr                     ! total length of axis / position step size
  allocate( init(d,SIZE) )             ! Row: Dimensionality, Column: Size of position lattice
  init = 0                             ! Initialize initial condition lattice to 0


!! Initial Condition (Gaussian)

  r = R0                                ! set initial position of lattice at the origin

  do i=1,d                              ! initialize each dimension of lattice
     do j=1,SIZE                        ! run through each position step in lattice
        init(i,j) = MU*exp( -( (r(i)-R0(i))**2 )/(2*SIGMA**2) ) * exp( complex*k(i)*r(i) )
        r(i)=r(i)+dr                    ! increment through all positions of lattice
     enddo
  enddo



!! Perform Numerical Analysis

  call opentextfiles

!! call subroutines

  call closetextfiles



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! How do I include normalization??



contains

  subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="Euler.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, Euler.dat file not opened properly------------"
    endif
    OPEN(UNIT=20,FILE="RungeKutta.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, RungeKutta.dat file not opened properly------------"
    endif
  endsubroutine opentextfiles

  subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=20)
  endsubroutine closetextfiles

!!=================================================================================

  subroutine euler(y_0,a,b,h)
    implicit none
    real(8), intent(in):: a, b      ! Initial and final times
    real(8), intent(in):: h         ! Timestep
    real(8), intent(in):: y_0       ! Initial condition
    real(8) :: y,t                  ! Solution (wavefunction) from previous timestep, time
    
    t = a          
    y = y_0     
    
    do while(t <= b)
       y = y + h*f(y,t)
       
       write(15,*)t, y
       t=t+h
    end do
  end subroutine euler

!!=================================================================================

  subroutine RungeKutta(y_0,a,b,h)
    implicit none
    
    real(8), intent(in):: a, b      ! Initial and final times
    real(8), intent(in):: h         ! Timestep
    real(8), intent(in):: y_0       ! Initial condition
    real(8) :: y,t                  ! Solution (wavefunction) from previous timestep, time
    real(8) :: k1, k2, k3, k4 
    
    t = a
    y = y_0
    
    do while(t <= b)
       k1 = h*f(y,t)
       k2 = h*f(y+(1d0/2d0)*k1, t+h*(1d0/2d0))
       k3 = h*f(y+(1d0/2d0)*k2, t+h*(1d0/2d0))
       k4 = h*f(y+k3,t+h)
       
       y = y + (k1 + 2d0*k2 + 2d0*k3 + k4) / 6d0
       
       write(20,*)t, y
       t=t+h
    end do
  end subroutine RungeKutta

!!=================================================================================  
  
  function f(y,t)           !! Contains the differential equation
    implicit none
    
    real(8)::f,y,t,pi

    pi=acos(-1d0)
    f=-y+sin(pi*2d0*t)       !! The differential equation
  end function f
  

end program wavepacket
