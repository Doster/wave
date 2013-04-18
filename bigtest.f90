program test2

  implicit none

  real(8)::a,b,h,y_0,t,y, pi

  write(*,*)"Enter the interval a,b, the value of the step-size h and the value of y_0"
  read(*,*)a,b,h,y_0

  call opentextfiles()


  pi = acos(-1d0)
  t = a

  do while(t <= b)
     y = (1d0+(2d0*pi)/(1d0+4d0*pi**2d0))*exp(-t)+(sin(2d0*pi*t)-2d0*pi*cos(2d0*pi*t))/(1d0+4d0*pi**2d0)
     write(15,*)t, y
     t = t+h
  end do

  call euler(y_0,a,b,h)

  call RungeKutta(y_0,a,b,h)


  call closetextfiles()



contains

 subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="Exact_values.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, Euler.dat file not opened properly------------"
    endif
    OPEN(UNIT=16,FILE="Euler.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, Euler.dat file not opened properly------------"
    endif
    OPEN(UNIT=17,FILE="Euler_error.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, Euler_error.dat file not opened properly------------"
    endif
    OPEN(UNIT=18,FILE="RungeKutta.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, RungeKutta.dat file not opened properly------------"
    endif
    OPEN(UNIT=19,FILE="RungeKutta_error.dat",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, RungeKutta_error.dat file not opened properly------------"
    endif
  endsubroutine opentextfiles

  subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=16)
    CLOSE(UNIT=17)
    CLOSE(UNIT=18)
    CLOSE(UNIT=19)
  endsubroutine closetextfiles


  subroutine euler(y_0,a,b,h)
    implicit none
    
    real(8), intent(inout)::a,h,b, y_0
    real(8):: y,t, valor_exacto, pi, error_euler
    
    pi = acos(-1d0)  !! Maximum precision for PI possible on any architecture
    t = a
    y = y_0
    
    do while(t<=b)
       y = y + h*f(y,t)
       valor_exacto = (1d0+(2d0*pi)/(1d0+4d0*pi**2d0))*exp(-t)+(sin(2d0*pi*t)-2d0*pi*cos(2d0*pi*t))/(1d0+4d0*pi**2d0)
       error_euler = abs(y-valor_exacto)
       
       write(16,*)t, y, valor_exacto
       write(17,*)error_euler
       t = t+h
    end do
  end subroutine euler
  

  subroutine RungeKutta(y_0,a,b,h)
    implicit none
  
    real(8), intent(inout)::a,h,b, y_0
    real(8):: y,t, F_1, F_2, F_3, F_4, valor_perfecto, pi, error_kutta

    pi = acos(-1d0)
    t = a
    y = y_0

    do while(t<=b)
       F_1=h*f(y,t)
       F_2=h*f(y+(1d0/2d0)*F_1, t+h*(1d0/2d0))
       F_3=h*f(y+(1d0/2d0)*F_2, t+h*(1d0/2d0))
       F_4=h*f(y+F_3,t+h)

       y = y+(F_1+2d0*F_2+2d0*F_3+F_4)/6d0

       valor_perfecto = (1d0+(2d0*pi)/(1d0+4d0*pi**2d0))*exp(-t)+(sin(2d0*pi*t)-2d0*pi*cos(2d0*pi*t))/(1d0+4d0*pi**2d0)
       error_kutta = abs(y-valor_perfecto)

       write(18,*)t, y, valor_perfecto
       write(19,*)error_kutta
       t=t+h
    end do
  end subroutine RungeKutta


  function f(y,t)
    implicit none
    real(8)::f,y,t,pi
    pi=acos(-1d0)
    f=-y+sin(pi*2d0*t) 
  end function f

end program test2
