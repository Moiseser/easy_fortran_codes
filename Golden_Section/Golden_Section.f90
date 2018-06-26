!=============================================================================!
!                               Golden Section                        !
!==============================================================================!
!    Discussion:
!    Golden Section method is used to find a local minima between an interval
!    by taking a difference of the interval and multiplying by a the 'golden
!    ratio'. This value is then added/subtracted to the selected intervals.
!    The function is then reevaluated and the boundries are readjusted. This
!    is repeated until a local minima is established. Local minima is selected
!    by taking an averages and comparing to a small#. When fluctuations are small
!    the local minima within your selected intervals is selected.
!
!    Modified:
!       24 May 2018
!    Author:
!       Moises Romero
!==============================================================================!
!Equation used in this example is x^2
module equation
  implicit none
  double precision :: y
contains
  subroutine fnx(x)
    implicit none
    double precision :: x
    y = x**2
  end subroutine fnx
end module equation
program Golden_Section
use equation
implicit none
double precision :: a,b,d,phi,x1,x2,avg,avg1,avg2,cut_off,fx1,fx2
  integer :: i
  logical :: H
  read(*,*) a
  read(*,*) b
!cut_off is a small # to compare when avg fluctuations are small
  cut_off = .0001
  H = .True.
!'Golden' ratio
  phi = (SQRT(5.0) -1.0) / 2.0
  i = 0
  avg = 0d0
  avg1 = 0d0
  avg2 = 0d0
!An infinite loop is needed until minimum
do while (H)
  i = i +1
  !Points on the interval are found
  d = phi * (b-a)
  x1 = a + d
  x2 = b - d
!Function is evaluated
 call fnx(x1)
 fx1 = y
 call fnx(x2)
 fx2 = y
!Interval is shifted
if (fx1 < fx2) then
  a = x2
else if (fx2<fx1) then
  b = x1
end if
!fluctuations are checked 
avg1 = avg
avg = (a+b)/2
avg2 = (avg - avg1)
avg2 = ABS(avg2)
if (avg2<cut_off) then
  H = .FALSE.
end if
end do
  write (*,*) 'This took', i , 'iterations'
  write (*,*) 'Minimum is' , avg
end program Golden_Section
