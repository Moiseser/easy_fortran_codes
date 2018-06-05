!=============================================================================!
!                               Newton Raphson Method                        !
!==============================================================================!
!    Discussion:
!
!
!    Modified:
!       30 May 2018
!    Author:
!       Moises Romero
!==============================================================================!
!Equation used in this example is x^3 = 1
module equation
  implicit none
  double precision :: y,y_prime
contains
  subroutine fnx(x)
    implicit none
    double precision :: x
    y = x**3 - 1
  end subroutine fnx
  subroutine deriv(x)
    implicit none
    double precision :: x
    y_prime = 3*x**2
  end subroutine
end module equation
program Newton_Raphson
  use equation
  implicit none
  double precision :: a1,a2, a,a0 , delta,fna,fna_prime, fna1 , fna2, cut_off,diff , error
  integer :: i
  logical :: H
  read (*,*) a
  cut_off = .0000000001
  i = 0d0
  H = .True.
  open(40,file = 'plot.dat')
  close(40)
  do while (H)
    i = i +1
call fnx(a)
fna = y
call deriv(a)
fna_prime = y_prime
a0 = a
error = ABS(1 -a0 )
open(40,position = 'append' , file = 'plot.dat')
write (40,*) i , error
close(40)
a = a - (fna/fna_prime)
diff = ABS(a-a0)
if (diff < cut_off) then
  write(*,*) 'Root is ' , a
  write(*,*) 'This took ' , i , 'iterations'
  H = .False.
end if
end do
end program Newton_Raphson
