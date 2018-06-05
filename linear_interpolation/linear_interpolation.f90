!=============================================================================!
!                               Linear Interpolation                        !
!==============================================================================!
!    Discussion:
!
!
!    Modified:
!       24 May 2018
!    Author:
!       Moises Romero
!==============================================================================!
!Equation used in this example is x^3 = 1
module equation
  implicit none
  double precision :: y
contains
  subroutine fnx(x)
    implicit none
    double precision :: x
    y = x**3 - 1
  end subroutine fnx
end module equation
program Linear_Interpolation
  use equation
  implicit none
  double precision :: a,b,fna,fnb,c,fnc,fna_fnb_prod,fluct,root, fna_fnc_prod , cut_off, error
  integer :: i
  logical :: H
  read(*,*) a
  read(*,*) b
  !cut_off is a small # to compare when avg fluctuations are small
    cut_off = .0000000001
    i = 0d0
    H = .True.
    open(40,file = 'plot.dat')
    close(40)
call fnx(a)
fna = y
call fnx(b)
fnb = y
fna_fnb_prod = fna*fnb
If (fna == 0) then
  write (*,*) 'Root is' , a
else if (fnb == 0) then
  write (*,*) 'Root is ' , b
else if (fna_fnb_prod > 0) then
  write (*,*) 'No Root in Interval'
else
  !An infinite loop is needed until minimum
do while (H)
  i = i +1
  call fnx(a)
  fna = y
  call fnx(b)
  fnb = y

c = (((((fnb-fna)/(b-a))*a)-fna)*(b-a))/(fnb-fna)
error = ABS(1-c)
open(40,position = 'append', file = 'plot.dat')
write(40,*) i,error
close(40)
if (fna < 0) then
  a = c
else
  b = c
end if
if (ABS(fna) < cut_off) then
  root = a
  write (*,*) 'The root is ' , root
  H = .False.
else if (ABS(fnb) < cut_off) then
  root = b
  write (*,*) 'The root is' , root
  H = .False.
  end if
end do
end if

end program Linear_Interpolation
