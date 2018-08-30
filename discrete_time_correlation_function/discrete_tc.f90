!=============================================================================80
!                 Discrete Time Correlation Functions
!=============================================================================80
!    This program takes in a set of data and calculates the normalized
!    timecorrelation functions for a single set of data.
!==============================================================================!
!    Modified:
!       16 August 2018
!    Author:
!       Moises Romero
!==============================================================================!
!Future Updates : I want have this code read in a matrix of data and ouput a
! a matrix of timecorrelation functions.
! I also would like to give it a cut off that matches my trajectory legnth.
!==============================================================================!
program time_correlation_function
  implicit none
  integer :: i , j , k , l , m , N
  double precision, allocatable :: x(:) , t_correlation(:) , final_tc(:)
  double precision :: shifted_sum , avgt , total , sum_squares
  character(50) :: distribution
!==============================================================================!
!  Read in Initial Parameters
!==============================================================================!
read(*,*) N
read(*,*) distribution
allocate (x(N),final_tc(N),t_correlation(N-2) )
!==============================================================================!
! Reading in Distributions
!==============================================================================!
open (13,file = distribution)
read (13,*) x(:)
close(13)
m = 0d0
sum_squares = 0d0
!==============================================================================!
! This part of the code calculates all the shifted terms for tc
! so everything but the first and last term (variance terms)
!==============================================================================!
do l = 1 , N - 2
do i = 1 , N - m - 1
   j = i + 1 + m
    shifted_sum = shifted_sum + x(i)*x(j)
  !  write(*,*)  shifted_sum
!  end do
end do
total = N - m - 1
  t_correlation(l) = (shifted_sum/(total))
!  write (*,*) t_correlation(l)
  m = m + 1
  shifted_sum = 0d0
end do
!==============================================================================!
! This part of the code calculates the first term of tc <x^2> and
! last term of the tc <x>^2 and appends them to a new array
!==============================================================================!
do i = 1 , N
  sum_squares = sum_squares + x(i)*x(i)
end do
final_tc(1) = sum_squares/N
do i = 2 , N -1
  final_tc(i) = t_correlation(i-1)
end do
final_tc(N) = (sum(x) /N)**2
!==============================================================================!
! Normalize TC and print to a file. 
!==============================================================================!
do i = 1 , N
  final_tc(i) = final_tc(i) / final_tc(1)
end do
!write(*,*) final_tc
open (14, file = 'time_correlation_function.dat')
do i = 1 , N
  write(14,*) final_tc(i)
end do
close(14)
end program
