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
do l = 1 , N - 2
do i = 1 , N - m - 1
   j = i + 1 + m
!    write (*,*) "values will be multiplied" , x(i) , x(j)
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
do i = 1 , N
  sum_squares = sum_squares + x(i)*x(i)
end do
final_tc(1) = sum_squares/N
do i = 2 , N -1
  final_tc(i) = t_correlation(i-1)
end do
final_tc(N) = (sum(x) /N)**2
!write(*,*) final_tc
open (14, file = 'time_correlation_function.dat')
do i = 1 , N
  write(14,*) final_tc(i)
end do
close(14)
end program
