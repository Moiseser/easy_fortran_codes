!=============================================================================80
!                 Standard Deviation of Monte_Carlo
!=============================================================================80
!    Discussion: Python 2 Implementation of Gaussian Elimination.
!       The code takes a square matrix and transforms it into a diagonal matrix
!       using basic row operations.
!       The code requires your first matrix to be square, and non-singular.
!       The second matrix must have the same number of rows as matrix 1, but not
!       necessairly square (any number of columns).
!       If the second matrix is a column vector, than the first matrix will be
!       interpreted as a system of equations, Ax = b, and it will solve for the
!       x vector.
!       Finally the inverse of the first matrix is computed.
!==============================================================================!
!       We assume that the matrix is non-singular and has no accidental 0s
!       If these occur in your matrix all bets are off!!!
!==============================================================================!
!    Modified:
!       15 July 2018
!    Author:
!       Moises Romero
!==============================================================================!
!==============================================================================!
!This code will take a list indepdent data and output Standard Deviations
program Standard_Deviation_MC
implicit none
double precision, allocatable :: mc_data(:), partial_avg(:) , variance(:)
double precision :: std_deviation, avg , old_avg
integer :: i, j , k , l ,N, M, partial_avg_size , skip , adjust , freq , count,x
character(50) :: datafile
!==============================================================================!
!Read in Initial Parameters
!==============================================================================!
read (*,*) skip
read (*,*) freq
read (*,*) N
read (*,*) M
read(*,*) datafile
allocate(mc_data(N))
write (*,*) 'Sucessfully Read in paramters from input file'
!==============================================================================!
!Reading in Data file
!==============================================================================!
!This will be turned into a matrix loop eventually for matrix sets of data
!Problem for another day
open (13,file = datafile)
read (13,*) mc_data(:)
close(13)
write(*,*) mc_data
!==============================================================================!
!Setting Varaibles for Calc Partial Averages
!==============================================================================!
! adjust = mod((N-skip),freq)
! partial_avg_size = (N - skip - adjust ) / freq
! allocate (partial_avg(partial_avg_size),variance(partial_avg_size))
! avg = 0d0
! count = 0d0
! partial_avg = 0D0
! x = 0d0
! do i = (1 + skip + adjust) , N
!   avg = avg + mc_data(i)
!   write (*,*) avg
!   count = count + 1
!   if (mod(count , freq) == 0 ) then
!     avg = avg / freq
!     x = x + 1
!    partial_avg(x) =  avg
!     avg = 0d0
!     count = 0d0
!   end if
! end do
! write(*,*) 'Avg Array is' , partial_avg
x = N - skip
allocate(variance(x))
avg = 0d0
do i = (1 + skip ) , N
  avg = avg + mc_data(i)
end do
avg = avg / (N - (skip))
do i = 1 , x 
  variance(i) = (mc_data(i+skip) - avg)**2
  write(*,*) variance
end do
std_deviation = sqrt(sum(variance)/(N-(1+skip)))
write(*,*) 'Standard Deviation' , std_deviation

end program
