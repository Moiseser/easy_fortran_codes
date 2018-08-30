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
double precision, allocatable :: mc_data(:,:), partial_avg(:,:) , variance(:,:)
double precision ,allocatable:: std_deviation(:,:), avg(:,:) ,  old_avg(:,:)
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
allocate(mc_data(N,M))
allocate(std_deviation(1,M),old_avg(1,M),avg(1,M))
write (*,*) 'Sucessfully Read in paramters from input file'
!==============================================================================!
!Reading in Data file
!==============================================================================!

open (13,file = datafile)
do i = 1 , N
read (13,*) mc_data(i,:)
write(*,*) mc_data(i,:)
end do
close(13)
!write(*,*) mc_data
!mc_data(:,1) = (/1.0,-3.0,1.0/)
!mc_data(:,2) = (/2.0,-8.0,8.0/)
!mc_data(:,3)= (/-6.0,3.0,-15.0/)
!mc_data(:,1) = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/)
!mc_data(:,2) = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/)
!mc_data(:,3) = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/)
!write (*,*) mc_data(5,:)
!==============================================================================!
!Setting Varaibles for Calc Partial Averages
!==============================================================================!
adjust = mod((N-skip),freq)
!write(*,*) 'Adjust is ' , adjust
partial_avg_size = (N - skip - adjust ) / freq
!write(*,*) 'partial_avg_size is' , partial_avg_size
allocate (partial_avg(partial_avg_size,M),variance(partial_avg_size,M))
avg(:,:) = 0d0
count = 0d0
partial_avg = 0D0
x = 0d0
do i = (1 + skip + adjust) , N
  write (*,*) 'avg is ' , avg
  write(*,*) 'mc_data is', mc_data(i,:)
  avg(1,:) = avg(1,:) + mc_data(i,:)
  write (*,*) avg
  write(*,*) 'some numbers' , avg
  count = count + 1
  write (*,*) count
  if (mod(count , freq) == 0 ) then
    avg(1,:) = avg(1,:) / freq
    write (*,*) 'avg/freq is' , avg
    x = x + 1
   partial_avg(x,:) =  avg(1,:)
    avg = 0d0
    count = 0d0
  end if
end do
write(*,*) 'Avg Array is:'
do i = 1 , partial_avg_size
  write (*,*) partial_avg(i,:)
end do

old_avg = 0d0
do i = (1 + skip + adjust) , N
  old_avg(1,:) = old_avg(1,:) + mc_data(i,:)
end do
old_avg(1,:) = old_avg(1,:) / (N - (skip+adjust))
write(*,*) 'old_avg is :' , old_avg
do i = 1 , partial_avg_size
  variance(i,:) = (partial_avg(i,:) - old_avg(1,:))**2
  write(*,*) variance(i,:)
end do
std_deviation(1,:) = sqrt(sum(variance)/(partial_avg_size-1))
write(*,*) 'Standard Deviation' , std_deviation

end program
