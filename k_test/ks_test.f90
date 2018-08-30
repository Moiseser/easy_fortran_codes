!=============================================================================80
!                 Kolomogorov - Smirnov Test
!=============================================================================80
!    Discussion: This code is adapated from Numerical Reciepes Second Edition.
!    It follows the outline presented in the book with some changes better fit
!    for a more contemporary f90 program.
!    The k-test algorithim compares two distributions of data
!    (or a refrence probability but this code uses the two distb method)
!    It takes two SORTED sets of data and creates Cumaltive Distribution Functions(CDF)
!    NOTES: The sets of data can be different lengths
!    From the CDF's you calculate the MAX distance(d).  d is the distance
!    between the two CDF's and will tell us if they are from different distributions.
!    d will vry from 0 to 1.
!    The closer to 1 we can confidentaly say the two distributions are from different PDFs
!    The closer to 0 we can not say they are different.
!==============================================================================!
!    Modified:
!       08 August 2018
!    Author:
!       Moises Romero
!==============================================================================!
module functions_globalvariables
  implicit none
  integer :: i , j , N
  double precision :: ks_prob
contains
!==============================================================================!
!The following Function is a work in progress and is not needed for the codes
! current function but might need it later
!==============================================================================!
!   function ks_probability(lambda)
!   implicit none
!   double precision :: lambda,a2,fac,term,termbf, ks_probability
!   double precision , parameter:: ESP1 = .001 , ESP2 = 1.0d-8
!   !logical :: H
!   a2 = -2d0 * lambda**2d0
!   fac = 2
!   ks_probability = 0d0
!   termbf = 0d0
! !  H = .FALSE.
!   do i = 1 , 100
! !    if (H .eqv. .TRUE.) then
! !      write(*,*)'found it '
! !    else
!     term = fac * exp(a2*(j**2))
!     ks_probability = ks_probability + term
!     write(*,*) ks_probability
!     if (abs(term) .le. (ESP1*termbf) .or. abs(term) .le. (ESP2*ks_prob)) then
!     write(*,*) 'found it : ' , ks_probability
! return
! end if
!   !    H = .TRUE.
!
! !    end if
!       fac = -fac
!       termbf = abs(term)
! !    end if
!     end do
! !    ks_probability = 1
!   end
end module functions_globalvariables
program KS_test
  use functions_globalvariables
  implicit none
integer :: N1 , N2 , j1 , j2
double precision :: d , d1 , d2, dt ,en, en1 , en2 , fn1, fn2
double precision , allocatable :: data1(:) , data2(:)
character(50) :: dist_langevindynamics , dist_randomvariate_pdf

!==============================================================================!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!IMPORTANT DATA SHOULD BE IN ASCENDING ORDER
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!Ill write in a flag later to check using a random integer
!==============================================================================!
!  Read in Initial Parameters
!==============================================================================!
read(*,*) N1
read(*,*) N2
read(*,*) dist_langevindynamics
read(*,*) dist_randomvariate_pdf
allocate (data1(N1) , data2(N2) )
!==============================================================================!
! Reading in Distributions
!==============================================================================!
open (13,file = dist_langevindynamics)
read (13,*) data1
close(13)
open (13,file = dist_randomvariate_pdf )
read (13,*) data2
close(13)
!==============================================================================!
! Initialize some Parameters
!==============================================================================!
en1 = N1
en2 = N2
j1 = 1d0
j2 = 1d0
fn1 = 0d0
fn2 = 0d0
!==============================================================================!
! Generates CDF and calculates "d" through the CDF only saving the maximum value of d 
!==============================================================================!
1 if (j1 .le. n1 .and. j2 .le. n2) then
  d1 = data1(j1)
  d2 = data2(j2)
  if (d1 .le. d2 ) then
    fn1 = j1 / en1
  end if
  if (d2 .le. d1) then
    fn2 = j2 / en2
    j2 = j2 + 1
  end if
  dt = abs(fn2-fn1)
  if (dt .gt. d) then
    d = dt
    goto 1
  end if
end if
  write (*,*) " Distance D is " , d
!==============================================================================!
!Pieces for function above
!==============================================================================!
! en = sqrt(en1*en2/(en1+en2))
! !call ks_probability((en+0.12d0 + 0.11d0 /en)*d)
! ks_prob = ks_probability((en+0.12d0 + 0.11d0 /en)*d)
! write(*,*) 'K-Test Value is:' , ks_prob


end program KS_test
