!=============================================================================80
!                 Kolomogorov - Smirnov Test
!=============================================================================80
!    Discussion: This code is adapated from Numerical Reciepes Second Edition.
!    The k-test algorithim compares a sample of data with a refrence probability
!    distribution(or a second sample of data but this code will use a ref prob.)
!    It does this by quantifying a distance(d) between the two probability
!    distributions. d is then used to calculate the significane level probability.
!    Small values of this prob. tells us that the two distributions are
!    significantly  different.
!==============================================================================!
!    Modified:
!       23 July 2018
!    Author:
!       Moises Romero
!==============================================================================!
module functions_n_stuff
  implicit none
  integer :: i , j , N
  double precision :: ks_prob , k , beta
  double precision , allocatable :: V(:),boltz_energy(:)
contains
  subroutine SHO_potential(position)
    implicit none
    double precision  :: position(N)
    do i = 1 , N
      V(i) = .5 * k * position(i) **2
    end do
  end subroutine

  subroutine Boltzman_Dist(energy)
    implicit none
    double precision :: energy(N)
    call SHO_potential(energy)
    do i = 1 , N
      boltz_energy(i) = exp(-beta * V(i))
    end do
  end subroutine



   subroutine ksone(energies)
     implicit none
   double precision ::  positions(N) ,energies(N),  D
   double precision ::  dt , en , ff , fn , fo
   call Boltzman_Dist(energies)
   en = N
   D = 0d0
   fo = 0d0
   do i = 1  , N
     fn = i / en
     ff = boltz_energy(i)
     dt = max(abs(fo-ff),abs(fn-ff))
     if (dt.gt.d) then
       d = dt
       write(*,*) 'New greater distance!:' , d
     end if
        fo = fn
     end do
     en = sqrt(en)
     call ks_probability((en+0.12 + 0.11/en)*d)
     write(*,*) 'K-Test Value is:' , ks_prob
   end subroutine


  subroutine ks_probability(lambda)
    implicit none
    double precision ::lambda,a2,fac,term,termbf
    double precision , parameter:: ESP1 = .001 , ESP2 = 1.0d-8
    logical :: H
    a2 = -2d0 * lambda**2d0
    fac = 2
    ks_prob = 0d0
    termbf = 0d0
    H = .FALSE.
    do i = 1 , 100
      if (H .eqv. .TRUE.) then
        write(*,*)'found it '
      else
      term = fac * exp(a2*(j**2))
      ks_prob = ks_prob + term
      write(*,*) ks_prob
      if (abs(term) .le. ESP1*termbf .or. abs(term) .le. ESP2*ks_prob) then
        H = .TRUE.
        write(*,*) 'found it : ' , ks_prob
      end if
        fac = -fac
        termbf = abs(term)
      end if
      end do
      ks_prob = 1.0
    end subroutine ks_probability
  end module functions_n_stuff
!==============================================================================!
program ktest
use functions_n_stuff
implicit none
double precision, allocatable :: x(:)
double precision :: test , y
character(50) :: traj_file
!==============================================================================!
!Read in Initial Parameters
!==============================================================================!
read (*,*) N
read (*,*) k
read (*,*) beta
Allocate (x(N),V(N),boltz_energy(N))
read (*,*) traj_file
!==============================================================================!
! Trajectories
!==============================================================================!
open (13,file = traj_file)
read (13,*) x(:)
close(13)

call ksone(x(:))





end program ktest
