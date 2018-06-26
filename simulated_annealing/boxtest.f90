module bits_and_pieces
implicit none
integer :: Natoms
double precision , allocatable :: Q(:,:)
double precision :: LJ
!Calc Energy using Lennard Jones Potential
contains

subroutine LennardJones(Q1,Q2)
  implicit none
  double precision , parameter :: sigma = 1d0 , epsilon = 1d0
  double precision :: Q1(3) , Q2(3)
  LJ = sigma / sum((Q2(:)-Q1(:))**2)**3
!write (*,*) 'difference' , Q2(:) - Q1(:)
  LJ = 4 * epsilon * (LJ**2 - LJ)
end subroutine

subroutine LennardJones_Energy(coord,ene)
implicit none
integer :: i,j
double precision :: coord(3,Natoms), ene
do i = 1 , Natoms
  do j = 1 , i - 1
    call LennardJones(coord(:,i),coord(:,j))
    ene = ene + LJ
  end do
end do
end subroutine
subroutine random_move(some_number)
  implicit none
  double precision some_number , C , v
  double precision , parameter :: A = -1d0 , B = 1d0
   Call random_number(C)
  some_number = A + (B-A) * C
end subroutine



end module bits_and_pieces

program boxtest
use bits_and_pieces
implicit none
integer :: i,j,m,n,k,l , NMC_steps , count
double precision :: energy,old_energy,new_energy, energy_min, deltaE
double precision :: acceptance_ratio , rando ,lambda_initial, lambda, rando2
double precision :: T, Tmin , Tmax , Tstep, box
double precision,allocatable :: trial(:,:) , min_config(:,:)
character*30 filename

!Reading in XYZ file and Coordinates
 open(13, file = 'new_jones.xyz')
  read(13,*) Natoms
allocate (Q(3,Natoms),trial(3,Natoms),min_config(3,Natoms))
  read(13,*)
do i = 1 , Natoms
  read(13,*) Q(:,i)
  end do

box = 2
write(*,*) box
lambda = .01
trial(:,:) = Q(:,:)
do i = 1 ,3
  do j = 1 , Natoms
    call random_move(rando)
    Trial (i,j) = Trial (i,j) + lambda * rando
    do while (ABS(Trial(i,j)) > box)
      write(*,*) 'Get back in the box!:' , trial(i,j)
      call random_move(rando)
      Trial (i,j) = Trial (i,j) + lambda * rando
    end do
    write(*,*) 'We in the box:' , trial(i,j)
  end do
end do
end program boxtest
