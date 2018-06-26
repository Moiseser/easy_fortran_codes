module functions
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



end module functions

program simulated_annealing
use functions
implicit none
integer :: i,j,m,n,k,l , NMC_steps , count
double precision :: energy,old_energy,new_energy, energy_min, deltaE
double precision :: acceptance_ratio , rando ,lambda_initial, lambda, rando2
double precision :: T, Tmin , Tmax , Tstep, box
double precision,allocatable :: trial(:,:) , min_config(:,:)
character*30 filename

!Reading in XYZ file and Coordinates
 open(13, file = 'LJ13_random.xyz')
  read(13,*) Natoms
allocate (Q(3,Natoms),trial(3,Natoms),min_config(3,Natoms))
  read(13,*)
do i = 1 , Natoms
  read(13,*) Q(:,i)
  end do

!Contraining box
box = 2.0
!Choose amount of Monte Carlo steps
NMC_steps = 10000
lambda_initial = .01
energy_min = 10E8
acceptance_ratio = .45
!Selecting Temperature Range Loop
Tmin = .00001
T = 1
Tstep = -.01
do while (T > (Tmin+Tstep))
  T = T+Tstep
!Setting initial parameters
count = 0d0
lambda = lambda_initial
!Beginning Monte Carlo Steps
  do n = 1 , NMC_steps
    old_energy = 0d0
    call LennardJones_Energy(Q,old_energy)
    trial(:,:) = Q(:,:)
    do i = 1 ,3
      do j = 1 , Natoms
        call random_move(rando)
        Trial (i,j) = Trial (i,j) + lambda * rando
        do while (ABS(Trial(i,j)) > box)
          call random_move(rando)
          Trial (i,j) = Trial (i,j) + lambda * rando
        end do
      end do
    end do
    !  write (*,'(3(F22.16))') Trial
      new_energy= 0d0
    call LennardJones_Energy(Trial,new_energy)
  !  write (*,*) 'old_energy is :' , old_energy
  !  write (*,*) 'new_energy is :' , new_energy
    deltaE = new_energy - old_energy
  !  write(*,*) 'deltaE is :' ,deltaE
    call random_number(rando2)
    if (Exp(-deltaE/T) > rando2) then
      Q(:,:) = Trial (:,:)
      count = count + 1
    end if
    if (mod(n,1000) == 0) then
        if (count/1000 > acceptance_ratio) then
          lambda = lambda * .99
        else
          lambda = lambda * 1.01
            count = 0d0
        end if
      !  write(*,*) 'lambda is now' , lambda
    end if
    if (new_energy < energy_min) then
      min_config(:,:) = Q(:,:)
      energy_min = new_energy
  !    write(*,*)'New lowest energy' , energy_min
  !    write(*,*) 'We found this at Temp:' , T
    !end if
    end if
end do
end do
write (*,*) 'Minimum energy found is:', energy_min
end program simulated_annealing
