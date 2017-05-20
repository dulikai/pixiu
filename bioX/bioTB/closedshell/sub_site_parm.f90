! fortran

! calculate the site energy and transfer integral with a specific pair of molecular orbitals..


function site_energy(wf, i_mo, frag_id)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  real(kind=8) :: site_energy
  integer :: i_mo, frag_id
  integer :: i, j
  integer :: n_bas, n_bas1, n_bas2
  real(kind=8),allocatable :: vtemp(:), vec1(:), vec2(:)

  n_bas = wf%dimer%n_basis
  n_bas1 = wf%frag1%n_basis
  n_bas2 = wf%frag2%n_basis

  allocate(vtemp(n_bas))
  allocate(vec1(n_bas))
  allocate(vec2(n_bas))

  vec1 = 0.0d0
  vec2 = 0.0d0
  if ( frag_id == 1 ) then
     vec1(1:n_bas1) = wf%frag1%mo_coefficient(:,i_mo)
     vec2(1:n_bas1) = wf%frag1%mo_coefficient(:,i_mo)
  else if (frag_id == 2) then
     vec1(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(:,i_mo)
     vec2(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(:,i_mo)
  else
     write(*,*) 'ERROR: ONLY FRAG_ID = 1/2 ARE POSSIBLE'
     STOP
  end if

  ! calculate site energy
  vtemp = 0.0d0
  do i = 1, n_bas
     do j = 1, n_bas
        vtemp(i) = vtemp(i) + vec1(j) * wf%dimer%fock(j,i)
     end do
  end do
  site_energy = 0.0d0
  do i = 1, n_bas
     site_energy = site_energy + vtemp(i) * vec2(i)
  end do

  deallocate(vtemp, vec1, vec2)

  return

end


function site_polarizable(wf, i_mo, j_mo, frag_id)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  real(kind=8) :: site_polarizable
  integer :: i_mo, j_mo, frag_id
  integer :: i, j
  integer :: n_bas, n_bas1, n_bas2
  real(kind=8),allocatable :: vtemp(:), vec1(:), vec2(:)

  n_bas = wf%dimer%n_basis
  n_bas1 = wf%frag1%n_basis
  n_bas2 = wf%frag2%n_basis

  allocate(vtemp(n_bas))
  allocate(vec1(n_bas))
  allocate(vec2(n_bas))

  vec1 = 0.0d0
  vec2 = 0.0d0
  if ( frag_id == 1 ) then
     vec1(1:n_bas1) = wf%frag1%mo_coefficient(:,i_mo)
     vec2(1:n_bas1) = wf%frag1%mo_coefficient(:,j_mo)
  else if (frag_id == 2) then
     vec1(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(:,i_mo)
     vec2(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(:,j_mo)
  else
     write(*,*) 'ERROR: ONLY FRAG_ID = 1/2 ARE POSSIBLE'
     STOP
  end if

  ! calculate transfer integral
  vtemp = 0.0d0
  do i = 1, n_bas
     do j = 1, n_bas
        vtemp(i) = vtemp(i) + vec1(j) * wf%dimer%fock(j,i)
     end do
  end do
  site_polarizable = 0.0d0
  do i = 1, n_bas
     site_polarizable = site_polarizable + vtemp(i) * vec2(i)
  end do

  deallocate(vtemp, vec1, vec2)

  return

end




function transfer_integral(wf, i_mo, j_mo)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  real(kind=8) :: transfer_integral 
  integer :: i_mo, j_mo
  integer :: i, j
  integer :: n_bas, n_bas1, n_bas2
  real(kind=8),allocatable :: vtemp(:), vec1(:), vec2(:)

  n_bas = wf%dimer%n_basis
  n_bas1 = wf%frag1%n_basis
  n_bas2 = wf%frag2%n_basis

  allocate(vtemp(n_bas))
  allocate(vec1(n_bas))
  allocate(vec2(n_bas))

  vec1 = 0.0d0
  vec2 = 0.0d0
  vec1(1:n_bas1) = wf%frag1%mo_coefficient(:,i_mo)
  vec2(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(:,j_mo)

  ! calculate transfer integral
  vtemp = 0.0d0
  do i = 1, n_bas
     do j = 1, n_bas
        vtemp(i) = vtemp(i) + vec1(j) * wf%dimer%fock(j,i)
     end do
  end do
  transfer_integral = 0.0d0
  do i = 1, n_bas
     transfer_integral = transfer_integral + vtemp(i) * vec2(i)
  end do

  deallocate(vtemp, vec1, vec2)
  
  return

end 


! overlap between two orbital between frag1 and frag2
function mo_overlap(wf, i_mo, j_mo)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  integer :: i_mo, j_mo, i, j
  real(kind=8) :: mo_overlap

  integer :: n_bas, n_bas1, n_bas2
  ! mo coefficients for each fragment
  real(kind=8), allocatable :: mo1(:), mo2(:), vtemp(:), ao_overlap(:,:)

  n_bas = wf%dimer%n_basis
  ao_overlap = wf%dimer%ao_overlap
  n_bas1 = wf%frag1%n_basis
  n_bas2 = wf%frag2%n_basis

  allocate(mo1(n_bas), mo2(n_bas), vtemp(n_bas))
  
  ! assign molecule orbital values
  mo1 = 0.0d0
  mo2 = 0.0d0
  mo1(1:n_bas1) = wf%frag1%mo_coefficient(1:n_bas1, i_mo)
  mo2(n_bas1+1:n_bas) = wf%frag2%mo_coefficient(1:n_bas2, j_mo)

  ! @ CALCULATE C_{MO1} S_{AO} C_{MO2} : <phi_1|V|phi_j>
  ! vector.matrix product
  vtemp = 0.0d0
  mo_overlap = 0.0d0
  do i = 1, n_bas
     do j = 1, n_bas
        vtemp(i) = vtemp(i) + mo1(j) * ao_overlap(j,i)
     end do
  end do
  do i = 1, n_bas
     mo_overlap = mo_overlap + vtemp(i) * mo2(i)
  end do

  deallocate(mo1, mo2, vtemp)

  return

end


subroutine sub_site_parm(wf, i_mo, j_mo)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  integer :: i_mo, j_mo
  real(kind = 8),external :: site_energy, transfer_integral, mo_overlap

  wf%parm(1,1) = site_energy(wf, i_mo, 1)
  wf%parm(2,2) = site_energy(wf, j_mo, 2)
  wf%parm(2,1) = transfer_integral(wf, i_mo, j_mo)
  wf%parm(1,2) = wf%parm(2,1)

  wf%overlap = mo_overlap(wf, i_mo, j_mo)

  wf%i_mo = i_mo
  wf%j_mo = j_mo

  return

end subroutine



