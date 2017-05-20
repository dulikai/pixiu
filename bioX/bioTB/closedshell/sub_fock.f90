! fortran 95

! get fock matrix info.
! the fock matrix can be routinely dumped by any QC package.
! HOWEVER, this can only possible when the debug option is applied.
! For a flexible implementation, SAO MO E to obtain FOCK matrix
!
! Lowdin orthogonalization procedure
! Lowdin, P. O. J. Chem. Phys. 1950, 18, 365
! JACS, 2006, 128, 9882-9886
! J. Chem. Phys. 2003, 119, 9809

! calculate the transfer integral from QC results
! fock matrix is got by the formula F=SCE(C^T)S

! http://blog.163.com/jey_df/blog/static/182550161201302345930930/ 

! build fock matrix under AO basis
subroutine sub_fock(wf)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  integer :: i, j, n_basis
  real(kind=8), allocatable :: ao(:,:), moc(:,:)
  real(kind=8), allocatable :: moe(:,:), vtemp(:,:)

  n_basis = wf%dimer%n_basis
  moc = wf%dimer%mo_coefficient
  ao  = wf%dimer%ao_overlap

  allocate(moe(n_basis, n_basis))
  allocate(vtemp(n_basis, n_basis))
  allocate(wf%dimer%fock(n_basis, n_basis))

  do i = 1, n_basis
    moe(i,i) = wf%dimer%mo_energy(i)
  end do
  
  ! the procedure to build fock matrix
  ! calculate S_AO C_MO E C_MO^T S_AO^T
  vtemp = matmul(ao, moc)
  vtemp = matmul(vtemp, moe)
  vtemp = matmul(vtemp, transpose(moc))
  wf%dimer%fock = matmul(vtemp, transpose(ao))

  deallocate(moe, vtemp)

end subroutine



