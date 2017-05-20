
module wavefunction
! molecular orbital and related wavefunction representation
type wfrep
    character(len = 30) :: name
    ! num-of-basis for a model
    integer :: n_basis
    ! atomic basis overlap integral matrix
    real(kind = 8), allocatable :: ao_overlap(:,:)

    ! for closed shell only >
    integer :: n_ele
    real(kind = 8), allocatable :: mo_coefficient(:,:)
    real(kind = 8), allocatable :: fock(:,:)
    real(kind = 8), allocatable :: mo_energy(:)

    ! for open shell >
    integer :: n_ele_alpha, n_ele_beta
    real(kind = 8), allocatable :: mo_coefficient_alpha(:,:), mo_coefficient_beta(:,:)
    real(kind = 8), allocatable :: fock_alpha(:,:), fock_beta(:,:)
    real(kind = 8), allocatable :: mo_energy_alpha(:), mo_energy_beta(:,:)
    
end type


! wavefunction information for the tight binding model
type wfinfo
    character(len = 30) :: name
    type(wfrep) :: dimer
    type(wfrep) :: frag1
    type(wfrep) :: frag2
    ! parameters & results
    integer :: i_mo, j_mo
    real(kind = 8) :: overlap
    real(kind = 8) :: parm(2,2), Heff(2,2)
end type 

end module

