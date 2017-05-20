! i/o routine


subroutine sub_read_dimer(wf)
    ! read in input filename list
    ! hard code: dimer.fchk, dimer.log
    use wavefunction
    implicit none
    type (wfinfo) :: wf
    character(len = 15) :: start, tmpstr
    character(len = 35) :: headline
    integer :: i, j, k, m
    integer :: n_ele, n_basis, n_block
    integer :: dimerfchk = 101, dimerlog = 102
    
    write(*, *) "READ FCHK FILE OF DIMER."
51  format(a15)
    open(dimerfchk, file = 'dimer.fchk')
    do
        read(dimerfchk, 51) start
        select case(start)
            case('Number of alpha')
                backspace(dimerfchk)
                read(dimerfchk,*) (tmpstr, i=1,5), n_ele
            case('Number of basis')
                backspace(dimerfchk)
                read(dimerfchk,*) (tmpstr, i=1,5), n_basis
                exit
        end select
    end do

    write(*,*) n_ele, n_basis
    wf%dimer%n_ele = n_ele
    wf%dimer%n_basis= n_basis
    allocate(wf%dimer%ao_overlap(n_basis,n_basis))
    allocate(wf%dimer%mo_coefficient(n_basis,n_basis))
    allocate(wf%dimer%mo_energy(n_basis))

    do
       read(dimerfchk,51) start
       select case(start)
          case('Alpha Orbital E')
             read(dimerfchk,*) wf%dimer%mo_energy
          case('Alpha MO coeffi')
             read(dimerfchk,*) wf%dimer%mo_coefficient
             exit
          end select
       end do

    close(dimerfchk)

    write(*,*) "READ LOG FILE OF DIMER"
52  format(1x,a35)
53  format(1x,a7,5e14.6)

    n_block = (n_basis + 4) / 5

    open(dimerlog, file = 'dimer.log')
    do
        read(dimerlog,52) headline
        headline = adjustl(headline)
        headline = trim(headline)
        if(headline .eq. '*** Overlap ***') then
            do i = 1, n_block
                read(dimerlog, *)
                do j = (i-1)*5 + 1, n_basis
                    m = j - (i-1)*5
                    if (m > 5) m = 5
                    read(dimerlog, 53) tmpstr, (wf%dimer%ao_overlap(j,k),k = (i-1)*5+1,(i-1)*5+m)
                end do
            end do
            exit
        end if
    end do
    close(dimerlog)
    ! fill the symmetric matrix
    do i = 1, n_basis-1
        do j = i+1, n_basis
            wf%dimer%ao_overlap(i,j) = wf%dimer%ao_overlap(j,i)
        end do
    end do

end subroutine



subroutine sub_read_frag(frag, file_id)
    use wavefunction
    implicit none
    type (wfrep) :: frag
    integer :: file_id
    character(len=30) :: tmpstr, start, myappdix
    integer :: n_ele, n_basis
    integer :: i
    integer :: frag1file = 1001
    
    write(*,*) "READ FCHK FILE OF ONE FRAG"
51  format(a15)
    ! filename list frag*.fchk
    write(myappdix, '(I10)') file_id
    open(frag1file, file = 'frag'//trim(adjustl(myappdix))//'.fchk')
    do
       read(frag1file, 51) start
       select case(start)
          case('Number of alpha')
             backspace(frag1file)
             read(frag1file,*) (tmpstr,i=1,5), n_ele
          case('Number of basis')
             backspace(frag1file)
             read(frag1file,*) (tmpstr,i=1,5), n_basis
             exit
          end select
       end do

    write(*,*) n_ele, n_basis

    frag%n_ele = n_ele
    frag%n_basis = n_basis
    allocate(frag%mo_coefficient(n_basis, n_basis))
    allocate(frag%mo_energy(n_basis))
    
    do
       read(frag1file,51) start
       select case(start)
          case('Alpha Orbital E')
             read(frag1file,*) frag%mo_energy
          case('Alpha MO coeffi')
             read(frag1file,*) frag%mo_coefficient
             exit
          end select
       end do

     end subroutine


subroutine sub_read(wf)
    use wavefunction
    implicit none
    type(wfinfo) :: wf
    
    call sub_read_dimer(wf)
    ! frag1, frag2 ..
    call sub_read_frag(wf%frag1, 1)
    call sub_read_frag(wf%frag2, 2)

end subroutine


!program main
!use wavefunction
!implicit none
!type (wfinfo) :: wf
!call sub_read(wf)
!end



