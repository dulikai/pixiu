! fortran 90

! only calculate HOMO LUMO contributions


subroutine sub_frontier(wf)
    use constant
    use wavefunction
    type(wfinfo) :: wf
    integer :: frag1_homo, frag1_lumo, frag2_homo, frag2_lumo
    real(kind = 8),external :: site_polarizable
    real(kind = 8) :: sp1, sp2
    integer :: file_id = 1001

    ! prepare input & data
    call sub_read(wf)
    call sub_fock(wf)

    open(file_id, file='frontier.dat')

    frag1_homo = wf%frag1%n_ele
    frag2_homo = wf%frag2%n_ele
    frag1_lumo = frag1_homo + 1
    frag2_lumo = frag2_homo + 1

    ! for HOMOs between fragment 1 and 2
    call sub_site_parm(wf, frag1_homo, frag2_homo)
    call sub_Lowdin2d(wf)
    call sub_write_frontier(file_id, wf)

    ! for LUMOs between fragment 1 and 2
    call sub_site_parm(wf, frag1_lumo, frag2_lumo)
    call sub_Lowdin2d(wf)
    call sub_write_frontier(file_id, wf)

    ! for HOMO-LUMO between fragment 1 and 2
    call sub_site_parm(wf, frag1_homo, frag2_lumo)
    call sub_Lowdin2d(wf)
    call sub_write_frontier(file_id, wf)

    ! for LUMO-HOMO between fragment 1 and 2
    call sub_site_parm(wf, frag1_lumo, frag2_homo)
    call sub_Lowdin2d(wf)
    call sub_write_frontier(file_id, wf)

    ! for HOMO-LUMO of fragment 1 and 2
    sp1 = site_polarizable(wf, frag1_homo, frag1_lumo, 1)
    sp2 = site_polarizable(wf, frag2_homo, frag2_lumo, 2)
51  format(2F15.6)
    write(file_id, '(1X,A30,2F15.6)') 'SITE POLARIZABLE ENERGY:', sp1*AUTOEV, sp2*AUTOEV

    close(file_id)

end subroutine


subroutine sub_write_frontier(file_id, wf)
    use constant
    use wavefunction
    type(wfinfo) :: wf
    integer :: file_id

101 format(2I10)
102 format(F15.6)
103 format(4F15.6)
    write(file_id,101) wf%i_mo, wf%j_mo
    write(file_id,102) wf%overlap
    write(file_id,103) wf%parm * AUTOEV
    write(file_id,103) wf%Heff * AUTOEV

end subroutine




