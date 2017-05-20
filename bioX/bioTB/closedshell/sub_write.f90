! fortran


subroutine sub_write(wf)
  use wavefunction
  type(wfinfo) :: wf
  
  integer :: outfileid = 2000

51 format(1x,2I10)
52 format(1x,A20,E15.6)
53 format(1x,A20,E15.6)
  open(outfileid, file = 'parm.dat')
  write(outfileid, 51) wf%i_mo, wf%j_mo

  write(outfileid, 52) 'OVERLAP:', wf%overlap

  write(outfileid, 53) 'H11 (eV):', wf%parm(1,1)*27.2116
  write(outfileid, 53) 'H22 (eV):', wf%parm(2,2)*27.2116
  write(outfileid, 53) 'H12 (eV):', wf%parm(2,1)*27.2116

  write(outfileid, 53) 'T11 (eV):', wf%Heff(1,1)*27.2116
  write(outfileid, 53) 'T22 (eV):', wf%Heff(2,2)*27.2116
  write(outfileid, 53) 'T12 (eV):', wf%Heff(2,1)*27.2116

  close(outfileid)

end subroutine

