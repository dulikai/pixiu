! fortran
!
! Lowdin orthogonalization procedure
! Lowdin, P. O. J. Chem. Phys. 1950, 18, 365
! JACS, 2006, 128, 9882-9886
! J. Chem. Phys. 2003, 119, 9809

!

subroutine sub_Lowdin2d(wf)
  use wavefunction
  implicit none
  type(wfinfo) :: wf
  real(kind = 8) :: h11, h22, h12, s12
  real(kind = 8) :: h12_eff, h11_eff, h22_eff
  real(kind = 8) :: vs, e12p, e12m

  h11 = wf%parm(1,1)
  h22 = wf%parm(2,2)
  h12 = wf%parm(2,1)
  s12 = wf%overlap

  vs = 1.0d0 - s12 * s12
  e12p = h11 + h22
  e12m = h11 - h22

  h12_eff = (h12-0.5d0*e12p*s12) / vs
  h11_eff = 0.5 * (e12p-2.0d0*h12*s12+e12m*DSQRT(vs)) / vs
  h22_eff = 0.5 * (e12p-2.0d0*h12*s12-e12m*DSQRT(vs)) / vs
  
  wf%Heff(1,1) = h11_eff
  wf%Heff(2,2) = h22_eff
  wf%Heff(2,1) = h12_eff
  wf%Heff(1,2) = wf%Heff(2,1)

  return

end subroutine

