!
! physical constants

module constant
    implicit none
    ! mass unit
    real(kind = 8) :: ATOMMASS, AMU
    ! energy unit
    real(kind = 8) :: AUTOEV, AUTOCM, KCALTOEV
    ! time unit
    real(kind = 8) :: AUTOFS
    ! length unit
    real(kind = 8) :: BOHRTOANG, ANGTOBOHR
    ! temperature
    real(kind = 8) :: AUTOTEMP
    ! other constant
    real(kind = 8) :: PI

    parameter (ATOMMASS=1822.88851633D+00, AMU=1.0D+0/ATOMMASS, &
        KCALTOEV=0.0433641D+00, &
        AUTOEV= 27.2113961D+00, AUTOCM=219474.0D+00, &
        AUTOFS=2.418918299D-02, &
        BOHRTOANG=0.52917720859D+00, ANGTOBOHR=1.0/BOHRTOANG, &
        AUTOTEMP = 3.1577464D+5, &
        PI=3.14159265358979323846 &
                 )

end module constant
