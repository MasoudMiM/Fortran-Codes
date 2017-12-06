!-----------------------------------------------------
! boundary condition function for ELMER:
! 
!-----------------------------------------------------
FUNCTION myBC_z( model, n, tim) RESULT(Displ)
! modules needed
USE DefUtils
IMPLICIT NONE

!------- variables
TYPE(Model_t) :: model
INTEGER :: n
REAL(KIND=dp):: tim, Rad, DispR, Displ

!------- variables for reading the file
LOGICAL :: FirstVisited = .TRUE.
integer, parameter :: nxd = 501, nyd = 1000
real(kind =8), dimension(nyd,nxd) :: Force_r = 0.0
real(kind =8 ), dimension(nxd) :: Radial
real(kind =8), dimension(nyd) :: TimeVector

integer :: row, col


!------- variables for the BC load
REAL(KIND=dp), PARAMETER :: PI_8  = 4 * atan (1.0_8)
REAL(KIND=dp) :: xcor, ycor, zcor, angle

!------- remember these variables
SAVE Force_r, Radial, TimeVector, FirstVisited

! get nodal coordinates and calculating angle
! -----------------
xcor = model % Nodes % x(n)
ycor = model % Nodes % y(n)
zcor = model % Nodes % z(n)
angle = ATAN2(ycor,xcor)

!-----------------------------------------------------------------------
!----------- Start: Reading the excitation signal file -----------------
If (FirstVisited) THEN 

!-------------------- Reading Force values (r-dir.)
open(10, file="outfile_Force_z.dat")
!READ(10,*) ((Force_r(row,col), row = 1, nyd), col = 1, nxd) ! Column-Wise
READ(10,*) ((Force_r(row,col), col = 1, nxd), row = 1, nyd) ! Row-Wise (Always use this one!)
close(10)
print *, "The size of Force Matrix is", SHAPE(Force_r)

!-------------------- Reading r-coordinate values
open(20, file="outfile_Loc_r.dat")
read(20,*) Radial
close(20)
print *, "The size of Radial Vector is", SIZE(Radial)

!-------------------- Reading time values 
open(30, file="outfile_Time.dat")
read(30,*) TimeVector
close(30)
print *, "The size of Time Vector is", SIZE(TimeVector)


FirstVisited = .FALSE.
END IF

! Use the values of xcor and ycor to find the location in cylindrical coordinate 
! and then use that value for radius and interpolation to find the displacement
Rad = SQRT(xcor**2+ycor**2)

!######################################################## Interpolating
IF (tim < TimeVector(nyd)) THEN
  call pwl_interp_2d ( nxd, nyd, Radial, TimeVector, transpose(Force_r), Rad, tim, DispR )
ELSE
  DispR = 0.0
END IF


!------ load calculation
!CALL nearest_interp_1d ( ndata, t, amp, 1, tim, Disp0 )
Displ = DispR

! ------------------------------------------------------ Module

contains
! Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    John Burkardt (slightly modified by Masoud Masoumi for one signle data at a time)
  
  subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, xii, yii, zii )
    implicit none
    integer ( kind = 4 ) nxd
    integer ( kind = 4 ) nyd
    real ( kind = 8 ) alpha
    real ( kind = 8 ) beta
    real ( kind = 8 ) det
    real ( kind = 8 ) dxa
    real ( kind = 8 ) dxb
    real ( kind = 8 ) dxi
    real ( kind = 8 ) dya
    real ( kind = 8 ) dyb
    real ( kind = 8 ) dyi
    real ( kind = 8 ) gamma
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    !real ( kind = 8 ) r8_huge
    !integer ( kind = 4 ) r8vec_bracket5_func
    real ( kind = 8 ) xd(nxd)
    real ( kind = 8 ) xi(1), xii
    real ( kind = 8 ) yd(nyd)
    real ( kind = 8 ) yi(1), yii
    real ( kind = 8 ) zd(nxd,nyd)
    real ( kind = 8 ) zi(1), zii

    !do k = 1, ni
    xi(1) = xii
    yi(1) = yii
    k = 1

      i = r8vec_bracket5 ( nxd, xd, xi(k) )
      if ( i == -1 ) then
        zi(k) = r8_huge ( )
        !cycle
      end if

      j = r8vec_bracket5 ( nyd, yd, yi(k) )
      if ( j == -1 ) then
        zi(k) = r8_huge ( )
        !cycle
      end if

      if ( yi(k) < yd(j+1) &
        + ( yd(j) - yd(j+1) ) * ( xi(i) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

        dxa = xd(i+1) - xd(i)
        dya = yd(j)   - yd(j)

        dxb = xd(i)   - xd(i)
        dyb = yd(j+1) - yd(j)

        dxi = xi(k)   - xd(i)
        dyi = yi(k)   - yd(j)

        det = dxa * dyb - dya * dxb

        alpha = ( dxi * dyb - dyi * dxb ) / det
        beta =  ( dxa * dyi - dya * dxi ) / det
        gamma = 1.0D+00 - alpha - beta

        zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

      else

        dxa = xd(i)   - xd(i+1)
        dya = yd(j+1) - yd(j+1)

        dxb = xd(i+1) - xd(i+1)
        dyb = yd(j)   - yd(j+1)

        dxi = xi(k)   - xd(i+1)
        dyi = yi(k)   - yd(j+1)

        det = dxa * dyb - dya * dxb

        alpha = ( dxi * dyb - dyi * dxb ) / det
        beta =  ( dxa * dyi - dya * dxi ) / det
        gamma = 1.0D+00 - alpha - beta

        zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

      end if

    !end do
    zii = zi(1)

    return
  end subroutine pwl_interp_2d
  
  function r8vec_bracket5 ( nd, xd, xi )

    implicit none

    integer ( kind = 4 ) nd

    integer ( kind = 4 ) b
    integer ( kind = 4 ) l
    integer ( kind = 4 ) m
    integer ( kind = 4 ) r
    integer ( kind = 4 ) r8vec_bracket5
    real ( kind = 8 ) xd(nd)
    real ( kind = 8 ) xi

    if ( xi < xd(1) .or. xd(nd) < xi ) then

      b = -1

    else

      l = 1
      r = nd

      do while ( l + 1 < r )
        m = ( l + r ) / 2
        if ( xi < xd(m) ) then
          r = m
        else
          l = m
        end if
      end do

      b = l

    end if

    r8vec_bracket5 = b

    return
  end function r8vec_bracket5

  function r8_huge ( )
    !implicit none

    real ( kind = 8 ) r8_huge

    r8_huge = 0.0 !1.79769313486231571D+308

    return
  end function r8_huge


END FUNCTION myBC_z

