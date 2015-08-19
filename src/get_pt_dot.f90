!
! Function
!   estimate D(pt)/Dt (i.e. diabatic heating) from total differential
!
! Arguements (in)
!   dt        : timestep of the input date
!   u         : x-wind
!   v         : y-wind
!   omega     : omega velocity
!   pt        : potential temperature at the pressure levels
!   pt_before : pt 1-step before 
!
! Arguements (out)
!   pt_dot    : D(pt)/Dt
!
! Note
!   -If omega=0 (no input), vertical advection of pt will be neglected.
!   -If pt=pt_before (e.g. first time step), local pt change will be neglected.
!
subroutine get_pt_dot_omega( dt, u, v, omega, pt, pt_before, &
     &                       pt_dot )
  use parameter, only: radius, pai
  use com_var, only : im, jm, km, alat, costbl, pin
  implicit none
  real(4),intent(in)  :: dt
  real(4),intent(in)  :: u(im, jm, km)
  real(4),intent(in)  :: v(im, jm, km)
  real(4),intent(in)  :: omega(im, jm, km)
  real(4),intent(in)  :: pt(im, jm, km)
  real(4),intent(in)  :: pt_before(im, jm, km)
  real(4),intent(out) :: pt_dot(im, jm, km)

  integer :: i, j, k

  do k=1, km                    
     do j=1, jm
        do i=1, im

           ! d(pt)/dt
           pt_dot(i,j,k) = ( pt(i,j,k) - pt_before(i,j,k) ) / dt 

           ! u d(pt)/dx
           if( i == 1 ) then
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + u(i,j,k) * ( pt(i+1,j,k) - pt(im,j,k) ) &
                   &          / ( 4 * pai * radius * costbl(j) / im ) 
           else if( i == im ) then
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + u(i,j,k) * ( pt(1,j,k) - pt(i-1,j,k) ) &
                   &          / ( 4 * pai * radius * costbl(j) / im )
           else
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + u(i,j,k) * ( pt(i+1,j,k) - pt(i-1,j,k) ) &
                   &          / ( 4 * pai * radius * costbl(j) / im )
           endif

           ! v d(pt)/dy
           if( j == 1 ) then
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + v(i,j,k) * ( pt(i,j+1,k) - pt(i,1,k) ) &
                   &          / ( alat(j+1) - alat(1) ) / radius 
           else if( j == jm ) then
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + v(i,j,k) * ( pt(i,jm,k) - pt(i,j-1,k) ) &
                   &          / ( alat(jm) - alat(j-1) ) / radius 
           else
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + v(i,j,k) * ( pt(i,j+1,k) - pt(i,j-1,k) ) &
                   &          / ( alat(j+1) - alat(j-1) ) / radius 
           endif

           ! omega d(pt)/dp
           if( k == 1 ) then
              continue
           else if( k == km ) then
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &           + omega(i,j,k) * (pt(i,j,k) - pt(i,j,k-1) ) &
                   &             / ( ( pin(k) - pin(k-1) ) * 100 ) 
           else
              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &          + omega(i,j,k) * ( pt(i,j,k+1) - pt(i,j,k-1) ) &
                   &            / ( ( pin(k+1) - pin(k-1) ) * 100 )
           end if

        end do
     end do
  end do

end subroutine get_pt_dot_omega



!
! Function
!   estimate D(pt)/Dt (i.e. diabatic heating) from diabatic heating input
!
! Arguements (in)
!   q_3d      : input diabatic heating data
!
! Arguements (out)
!   pt_dot    : D(pt)/Dt
!
subroutine get_pt_dot_q( q_3d, pt_dot )
  use parameter, only: cp, rkappa
  use com_var, only : im, jm, km, pin
  implicit none
  real(4),intent(in)  :: q_3d(im, jm, km)  ! [J/(kg s)]
  real(4),intent(out) :: pt_dot(im, jm, km)
  
  integer :: i, j, k

  do k=1, km
     do j=1, jm
        do i=1, im
           pt_dot(i,j,k) = q_3d(i,j,k) / ( cp * (pin(k)/1000.0)**rkappa )
        end do
     end do
  end do

end subroutine get_pt_dot_q


subroutine get_pt_dot_q_inv( pt_dot, q_3d )
  use parameter, only: cp, rkappa
  use com_var, only : im, jm, km, pin
  implicit none
  real(4),intent(in)  :: pt_dot(im, jm, km)
  real(4),intent(out) :: q_3d(im, jm, km)  ! [J/(kg s)]
  
  integer :: i, j, k

  do k=1, km
     do j=1, jm
        do i=1, im
           q_3d(i,j,k) =  pt_dot(i,j,k) * ( cp * (pin(k)/1000.0)**rkappa )
        end do
     end do
  end do

end subroutine get_pt_dot_q_inv


