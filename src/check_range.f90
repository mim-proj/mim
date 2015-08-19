!
! Function
!   check range of the variables
!
! Arguements (in)
!   im   : number of x-direction grid
!   jm   : number of y-direction grid
!   km   : number of z-direction grid
!   var  : variable to be checked
!   min  : possible minimum value of var
!   max  : possible maximum value of var
!   from : subroutine name (for message)
!   name : variable name (for message)
!
subroutine check_range( im, jm, km, var, min, max, from, name )
  implicit none
  integer,intent(in) :: im, jm, km
  real(4),intent(in) :: var(im, jm, km)
  real(4),intent(in) :: min, max
  character(*),intent(in) :: from, name
  integer,save :: warn_count = 1
  integer,parameter :: max_warn_count = 100
  
  integer :: i, j, k
  
  do k=1, km
     do j=1, jm
        do i=1, im

           if( var(i,j,k) < min .or. var(i,j,k) > max ) then
              write(*,*) 'warning in check_range() : ' 
              write(*,*) trim(name) // ' should be between ', min, 'and', max
              write(*,*) '(i,j,k) = ', i, j, k
              write(*,*) 'var = ', var(i,j,k)
              write(*,*) 'called from ' // trim(from)
              write(*,*) 'warning count = ', warn_count
              
              if( warn_count > max_warn_count ) then
                 write(*,*) 'warning count exceeds maximum number'
                 stop
              end if
              warn_count = warn_count + 1

           end if

        end do
     end do
  end do
  
end subroutine check_range
