! by C.Kodama 2008.08.04
program main
  implicit none

  ! dimension
  integer,parameter :: xdef=144, ydef=73
  integer,parameter :: zdef=12
  integer,parameter :: tdef=1460

  character(1024) :: fin  = 'omega.2007.nc'  ! input filename
  character(1024) :: fout = 'omega2007.grd'     ! output filename
  character(1024) :: netcdf_name = 'omega'   ! NetCDF variable name

  ! variable
  real(4) :: value(xdef,ydef,zdef), temp(xdef,ydef,zdef)

  ! for NetCDF open
  include '/usr/local/include/netcdf.inc'  ! include NetCDF
  integer,parameter :: omode=0
  integer :: ncID, varID, status
  real(4) :: scale, offset
  integer :: start(4) = (/1, 1, 1, 1/)
  integer :: cnt(4) = (/xdef, ydef, zdef, 1/)

  integer :: t, irecl

  ! open input NetCDF
  status = nf_open( fin, omode, ncID ) ! open
  status = nf_inq_varid( ncID, TRIM(netcdf_name), varID )
  status = nf_get_att_real( ncID, varID, 'scale_factor', scale )
  status = nf_get_att_real( ncID, varID, 'add_offset', offset )
    
  ! open output GrADS file
  open( 50, file=fout,  status='unknown',        &
       &  form='unformatted', access='direct', recl=xdef*ydef*zdef*4 )

  irecl = 1
  do t=1, tdef
     write(*,*) t, "/", tdef

     ! read
     start(4) = t
     status = nf_get_vara_real( ncID, varID, start, cnt, temp ) 
     if( status .ne. NF_NOERR ) then
        temp(:,:,:) = -9.99e+33
     end if

     value(:,:,:) = temp(:,:,:) * scale + offset 

     ! write
     write(50,rec=irecl) value(:,:,:)
     irecl = irecl + 1

  end do

  close(50)

end program main
