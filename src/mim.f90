program mim
!############################################################
!
!  Program MIM  
!
!
!    Please see Readme_namelist_en.txt for details
!
!
!############################################################
  use namelist
  use parameter
  use com_var, only : im, jm, km, ko, wmax, &
       &              pin, pout, rho, alat, sintbl, costbl, &
       &              com_var_ini, com_var_end, &
       &              com_var_pin, com_var_pout, com_var_zd, com_var_rho, &
       &              com_var_alat, &
       &              com_var_sintbl, com_var_costbl, com_var_tantbl
  use mim_var
  use grads, only : grads_info, &
       &            grads_open, grads_close, &
       &            grads_read, grads_write
  use biseki, only : biseki_ini, biseki_biseki, biseki_sekibun, biseki_bibun
  use biseki_y, only : biseki_y_ini, biseki_y_biseki_y, &
       &               biseki_y_sekibun_y, biseki_y_bibun_y
  implicit none

  ! temporal variables
  integer :: tstep
  integer :: i, j, k, w
  integer :: icount
  real(4) :: cons
  
  ! input information
  type(grads_info) :: ginfo_topo  ! for INPUT_TOPO
  type(grads_info) :: ginfo_uvt   ! for INPUT_UVT
  type(grads_info) :: ginfo_u     ! for INPUT_U
  type(grads_info) :: ginfo_v     ! for INPUT_V
  type(grads_info) :: ginfo_t     ! for INPUT_T
  type(grads_info) :: ginfo_z     ! for INPUT_Z
  type(grads_info) :: ginfo_omega ! for INPUT_OMEGA
  type(grads_info) :: ginfo_ps    ! for INPUT_PS
  type(grads_info) :: ginfo_msl   ! for INPUT_MSL
  type(grads_info) :: ginfo_ts    ! for INPUT_TS
  type(grads_info) :: ginfo_q     ! for INPUT_Q

  ! output information
  type(grads_info) :: ginfo_zonal  ! for OUTPUT_ZONAL
  type(grads_info) :: ginfo_vint   ! for OUTPUT_VINT
  type(grads_info) :: ginfo_gmean  ! for OUTPUT_GMEAN
  type(grads_info) :: ginfo_wave   ! for OUTPUT_WAVE


!  write(*,*) "MIM Version dev"
  write(*,*) "MIM Version 0.33r1" ! Please rewrite Makefile


  !********** initiale & allocate memory **********

  ! initialize namelist
  call namelist_init()

  ! set common variables
  call com_var_ini( INPUT_XDEF_NUM, INPUT_YDEF_NUM, INPUT_ZDEF_NUM, &
       &            OUTPUT_ZDEF_NUM, &
       &            WAVE_MAX_NUMBER )

  ! allocate arrays for mim.f90
  call mim_var_ini( INPUT_XDEF_NUM, INPUT_YDEF_NUM, INPUT_ZDEF_NUM, &
       &            OUTPUT_ZDEF_NUM, &
       &            WAVE_MAX_NUMBER )

  ! standard pressure and density
  call com_var_pin()
  call com_var_pout()
  call com_var_zd()
  call com_var_rho()

  ! latitude ( alat [rad] ) and sin/cos/tan table 
  !   alat: YREV (i.e. north -> south)
  call com_var_alat()
  call com_var_sintbl()
  call com_var_costbl()
  call com_var_tantbl()

  ! calendar
  call calendar_tdef( tstep )
  write(*,*) 'tstep =', tstep


  !********** open input file **********!
  ! - all variables should be YREV (north->south) 
  !   and ZREV (upper->lower) in mim.f90.

  !***** input *****!
  ! - convert input data to YREV IF input data is not YREV.
  ! - convert input data to ZREV, ASSUMING input data is not YREV.

  if( INPUT_UVT_FILENAME /= '' ) then
     call grads_open( 10, INPUT_UVT_FILENAME, im, jm, 1, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_UVT_INT, &
          &           ginfo_uvt )
  else
     call grads_open( 11, INPUT_U_FILENAME, im, jm, km, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_U_INT, &
          &           ginfo_u )
     call grads_open( 12, INPUT_V_FILENAME, im, jm, km, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_V_INT, &
          &           ginfo_v )
     call grads_open( 13, INPUT_T_FILENAME, im, jm, km, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_T_INT, &
          &           ginfo_t )
  end if

  if( INPUT_PS_FILENAME /= '' ) then
     call grads_open( 20, INPUT_PS_FILENAME, im, jm, 1, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1, &
          &           INPUT_ENDIAN_PS_INT, &
          &           ginfo_ps )
  else
     call grads_open( 21, INPUT_MSL_FILENAME, im, jm, 1, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1, &
          &           INPUT_ENDIAN_MSL_INT, &
          &           ginfo_msl )
     call grads_open( 22, INPUT_TS_FILENAME, im, jm, 1, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1, &
          &           INPUT_ENDIAN_TS_INT, &
          &           ginfo_ts )
  end if

  call grads_open( 30, INPUT_TOPO_FILENAME, im, jm, 1, &
       &           0, 1-INPUT_YDEF_YREV_TOPO, 1, &
       &           INPUT_ENDIAN_TOPO_INT, &
       &           ginfo_topo )

  call grads_open( 40, INPUT_Z_FILENAME, im, jm, km, &
       &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
       &           INPUT_ENDIAN_Z_INT, &
       &           ginfo_z )

  if( INPUT_OMEGA_FILENAME /= '' ) then
     call grads_open( 50, INPUT_OMEGA_FILENAME, im, jm, 1, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_OMEGA_INT, &
          &           ginfo_omega )
  end if

  if( INPUT_Q_FILENAME /= '' ) then
     call grads_open( 60, INPUT_Q_FILENAME, im, jm, km, &
          &           0, 1-INPUT_YDEF_YREV_DEFAULT, 1-INPUT_ZDEF_ZREV, &
          &           INPUT_ENDIAN_Q_INT, &
          &           ginfo_q )
  end if


  !***** output *****!
  ! - output data is YREV, but NOT ZREV.
  ! - output data is LITTLE ENDIAN.

  call grads_open( 100, OUTPUT_ZONAL_FILENAME, 1, jm, ko, &
       &           0, 0, 1, 1, &
       &           ginfo_zonal )

  if( OUTPUT_VINT_FILENAME /= '' ) then
     call grads_open( 110, OUTPUT_VINT_FILENAME, 1, jm, 1, &
          &           0, 0, 1, 1, &
          &           ginfo_vint )
  end if

  if( OUTPUT_GMEAN_FILENAME /= '' ) then
     call grads_open( 120, OUTPUT_GMEAN_FILENAME, 1, 1, 1, &
          &           0, 0, 1, 1, &
          &           ginfo_gmean )
  end if

  if( OUTPUT_WAVE_FILENAME /= '' .and. WAVE_MAX_NUMBER > 0 ) then
     call grads_open( 130, OUTPUT_WAVE_FILENAME, 1, jm, ko, &
          &           0, 0, 1, 1, &
          &           ginfo_wave )
  end if


  ! write iteration error on log_file
  open(65, file=OUTPUT_ERROR_FILENAME)
  write(65,*) 'OUTPUT FILE', trim(OUTPUT_ZONAL_FILENAME)
  write(65,*) "warning(getpt):"
  write(65,*) "iteration was finished not approriately below:"


  !***** load topography ( alt [m] ) *****!
  ! - convert topography data to YREV IF it is not YREV.
  call grads_read( ginfo_topo, alt )

  if( INPUT_UNIT_TOPO == 'm^2/s^2' ) then
     alt(:,:) = alt(:,:) / grav  ! for geopotential [m^2/s^2] -> [m]

  else if( INPUT_UNIT_TOPO == 'm' ) then
     alt(:,:) = alt(:,:)  ! for geopotential height [m]

  else
     write(0,*) 'error in lagmain : INPUT_UNIT_TOPO = ' &
          &     // INPUT_UNIT_TOPO // 'is invalid'
     stop
  end if

  call check_range( im, jm, 1, alt, alt_min, alt_max, 'lagmain()', 'alt' )

  if( sum( alt(:,jm/4) ) < sum( alt(:,jm-jm/4) ) ) then  ! check YREV
     write(0,*) 'error in lagmain.f90: Altitude seems to be wrong'
     write(0,*) 'Probably, INPUT_YDEF_YREV_TOPO in namelist should be changed'
     stop
  end if


  !**************************************************!
  !                                                  !
  !                    time loop                     !
  !                                                  !
  !**************************************************!
  do icount=1, tstep
     write(*,*) 'icount=', icount, '/' , tstep


     !********** load data **********!

     !***** UVT(RH) *****!
     if( INPUT_UVT_FILENAME /= '' ) then
        do k=km, 1, -1
           call grads_read( ginfo_uvt, u(1:im,1:jm,k) )
        end do

        do k=km, 1, -1
           call grads_read( ginfo_uvt, v(1:im,1:jm,k) )
        end do
     
        do k=km, 1, -1
           call grads_read( ginfo_uvt, t(1:im,1:jm,k) )
        end do

        ! for classic NCEP-tohoku style (incluing RH)
        if( INPUT_TYPE == 'include RH8' ) then
           ginfo_uvt%record = ginfo_uvt%record + 8
        end if
     else
        call grads_read( ginfo_u, u )
        call grads_read( ginfo_v, v )
        call grads_read( ginfo_t, t )
     end if


     !***** surface pressure [hPa] *****!     
     if( INPUT_PS_FILENAME /= '' ) then
        call grads_read( ginfo_ps, p_sfc )
     
        if( INPUT_UNIT_PS == 'Pa' ) then
           p_sfc(:,:) = p_sfc(:,:) / 100.0  ! [Pa] -> [hPa]
        end if
     else
        call grads_read( ginfo_msl, msl )
        call grads_read( ginfo_ts, t_sfc )
        if( INPUT_UNIT_MSL == 'Pa' ) then
           msl(:,:) = msl(:,:) / 100.0  ! [Pa] -> [hPa]
        end if
        ! mean sea level -> surface
        p_sfc(:,:) = msl(:,:) * ( 1 + gamma * alt(:,:) / t_sfc(:,:) ) &
             &                  ** ( -grav / ( gasr * gamma ) )
     end if

     !***** geopotential height [m] *****!
     call grads_read( ginfo_z, z )

     if( INPUT_UNIT_Z == 'm^2/s^2' ) then
        z(:,:,:) = z(:,:,:) / grav  ! for geopotential [m^2/s^2] -> [m]

     else if( INPUT_UNIT_Z == 'm' ) then
        z(:,:,:) = z(:,:,:)  ! for geopotential height [m]

     else
        write(0,*) 'error in lagmain : INPUT_UNIT_Z = ' &
          &     // INPUT_UNIT_Z // 'is invalid'
        stop
     end if


     !***** omega *****!
     ! omega is assumed to be 0 if (part of) omega is not prepared.
     if( INPUT_OMEGA_FILENAME /= "" ) then
        
        do k=km, km-INPUT_ZDEF_NUM_OMEGA+1, -1
           call grads_read( ginfo_omega, omega(1:im,1:jm,k) )
        end do

        do k=km-INPUT_ZDEF_NUM_OMEGA, 1, -1
           omega(1:im,1:jm,k) = 0
        end do
        
     else
        omega(:,:,:) = 0
     end if


     !***** diabatic Heating *****!
     if( INPUT_Q_FILENAME /= "" ) then
        call grads_read( ginfo_q, q_3d )
     else
        q_3d(:,:,:) = 0
     end if

     where( q_3d /= INPUT_UNDEF_Q_REAL )
        q_3d = q_3d * cp        ! [K/s] (dT/dt) -> [J/(kg s)]
     end where


     !********** interpolate/extrapolate to undef data **********!
     call undef_fill( im, jm, km, INPUT_UNDEF_U_REAL    , pin, u )
     call undef_fill( im, jm, km, INPUT_UNDEF_V_REAL    , pin, v )
     call undef_fill( im, jm, km, INPUT_UNDEF_T_REAL    , pin, t )
     call undef_fill( im, jm, km, INPUT_UNDEF_Z_REAL    , pin, z )
     call undef_fill( im, jm, km, INPUT_UNDEF_OMEGA_REAL, pin, omega )
     call undef_fill( im, jm, km, INPUT_UNDEF_Q_REAL    , pin, q_3d )
     

     !********** check value **********!
     call check_range( im, jm, km, u, wind_min, wind_max, 'mim()', 'u' )
     call check_range( im, jm, km, v, wind_min, wind_max, 'mim()', 'v' )
     call check_range( im, jm, km, t, t_min, t_max, 'mim()', 't' )
     call check_range( im, jm, 1, p_sfc, p_min, p_max, 'mim()', 'p_sfc' )
     call check_range( im, jm, km, z, z_min, z_max, 'mim()', 'z' )
     call check_range( im, jm, km, omega, omega_min, omega_max, &
          &            'mim()', 'omega' )


     !********** start diagnosis **********!
  
     !***** surface pressure *****!
     p_pds = sum( p_sfc, dim=1 ) / im   ! p+s (p-dagger at the surface)
                                        ! i.e. zonal mean p_sfc
     call integral_meridional( 1, jm, 1, alat, p_pds, &
          &                    p_pdds ) ! p++s (p-dagger-dagger at the surface)
                                        ! i.e. global mean p_sfc


     !***** prepare for p -> p+ transformation *****!
     ! p+ or pd : zonal mean pressure on the isentropic surface

     ! pt     : 3-dimensional potential temperature at the standart p levels
     ! pt_sfc : 2-dimensional porential temperature at the surface
     call setpt0( t, p_sfc, &
          &       pt, pt_sfc )
     if( icount == 1 ) pt_past(:,:,:) = pt(:,:,:)

     ! pt_pds : potential temperature at the surface in the p+ coordinate
     !          i.e. min(pt_sfc)
     pt_pds(:) = minval( pt_sfc, dim=1 )

     ! pt_zm : (zonal mean) pt at the standard p+ levels
     ! p_pd  : 3-dimensional pressure at the standard p+ levels
     !         p_pd = p_pds under the ground (i.e. pout > p_pds)
     ! p_zm  : zonal mean pressure at the standard p+ levels
     !         p_zm must be almost equal to pout (standard p+ levels)
     !         except under the ground (p_zm = p_pds)
     do j=1, jm
        call getpt( im, km, ko, icount, pin, pout, p_sfc(:,j), p_pds(j), &
             &      pt(:,j,:), pt_sfc(:,j), pt_pds(j), &
             &      dlev(:,j,:), nlev(:,j,:), &
             &      p_pd(:,j,:), p_zm(j,:), pt_zm(j,:) )
     end do
     if( icount == 1 ) p_pd_past(:,:,:) = p_pd(:,:,:)

     ! T+ : temperature on p+ (different from t_zm)
     t_dagger(:,:) = pt_zm(:,:) * ( spread(pout,1,jm) / 1000.0 )**rkappa
     call check_range( 1, jm, ko, t_dagger, t_min, t_max, &
          &            'lagmain()', 't_dagger' )

     ! pd_p : 3-dimensional p+ on the standard pressure levels
     call intpl_pd_p( im, jm, km, ko, p_zm, pt, pt_zm, &
          &           pd_p )


     !***** prepare for p+ -> p++ transformation *****!
     ! p++ : global mean pressure on the isentropic surface

     ! pt_pdds : potential temperature at the surface in the p++ coordinate
     pt_pdds = minval( pt_pds, dim=1 )

     ! pt_ym  : (global mean) pt at the standard p++ levels
     ! pd_pdd : 2-dimensional p+ at the standard p++ levels
     !          pd_pt = p_pdds under the ground (i.e. pout > p_pdds)
     ! pd_ym  : global mean pressure at the standard p++ levels
     !          pd_ym must be equal to pout (standard p++ levels itself)
     !          provided the accuracy is perfect.
     !          pd_ym = p_pdds under the ground (i.e. pout > p_pdds).
     call getpt_y( jm, ko, ko, icount, &
          &        pout, pout, alat, p_pds, p_pdds, pt_zm, pt_pds, pt_pdds, &
          &        dlev_y, nlev_y, pd_pdd, pd_ym, pt_ym )

     ! pdd_pd : 2-dimensional p++ on the standard p+ levels
     call intpl_pdd_pd( jm, ko, pd_ym, pt_zm, pt_ym, &
          &             pdd_pd )

!     do k=1, ko
!        do j=1, jm
!           write(*,*) j, k, pt_zm(j,nlev_y(j,k)), pt_ym(k), pt_zm(j,nlev_y(j,k)+1)
!        end do
!        write(*,*) k, pt_ym(k)
!        call integral_meridional(1, jm, 1, alat, pd_pdd(:,k), temp)
!        write(*,*) temp
!     end do
!     stop


     !***** initialize biseki *****!
     call biseki_ini( p_pd, p_sfc, p_pds, nlev, dlev )
     call biseki_y_ini( pd_pdd, p_pds, p_pdds, nlev_y, dlev_y )


     !***** zonal wind *****
     call biseki_biseki( u, u_zm )
     call check_range( 1, jm, ko, u_zm, wind_min, wind_max, 'mim()', 'u_zm' )


!     !***** temperature (not temperature dagger) *****
!     call biseki_biseki( t, t_zm )


     !***** meridional wind & streamfunction *****!
     call biseki_sekibun( v, x_pd, xint_zm )  ! xint_zm���g�p
     call biseki_bibun( xint_zm, v_zm )
     call check_range( 1, jm, ko, v_zm, wind_min, wind_max, 'mim()', 'v_zm' )
     
     do j=1, jm
        cons = 2.0 * pai * radius * 100.0 / grav * costbl(j)        
        do k=1, ko
           st_zm(j,k) = xint_zm(j,k) * cons
        end do
     end do
     call check_range( 1, jm, ko, st_zm, st_min, st_max, 'mim()', 'st_zm' )


     !***** vertical velocity derived from streamfunction *****!
     call w_from_st( jm, ko, sintbl, pout, p_pds, st_zm, &
          &          w_zm )
     call check_range( 1, jm, ko, w_zm, w_min, w_max, 'mim()', 'w_zm' )
     

     !***** D(pt)/Dt (D: Lagrangian) *****!
     if( INPUT_Q_FILENAME == "" ) then
        call get_pt_dot_omega( INPUT_TDEF_DT, u, v, omega, pt, pt_past, &
             &                 pt_dot )        ! pt, u, v, w -> D(pt)/dt
        call get_pt_dot_q_inv( pt_dot, q_3d )  ! D(pt)/dt -> Q

     else
        call get_pt_dot_q( q_3d, pt_dot )  ! Q -> D(pt)/dt
     end if
     call biseki_biseki( pt_dot, pt_dot_zm )


     !***** correlation (1) : wind and/or pt_dot *****!
     ! u_v_x_zm : (u' v')_zm 
     work(:,:,:) = u(:,:,:) * v(:,:,:)
     call biseki_biseki( work, u_v_zm )
     u_v_x_zm(:,:) = u_v_zm(:,:) - u_zm(:,:) * v_zm(:,:)

     ! u_u_x_zm : (u'^2)_zm 
     work(:,:,:) = u(:,:,:)**2
     call biseki_biseki( work, u_u_zm )
     u_u_x_zm(:,:) = u_u_zm(:,:) - u_zm(:,:)**2

     ! v_v_x_zm : (v'^2)_zm 
     work(:,:,:) = v(:,:,:)**2
     call biseki_biseki( work, v_v_zm )
     v_v_x_zm(:,:) = v_v_zm(:,:) - v_zm(:,:)**2

     ! u_u_v_zm : (u^2 v)_zm
     work(:,:,:) = u(:,:,:)**2 * v(:,:,:)
     call biseki_biseki( work, u_u_v_zm )

     ! v_v_v_zm : (v^3)_zm
     work(:,:,:) = v(:,:,:)**3
     call biseki_biseki( work, v_v_v_zm )

     ! u_pt_dot_x_zm : (u' pt_dot')_zm 
     work(:,:,:) = u(:,:,:) * pt_dot(:,:,:)
     call biseki_biseki( work, u_pt_dot_zm )
     u_pt_dot_x_zm(:,:) = u_pt_dot_zm(:,:) - u_zm(:,:) * pt_dot_zm(:,:)

     ! v_pt_dot_x_zm : (v' pt_dot')_zm
     work(:,:,:) = v(:,:,:) * pt_dot(:,:,:)
     call biseki_biseki( work, v_pt_dot_zm )
     v_pt_dot_x_zm(:,:) = v_pt_dot_zm(:,:) - v_zm(:,:) * pt_dot_zm(:,:)

     ! u_u_pt_dot_zm : (u^2 pt_dot)_zm
     work(:,:,:) = u(:,:,:)**2 * pt_dot(:,:,:)
     call biseki_biseki( work, u_u_pt_dot_zm )

     ! v_v_pt_dot_zm : (v^2 pt_dot)_zm 
     work(:,:,:) = v(:,:,:)**2 * pt_dot(:,:,:)
     call biseki_biseki( work, v_v_pt_dot_zm )


     !***** geopotential height [m] *****!
     ! z_pd : 3-dimensional geopotential height on the standard p+ levels
     ! z_zm : zonal mean geopotential height on the standard p+ levels
     ! 
     ! for form drag calculation
     call get_z_pt( alt, t, z, p_pd, p_pds, &
          &         z_pd, z_zm )
     if( icount == 1 ) z_pd_past(:,:,:) = z_pd(:,:,:)


     !***** epy : meridional component of EP flux *****!
     call epflux_y( u_v_x_zm, epy )
     call epflux_div_y( epy, depy )


     !***** epz_form : vertical component of EP flux (form drag) *****
     call epflux_z_form( p_pd, z_pd, p_sfc, &
          &              epz_form )
     call epflux_div_z( epz_form, depz_form )


     !***** wavenumber decomposition of form drag *****
     if( OUTPUT_WAVE_FILENAME /= '' .and. WAVE_MAX_NUMBER > 0 ) then
        call epflux_z_form_wave( p_pd, z_pd, p_sfc, &
             &                   wmax, p_pt_wave, z_pt_wave, &
             &                   epz_wave )
     end if


     !***** epz_uw : vertical component of EP flux (except form drag) *****
     ! epz_uw = epz_uv + epz_ut
     call epflux_z_uw( pt_zm, u_v_x_zm, u_pt_dot_x_zm, &
          &            epz_uv, epz_ut, epz_uw )
     call epflux_div_z( epz_uv, depz_uv )
     call epflux_div_z( epz_ut, depz_ut )
     call epflux_div_z( epz_uw, depz_uw )


     !***** vertical component of EP flux & its divergence *****!
     epz(:,:) = epz_form(:,:) + epz_uw(:,:)
     call epflux_div_z( epz, depz ) 


     !***** divF *****!
     divf(:,:) = depy(:,:) + depz(:,:)


     !***** G Flux (meridional momentum equation) *****!
     call gflux_y( v_zm, v_v_zm, &
          &        gy )
     call gflux_div_y( gy, dgy )
     call gflux_z( pt_zm, v_zm, v_v_zm, v_pt_dot_x_zm, &
          &        gz )
     call gflux_div_z( gz, dgz )


     !***** phi_dagger *****!
     ! different from z_zm
     call get_phi_dagger( alt, p_pds, pt_pds, t_dagger, &
          &               phi_dagger )
     if( icount == 1 ) phi_dagger_past(:,:) = phi_dagger(:,:)
     

     !***** (dz/dy)_zm & v * dz/dlat *****!
     call derivative_y( im, jm, km, alat, z, &
          &             work )
     call biseki_biseki( work, dz_dlat_zm )
     work(:,:,:) = v(:,:,:) * work(:,:,:)
     call biseki_biseki( work, v_dz_dlat_zm )

     
     !***** u * dz/dx *****!
     call derivative_x( im, jm, km, z, &
          &             work )
     work(:,:,:) = u(:,:,:) * work(:,:,:)
     call biseki_biseki( work, u_dz_dlon_zm )


     !***** p dz/dt *****!
     p_dz_dt(:,:,:) = ( z_pd(:,:,:) - z_pd_past(:,:,:) ) &
          &         * ( p_pd(:,:,:) + p_pd_past(:,:,:) ) &
          &         / 2 * 100 / INPUT_TDEF_DT
     p_dz_dt_zm = sum( p_dz_dt, dim=1 ) / real(im)


     !***** p+ d(phi_dagger)/dt *****!
     p_dphi_dt(:,:) = ( phi_dagger(:,:) - phi_dagger_past(:,:) ) &
          &         * spread( pout, 1, jm ) * 100 / INPUT_TDEF_DT


     !***** energy *****!

     ! kz : zonal kinetic energy
     kz_zm(:,:) = 0.5 * ( u_zm(:,:)**2 + v_zm(:,:)**2 )

     ! ke : eddy kinetic energy
     ke_zm(:,:) = 0.5 * ( u_u_x_zm(:,:) + v_v_x_zm(:,:) )

     ! pz : zonal potential energy (NOT available potential energy)
     pz_zm(:,:) = ( cp - gasr ) * t_dagger(:,:) + phi_dagger(:,:)

     ! ae_total_zm : eddy available potential energy (including surface term)
     !               (recommend not to use)
     call energy_ae_total( p_pd, p_zm, p_pds, pt_zm, &
          &                ae_total_zm )

     ! ae_zm_vint : vertically integrated eddy available potential energy
     !              (not equals to vertically integrated ae_total_zm)
     call energy_ae_vint( p_pd, p_zm, p_sfc, p_pds, pt_zm, pt_pds, &
          &               ae_zm_vint )

     ! az_zm_vint : vertically integrated zonal available potential energy
     call energy_az_vint( p_pds, p_pdds, pd_pdd, pt_pdds, pt_ym, &
          &               az_zm_vint )

     ! az_gmean   : global mean zonal available potential energy
     !              (for accuracy, az_zm_vint is not used here)
     call energy_az_gmean( p_pds, p_pdds, pd_pdd, pd_ym, &
          &                pt_pdds, pt_ym, &
          &                az_gmean )


     !***** correlation (2) : Ke and something *****!
     v_ke_zm(:,:) = 0.5 * ( u_u_v_zm(:,:) - 2 * u_v_zm(:,:) * u_zm(:,:) &
          &               + u_zm(:,:)**2 * v_zm(:,:) &
          &               + v_v_v_zm(:,:) - 2 * v_v_zm(:,:) * v_zm(:,:) &
          &               + v_zm(:,:)**3 )
     
     pt_dot_ke_zm(:,:) = 0.5 * ( u_u_pt_dot_zm(:,:) &
          &                    - 2 * u_pt_dot_zm(:,:) * u_zm(:,:) &
          &                    + u_zm(:,:)**2 * pt_dot_zm(:,:) &
          &                    + v_v_pt_dot_zm(:,:) &
          &                    - 2 * v_pt_dot_zm(:,:) * v_zm(:,:) &
          &                    + v_zm(:,:)**2 * pt_dot_zm(:,:) )


     !***** energy conversion (global mean) *****!

     ! C(Az,Kz)
     call energy_conv_az_kz( alt, p_pds, v_zm, pt_zm, pt_sfc, phi_dagger, &
     &                       c_az_kz )
     
     ! C(Ae,Kz)
     call energy_conv_kz_ae( u_zm, v_zm, depz_form, dz_dlat_zm, c_az_kz, &
          &                  c_kz_ae_u, c_kz_ae_v, c_kz_ae )
     
     ! C(Ae,Ke)
     call energy_conv_ae_ke( v_zm, c_kz_ae_u, &
          &                  dz_dlat_zm, v_dz_dlat_zm, u_dz_dlon_zm, &
          &                  c_ae_ke_u, c_ae_ke_v, c_ae_ke )
     
     ! C(Kz,Ke)
     call energy_conv_kz_ke( u_zm, v_zm, u_u_x_zm,  &
          &                  depy, depz_uw, dgy, dgz, &
          &                  c_kz_ke_uy, c_kz_ke_uz, &
          &                  c_kz_ke_vy, c_kz_ke_vz, c_kz_ke_tan, &
          &                  c_kz_ke )

     ! C(Kz,W)
     c_kz_w(:,:) = c_kz_ke(:,:) + c_kz_ae(:,:)



     !***** 2-dimensional energy conversion *****!

     ! dKz/dt
     call energy_tendency_dkzdt_vkz( v_zm, kz_zm, &
          &                          dkzdt_vkz )
     call energy_tendency_dkzdt_wkz( w_zm, kz_zm, &
          &                          dkzdt_wkz )

     ! dKe/dt
     call energy_tendency_dkedt_uy( u_zm, epy, &
          &                         dkedt_uy )
     call energy_tendency_dkedt_vy( v_zm, gy, &
          &                         dkedt_vy )
     call energy_tendency_dkedt_uz( u_zm, epz_uw, &
          &                         dkedt_uz )
     call energy_tendency_dkedt_vz( v_zm, gz, &
          &                         dkedt_vz )

     call energy_tendency_dkedt_vke( v_ke_zm, dkedt_vke )

     call energy_tendency_dkedt_wke( pt_zm, v_ke_zm, pt_dot_ke_zm, &
          &                          dkedt_wke )


     !***** Diabatic Heating *****!
     ! q_zm : zonal mean diabatic heating
     call biseki_biseki( q_3d, q_zm )

     ! q_3d [J/(kg s)] -> q_ex_3d [K/s] ( = Q / Pi )
     do k=1, km
        do j=1, jm
           do i=1, im
              q_ex_3d(i,j,k) = q_3d(i,j,k) &
                   &         / ( cp * ( pin(k) / 1000.0 )**rkappa )
           end do
        end do
     end do
     call biseki_biseki( q_ex_3d, q_ex_zm )

     ! Qgz
     do k=1, ko
        do j=1, jm
           qgz_zm(j,k) = q_ex_zm(j,k) * ( pout(k) / 1000.0 )**rkappa * cp
        end do
     end do

     ! qe_zm : zonal mean eddy diabatic heating
     call diabatic_qe( im, jm, km, ko, q_3d, pin, pd_p, &
          &            qe_zm )

     ! qz_gmean : global mean zonal diabatic heating
     call diabatic_qz( jm, ko, q_ex_zm, pout, pdd_pd, &
          &            qz_pdd )
     call integral_p( 1, ko, pout, p_pdds, qz_pdd, &
          &           qz_gmean )


     ! ***** below not checked *****!
     ! for future reuse
     ! local____time differential
     
!     call derivative_p_nops(1, jm, ko, pout, p_dz_dt_zm, &
!          &                 divz_tzm)
     call derivative_p(1, jm, ko, pout, p_pds, p_dz_dt_zm, &
          &            divz_tzm)
     divz_tzm(:,:) = -divz_tzm(:,:) / 100.0 * grav
     
     call derivative_p(1, jm, ko, pout, p_pds, p_dphi_dt, &
          &            divphi_t)
!     divphi_t(:,:) = -divphi_t(:,:) / 100.0 * grav
     divphi_t(:,:) = -divphi_t(:,:) / 100.0

     do k=1, ko
        do j=1, jm
           uuv_tmp(j,k) = u_zm(j,k) * ( epz_uv(j,k) + epz_ut(j,k) ) &
                &         / radius / costbl(j)
        end do
     end do
     call derivative_p(1, jm, ko, pout, p_pds, uuv_tmp, &
          &            d_u_epz)
     d_u_epz(:,:) = -d_u_epz(:,:) / 100.0 * grav
     
     !    local energy budget about wave energy
     divz_tzm(:,:) = -divz_tzm(:,:)
     divphi_t(:,:) = divphi_t(:,:)
     dwdt(:,:) = divz_tzm(:,:) + divphi_t(:,:) + dkedt_uy(:,:) + d_u_epz(:,:)

     ! ***** above not checked *****!
    




     !********** output GrADS data **********

     ! lat-lev
     if( OUTPUT_ZONAL_FILENAME /= '' ) then
        call grads_write( ginfo_zonal, u_zm )
        call grads_write( ginfo_zonal, v_zm )
        call grads_write( ginfo_zonal, pt_zm )
        call grads_write( ginfo_zonal, t_dagger )
        call grads_write( ginfo_zonal, st_zm )
        
        call grads_write( ginfo_zonal, w_zm )
        call grads_write( ginfo_zonal, z_zm )
        call grads_write( ginfo_zonal, epy )
        call grads_write( ginfo_zonal, depy )
        call grads_write( ginfo_zonal, epz_form )
        
        call grads_write( ginfo_zonal, depz_form )
        call grads_write( ginfo_zonal, epz_uv )
        call grads_write( ginfo_zonal, depz_uv )
        call grads_write( ginfo_zonal, epz_ut )
        call grads_write( ginfo_zonal, depz_ut )
        
        call grads_write( ginfo_zonal, epz )
        call grads_write( ginfo_zonal, depz )
        call grads_write( ginfo_zonal, divf )
        call grads_write( ginfo_zonal, gy )
        call grads_write( ginfo_zonal, dgy )
                
        call grads_write( ginfo_zonal, gz )
        call grads_write( ginfo_zonal, dgz )
        call grads_write( ginfo_zonal, u_u_x_zm )
        call grads_write( ginfo_zonal, c_az_kz )
        call grads_write( ginfo_zonal, c_kz_ae_u )
        
        call grads_write( ginfo_zonal, c_kz_ae_v )
        call grads_write( ginfo_zonal, c_kz_ae )
        call grads_write( ginfo_zonal, c_ae_ke_u )
        call grads_write( ginfo_zonal, c_ae_ke_v )
        call grads_write( ginfo_zonal, c_ae_ke )
        
        call grads_write( ginfo_zonal, c_kz_ke_uy )
        call grads_write( ginfo_zonal, c_kz_ke_uz )
        call grads_write( ginfo_zonal, c_kz_ke_vy )
        call grads_write( ginfo_zonal, c_kz_ke_vz )
        call grads_write( ginfo_zonal, c_kz_ke_tan )
        
        call grads_write( ginfo_zonal, c_kz_ke )
        call grads_write( ginfo_zonal, c_kz_w )
        call grads_write( ginfo_zonal, q_zm )
        call grads_write( ginfo_zonal, qgz_zm )
        call grads_write( ginfo_zonal, qe_zm )
        
        call grads_write( ginfo_zonal, kz_zm )
        call grads_write( ginfo_zonal, ke_zm )
        call grads_write( ginfo_zonal, pz_zm )
        call grads_write( ginfo_zonal, ae_total_zm )
        call grads_write( ginfo_zonal, dkzdt_vkz )

        call grads_write( ginfo_zonal, dkzdt_wkz )
        call grads_write( ginfo_zonal, dkedt_uy )        
        call grads_write( ginfo_zonal, dkedt_vy )
        call grads_write( ginfo_zonal, dkedt_uz )
        call grads_write( ginfo_zonal, dkedt_vz )

        call grads_write( ginfo_zonal, dkedt_vke )
        call grads_write( ginfo_zonal, dkedt_wke )
        call grads_write( ginfo_zonal, dpedt_vt )  ! not used
        call grads_write( ginfo_zonal, d_u_epz )   ! not checked 
        call grads_write( ginfo_zonal, divz_tzm )  ! not checked

        call grads_write( ginfo_zonal, divphi_t )  ! not checked
        call grads_write( ginfo_zonal, dwdt )      ! not checked
     end if
     
     ! lat
     if( OUTPUT_VINT_FILENAME /= '' ) then
        call integral_p( jm, ko, pout, p_pds, kz_zm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, ke_zm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call grads_write( ginfo_vint, az_zm_vint )
        
        call grads_write( ginfo_vint, ae_zm_vint )

        call integral_p( jm, ko, pout, p_pds, c_az_kz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ae_u, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ae_v, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ae, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_ae_ke_u, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_ae_ke_v, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_ae_ke, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ke_uy, temp_vint )
        call grads_write( ginfo_vint, temp_vint )

        call integral_p( jm, ko, pout, p_pds, c_kz_ke_uz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ke_vy, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ke_vz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ke_tan, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_ke, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, c_kz_w, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, q_zm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, qgz_zm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, qe_zm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkzdt_vkz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )

        call integral_p( jm, ko, pout, p_pds, dkzdt_wkz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )

        call integral_p( jm, ko, pout, p_pds, dkedt_uy, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkedt_vy, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkedt_uz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkedt_vz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkedt_vke, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dkedt_wke, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dpedt_vt, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, d_u_epz, temp_vint )
        call grads_write( ginfo_vint, temp_vint )

        call integral_p( jm, ko, pout, p_pds, divz_tzm, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, divphi_t, temp_vint )
        call grads_write( ginfo_vint, temp_vint )
        
        call integral_p( jm, ko, pout, p_pds, dwdt, temp_vint )
        call grads_write( ginfo_vint, temp_vint )

     end if

     ! global mean
     if( OUTPUT_GMEAN_FILENAME /= '' ) then
        call grads_write( ginfo_gmean, az_gmean )
        call grads_write( ginfo_gmean, qz_gmean )
     end if

     ! wavenumber decomposition
     if( OUTPUT_WAVE_FILENAME /= '' ) then
        do w=1, wmax
           call grads_write( ginfo_wave, epz_wave(w,1:jm,1:ko) )
        end do
     end if


     ! save *_past
     pt_past(:,:,:) = pt(:,:,:)
     z_pd_past(:,:,:) = z_pd(:,:,:)
     p_pd_past(:,:,:) = p_pd(:,:,:)
     phi_dagger_past(:,:) = phi_dagger(:,:)
     
  end do

  ! deallocate memory
  call mim_var_end()
  call com_var_end()

end program mim
