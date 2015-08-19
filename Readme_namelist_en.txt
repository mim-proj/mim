MIM (namelist version) by C. Kodama


-4 bytes unformatted binary data (GrADS format) is acceptable.
-These explanation are based on GrADS control file style
-See src/namelist.f90 for more imformation, e.g. variable types.



##### input file setting #####

&INPUT : input filename setting
  INPUT_TYPE           : = "general" (4 bytes unformatted binary data, default)
                           "tohoku_ncep" (only for Tohoku Univ. legacy)

  INPUT_UVT_FILENAME   : U, V, T filename (default='')
  INPUT_U_FILENAME     : U filename 
                         (valid if INPUT_UVT_FILENAME is NOT specified)
  INPUT_V_FILENAME     : V filename
                         (valid if INPUT_UVT_FILENAME is NOT specified)
  INPUT_T_FILENAME     : T filename
                         (valid if INPUT_UVT_FILENAME is NOT specified)
  INPUT_PS_FILENAME    : PS filename (default='')
  INPUT_MSL_FILENAME   : MSL pressure filename
                         (valid if INPUT_PS_FILENAME is NOT specified)
  INPUT_TS_FILENAME    : surface temperature filename
                         (valid if INPUT_PS_FILENAME is NOT specified)
  INPUT_Z_FILENAME     : Z filename
  INPUT_OMEGA_FILENAME : OMEGA filename
                         (without specified, OMEGA is assumed to be 0)
  INPUT_TOPO_FILENAME  : topography filename
  INPUT_Q_FILENAME     : diabatic heat filename
                         (without specified, Q is assumed to be 0)


&INPUT_UNIT : unit setting
  INPUT_UNIT_Z        : = "m" (height, default)
                        = "m^2/s^2" (geopotential)
  INPUT_UNIT_PS       : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_MSL      : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_TOPO     : = "m" (height, default)
                        = "m^2/s^2" (geopotential)


&INPUT_UNDEF : undef setting
  INPUT_UNDEF_DEFAULT : common UNDEF value (default=9.999e+20)
  INPUT_UNDEF_UVT     : UVT file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_U       : U file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_V       : V file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_T       : T file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_PS      : PS file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_MSL     : MSL file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_TS      : TS file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Z       : Z file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_OMEGA   : OMEGA file UNDEF value (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Q       : Q file UNDEF value (default=INPUT_UNDEF_DEFAULT)


&INPUT_ENDIAN : endian setting
  INPUT_ENDIAN_DEFAULT : common endian setting
                         = "little" (little endian, default)
                         = "big"    (big endian)
  INPUT_ENDIAN_UVT   : UVT file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_U     : U file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_V     : V file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_T     : T file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_PS    : PS file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_MSL   : MSL file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TS    : TS file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Z     : Z file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_OMEGA : OMEGA file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Q     : Q file endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TOPO  : TOPO file endian (default=INPUT_ENDIAN_DEFAULT)


&INPUT_XDEF : X direction setting
              only for uniform grid
              output setting is same as input setting
  INPUT_XDEF_NUM  : number of grid point



&INPUT_YDEF : Y direction setting
              output setting is same as input setting
  INPUT_YDEF_TYPE   : level type
                      = "lat_degree"  (levels [degree], default)
                      = "lat_radian"  (levels [radian])
                      = "linear"      (uniform interval [degree])
  INPUT_YDEF_NUM    : number of grid point
  INPUT_YDEF_LEVEL  : latitude levels ( South->North or North->South )
                      valid if INPUT_YDEF_TYPE = "lat_degree"  or "lat_radian"
  INPUT_YDEF_SOUTH  : southern edge latitude [degree] (default=-90)
                      valid if INPUT_YDEF_TYPE = "linear" 
  INPUT_YDEF_NORTH  : northern edge latitude [degree] (default=90)
                      valid if INPUT_YDEF_TYPE = "linear" 
  INPUT_YDEF_YREV_DEFAULT
                    : common YREV setting
                       = 0 (NOT YREV i.e. south -> north, default)
                       = 1 (yrev i.e. north -> south)
  INPUT_YDEF_YREV_TOPO
                    : YREV setting for TOPO file
                       = 0 (NOT YREV i.e. south -> north)
                       = 1 (yrev i.e. north -> south)
                      without specified, INPUT_YDEF_YREV_TOPO equals 
                      INPUT_YDEF_YREV_DEFAULT.


&INPUT_ZDEF : Z direction setting
  INPUT_ZDEF_NUM       : number of grid point
  INPUT_ZDEF_NUM_OMEGA : number of OMEGA grid point (default=INPUT_ZDEF_NUM)
                         for NCEP/NCAR reanalysis, for example
  INPUT_ZDEF_LEVEL     : pressure levels ( Upper->Lower or Lower->Upper )


&INPUT_TDEF : time
  INPUT_TDEF_TYPE : time step type
                    = "tstep"   (directly specify time step)
                    = "monthly" (input file has 1 month time step)
                      "annual"  (input file has 1 year time step)
  INPUT_TDEF_DAYNUM : number of step per day
                      valid if INPUT_TDEF_TYPE = "monthly" or "annual"
  INPUT_TDEF_365DAY : 365 day/year switch
                      = 0 (default, leap year considered)
                      = 1 (fixed 365 day/year, i.e. no leap year)
                      valid if INPUT_TDEF_TYPE = "annual" or "monthly"
  INPUT_TDEF_YEAR   : year
                      valid if INPUT_TDEF_TYPE = "annual" or "monthly"
  INPUT_TDEF_MONTH  : month
                      valid if INPUT_TDEF_TYPE = "monthly"
  INPUT_TDEF_TSTEP  : number of time step
                      valid if INPUT_TDEF_TYPE = "tstep"



##### analysis #####

&WAVE : wavenumber decomposition of form-drag
  WAVE_MAX_NUMBER : maximum number of wavenumber
                    specify 0 if wavenumber decomposition is not needed
                    (default=0)


##### output file setting #####
-Output file setting is same as input file setting except below.

&OUTPUT : input filename setting
  OUTPUT_ZONAL_FILENAME  : zonal mean field (2D)
  OUTPUT_VINT_FILENAME   : vertically integrated value (optional)
  OUTPUT_GMEAN_FILENAME  : global mean value (optional)
  OUTPUT_WAVE_FILENAME   : wave decomposition of form drag (optional)
  OUTPUT_ERROR_FILENAME  : error log

&OUTPUT_ZDEF : Z direction setting
  OUTPUT_ZDEF_NUM     : number of grid point
  OUTPUT_ZDEF_LEVEL   : pressure dagger levels ( Upper->Lower or Lower->Upper )
