DSET ^vint2007.grd
UNDEF  -9.99E33
OPTIONS LITTLE_ENDIAN YREV
XDEF     1  LINEAR  0.0  2.5 
YDEF    73  LINEAR  -90.0 2.5
ZDEF     1  LEVELS  1000 
TDEF  1460  LINEAR  00z01jan2007 6HR
vars 34
kz           1 99  VI (=Vertically Integrated) Zonal Kinetic Energy [J/m^2]
ke           1 99  VI Eddy Kinetic Energy [J/m^2]
az           1 99  VI Zonal Available Potential Energy (not for gmean) [J/m^2]
ae           1 99  VI Eddy Available Potential Energy [J/m^2]
c_az_kz      1 99  VI C(Az->Kz) [W/m^2]
c_az_pe_u    1 99  VI C(Kz->Ae) by U (i.e. form drag) [W/m^2]
c_kz_ae_v    1 99  VI C(Kz->Ae) by V [W/m^2]
c_kz_ae      1 99  VI C(Kz->Ae) [W/m^2]
c_ae_ke_u    1 99  VI C(Ae->Ke) by U [W/m^2]
c_ae_ke_v    1 99  VI C(Ae->Ke) by V [W/m^2]
c_ae_ke      1 99  VI C(Ae->Ke) [W/m^2]
c_kz_ke_uy   1 99  VI C(Kz->Ke) by u * DFy [W/m^2]
c_kz_ke_uz   1 99  VI C(Kz->Ke) by u * DFz [W/m^2]
c_kz_ke_vy   1 99  VI C(Kz->Ke) by u * DGy [W/m^2]
c_kz_ke_vz   1 99  VI C(Kz->Ke) by u * DGz [W/m^2]
c_kz_ke_tan  1 99  VI C(Kz->Ke) by tan [W/m^2]
c_kz_ke      1 99  VI C(Kz->Ke) [W/m^2]
c_kz_w       1 99  VI C(Kz->W) = C(Kz->Ke) + C(Kz->Ae) [W/m^2]
q            1 99  Diabatic Heating [W/m^2]
qgz          1 99  Zonal Diabatic Heating (+Ground State) [W/m^2]
qe           1 99  Eddy Available Diabatic Heating [W/m^2]
dkzdt_vkz    1 99  VI Kz advection by v [kg/(m s^3)]
dkzdt_wkz    1 99  VI Kz advection by w+ [kg/(m s^3)]
dkedt_uy     1 99  VI Wave Energy Flux Div. d(u Fy)/dy [kg/(m s^3)]
dkedt_vy     1 99  VI Wave Energy Flux Div. d(v Gy)/dy [kg/(m s^3)]
dkedt_uz     1 99  VI Wave Energy Flux Div. d(u Fz^uw)/dz [kg/(m s^3)]
dkedt_vz     1 99  VI Wave Energy Flux Div. d(v Gz)/dz [kg/(m s^3)]
dkedt_vke    1 99  VI Ke advection by v [kg/(m s^3)]
dkedt_wke    1 99  VI Ke advection by w+ [kg/(m s^3)]
dpedt_vt     1 99  VI  [kg/(m s^3)] (not used)
duepz        1 99  duepz (not checked)
dzdtm        1 99  divz_tzm (not checked)
dzdt         1 99  divphi_t (not checked)
dwdt         1 99  dwdt (not checked)
ENDVARS
