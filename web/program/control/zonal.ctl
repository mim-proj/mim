DSET ^zonal2007.grd
UNDEF  -9.99E33
OPTIONS LITTLE_ENDIAN YREV
XDEF     1  LINEAR  0.0  2.5 
YDEF    73  LINEAR  -90.0 2.5
ZDEF    17 LEVELS
1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10
TDEF  1460    LINEAR 00z01jan2007 6hr
vars 57
u           17 99  Zonal Wind [m/s]
v           17 99  Meridional Wind [m/s]
pt          17 99  Potential Temperature [K]
t           17 99  Temperature Dagger derived from pt [K]
st          17 99  Mass Streamfunction [kg/s]
w	    17 99  Vertical Velocity derived from st [m/s]
z           17 99  Geopotential Height [m]
epy	    17 99  Meridional Component of EP Flux [kg/s^2]
depy        17 99  EP Flux Divergence dut to epy [m/s^2]
epz_form    17 99  Vertical Component of EP Flux (Form drag) [kg/s^2]
depz_form   17 99  EP Flux Divergence dut to epz_form [m/s^2]
epz_w       17 99  Vertical Component of EP Flux (W) [kg/s^2]
depz_w      17 99  EP Flux Divergence dut to epz_w [m/s^2]
epz_ut      17 99  Vertical Component of EP Flux (u'T') [kg/s^2]
depz_ut     17 99  EP Flux Divergence dut to epz_ut [m/s^2]
epz         17 99  Vertical Component of EP Flux (Total) [kg/s^2]
depz        17 99  EP Flux Divergence dut to epz [m/s^2]
divf        17 99  EP Flux Divergence [m/s^2]
gy          17 99  Meridional Comp. of G (Meridional Momentum) Flux [kg/s^2]
dgy         17 99  G Flux Divergence due to gy [m/s^2]
gz          17 99  Vertical Component of G (Meridional Momentum) Flux [kg/s^2]
dgz         17 99  G Flux Divergence due to gz [m/s^2]
uux         17 99  Zonal mean (u'u') [m^2/s^2]
c_az_kz     17 99  C(Az->Kz) [J/(kg s)] or [W/kg]
c_kz_ae_u   17 99  C(Kz->Ae) by U (i.e. form drag) [J/(kg s)] or [W/kg]
c_kz_ae_v   17 99  C(Kz->Ae) by V [J/(kg s)] or [W/kg]
c_kz_ae     17 99  C(Kz->Ae) [J/(kg s)] or [W/kg]
c_ae_ke_u   17 99  C(Ae->Ke) by U [J/(kg s)] or [W/kg]
c_ae_ke_v   17 99  C(Ae->Ke) by V [J/(kg s)] or [W/kg]
c_ae_ke     17 99  C(Ae->Ke) [J/(kg s)] or [W/kg]
c_kz_ke_uy  17 99  C(Kz->Ke) by u * DFy [J/(kg s)] or [W/kg]
c_kz_ke_uz  17 99  C(Kz->Ke) by u * DFz^uw [J/(kg s)] or [W/kg]
c_kz_ke_vy  17 99  C(Kz->Ke) by v * DGy [J/(kg s)] or [W/kg]
c_kz_ke_vz  17 99  C(Kz->Ke) by v * DGz [J/(kg s)] or [W/kg]
c_kz_ke_tan 17 99  C(Kz->Ke) by tan [J/(kg s)] or [W/kg]
c_kz_ke     17 99  C(Kz->Ke) [J/(kg s)] or [W/kg]
c_kz_w      17 99  C(Kz->W) = C(Kz->Ke) + C(Kz->Ae) [J/(kg s)] or [W/kg]
q           17 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)]
qgz         17 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)]
qe          17 99  Eddy Available Diabatic Heating [J/(kg s)]
kz          17 99  Zonal Kinetic Energy [m^2/s^2] or [J/kg]
ke          17 99  Eddy Kinetic Energy [m^2/s^2] or [J/kg]
pz          17 99  Potential Energy (including groud state) [m^2/s^2] or [J/kg]
ae_total    17 99  Eddy Available Potential Energy (including surface effect)
dkzdt_vkz   17 99  Kz advection by v [m/s^3]
dkzdt_wkz   17 99  Kz advection by w [m/s^3]
dkedt_uy    17 99  Wave Energy Flux Div. d(u Fy)/dy [m/s^3]
dkedt_vy    17 99  Wave Energy Flux Div. d(v Gy)/dy [m/s^3]
dkedt_uz    17 99  Wave Energy Flux Div. d(u Fz^uw)/dz [m/s^3]
dkedt_vz    17 99  Wave Energy Flux Div. d(v Gz)/dz [m/s^3]
dkedt_vke   17 99  Ke advection by v [m/s^3]
dkedt_wke   17 99  Ke advection by w [m/s^3]
dpedt_vt    17 99  [m/s^3] (not used)
duepz       17 99  div_uepz(local) (not checked)
dzdtm       17 99  divz_tzm(local) (not checked)
dzdt        17 99  divphi_t(local) (not checked)
dwdt        17 99  dwdt(local)     (not checked)
ENDVARS
