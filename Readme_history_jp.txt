��TODO���X�g
  ������
  Pz�̈ܓx���z�A�e�[���[�W�J�̕K�v����(Pz�̑S�����ς͂��̂܂܂̕���better)
  wave�𓝍�(aqua��ana�Ŏg���Ă�����)


��mochi�łƒl�������Ⴄ���R(���ɒn�ʕt�߂ő�B��ʂɂ͈Ⴂ�͖��ɂȂ�Ȃ����x)
    biseki_bibun() (��bibun()) �� ncalc=1(kodama), ncalc=2(mochi)
      c_kz_ke��4%���x����
    pt_dot��omega���g��(kodama)���g��Ȃ�(mochi)��
      Q���g���Ȃ��ꍇ:  epz_ut(kodama)  <-- �Ή� --> epzut(mochi)
                        c_kz_ke(kodama) <-- �Ή� --> kzke(mochi)
      Q���g����ꍇ:    epz_ut(kodama)  <-- �Ή� --> epzut2(mochi)
                        c_kz_ke(kodama) <-- �Ή� --> (kzkeuw2)(mochi)
    pe_1(kodama) <-- �Ή� --> pe11(mochi)
    mochi�ł�pe1�͌o�x���z������̂ɗL��

�����̃v���O�����̌n��(����)
Iwasaki �� Shimizugichi �� Tanaka �� Ujiie �� Miyazaki �� Kodama(0.10�n)
                                     ��                        ��
                                     ���� Uno �� Mochi �� Kodama(0.20�n)

���X�V����(�V�������̂���)
2020.10.23
  p+�ʒT���̔����v�Z�������B
  ��getpt_ptiter()��if����ł́A
  pout(k) > p_zm(k+1)�̂Ƃ��ɒn��f�[�^���g���Ă������߁A
  �C�^���[�V���������邲�Ƃ�p_zm���ڕW�lpout���牓�����邱�Ƃ��������B
  ����͓���CMIP5�̌v�Z�Ō����ł������B
  �����ŁAget_ptiter()�̃A���S���Y���̑S�ʓI�ȉ������s�����B
  �v�Z���O��p_zm��pout�̍����m�F����ƁA�V��@�ł͒n�\�t�߂݂̂ō�����〈����B
  �܂��Agetpt()���ŁA�t�]�w���o�̃T�u���[�`�����g���Ă��Ȃ������̂ŁA
  �g�����ƂɕύX�B
  p_zm�̋t�]�����o����getp_stable()���ǉ������B

2016.12.21
  omega�ǂݍ��݂̂Ƃ���zrev�ɑΉ��B

2015.08.19
  namelist��INPUT_ZDEV_ZREV��ǉ��B

2009.03.27
  Version 0.33 Release 1

2009.03.19
  getpt_p()��p_pd���}�ŁAp���`���}��log-p���`���}�ɕύX�B
    ����ɂ��A�n�\�ʂł̋C���̒l�̖������قډ����B�S�ϐ���1%���x�̂���B

2009.03.10
  ��f�M���Mq_3d�̓��͂��Ȃ��ꍇ�A�]����0�ɂ��Ă������AD(pt)/Dt��
    �S��������ߎ�����悤�ɕύX�B�덷���傫���̂Œ��ӁB
    �܂��A�ŏ��̃^�C���X�e�b�v�ł͕Δ���d(pt)/dt=0�Ƃ����̂����ӁB

2009.02.24
  Version 0.32 Release 2
  namelist��INPUT_UNIT_Z��ǉ��B

2009.02.13
  namelist_check()��OUTPUT_WAVE_FILENAME�Ɋւ���error�̊��ύX�B

2009.02.12
  Version 0.32 Release 1
  �R���p�C���I�v�V������-heap-array��ǉ��B����ɂ��A���𑜓x�f�[�^
    (e.g. JRA125)�������ۂ�stack�s����segmentation fault���N����̂�}�~�B
    �������Aifort version 10.0�ȍ~�̂ݑΉ��B
  get_phi_dagger(), mount_modify() �ŁApout(ko+1)�ɃA�N�Z�X���Ă��܂��o�O���C��
    �A���A���ʂւ̉e���͂��蓾�Ȃ��B
  �o�O: biseki_biseki() �� x_pt(im,jm,km) �� x_pt(im,jm,ko)�ɏC��
    �A���Akm=ko�̏ꍇ�͉e���Ȃ��B

2008.09.19
  Gy�̒�`��ύX�B
    �]��: -rho_0 (v'^2)_zm
    ����: -rho_0 a (v'^2)_zm

2008.09.07
  wave�v�Z��T63�f�[�^�Ŏ��s����s��������B
    �������A�s��̌����ȗ��R�͕s��(where���̐���?)

2008.09.01
  Version 0.31 Release 1
  �\�[�X�R�[�h�̉p�ꉻ����(fft.f90������)�B

2008.08.27
  mountain_modify()�̍Ō��c_pz_kz2�Z�o�����ŁA
    k=1�̂Ƃ���pout_new2���g���悤�ɕύX�B
    ����͒P�Ȃ�v���O�����̕�����₷���̖��ł����āA�{���I�ȕύX�ł͂Ȃ��B
    �Ⴂ�͂��������킸���B
  integral_pt_ym()��V�݁A�ϕ���Ԃ̏�[���s�K�؂ł������̂��C���B
    �������A���ʂւ̉e���͂����킸��(az_vint�ɉe��:0.01%���x)�B
  check_range()�ŁA���񐔈ȏ�warning�����m���ꂽ�痎����悤�Ɏd�l�ύX�B

2008.08.26
  get_phi_zm.f90��get_phi_dagger.f90�ɉ����B
  get_phi_dagger()�ŁAp+�n�ɂ�����pt�̒n�\�ʒl��
    �P�Ȃ�zonal mean�ɂȂ��Ă����o�O���C���B
    C(AZ->KZ)�ɉe���B�A���S�����ϒl��0.1%���x�Ƃ����킸���B

2008.08.24
  intpl()��p�~�B
  getpt()�܂��̕ϐ��𐮗��B
  integral_pt()�ŁA�ϕ���Ԃ̏�[���s�K�؂ł������̂��C���B
    �������A���ʂւ̉e���͂����킸��(ae_vint�ɉe��:0.01%���x)�B
    <-����1%���x�̊�^???

2008.08.21
  biseki_sekibun()�̓��}���ŁAd<0�ɂȂ����Ƃ��s�K�؂Ȓl�ɂȂ�o�O���C���B
    �������A���ʂւ̉e���͂����킸���B
  biseki_y.f90�� p_pt, x_pt �� p_pdd, x_pdd �ɕύX�B

2008.08.20
  Version 0.30 Release 1
  PS(�n�\�ʋC��)�̑��� MSL(�C�ʊr���C��), TS(�n�\�C��)�̓��͂ɑΉ��B
    �덷�͏�����(divf���͒n�\�t�߂ł��傾�����Ȃ����x)

2008.08.19
  biseki_sekibun()�̓��}���ŁAd<0�ɂȂ����Ƃ��s�K�؂Ȓl�ɂȂ�o�O���C���B
    �܂��Ax_pt(im,jm,km)��x_pt(im,jm,ko)�ɏC���B
    �������A���ʂւ̉e���͏�[�t�߂����S�ł���A�����킸���B
  biseki.f90�� p_pt, ps, x_pt �� p_pd, p_sfc, x_pd �ɕύX�B

2008.08.11
  get_pt_dot()��get_pt_dot_omega()�ɉ����Aget_pt_dot_q()��V�݁B
  pd_pt �� pd_pdd �ɕύX�B
  p_pt, z_pt �� p_pd, z_pd �ɕύX�B
  lagmain*.f90 �� mim*.f90 �ɕύX�B

2008.08.04
  Version 0.29 Release 2
  �h�L�������g�𐮗��B

2008.08.01
  Version 0.29 Release 1
  pt_dot�̌v�Z��Q���g���悤�ɕύX(Q���g�p�\�Ȃ�)�B

2008.07.12
  Version 0.28 Release 1

2008.07.10
  namelist��OUTPUT_DEBUG_FILENAME*��p�~�B
  namelist��WAVE_NUMBER=0�������悤�ɕύX�B
  namelist��INPUT_ZDEF_LEVEL, OUTPUT_ZDEF_LEVEL���A�ォ�牺�A�������A
    �ǂ���Ŏw�肵�Ă��󂯕t����悤�ɕύX�B
  namelist��INPUT_YDEF_SOUTH, INPUT_YDEF_NORTH��ǉ��B
  check_range()�Œl���O�ꂽ�Ƃ��A���s���~�������s����悤�ɕύX�B

2008.07.09
  zonal.ctl�̕ϐ����̃~�X(dgz�̂͂���dgy�ɂȂ��Ă���)���C���B

2008.07.04
  parameter.f90��econv, divf�̒l�͈̔͂��A���ꂼ��10�{�Ɋɘa�B

2008.06.02
  Version 0.27 Release 3
  parameter.f90��z_min��-1000����-3000�ɕύX�B
  �����getpt_lev()��ppp�v�Z������log(p)���}�ɖ߂��B
    z_zm�����ُ�l���݁B

2008.06.01
  Version 0.27 Release 2
  getpt_lev()��dlev�v�Z������log(p)���}�ɖ߂��B
    ���^����depz_w, dept_ut�v�Z���ɔ��U���邱�Ƃ�����(�����s��)

2008.05.30
  Version 0.27 Release 1
  dkzdt_vkz, dkzdt_wkz�̌v�Z��ǉ��B

2008.05.29
  Ke�̌v�Z�̏�����ύX�A�l�������ɂ����B
  dkedt_vke, dkedt_wke�̌v�Z��ǉ��B

2008.05.27
  Namelist��endian�𓱓��B
  get_z_pt()�̌v�Z������ύX�B���ʂɖ����ł�����x�̂���B

2008.05.26
  Version 0.26 Release 1
  Qz���C���B
  ���}����log(p)��p����̂�������߂�B
    log�v�Z�̌덷���l����Ƒf���Ȑ��^���}�œ��ꂵ�������悳�����B
    ������p���̂��̂���}����ꍇ�͏�[�ŕ��ɂȂ�\��������̂ŁA
    ��[�̂�log���g�p����ׂ��B
  pd_p�̓��}���@��log(p)�^�ɕύX�A�����ώZQE�̐�%�̂���B

2008.05.13
  pd_p (�W��p�ʂł�p+) �̌v�Z��ǉ�
  Ae�̉����ώZ�l���C���A�S�����ς�1%���x�̂���B

2008.05.12
  vint.ctl����ae_s, pe_total���폜�B
  Az�̉����ώZ����ymin�ɕύX�B�ϕ��͈͂̕s��v����������␳�����ǉ��B
    �����������ώZ�A�S�����ςւ̊�^�͋ɂ߂ď������B
  gmean.ctl��Az��ǉ��B

2008.04.30
  diabatic_aqz()��diabatic_qz_gmean()�ɉ����B

2008.04.29
  energy_ae()��energy_ae_total()�ɉ����B

2008.04.21
  P��A�̋�ʂ��������B(�\�[�X�R�[�h�A�R���g���[���t�@�C��)
  aqe��qe�ɉ����B
  qz_zm��qgz_zm�ɉ����AA_Z��P_Z�̐�������ʂ��邽�߁B

2008.03.28
  Version 0.25 Release 1

2008.03.27
  pt(p++s)�̒�`�̊ԈႢ���C���BPz���ꏊ�ɂ����1�������B
  aqz�̑S�����ς��o�́B
  getpt_*()��getpt_etc.f90�ɓ����A�����B
  getpt1_y()��getpt1()�ɓ����B
  ptiter.f90�̓��}���@���C���B

2008.03.25
  zonal.ctl�� G Flux ��ǉ��B
  EP Flux ���U�̒P�ʂ�m/s^2�ɕύX�B
  Namelist�̑啝�����B
  �\�[�X�R�[�h��src�ȉ��ֈړ��AMakefile�̑啝���������B

2008.01.20
  Version 0.24 Release 2
  �`���R�̔g���W�J�̌v�Z���ł��Ȃ�bug������
  Version 0.24 Release 1

2008.01.12
  pe_total_zm���v�Z�B����͕����ϕ����s���O��Ae�ł���B
  integral_pt()�̉��[���A��s�̓������ς������̂���min�ɕύX�B
  ����ɂ��Ape�̑S�����ϒl��2�����x�����B

2008.01.08
  integral_pt()��V�݁B
  merifional_integral.f90��integral.f90�ɉ����Azsekibun.f90���z��

2008.01.07
  lagmain_var.f90�̕ϐ��𐮗��B
  phi_zm��phi_dagger�ɉ����B���Ԃɍ�����B
  namelist��FILE_DEBUG1, FILE_DEBUG2��ǉ��B

2007.12.25
  bibun(), sekibun() �����W���[��biseki�Ƃ��ē����B
  epz, gz ��div�����߂�X�L�[����derivative_z()�ɔC����B
    ����ɂ��A���ʂ����������B
  derivative_z()��V�݁B����� half level ��p����z�����B

2007.12.24
  div_epflux.f90 �� epflux_div.f90 �ɉ����B
  epflux_*.f90 �� epflux.f90 �ɓ����B
  get_var_w_x_zm.f90 ��V�݁Aepflux_z_uw.f90 �̖w�ǂ̏������ڍs�B
  gflux.f90��V�݁A�ϐ�gy, dgy�Ȃǂ𓱓��B
  energy_conv_*()�̒l�`�F�b�N��������B
  energy_tendency.f90��V�݁Au_epflux.f90���z���B
  epy�v�Z�ɕϐ�rho���g�p����悤�ɂ���B
  c_kz_ke�̍��̕���������ύX�Bgrads�o�͕ϐ����ύX�B
  c_pe_kz* �� c_kz_pe*�ɕύX�B
  phi_zm�̒�`���Ageopotential height ���� geopotential �ɕύX(�{���̒�`)�B
  �ϐ� epzw or epz_w �� epz_uv�ɁAepzut �� epz_ut �ɕύX���Aepz_uw��V�݁B
  depz_w�������l�B
  namelist��FILE_TOPO_KIND�̎w��̎d����ύX�B
    geopotential�̏ꍇ�A"g"�Ǝw��ł���悤�ɂȂ�B
    �ߋ��o�[�W�����Ƃ̌݊�����ۂ��߁A"gh"���c���B

2007.12.17
  Version 0.23 Release 1
  JRA25����m�F�B
  namelist��FILE_U, FILE_V, FILE_T ��ǉ��B

2007.12.01
  com_var module �� tantbl ��ǉ��B
  parameter module �� sec_day ��ǉ��B
  mountain_modify.f90 ����s�v�ȃR�[�h���폜�B��1000�s����300�s�B
  mountain_modify.f90 �ŁA���}�v�Z�̃o�O���C��
  1hPa�ȉ���c_pz_kz��0�ɂȂ�ꍇ���������B���݉����ρB
  �����ώZ�l����c_pz_kz�ւ̉e����1%���x�B
  energy_conv.f90��V�݁B�G�l���M�[�ϊ��̃R�[�h���ړ��B
  div_z()��p���Ă���divz_tzm, divphi_t, d_u_epz��derivative_p()�ɕύX�B
  �n�\�ʂ��l���Ă���̂ŁAderivative_p_nops()���������͂��B
  �n�\�t�߂ł����B�����ώZ�������(10%���x����)
  �l�̊m�F�͂��Ă��Ȃ�(���Xno checked)
  div_z.f90��p�~���Aderivative_p_nops()��V�݁B���҂̈Ⴂ�͒P�ʂƏd�͂̂݁B
  bibun()��ncalc=1�����ɁAderivative_p()��V�݁Bbibun()�̎����������ړ��B
  phi_t �� p_dphi_dt* �ɉ���
  z_t* �� p_dz_dt* �ɉ����B
  get_phi_zm.f90 : com_var��use�B
  *_past(:,:)��������"�ϐ���_past(:,:)"�ƂȂ�悤�ɉ����B
  pt_before(:,:) �� pt_past(:,:)�ɉ����B
  pt_past(:,:)��p�~�B
  get_pt_dot.f90 : com_var��use�B
  �n�\�t�߂�p_pt�̒l���A�v���O�����O���ƌ㔼�ňقȂ�s�����������B
  �]����epflux_z_form_z_pt()�Œu��������Ă����̂�p�~�B
  p_pt�̒n���l��ps�Ɉ�v�����Ȃ��B�K�v�Ȃ炻�̓s�x���(e.g. form drag)
  ����ɂ��e���͒n�\�ʕt�߂ōő�1%���x�A�����ώZ�l�ł�1%���x�B
  ������ɂ���A�n�ʕt�߂̒l�͌���ׂ��ł͂Ȃ��B

2007.11.30
  epflux_z_uw.f90 : com_var��use�B
  epflux_z_form.f90 : com_var��use�B
  epflux_z.f90 �� epflux_z_form.f90 �ɉ���
  epflux_y() : com_var��use�B
  [�d�v] divf�̒l��Fy��Fz��form frag�����l������Ă��Ȃ������o�O�������B
    ���g���_����divF����������Ă���\������B����OK�B�v�����B
  div_epflux_y(), div_epflux_z() : com_var��use�B
  depy.f90 �� depz.f90 �� div_epflux.f90 �ɓ����B
  bibun() : com_var��use�B
  sekibun() : com_var��use�B
  setpt0() : com_var��use�B
  namelist: WAVE_OUTPUT ��FILE_WAVE�ɕύX�AWAVE_SWITCH��p�~�B
  FILE_P_Z�̏o�͂��x�~
  grads.f90�����W���[�����Bgrads_open()�Agrads_write()����lagmain.f90�ɓ����B
  ���W���[�����A�T�u���[�`�����ƃt�@�C������Ή�������

2007.11.29
  com_var.f90��V�݁A�v���O�����S�̂Ŏg�p���鋤�L�ϐ����`����B
  var.f90 -> lagmain_var.f90�ɉ����Alagmain.f90�݂̂�use���邱�Ƃ������B

2007.11.26
  Version 0.22 Release 2
  check_range()���쐬�A�l�̃`�F�b�N�Ɏg�p�B
  param.f90�̒萔��parameter������t����
  getlev()��dlev���}������ύX�B
    �S��log(p)���}�ł������̂��f�[�^Top�t�߂̂ݐ��`���}�ɂ���B
    ����ɂ��Adlev���������Ȃ肷����p_pt��0�ɂȂ�s�(�ő��ɋN���Ȃ���)��
    ����B

2007.11.02
  grads()��p�~(term_main.f90�͖��Ή�)
  grads_read_2d()�Agrads_read_3d()��p�~�Agrads_read()�ɒu������
  �J�E���g���̉�ʕ\����ύX
  nlist.f90��use��lagmain.f90�ȊO�ɂ�����
    calendar.f90��use

2007.10.22
  Version 0.22 Release 1
  qz, qe�̏o�͂�ǉ�
  namelist��ZONAL_Q��p�~�Aq_zm�̒P�ʕύX�Azonal.ctl�ɑg����
  c_kz_ke_uv��2�ɕ�����

2007.10.19
  dkedt_uvu_y(Ke�̐�������(u'v') d/dy(u_zm)�̍�)��ǉ�
  zonal.ctl��uepy��duepy�ɖ��̕ύX�A�o�͏��Ԃ̕ύX

2007.10.18
  Version 0.21 Release 1
  epflux_z_wave()�ɐVfft()���̗p
  fft()�̃f�[�^����C�ӂɐݒ�ł���悤�ɉ���
    ��fft()��fft_quick()�Ɩ��̕ύX
    �v�Z��@�̓}�j���A�����Q��

2007.10.09
  Version 0.20 Release 3
  bibun()�́`0��������
    ���ɃG�l���M�[�v�Z�ɉe���B���܂ł��܂������Ă������R�͕s���B
    ���g������cntl��2002.12.16�Ȃǂɔ���

2007.10.03
  Version 0.20 Release 2

2007.10.02
  mount_modify()�𔽉f
  getlev()��dlev��p����log(p)�Ō��肷��悤�ɕύX

2007.10.01
  getpt.f90�̔�����: 6 -> 8
  �n�`�f�[�^��yrev�L�����AFILE_TOPO_YREV�ł͂Ȃ�YDEF_YREV�Ŕ��f�����o�O������

2007.09.17
  Version 0.20 Release 1
  TERM, PHASE�ւ̑Ή�����U�ł��؂�
  ifort�ł̃R���p�C���ɐ���
  ptiter()�� 0/0 �ƂȂ�bug���C��(fujitsu f90�ł�0�ɂȂ邽�ߖ��Ȃ�)
  div_y.f90��u_epflux.f90�ɉ���
  2�����ւ̖����K��(e.g. u_d_v_d_zm)�𐧒�A�����̕ϐ����ɓK�p
  zsekibun()�𐮗�
  �ϐ����̕ύX�ƕ��בւ�
    epz -> epz_form
    depz -> depz_form
    ep_z -> epz (�v����)
    �G�l���M�[�ϊ� -> c_*
  zonal.ctl��t��ǉ�
  lag.ctl��zonal.ctl�ɉ���

2007.09.15
  mount_modify.f90��mochi�ł��瓱���B����������ł͕s�g�p�B
  energy_az()��ǉ�
  energy_ae()��ǉ�

2007.09.12
  epz_uw()�̈�������z_zm���폜
  c(pz,kz)�̌v�Z��derivative_y()���g���悤�ɕύX
  derivative_x()��p_dagger����ڐA�A�ǉ��Bdz_dlon�̌v�Z�ɗ��p
  derivative_y()��p_dagger����ڐA�A�ǉ��Bdz_dlat�̌v�Z�ɗ��p
  get_phi_zm()��ǉ��Alagmain���珈�����ړ�
  lagmain��z_zm��modify��get_z_pt()�Ɉړ�

2007.09.09
  get_pt_dot()��ǉ��Alagmain���珈�����ړ����A�v�Z���x���グ��(pt_dot�̌v�Z�@)
  namelist��TIMESTEP_DT��ǉ�
  array()��p�~�A�z�񎮂ɒu������
  sekibun(), bibun()�𐮌`
  getpt_y()��ǉ�( getpt()�̓�k�� )
  getlev()�𐮌`
  meridional_integral()��ǉ�

2007.09.08
  getpt()��3������
  stable()��3������
  getpzm()��3������
  intpl()��3�������A��ʂ�p->�Ɠ��}�ɂ��g�p�ł���悤�ɂ���
  getpt1()��3������
  zmean()��p�~���A�g�ݍ��݊֐�sum()�Œu������
  getpt()����call���Ă���setpt0()��lagmain�ֈړ�
  setpt0()��3������
  lagmain�� call getpt() �t�߂̏������ȗ���
  ps_zm�̌v�Z��g�ݍ��݊֐�sum()��p���Ċȗ����A�ꎞ�ϐ�pss���폜
    �v�Z�̏��Ԃ��ς�邽�߁A�ꕔ�̕ϐ��̌��ʂ������ɂ����
  namelist��FILE_Q_UNDEF��ǉ�
  undef.f90��ǉ��Alagmain����undef������undef_fill()�Ɉړ�
  omega�t�@�C����open�ʒu��ύX
  version 0.20�n�̊J���J�n(mochi��(ene-sam)�Ƃ̓���)

2007.04.08
  version 0.13 release 3
  Pe�̐ϕ��l�̏o�͂�ǉ�

2007.04.06
  version 0.13 release 2

2007.04.05
  namelist��FILE_VINT��ǉ�

2007.03.16
  version 0.13 release 1
  namelist��ZDEF_ZDEF_IN_OMEGA�𓱓�

2007.02.06
  version 0.13 beta 2
  getpt.f90�̔�����: 5 -> 6

2006.11.13
  version 0.13 beta 1
  Pe���o�͂���悤�ɂ���(���Ȃ������)
  Pz���C��(�܂�������)

2006.09.21
  version 0.12 release 1
  namelist��FILE_TOPO_YREV��ǉ�(���`�F�b�N)

2006.03.30
  version 0.11 release 3
  divepz�̃o�O�C��(z -> z+)

2006.03.14
  version 0.11 release 2
  FILE_ZONAL_Q���ȗ��\�ɂ���D
  version 0.11 release 1
  namelist�̒l�`�F�b�N������
  OMEGA���ȗ��\�ɂ���
  namelist��UNDEF��ǉ�

2006.03.01
  FILE_Q���w�莞�̃��b�Z�[�W������

2006.01.04
  getpt()�ŁApout(:)��intent��in�ł͂Ȃ�out�ɂȂ��Ă����o�O���C��
  bibun()�̈����ύX(name��p�~)�ɔ����ύX
  sekibun()�̈����ύX(name��p�~)�ɔ����ύX

2005.10.02
  namelist.f90��V��

2005.09.28
  epz�̔g���W�J���[�`����ǉ�

2005.08.05
  YDEF_YDEF �� YDEF_YREV �̍����ɂ��o�O������

2005.07.12
  parameter�����W���[����

2005.04.27
  XDEF, YDEF��g�ݍ���

2005.04.26
  kzkevv�Ƀo�O -> �œK���̂����Ō��ʂ����ꂽ(�����ł��邭�炢������)
  allocatable����

2005.04.12
  namelist��ZDEF�O���[�v��ǉ�

2005.04.11
  �N����wp��������calendar.f90��ǉ�
  namelist�𓱓�, �S�Ă̓��o�̓t�@�C�������w��\�ɂ���
