MIM namelist�� by ���� (���{)


-�f�[�^��4�o�C�g�����Ȃ��o�C�i���f�[�^ (GrADS�`��) �ɑΉ����Ă��܂��B
-������ GrADS �̃R���g���[���t�@�C���̃X�^�C���ɉ����Ă��܂�
-�e�X�̌^��m�肽���ꍇ�́Asrc/namelist.f90 ���Q�Ƃ��Ă��������B



##### ���̓t�@�C���Ɋւ���ݒ� #####

&INPUT : ���̓t�@�C���̎w��
  INPUT_TYPE           : = "general" (�����Ȃ�4�o�C�g�����f�[�^, default)
                           "tohoku_ncep" (�Ǝ��`����ncep1,ncep2�p)

  INPUT_UVT_FILENAME   : U, V, T �t�@�C���� (default='')
  INPUT_U_FILENAME     : U�t�@�C����
                         (INPUT_UVT_FILENAME���w�肳��Ă��Ȃ��ꍇ�ɗL��)
  INPUT_V_FILENAME     : V�t�@�C����
                         (INPUT_UVT_FILENAME���w�肳��Ă��Ȃ��ꍇ�ɗL��)
  INPUT_T_FILENAME     : T�t�@�C����
                         (INPUT_UVT_FILENAME���w�肳��Ă��Ȃ��ꍇ�ɗL��)
  INPUT_PS_FILENAME    : PS�̃t�@�C���� (default='')
  INPUT_MSL_FILENAME   : MSL�̃t�@�C����
                         (INPUT_PS_FILENAME���w�肳��Ă��Ȃ��ꍇ�ɗL��)
  INPUT_TS_FILENAME    : �n�\T�̃t�@�C����
                         (INPUT_PS_FILENAME���w�肳��Ă��Ȃ��ꍇ�ɗL��)
  INPUT_Z_FILENAME     : Z�t�@�C����
  INPUT_OMEGA_FILENAME : OMEGA�̃t�@�C����
                         (�ȗ��A�ȗ��� OMEGA=0 �Ƃ��Ď�舵��)
  INPUT_TOPO_FILENAME  : �n�`�t�@�C����
  INPUT_Q_FILENAME     : ��f�M���M�t�@�C����
                         (�ȗ��A�ȗ��� Q=0 �Ƃ��Ď�舵��)


&INPUT_UNIT : �P�ʂ̎w��
  INPUT_UNIT_Z        : = "m" (�n�`���x, default)
                        = "m^2/s^2" (�W�I�|�e���V����)
  INPUT_UNIT_PS       : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_MSL      : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_TOPO     : = "m" (�n�`���x, default)
                        = "m^2/s^2" (�W�I�|�e���V����)


&INPUT_UNDEF : undef�̎w��
  INPUT_UNDEF_DEFAULT : �e�t�@�C�����ʂ�UNDEF�l (default=9.999e+20)
  INPUT_UNDEF_UVT     : UVT�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_U       : U�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_V       : V�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_T       : T�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_PS      : PS�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_MSL     : MSL�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_TS      : TS�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Z       : Z�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_OMEGA   : OMEGA�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Q       : Q�t�@�C����UNDEF�l (default=INPUT_UNDEF_DEFAULT)


&INPUT_ENDIAN : endian�̐ݒ�
  INPUT_ENDIAN_DEFAULT : �e�t�@�C�����ʂ�endian
                         = "little" (little endian, default)
                         = "big"    (big endian)
  INPUT_ENDIAN_UVT   : UVT�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_U     : U�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_V     : V�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_T     : T�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_PS    : PS�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_MSL   : MSL�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TS    : TS�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Z     : Z�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_OMEGA : OMEGA�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Q     : Q�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TOPO  : TOPO�t�@�C����endian (default=INPUT_ENDIAN_DEFAULT)


&INPUT_XDEF : X�����̐ݒ�
              ��l�Ԋu�̂ݑΉ�
              output�ɂ����̂܂܌p�������
  INPUT_XDEF_NUM  : �i�q��



&INPUT_YDEF : Y�����̐ݒ�
              ��́Aoutput�ɂ����̂܂܌p�������
  INPUT_YDEF_TYPE   : level�̎w����@
                      = "lat_degree"   (�z��[�x], default)
                      = "lat_radian"   (�z��[radian])
                      = "linear"       (��l�Ԋu[�x])
  INPUT_YDEF_NUM    : �i�q��
  INPUT_YDEF_LEVEL  : �ܓx��\���z�� ( South->North or North->South )
                      INPUT_YDEF_TYPE �� "lat_degree" ���� "lat_radian" ��
                      �ꍇ�̂ݗL��
  INPUT_YDEF_SOUTH  : ��[�̈ܓx[��] (default=-90)
                      INPUT_YDEF_TYPE �� "linear" �̏ꍇ�̂ݗL��(���`�F�b�N)
  INPUT_YDEF_NORTH  : �k�[�̈ܓx[��] (default=90)
                      INPUT_YDEF_TYPE �� "linear" �̏ꍇ�̂ݗL��(���`�F�b�N)
  INPUT_YDEF_YREV_DEFAULT
                    : �e�t�@�C���̋��ʂ�YREV
                       = 0 (yrev�Ȃ� i.e. �쁨�k, default)
                       = 1 (yrev���� i.e. �k����)
  INPUT_YDEF_YREV_TOPO
                    : �n�`�t�@�C����YREV(�f�[�^�z��̌���)
                      = 0 (yrev�Ȃ� i.e. �쁨�k)
                      = 1 (yrev���� i.e. �k����)
                      �w�肪�Ȃ��ꍇ�AINPUT_YDEF_YREV_DEFAULT�ɏ]��


&INPUT_ZDEF : Z�����̐ݒ�
  INPUT_ZDEF_NUM       : �i�q��
  INPUT_ZDEF_NUM_OMEGA : OMEGA�̊i�q�� (default=INPUT_ZDEF_NUM)
                         NCEP/NCAR�ĉ�͂̂悤�ɁA
                         OMEGA�����w�����قȂ�ꍇ�̑΍�
  INPUT_ZDEF_LEVEL     : ���̓f�[�^�̋C���� ( Upper->Lower or Lower->Upper )


&INPUT_TDEF : ����
  INPUT_TDEF_TYPE : time�̎w����@
                    = "tstep"   (��̓X�e�b�v���𒼐ڎw��)
                    = "monthly" (1�t�@�C��/1���Ƃ��āA�N������X�e�b�v�����v�Z)
                      "annual"  (1�t�@�C��/1�N�Ƃ��āA�N����X�e�b�v�����v�Z)
  INPUT_TDEF_DAYNUM : 1��������̃X�e�b�v��
                      INPUT_TDEF_TYPE �� "monthly" ���� "annual" �̏ꍇ�̂ݗL��
  INPUT_TDEF_365DAY : 1�N��365���ɌŒ肷�邩�ǂ���
                      = 0 (default, �Œ肵�Ȃ� = �[�N����)
                      = 1 (�Œ� i.e. �[�N�Ȃ�)
                      INPUT_TDEF_TYPE �� "monthly" ���� "annual" �̏ꍇ�̂ݗL��
  INPUT_TDEF_YEAR   : �N
                      INPUT_TDEF_TYPE �� "monthly" ���� "annual" �̏ꍇ�̂ݗL��
  INPUT_TDEF_MONTH  : ��
                      INPUT_TDEF_TYPE �� "monthly" �̏ꍇ�̂ݗL��
  INPUT_TDEF_TSTEP  : �X�e�b�v�� (tstep)
                      INPUT_TDEF_TYPE �� "tstep" �̏ꍇ�̂ݗL��



##### �v�Z�ݒ� #####

&WAVE : form-drag�̔g���W�J�p�̐ݒ�
  WAVE_MAX_NUMBER : �g���W�J���s�Ȃ��ő�g���B
                    0���w�肷��Ɣg���W�J���Ȃ��B
                    (default=0)


##### �o�̓t�@�C���Ɋւ����� #####
-XDEF,YDEF�ȂǁC�Ȃ��ꍇ��INPUT_*�Ɠ���

&OUTPUT : �o�̓t�@�C���̎w��
  OUTPUT_ZONAL_FILENAME  : �я󕽋Ϗ� (2����)
  OUTPUT_VINT_FILENAME   : �����ώZ�l (�ȗ���)
  OUTPUT_GMEAN_FILENAME  : �S�����ϒl (�ȗ���)
  OUTPUT_WAVE_FILENAME   : form drag�̔g���W�J�̏o�͐� (�ȗ���)
  OUTPUT_ERROR_FILENAME  : �G���[���O

&OUTPUT_ZDEF : Z�����̐ݒ�
  OUTPUT_ZDEF_NUM     : Z�����̊i�q��
  OUTPUT_ZDEF_LEVEL   : p+�� ( Upper->Lower or Lower->Upper )
