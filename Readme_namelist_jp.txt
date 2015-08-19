MIM namelist版 by 小玉 (原本)


-データは4バイト書式なしバイナリデータ (GrADS形式) に対応しています。
-説明は GrADS のコントロールファイルのスタイルに沿っています
-各々の型を知りたい場合は、src/namelist.f90 を参照してください。



##### 入力ファイルに関する設定 #####

&INPUT : 入力ファイルの指定
  INPUT_TYPE           : = "general" (書式なし4バイト実数データ, default)
                           "tohoku_ncep" (独自形式のncep1,ncep2用)

  INPUT_UVT_FILENAME   : U, V, T ファイル名 (default='')
  INPUT_U_FILENAME     : Uファイル名
                         (INPUT_UVT_FILENAMEが指定されていない場合に有効)
  INPUT_V_FILENAME     : Vファイル名
                         (INPUT_UVT_FILENAMEが指定されていない場合に有効)
  INPUT_T_FILENAME     : Tファイル名
                         (INPUT_UVT_FILENAMEが指定されていない場合に有効)
  INPUT_PS_FILENAME    : PSのファイル名 (default='')
  INPUT_MSL_FILENAME   : MSLのファイル名
                         (INPUT_PS_FILENAMEが指定されていない場合に有効)
  INPUT_TS_FILENAME    : 地表Tのファイル名
                         (INPUT_PS_FILENAMEが指定されていない場合に有効)
  INPUT_Z_FILENAME     : Zファイル名
  INPUT_OMEGA_FILENAME : OMEGAのファイル名
                         (省略可、省略時 OMEGA=0 として取り扱う)
  INPUT_TOPO_FILENAME  : 地形ファイル名
  INPUT_Q_FILENAME     : 非断熱加熱ファイル名
                         (省略可、省略時 Q=0 として取り扱う)


&INPUT_UNIT : 単位の指定
  INPUT_UNIT_Z        : = "m" (地形高度, default)
                        = "m^2/s^2" (ジオポテンシャル)
  INPUT_UNIT_PS       : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_MSL      : = "hPa" (default)
                        = "Pa"
  INPUT_UNIT_TOPO     : = "m" (地形高度, default)
                        = "m^2/s^2" (ジオポテンシャル)


&INPUT_UNDEF : undefの指定
  INPUT_UNDEF_DEFAULT : 各ファイル共通のUNDEF値 (default=9.999e+20)
  INPUT_UNDEF_UVT     : UVTファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_U       : UファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_V       : VファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_T       : TファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_PS      : PSファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_MSL     : MSLファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_TS      : TSファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Z       : ZファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_OMEGA   : OMEGAファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)
  INPUT_UNDEF_Q       : QファイルのUNDEF値 (default=INPUT_UNDEF_DEFAULT)


&INPUT_ENDIAN : endianの設定
  INPUT_ENDIAN_DEFAULT : 各ファイル共通のendian
                         = "little" (little endian, default)
                         = "big"    (big endian)
  INPUT_ENDIAN_UVT   : UVTファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_U     : Uファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_V     : Vファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_T     : Tファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_PS    : PSファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_MSL   : MSLファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TS    : TSファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Z     : Zファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_OMEGA : OMEGAファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_Q     : Qファイルのendian (default=INPUT_ENDIAN_DEFAULT)
  INPUT_ENDIAN_TOPO  : TOPOファイルのendian (default=INPUT_ENDIAN_DEFAULT)


&INPUT_XDEF : X方向の設定
              一様間隔のみ対応
              outputにもそのまま継承される
  INPUT_XDEF_NUM  : 格子数



&INPUT_YDEF : Y方向の設定
              解析、outputにもそのまま継承される
  INPUT_YDEF_TYPE   : levelの指定方法
                      = "lat_degree"   (配列[度], default)
                      = "lat_radian"   (配列[radian])
                      = "linear"       (一様間隔[度])
  INPUT_YDEF_NUM    : 格子数
  INPUT_YDEF_LEVEL  : 緯度を表す配列 ( South->North or North->South )
                      INPUT_YDEF_TYPE が "lat_degree" 又は "lat_radian" の
                      場合のみ有効
  INPUT_YDEF_SOUTH  : 南端の緯度[°] (default=-90)
                      INPUT_YDEF_TYPE が "linear" の場合のみ有効(未チェック)
  INPUT_YDEF_NORTH  : 北端の緯度[°] (default=90)
                      INPUT_YDEF_TYPE が "linear" の場合のみ有効(未チェック)
  INPUT_YDEF_YREV_DEFAULT
                    : 各ファイルの共通のYREV
                       = 0 (yrevなし i.e. 南→北, default)
                       = 1 (yrevあり i.e. 北→南)
  INPUT_YDEF_YREV_TOPO
                    : 地形ファイルのYREV(データ配列の向き)
                      = 0 (yrevなし i.e. 南→北)
                      = 1 (yrevあり i.e. 北→南)
                      指定がない場合、INPUT_YDEF_YREV_DEFAULTに従う


&INPUT_ZDEF : Z方向の設定
  INPUT_ZDEF_NUM       : 格子数
  INPUT_ZDEF_NUM_OMEGA : OMEGAの格子数 (default=INPUT_ZDEF_NUM)
                         NCEP/NCAR再解析のように、
                         OMEGAだけ層数が異なる場合の対策
  INPUT_ZDEF_LEVEL     : 入力データの気圧面 ( Upper->Lower or Lower->Upper )


&INPUT_TDEF : 時刻
  INPUT_TDEF_TYPE : timeの指定方法
                    = "tstep"   (解析ステップ数を直接指定)
                    = "monthly" (1ファイル/1月として、年月からステップ数を計算)
                      "annual"  (1ファイル/1年として、年からステップ数を計算)
  INPUT_TDEF_DAYNUM : 1日あたりのステップ数
                      INPUT_TDEF_TYPE が "monthly" 又は "annual" の場合のみ有効
  INPUT_TDEF_365DAY : 1年を365日に固定するかどうか
                      = 0 (default, 固定しない = 閏年あり)
                      = 1 (固定 i.e. 閏年なし)
                      INPUT_TDEF_TYPE が "monthly" 又は "annual" の場合のみ有効
  INPUT_TDEF_YEAR   : 年
                      INPUT_TDEF_TYPE が "monthly" 又は "annual" の場合のみ有効
  INPUT_TDEF_MONTH  : 月
                      INPUT_TDEF_TYPE が "monthly" の場合のみ有効
  INPUT_TDEF_TSTEP  : ステップ数 (tstep)
                      INPUT_TDEF_TYPE が "tstep" の場合のみ有効



##### 計算設定 #####

&WAVE : form-dragの波数展開用の設定
  WAVE_MAX_NUMBER : 波数展開を行なう最大波数。
                    0を指定すると波数展開しない。
                    (default=0)


##### 出力ファイルに関する情報 #####
-XDEF,YDEFなど，ない場合はINPUT_*と同じ

&OUTPUT : 出力ファイルの指定
  OUTPUT_ZONAL_FILENAME  : 帯状平均場 (2次元)
  OUTPUT_VINT_FILENAME   : 鉛直積算値 (省略可)
  OUTPUT_GMEAN_FILENAME  : 全球平均値 (省略可)
  OUTPUT_WAVE_FILENAME   : form dragの波数展開の出力先 (省略可)
  OUTPUT_ERROR_FILENAME  : エラーログ

&OUTPUT_ZDEF : Z方向の設定
  OUTPUT_ZDEF_NUM     : Z方向の格子数
  OUTPUT_ZDEF_LEVEL   : p+面 ( Upper->Lower or Lower->Upper )
