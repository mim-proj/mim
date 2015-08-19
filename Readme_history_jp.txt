□TODOリスト
  高速化
  Pzの緯度分布、テーラー展開の必要あり(Pzの全球平均はそのままの方がbetter)
  waveを統合(aquaのanaで使っているやつ)


□mochi版と値が少し違う理由(特に地面付近で大。一般には違いは問題にならない程度)
    biseki_bibun() (旧bibun()) で ncalc=1(kodama), ncalc=2(mochi)
      c_kz_keが4%程度動く
    pt_dotでomegaを使う(kodama)か使わない(mochi)か
      Qが使えない場合:  epz_ut(kodama)  <-- 対応 --> epzut(mochi)
                        c_kz_ke(kodama) <-- 対応 --> kzke(mochi)
      Qが使える場合:    epz_ut(kodama)  <-- 対応 --> epzut2(mochi)
                        c_kz_ke(kodama) <-- 対応 --> (kzkeuw2)(mochi)
    pe_1(kodama) <-- 対応 --> pe11(mochi)
    mochi版のpe1は経度分布を見るのに有効

□このプログラムの系譜(推定)
Iwasaki → Shimizugichi → Tanaka → Ujiie → Miyazaki → Kodama(0.10系)
                                     │                        ↓
                                     └→ Uno → Mochi → Kodama(0.20系)

□更新履歴(新しいものが上)

2009.03.27
  Version 0.33 Release 1

2009.03.19
  getpt_p()のp_pd内挿で、p線形内挿をlog-p線形内挿に変更。
    これにより、地表面での気圧の値の矛盾をほぼ解消。全変数に1%程度のずれ。

2009.03.10
  非断熱加熱q_3dの入力がない場合、従来は0にしていたが、D(pt)/Dtの
    全微分から近似するように変更。誤差が大きいので注意。
    また、最初のタイムステップでは偏微分d(pt)/dt=0とされるのも注意。

2009.02.24
  Version 0.32 Release 2
  namelistにINPUT_UNIT_Zを追加。

2009.02.13
  namelist_check()でOUTPUT_WAVE_FILENAMEに関するerrorの基準を変更。

2009.02.12
  Version 0.32 Release 1
  コンパイルオプションに-heap-arrayを追加。これにより、高解像度データ
    (e.g. JRA125)を扱う際にstack不足でsegmentation faultが起きるのを抑止。
    ただし、ifort version 10.0以降のみ対応。
  get_phi_dagger(), mount_modify() で、pout(ko+1)にアクセスしてしまうバグを修正
    但し、結果への影響はあり得ない。
  バグ: biseki_biseki() の x_pt(im,jm,km) を x_pt(im,jm,ko)に修正
    但し、km=koの場合は影響なし。

2008.09.19
  Gyの定義を変更。
    従来: -rho_0 (v'^2)_zm
    今回: -rho_0 a (v'^2)_zm

2008.09.07
  wave計算がT63データで失敗する不具合を解決。
    ただし、不具合の厳密な理由は不明(where文の制限?)

2008.09.01
  Version 0.31 Release 1
  ソースコードの英語化完了(fft.f90を除く)。

2008.08.27
  mountain_modify()の最後のc_pz_kz2算出部分で、
    k=1のときもpout_new2を使うように変更。
    これは単なるプログラムの分かりやすさの問題であって、本質的な変更ではない。
    違いはごくごくわずか。
  integral_pt_ym()を新設、積分区間の上端が不適切であったのを修正。
    ただし、結果への影響はごくわずか(az_vintに影響:0.01%程度)。
  check_range()で、一定回数以上warningが検知されたら落ちるように仕様変更。

2008.08.26
  get_phi_zm.f90をget_phi_dagger.f90に改名。
  get_phi_dagger()で、p+系におけるptの地表面値が
    単なるzonal meanになっていたバグを修正。
    C(AZ->KZ)に影響。但し全球平均値で0.1%程度とごくわずか。

2008.08.24
  intpl()を廃止。
  getpt()まわりの変数を整理。
  integral_pt()で、積分区間の上端が不適切であったのを修正。
    ただし、結果への影響はごくわずか(ae_vintに影響:0.01%程度)。
    <-実は1%程度の寄与???

2008.08.21
  biseki_sekibun()の内挿式で、d<0になったとき不適切な値になるバグを修正。
    ただし、結果への影響はごくわずか。
  biseki_y.f90の p_pt, x_pt を p_pdd, x_pdd に変更。

2008.08.20
  Version 0.30 Release 1
  PS(地表面気圧)の代りに MSL(海面較正気圧), TS(地表気温)の入力に対応。
    誤差は小さい(divf等は地表付近でやや大だが問題ない程度)

2008.08.19
  biseki_sekibun()の内挿式で、d<0になったとき不適切な値になるバグを修正。
    また、x_pt(im,jm,km)をx_pt(im,jm,ko)に修正。
    ただし、結果への影響は上端付近が中心であり、ごくわずか。
  biseki.f90の p_pt, ps, x_pt を p_pd, p_sfc, x_pd に変更。

2008.08.11
  get_pt_dot()をget_pt_dot_omega()に改名、get_pt_dot_q()を新設。
  pd_pt を pd_pdd に変更。
  p_pt, z_pt を p_pd, z_pd に変更。
  lagmain*.f90 を mim*.f90 に変更。

2008.08.04
  Version 0.29 Release 2
  ドキュメントを整理。

2008.08.01
  Version 0.29 Release 1
  pt_dotの計算でQを使うように変更(Qが使用可能なら)。

2008.07.12
  Version 0.28 Release 1

2008.07.10
  namelistのOUTPUT_DEBUG_FILENAME*を廃止。
  namelistでWAVE_NUMBER=0を許すように変更。
  namelistのINPUT_ZDEF_LEVEL, OUTPUT_ZDEF_LEVELを、上から下、下から上、
    どちらで指定しても受け付けるように変更。
  namelistにINPUT_YDEF_SOUTH, INPUT_YDEF_NORTHを追加。
  check_range()で値を外れたとき、実行を停止せず続行するように変更。

2008.07.09
  zonal.ctlの変数名のミス(dgzのはずがdgyになっていた)を修正。

2008.07.04
  parameter.f90のeconv, divfの値の範囲を、それぞれ10倍に緩和。

2008.06.02
  Version 0.27 Release 3
  parameter.f90のz_minを-1000から-3000に変更。
  さらにgetpt_lev()のppp計算部分をlog(p)内挿に戻す。
    z_zmがやや異常値ぎみ。

2008.06.01
  Version 0.27 Release 2
  getpt_lev()のdlev計算部分をlog(p)内挿に戻す。
    線型だとdepz_w, dept_ut計算時に発散することがある(原因不明)

2008.05.30
  Version 0.27 Release 1
  dkzdt_vkz, dkzdt_wkzの計算を追加。

2008.05.29
  Keの計算の順序を変更、値が微妙にずれる。
  dkedt_vke, dkedt_wkeの計算を追加。

2008.05.27
  Namelistにendianを導入。
  get_z_pt()の計算順序を変更。結果に無視できる程度のずれ。

2008.05.26
  Version 0.26 Release 1
  Qzを修正。
  内挿時にlog(p)を用いるのを原則やめる。
    log計算の誤差を考えると素直な線型内挿で統一した方がよさそう。
    ただしpそのものを内挿する場合は上端で負になる可能性があるので、
    上端のみlogを使用するべき。
  pd_pの内挿方法をlog(p)型に変更、鉛直積算QEの数%のずれ。

2008.05.13
  pd_p (標準p面でのp+) の計算を追加
  Aeの鉛直積算値も修正、全球平均で1%程度のずれ。

2008.05.12
  vint.ctlからae_s, pe_totalを削除。
  Azの鉛直積算をθyminに変更。積分範囲の不一致を解消する補正項も追加。
    ただし鉛直積算、全球平均への寄与は極めて小さい。
  gmean.ctlにAzを追加。

2008.04.30
  diabatic_aqz()をdiabatic_qz_gmean()に改名。

2008.04.29
  energy_ae()をenergy_ae_total()に改名。

2008.04.21
  PとAの区別を厳密化。(ソースコード、コントロールファイル)
  aqeをqeに改名。
  qz_zmをqgz_zmに改名、A_ZとP_Zの生成を区別するため。

2008.03.28
  Version 0.25 Release 1

2008.03.27
  pt(p++s)の定義の間違いを修正。Pzが場所によって1割ずれる。
  aqzの全球平均を出力。
  getpt_*()をgetpt_etc.f90に統合、整理。
  getpt1_y()をgetpt1()に統合。
  ptiter.f90の内挿方法を修正。

2008.03.25
  zonal.ctlに G Flux を追加。
  EP Flux 発散の単位をm/s^2に変更。
  Namelistの大幅改訂。
  ソースコードをsrc以下へ移動、Makefileの大幅書き換え。

2008.01.20
  Version 0.24 Release 2
  形状抵抗の波数展開の計算ができないbugを解決
  Version 0.24 Release 1  

2008.01.12
  pe_total_zmを計算。これは部分積分を行う前のAeである。
  integral_pt()の下端が、θsの東西平均だったのをθminに変更。
  これにより、peの全球平均値が2割程度増加。

2008.01.08
  integral_pt()を新設。
  merifional_integral.f90をintegral.f90に改名、zsekibun.f90を吸収

2008.01.07
  lagmain_var.f90の変数を整理。
  phi_zmをphi_daggerに改名。実態に合せる。
  namelistにFILE_DEBUG1, FILE_DEBUG2を追加。

2007.12.25
  bibun(), sekibun() をモジュールbisekiとして統合。
  epz, gz のdivを求めるスキームをderivative_z()に任せる。
    これにより、結果が少しずれる。
  derivative_z()を新設。これは half level を用いたz微分。

2007.12.24
  div_epflux.f90 を epflux_div.f90 に改名。
  epflux_*.f90 を epflux.f90 に統合。
  get_var_w_x_zm.f90 を新設、epflux_z_uw.f90 の殆どの処理を移行。
  gflux.f90を新設、変数gy, dgyなどを導入。
  energy_conv_*()の値チェックを加える。
  energy_tendency.f90を新設、u_epflux.f90を吸収。
  epy計算に変数rhoを使用するようにする。
  c_kz_keの項の分けかたを変更。grads出力変数も変更。
  c_pe_kz* を c_kz_pe*に変更。
  phi_zmの定義を、geopotential height から geopotential に変更(本来の定義)。
  変数 epzw or epz_w を epz_uvに、epzut を epz_ut に変更し、epz_uwを新設。
  depz_w等も同様。
  namelistのFILE_TOPO_KINDの指定の仕方を変更。
    geopotentialの場合、"g"と指定できるようになる。
    過去バージョンとの互換性を保つため、"gh"も残す。

2007.12.17
  Version 0.23 Release 1  
  JRA25動作確認。
  namelistにFILE_U, FILE_V, FILE_T を追加。

2007.12.01
  com_var module に tantbl を追加。
  parameter module に sec_day を追加。
  mountain_modify.f90 から不要なコードを削除。約1000行→約300行。
  mountain_modify.f90 で、内挿計算のバグを修正
  1hPa以下でc_pz_kzが0になる場合があった。現在解消済。
  鉛直積算値したc_pz_kzへの影響は1%程度。
  energy_conv.f90を新設。エネルギー変換のコードを移動。  
  div_z()を用いていたdivz_tzm, divphi_t, d_u_epzをderivative_p()に変更。
  地表面を考えているので、derivative_p_nops()よりもいいはず。
  地表付近でずれる。鉛直積算もずれる(10%程度かと)
  値の確認はしていない(元々no checked)
  div_z.f90を廃止し、derivative_p_nops()を新設。両者の違いは単位と重力のみ。
  bibun()のncalc=1を元に、derivative_p()を新設。bibun()の実質処理を移動。
  phi_t を p_dphi_dt* に改名
  z_t* を p_dz_dt* に改名。
  get_phi_zm.f90 : com_varをuse。
  *_past(:,:)を厳密に"変数名_past(:,:)"となるように改名。
  pt_before(:,:) を pt_past(:,:)に改名。
  pt_past(:,:)を廃止。
  get_pt_dot.f90 : com_varをuse。
  地表付近のp_ptの値が、プログラム前半と後半で異なる不整合を解消。
  従来はepflux_z_form_z_pt()で置き換わっていたのを廃止。
  p_ptの地中値はpsに一致させない。必要ならその都度作る(e.g. form drag)
  これによる影響は地表面付近で最大1%程度、鉛直積算値でも1%程度。
  いずれにせよ、地面付近の値は見るべきではない。

2007.11.30
  epflux_z_uw.f90 : com_varをuse。
  epflux_z_form.f90 : com_varをuse。
  epflux_z.f90 を epflux_z_form.f90 に改名
  epflux_y() : com_varをuse。
  [重要] divfの値にFyとFzのform fragしか考慮されていなかったバグを解決。
    温暖化論文のdivFも少しずれている可能性あり。矢印はOK。要調査。
  div_epflux_y(), div_epflux_z() : com_varをuse。
  depy.f90 と depz.f90 を div_epflux.f90 に統合。
  bibun() : com_varをuse。
  sekibun() : com_varをuse。
  setpt0() : com_varをuse。
  namelist: WAVE_OUTPUT をFILE_WAVEに変更、WAVE_SWITCHを廃止。
  FILE_P_Zの出力を休止
  grads.f90をモジュール化。grads_open()、grads_write()等をlagmain.f90に導入。
  モジュール名、サブルーチン名とファイル名を対応させる

2007.11.29
  com_var.f90を新設、プログラム全体で使用する共有変数を定義する。
  var.f90 -> lagmain_var.f90に改名、lagmain.f90のみでuseすることを強調。

2007.11.26
  Version 0.22 Release 2
  check_range()を作成、値のチェックに使用。
  param.f90の定数にparameter属性を付ける
  getlev()のdlev内挿部分を変更。
    全てlog(p)内挿であったのをデータTop付近のみ線形内挿にする。
    これにより、dlevが小さくなりすぎてp_ptが0になる不具合(滅多に起きないが)を
    回避。

2007.11.02
  grads()を廃止(term_main.f90は未対応)
  grads_read_2d()、grads_read_3d()を廃止、grads_read()に置き換え
  カウント数の画面表示を変更
  nlist.f90のuseをlagmain.f90以外にも解禁
    calendar.f90でuse

2007.10.22
  Version 0.22 Release 1
  qz, qeの出力を追加
  namelistのZONAL_Qを廃止、q_zmの単位変更、zonal.ctlに組込む
  c_kz_ke_uvを2つに分ける

2007.10.19
  dkedt_uvu_y(Keの生成率の(u'v') d/dy(u_zm)の項)を追加
  zonal.ctlのuepyをduepyに名称変更、出力順番の変更

2007.10.18
  Version 0.21 Release 1
  epflux_z_wave()に新fft()を採用
  fft()のデータ数を任意に設定できるように改良
    旧fft()はfft_quick()と名称変更
    計算手法はマニュアルを参照

2007.10.09
  Version 0.20 Release 3
  bibun()の～0割りを回避
    特にエネルギー計算に影響。今までうまくいっていた理由は不明。
    温暖化実験cntlの2002.12.16などに発生

2007.10.03
  Version 0.20 Release 2

2007.10.02
  mount_modify()を反映
  getlev()でdlevをpからlog(p)で決定するように変更

2007.10.01
  getpt.f90の反復回数: 6 -> 8
  地形データのyrev有無が、FILE_TOPO_YREVではなくYDEF_YREVで判断されるバグを解決

2007.09.17
  Version 0.20 Release 1
  TERM, PHASEへの対応を一旦打ち切り
  ifortでのコンパイルに成功
  ptiter()で 0/0 となるbugを修正(fujitsu f90では0になるため問題なし)
  div_y.f90をu_epflux.f90に改名
  2次相関の命名規則(e.g. u_d_v_d_zm)を制定、既存の変数名に適用
  zsekibun()を整理
  変数名の変更と並べ替え
    epz -> epz_form
    depz -> depz_form
    ep_z -> epz (要注意)
    エネルギー変換 -> c_*
  zonal.ctlにtを追加
  lag.ctlをzonal.ctlに改名

2007.09.15
  mount_modify.f90をmochi版から導入。ただし現状では不使用。
  energy_az()を追加
  energy_ae()を追加

2007.09.12
  epz_uw()の引数からz_zmを削除
  c(pz,kz)の計算でderivative_y()を使うように変更
  derivative_x()をp_daggerから移植、追加。dz_dlonの計算に利用
  derivative_y()をp_daggerから移植、追加。dz_dlatの計算に利用
  get_phi_zm()を追加、lagmainから処理を移動
  lagmainのz_zmのmodifyをget_z_pt()に移動

2007.09.09
  get_pt_dot()を追加、lagmainから処理を移動し、計算精度を上げる(pt_dotの計算法)
  namelistにTIMESTEP_DTを追加
  array()を廃止、配列式に置き換え
  sekibun(), bibun()を整形
  getpt_y()を追加( getpt()の南北版 )
  getlev()を整形
  meridional_integral()を追加

2007.09.08  
  getpt()を3次元化
  stable()を3次元化
  getpzm()を3次元化
  intpl()を3次元化、一般のp->θ内挿にも使用できるようにする
  getpt1()を3次元化
  zmean()を廃止し、組み込み関数sum()で置き換え
  getpt()内でcallしていたsetpt0()をlagmainへ移動
  setpt0()を3次元化
  lagmainの call getpt() 付近の処理を簡略化
  ps_zmの計算を組み込み関数sum()を用いて簡略化、一時変数pssを削除
    計算の順番が変わるため、一部の変数の結果が微妙にずれる
  namelistにFILE_Q_UNDEFを追加
  undef.f90を追加、lagmain中のundef処理をundef_fill()に移動
  omegaファイルのopen位置を変更
  version 0.20系の開発開始(mochi版(ene-sam)との統合)

2007.04.08
  version 0.13 release 3
  Peの積分値の出力を追加

2007.04.06
  version 0.13 release 2

2007.04.05
  namelistにFILE_VINTを追加

2007.03.16
  version 0.13 release 1
  namelistにZDEF_ZDEF_IN_OMEGAを導入

2007.02.06
  version 0.13 beta 2
  getpt.f90の反復回数: 5 -> 6

2006.11.13
  version 0.13 beta 1
  Peを出力するようにする(かなり怪しい)
  Pzを修正(まだ怪しい)

2006.09.21
  version 0.12 release 1
  namelistにFILE_TOPO_YREVを追加(未チェック)

2006.03.30
  version 0.11 release 3
  divepzのバグ修正(z -> z+)

2006.03.14
  version 0.11 release 2
  FILE_ZONAL_Qを省略可能にする．
  version 0.11 release 1
  namelistの値チェックを強化
  OMEGAを省略可能にする
  namelistにUNDEFを追加

2006.03.01
  FILE_Q未指定時のメッセージを消去

2006.01.04
  getpt()で、pout(:)のintentがinではなくoutになっていたバグを修正
  bibun()の引数変更(nameを廃止)に伴う変更
  sekibun()の引数変更(nameを廃止)に伴う変更

2005.10.02
  namelist.f90を新設

2005.09.28
  epzの波数展開ルーチンを追加

2005.08.05
  YDEF_YDEF と YDEF_YREV の混同によるバグを解消

2005.07.12
  parameterをモジュール化

2005.04.27
  XDEF, YDEFを組み込む

2005.04.26
  kzkevvにバグ -> 最適化のせいで結果がずれた(無視できるくらい小さい)
  allocatable導入

2005.04.12
  namelistにZDEFグループを追加

2005.04.11
  年月日wp処理するcalendar.f90を追加
  namelistを導入, 全ての入出力ファイル名を指定可能にする
