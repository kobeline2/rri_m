# rri_m
RRI model written in MATLAB

## パラメタに関するメモ: RRI_Read
1. utm: 座標系(1:緯度経度，0:メートル)
1. eightFlowDir: 流向パターン(1:8方向，0:4方向)
1. lasth: 計算時間 [hour]
1. dt: time step [sec]
1. dt_riv: time step for river[sec]

## パラメタに関するメモ: RRI_Sub
1. maxt: 計算回数 [-]
1. NX, NY: 横, 縦のセル数
1. XLLCORNER, YLLCORNER: 矩形領域の左下の経度，緯度
1. CELLSIZE: セルの一辺の長さ [度]（単位はUTMを用いて変更可能）
1. utm: 0: 緯度経度, 1: メートル
1. NODATA: NODATAを示す値
1. zs: 地表の標高
1. flowAcc: 集水面積
1. flowDir: 流向
1. land: 土地利用
1. dx: セルの横方向の長さ[m]
1. dy: セルの縦方向の長さ[m]
1. len: セルの対角線長( =sqrt(dx*dy) )
1. area: セルの面積( = dx*dy )
1. riv: 1：riverありのセル(flowAcc>rivThresh)，0：riverなしのセル
1. width: 川幅
1. depth: 河道深さ
1. leveeHeight: 堤防高
1. len_riv: セルの対角線長を河川領域にのみ定義したもの
1. area_ratio: セルにおけるriverの占める面積の割合
1. zb: 不透水層の標高
1. zb_riv: 河川の不透水層の標高
1. domain: 0:解析範囲外，1:範囲内, 2:河口，端
1. numOfcell: 解析対象のセル数( =nnz(domain)　)

#### riverに関するもの（idx: 2次元のデータを１次元(rivありセルのみ)にしたもの）
1. riv_count: riverありのセル数
1. domAndRiv: 1:rivセル， 0:範囲外 or sloのみ                        
1. riv_idx2i: rivセルのy座標
1. riv_idx2j: rivセルのx座標
1. riv_ij2idx: rivセルの番号  
1. dis_riv_idx: 次のセルとの距離 
1. down_riv_idx　次のセル( riv_ij2idx of next cell )

## 名前変更した変数
1. area -> cellarea (230413, RRI_Sub)
1. height -> leveeHeight (230413, RRI_Sub)

## 課題
1. area_ratio = width .* len_riv / cellareaの意味 (230413, RRI_Sub)
1. ddt_chk_riv = dt_rivの意味 (230413, RRI_Sub) 多分使われない．

