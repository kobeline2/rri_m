# rri_m
RRI model written in MATLAB

## パラメタに関するメモ: RRI_Sub
1. lasth: 計算時間 [hour]
1. dt: time step [sec]
1. maxt: 計算回数 [-]
1. NX, NY: 横, 縦のセル数
1. XLLCORNER, YLLCORNER: 矩形領域の左下の経度，緯度
1. CELLSIZE: セルの一辺の長さ [度]（単位はUTMを用いて変更可能）
1. utm: 0: 緯度経度, 1: メートル
1. NODATA: NODATAを示す値
1. NODATA: NODATAを示す値
1. land: 土地利用

## 名前変更した変数
1. area -> cellarea (230413, RRI_Sub)
1. height -> leveeHeight (230413, RRI_Sub)

## 課題
1. area_ratio = width .* len_riv / cellareaの意味 (230413, RRI_Sub)
1. ddt_chk_riv = dt_rivの意味 (230413, RRI_Sub) 多分使われない．

