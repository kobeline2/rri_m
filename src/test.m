%% test readGisFile

tic
% input files
demfile = "topo\adem_trim.txt";
accfile = "topo\acc_trim.txt";
dirfile = "topo\adir_trim.txt";

% read parameters
params = readvars(demfile, 'Range', 'A1:C6');
params = cellfun(@split, params, 'UniformOutput', false);

NX = int32(str2double(cell2mat(params{1}(2))));
NY = int32(str2double(cell2mat(params{2}(2))));
XLLCORNER = str2double(cell2mat(params{3}(2)));
YLLCORNER = str2double(cell2mat(params{4}(2)));
CELLSIZE  = str2double(cell2mat(params{5}(2)));
NODATA    = str2double(cell2mat(params{6}(2))); 

% 地表の標高
zs      = readGisFile(demfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE);
% 集水面積
flowAcc = readGisFile(accfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE);
% 流向
flowDir = readGisFile(dirfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE);

toc

% error check
% flowAcc      = readGisFile(demfile,NX-1,NY,XLLCORNER,YLLCORNER,CELLSIZE);


%% test riv_idx_setting

riv_count = sum(domain > 0 & riv == 1,'all'); % number of river cell
  
down_riv_idx = zeros(riv_count,1);
dis_riv_idx = zeros(riv_count,1); % rivセルの距離


domAndRiv = (domain>0 & riv==1);                   % 1:rivセル， 0:範囲外 or sloのみ                        
[riv_idx2i,riv_idx2j] = find(domAndRiv');          % rivセルのx,y座標
riv_ij2idx = zeros(size(domAndRiv));              % rivセルの座標
riv_ij2idx(domAndRiv) = 1:sum(domAndRiv, 'all');  % rivセルの番号を振る  
domain_riv_idx         = domain(domAndRiv);  % domain
width_idx              = width(domAndRiv);   % 川幅
depth_idx              = depth(domAndRiv);   % 河道深さ
% height_idx = height(domAndRiv);   % 堤防高
% area_ratio_idx = area_ratio(domAndRiv);   % セルにおける河川の割合
% zb_riv_idx = zb_riv(domAndRiv);   % 不透水層の標高
% dif_riv_idx = dif(domAndRiv);   %%%%%　修正　%%%%%%%
% sec_map_idx = sec_map(domAndRiv);   % 
% len_riv_idx = len_riv(domAndRiv);   % 

[dis_riv_idx,down_riv_idx] = rivIdxSetting(riv_count,NX,NY,domain,riv,flowDir,dx,dy,riv_ij2idx);
disp(size(dis_riv_idx))
disp(size(down_riv_idx))


%% test slo_idx_setting

