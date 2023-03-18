%% RRI_Sub.f90

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 0 : FILE NAME AND PARAMETER SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RRI_Read.m を実行する

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1 : FILE READING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxt = lasth * 3600 / dt;

% read parameters
params = readvars(demfile, 'Range', 'A1:C6');
params = cellfun(@split, params, 'UniformOutput', false);

NX = int32(str2double(cell2mat(params{1}(2))));
NY = int32(str2double(cell2mat(params{2}(2))));
XLLCORNER = str2double(cell2mat(params{3}(2)));
YLLCORNER = str2double(cell2mat(params{4}(2)));
CELLSIZE  = str2double(cell2mat(params{5}(2)));
NODATA    = str2double(cell2mat(params{6}(2))); 


zs = zeros(NX,NY);        % 地表の標高
zb = zeros(NX,NY);        % 不透水層の標高
zb_riv = zeros(NX,NY);    % 河川セルの標高
domain = zeros(NX,NY);    % 0:範囲外，1:範囲内,2:河口，端
flowAcc = zeros(NX,NY);   % 集水面積
flowDir = zeros(NX,NY);   % 流向

zs      = readGisFile(demfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)'; % 転置！
flowAcc = readGisFile(accfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)';
flowDir = readGisFile(dirfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)';

land = ones(NX,NY);   %土地利用
if land_switch == 1
    land = readGisFile(landfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)';
end
disp(['num_of_landuse : ' , num2str(num_of_landuse)]);
land(land <= 0 | land > num_of_landuse) = num_of_landuse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2 : CALC PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % d1:south side length (x1,y1):左下，(x2,y2):右下
% x1 = XLLCOTNER;
% y1 = YLLCOTNER;
% x2 = XLLCOTNER + NX * CELLSIZE;
% y2 = YLLCOTNER;
% % if( utm == 0 ); d1= hubeNY_sub( x1, y1, x2, y2 ); end 
% 
% % d2:north side length (x1,y1):左上，(x2,y2):右上
% x1 = XLLCOTNER;
% y1 = YLLCOTNER + NY * CELLSIZE;
% x2 = XLLCOTNER + NX * CELLSIZE;
% y2 = YLLCOTNER + NY * CELLSIZE;
% % if( utm == 0 ); d2= hubeNY_sub( x1, y1, x2, y2 ); end 
% 
% % d3:west side length (x1,y1):左下，(x2,y2):左上
% x1 = XLLCOTNER;
% y1 = YLLCOTNER;
% x2 = XLLCOTNER;
% y2 = YLLCOTNER + NY * CELLSIZE;
% % if( utm == 0 ); d3= hubeNY_sub( x1, y1, x2, y2 ); end 
% 
% % d1:east side length (x1,y1):右下，(x2,y2):右上
% x1 = XLLCOTNER + NX * CELLSIZE;
% y1 = YLLCOTNER;
% x2 = XLLCOTNER + NX * CELLSIZE;
% y2 = YLLCOTNER + NY * CELLSIZE;
% % if( utm == 0 ); d4= hubeNY_sub( x1, y1, x2, y2 ); end 

% if utm == 1
    dx = CELLSIZE;
    dy = CELLSIZE;
% else
%     dx = (d1 + d2) / 2 / NX;
%     dy = (d3 + d4) / 2 / NY;
% end
disp(['dx [m] : ', num2str(dx),'    dy [m] : ', num2str(dy)] )

height = zeros(NX,NY); % width, depth, area_ratio は事前割り当て不要
len_riv = zeros(NX,NY);

len = sqrt(dx * dy);
area = dx * dy;

riv = (flowAcc > rivThresh);
width = (width_param_c * (flowAcc*dx*dy*1d-6) .^ width_param_s) .* riv; % 川幅
depth = (depth_param_c * (flowAcc*dx*dy*1d-6) .^ depth_param_s) .* riv; % 河道深さ
height(riv == 1 & flowAcc > height_limit_param) = height_param;         % 堤防高

% if rivfile_swtch >= 1
%     
% end
len_riv(riv == 1) = len;

% sec_switch, sec_length_switch に関するオプション

area_ratio = width .* len / area;

zb_riv = zs;
for I = 1:NX
    for J = 1:NY
        zb(I,J) = zs(I,J) - soildepth(land(I,J));
        if riv(I,J)==1; zb_riv(I,J) = zs(I,J) - depth(I,J);end
    end
end

domain(zs>-100) = 1;
domain(flowDir == 0 | flowDir == -1) = 2;
numOfCell = sum(domain >= 1,'all');
disp(['num_of_cell : ', num2str(numOfCell)])
disp(['total area [km2] : ', num2str(numOfCell * area / (10 ^ 6))])

%%% ------------------------ call riv_idx_setting ---------------------------- %%%
riv_count = sum(domain > 0 & riv == 1,'all'); % number of river cell

domAndRiv             = (domain>0 & riv==1);      % 1:rivセル， 0:範囲外 or sloのみ                        
[riv_idx2i,riv_idx2j] = find(domAndRiv');         % rivセルのx,y座標
riv_ij2idx            = zeros(size(domAndRiv));   % rivセルの座標
riv_ij2idx(domAndRiv) = 1:sum(domAndRiv, 'all');  % rivセルの番号を振る  
domain_riv_idx        = domain(domAndRiv);        % domain
width_idx             = width(domAndRiv);         % 川幅
depth_idx             = depth(domAndRiv);         % 河道深さ
% height_idx            = height(domAndRiv);        % 堤防高
area_ratio_idx        = area_ratio(domAndRiv);    % セルにおける河川の割合
zb_riv_idx            = zb_riv(domAndRiv);        % 不透水層の標高
dif_riv_idx           = dif(land(domAndRiv))';    % 1:拡散波近似，2:kinematic
% RRI_Input.txtで指定する土地利用ごとのパラメータを縦ベクトルに統一するため転置する(Sloも同様)
% sec_map_idx           = sec_map(domAndRiv);       % 
len_riv_idx           = len_riv(domAndRiv);       % 河川の長さ

[dis_riv_idx,down_riv_idx] = rivIdxSetting(riv_count,NX,NY,domain,riv,flowDir,dx,dy,riv_ij2idx);

%%% ------------------------ call slo_idx_setting ---------------------------- %%%
i4 = 4;

slo_count = sum(domain > 0,'all'); % number of slope cell
domAndSlo             = (domain>0);               % 1:sloセル， 0:範囲外                       
[slo_idx2i,slo_idx2j] = find(domAndSlo');         % sloセルのx,y座標
slo_ij2idx            = zeros(size(domAndSlo));   % sloセルの座標
slo_ij2idx(domAndSlo) = 1:sum(domAndSlo, 'all');  % sloセルの番号を振る  
domain_slo_idx        = domain(domAndSlo);        % domain
zb_slo_idx            = zb(domAndSlo);            % 不透水層の標高
acc_slo_idx           = flowAcc(domAndSlo);       % 集水面積
land_idx              = land(domAndSlo);          % 土地利用

isLand = land(domAndSlo);
dif_slo_idx      = dif(isLand)';         % 1:拡散波近似，2:kinematic
ns_slo_idx       = ns_slope(isLand)';    % ns_slope
soildepth_idx    = soildepth(isLand)';   %
gammaa_idx       = gammaa(isLand)';      %
ksv_idx          = ksv(isLand)';         %
faif_idx         = faif(isLand)';        %
infilt_limit_idx = infilt_limit(isLand)';%

ka_idx           = ka(isLand)';          %
gammam_idx       = gammam(isLand)';      %
beta_idx         = beta(isLand)';        %
da_idx           = da(isLand)';          %
dm_idx           = dm(isLand)';          %
ksg_idx          = ksg(isLand)';         %
gammag_idx       = gammag(isLand)';      %
kg0_idx          = kg0(isLand)';         %
fpg_idx          = fpg(isLand)';         %
rgl_idx          = rgl(isLand)';         %

[down_slo_idx,dis_slo_idx,len_slo_idx,down_slo_1d_idx,dis_slo_1d_idx,len_slo_1d_idx] ...
    = sloIdxSetting(slo_count,NX,NY,domain,flowDir,dx,dy,slo_ij2idx,i4,eightFlowDir);

%%%-------------------------- call dam_read ----------------------------%%%

hs = zeros(NX,NY);
hr = zeros(NX,NY);
hg = zeros(NX,NY);
gampt_ff = zeros(NX,NY);
gampt_f = zeros(NX,NY);
qrs = zeros(NX,NY);


% where(riv ~= 1) hr = -0.1
% where(domain == 0) hs = -0.1

%%%-------------------------- initial condition ----------------------------%%%

%%%-------------------------- boundary condition ----------------------------%%%

%%%-------------------------- div file ----------------------------%%%

df = readmatrix(location_file);
for I = 1:size(df,1)
    hydro_i(I) = df(I, 2); % Y座標(列番号)
    hydro_j(I) = df(I, 3); % X座標(行番号)
end

%%%-------------------------- array initialization ----------------------------%%%

rain_i = zeros(NY, 1);
rain_j = zeros(NX, 1);


%%%-------------------------- gw initial setting ----------------------------%%%
if init_gw_switch ~= 1
    hg_idx = zeros(slo_count, 1); %  call hg_init( hg_idx )
%  call sub_slo_idx2ij( hg_idx, hg )
end

%%%-------------------------- initial srtorage calculation ----------------------------%%%
rain_sum = 0;
aevp_sum = 0;
pevp_sum = 0;
sout = 0;
si = 0;
sg = 0;
% open( 1000, file = outfile_storage )
% call storage_calc(hs, hr, hg, ss, sr, si, sg)
% sinit = ss + sr + si + sg
% write(1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
%   (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg

%%%-------------------------- reading rainfall data ----------------------------%%% 
df = readmatrix(rainfile,'NumHeaderLines',0);
nx_rain = df(1,2);
ny_rain = df(1,3);
tt_max_rain = size(df,1) / (ny_rain + 1);   % fortran - 1
disp([tt_max_rain,nx_rain, ny_rain])        

qp = zeros(nx_rain, ny_rain, tt_max_rain);  % 順番変更
for TT = 1:tt_max_rain
    t_rain(TT) = df((ny_rain+1)*(TT-1)+1,1);
    qp(:,:,TT) = df((ny_rain+1)*(TT-1)+2:(ny_rain+1)*TT,1:nx_rain)'; % 転置して地形データに合わせる
end

qp = qp / 3600/ 1000;  % mm/h -> m/s

for I = 1:NX
    rain_j(I) = fix((XLLCORNER + (double(I) - 0.5) * CELLSIZE - XLLCORNER_RAIN) / CELLSIZE_RAIN_X)  + 1;
end
for I = 1:NY
    rain_i(I) = ny_rain - fix((YLLCORNER + (double(NY) - double(I) + 0.5)*CELLSIZE - YLLCORNER_RAIN) / CELLSIZE_RAIN_Y );
end
disp("done: reading rain file")

%%%-------------------------- reading evp data ----------------------------%%%


% call sub_slo_ij2idx( hs, hs_idx ) 
% call sub_riv_ij2idx( hr, hr_idx ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3 : CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rain_sum = 0;
aevp_sum = 0;
sout = 0;

out_dt = maxt / outnum;
out_dt = max(1, out_dt);
out_next = round(out_dt);
TT = 0



toc


