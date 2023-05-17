%% RRI.f90

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 0 : FILE NAME AND PARAMETER SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RRI_Read

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1 : FILE READING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxt = lasth * 3600 / dt;

% read parameters
params = readvars(demfile, 'Range', 'A1:C6');
params = cellfun(@split, params, 'UniformOutput', false);

NX = str2double(cell2mat(params{1}(2)));
NY = str2double(cell2mat(params{2}(2)));
XLLCORNER = str2double(cell2mat(params{3}(2)));
YLLCORNER = str2double(cell2mat(params{4}(2)));
CELLSIZE  = str2double(cell2mat(params{5}(2)));

% 注意：以下は転値している
zs      = readGisFile(demfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)'; % 地表の標高
flowAcc = readGisFile(accfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)'; % 集水面積
flowDir = readGisFile(dirfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)'; % 流向

land = ones(NX, NY);   %土地利用
if land_switch == 1
    land = readGisFile(landfile,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)';
end
disp(['num_of_landuse : ' , num2str(num_of_landuse)]);
land(land <= 0) = num_of_landuse;
if sum(unique(land) > num_of_landuse) > 0
    error("Invalid landuse number included.")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2 : CALC PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dx, dy] = metricCellsize(XLLCORNER, YLLCORNER, NX, NY, CELLSIZE, utm);

leveeHeight  = zeros(NX,NY); % 堤防の高さ
len_riv = zeros(NX,NY);      % セルの対角線長を河川領域にのみ定義したもの

len = sqrt(dx * dy);         % セルの対角線長
cellarea = dx * dy;

riv = (flowAcc > rivThresh); % 1:river, 0:others
width = (width_param_c * (flowAcc*dx*dy*1d-6) .^ width_param_s) .* riv; % 川幅
depth = (depth_param_c * (flowAcc*dx*dy*1d-6) .^ depth_param_s) .* riv; % 河道深さ
leveeHeight(riv & (flowAcc > height_limit_param)) = height_param;       % 堤防高

% % 川幅, 河道深さ, 堤防高が与えられたとき
% if rivfile_swtch >= 1
%     
% end
len_riv(riv) = len;

% sec_switch, sec_length_switch に関するオプション

area_ratio = width .* len_riv / cellarea;
zb_riv = zeros(size(zs)); % 河川の不透水層の標高
zb_riv(riv) = zs(riv) - depth(riv);
zb = zs - soildepth(land); % 不透水層の標高

domain = zeros(NX,NY);    % 0:解析範囲外，1:範囲内, 2:河口，端
domain(zs>-100) = 1;
domain(flowDir == 0 | flowDir == -1) = 2; % 0 or -1 is outlet point
numOfCell = nnz(domain); % count non-zero cell number
disp(['num_of_cell : ', num2str(numOfCell)])
disp(['total area [km^2] : ', num2str(numOfCell * cellarea / (10^6))])

%%% ----------------------- call riv_idx_setting ---------------------- %%%
riv_count = nnz((domain > 0) & riv); % number of river cell

domAndRiv             = (domain>0) & riv;         % 1:rivセル， 0:範囲外 or sloのみ                        
[riv_idx2i,riv_idx2j] = find(domAndRiv');         % rivセルのx,y座標
riv_ij2idx            = zeros(size(domAndRiv));   % rivセルの座標
riv_ij2idx(domAndRiv) = 1:sum(domAndRiv, 'all');  % rivセルの番号を振る  
domain_riv_idx        = domain(domAndRiv);        % domain
width_idx             = width(domAndRiv);         % 川幅
depth_idx             = depth(domAndRiv);         % 河道深さ
height_idx            = leveeHeight(domAndRiv);   % 堤防高
area_ratio_idx        = area_ratio(domAndRiv);    % セルにおける河川の割合
zb_riv_idx            = zb_riv(domAndRiv);        % 不透水層の標高
dif_riv_idx           = dif(land(domAndRiv))';    % 1:拡散波近似，0:kinematic
% RRI_Input.txtで指定する土地利用ごとのパラメータを縦ベクトルに統一するため転置する(Sloも同様)
% sec_map_idx           = sec_map(domAndRiv);     % 
len_riv_idx           = len_riv(domAndRiv);       % 河川の長さ

[dis_riv_idx, down_riv_idx] = rivIdxSetting(riv_count,NX,NY,domain,riv,flowDir,dx,dy,riv_ij2idx);

%%% ------------------------ call slo_idx_setting --------------------- %%%
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

dif_slo_idx      = dif(land(domAndSlo))';         % 1:拡散波近似，0:kinematic
ns_slo_idx       = ns_slope(land(domAndSlo))';    % 粗度係数
soildepth_idx    = soildepth(land(domAndSlo))';   % 土層厚
gammaa_idx       = gammaa(land(domAndSlo))';      % 間隙率
ksv_idx          = ksv(land(domAndSlo))';         % 鉛直方向の透水係数
faif_idx         = faif(land(domAndSlo))';        % 湿潤前線における吸引圧
infilt_limit_idx = infilt_limit(land(domAndSlo))';% 浸透する水の上限　(= soildepth * gammaa)

ka_idx           = ka(land(domAndSlo))';          % 水平方向の透水係数
gammam_idx       = gammam(land(domAndSlo))';      % 不飽和間隙率
beta_idx         = beta(land(domAndSlo))';        % パラメータ
da_idx           = da(land(domAndSlo))';          % 間隙 (= soildepth * gammaa)
dm_idx           = dm(land(domAndSlo))';          % 不飽和間隙 (= soildepth * gammam)
% ksg_idx          = ksg(land(domAndSlo))';         %
% gammag_idx       = gammag(land(domAndSlo))';      %
% kg0_idx          = kg0(land(domAndSlo))';         %
% fpg_idx          = fpg(land(domAndSlo))';         %
% rgl_idx          = rgl(land(domAndSlo))';         %

[down_slo_idx,dis_slo_idx,len_slo_idx,down_slo_1d_idx,dis_slo_1d_idx,len_slo_1d_idx] ...
    = sloIdxSetting(slo_count,NX,NY,domain,flowDir,dx,dy,slo_ij2idx,i4,eightFlowDir);

if eightFlowDir == 0; lmax = 2; else; lmax = 4; end

%%%-------------------------- call dam_read ----------------------------%%%

%%%-------------------------- initial condition ------------------------%%%

hs = zeros(NX,NY); % slope(not river)の水深
hr = zeros(NX,NY); % river の水深
hg = zeros(NX,NY); % 地下水位（オプション）
gampt_ff = zeros(NX,NY); % accumulated filtration depth(今浸透している量)[m]
gampt_f  = zeros(NX,NY); % 単位時間あたりに浸透できる流量 (infiltration capacity)[m]
qrs = zeros(NX,NY);      % slopeとriver間の水のやりとりの量

hr(~riv) = -0.1;
hs(~domAndSlo) = -0.1;
hg(:,:) = -0.1;

%%%-------------------------- boundary condition -----------------------%%%

%%%-------------------------- div file ---------------------------------%%%

df = readmatrix(location_file);
for I = 1:size(df,1)
    hydro_i(I) = df(I, 2); % Y座標(列番号)
    hydro_j(I) = df(I, 3); % X座標(行番号)
end

%%%-------------------------- array initialization ---------------------%%%

qr_idx = zeros(riv_count, 1);        % 

hr_idx = zeros(riv_count, 1);        % rivセルの水深
vr_idx = zeros(riv_count, 1);        % rivセルの流量

qs_ave_temp_idx = zeros(i4, slo_count);
qs_ave_idx = zeros(i4, slo_count);
qs_idx = zeros(i4, slo_count);
hs_idx = zeros(slo_count, 1);        % sloセルの水深
qp_t_idx = zeros(slo_count, 1);      % sloセルの降水量

gampt_ff_idx = zeros(slo_count, 1);  % accumulated infiltration depth [m]
gampt_f_idx = zeros(slo_count, 1);   % infiltration capacity [m/s]

rain_i = zeros(NY, 1);               % 該当するrainのy座標（横方向）
rain_j = zeros(NX, 1);               % 該当するrainのx座標（縦方向）


%%%-------------------------- gw initial setting -----------------------%%%
% if init_gw_switch ~= 1
%     hg_idx = zeros(slo_count, 1); %  call hg_init( hg_idx )
%  call sub_slo_idx2ij( hg_idx, hg )
% end

%%%-------------------------- initial srtorage calculation -------------%%%
% rain_sum = 0;
% aevp_sum = 0;
% pevp_sum = 0;
% sout = 0;
% si = 0;
% sg = 0;
% open( 1000, file = outfile_storage )
% call storage_calc(hs, hr, hg, ss, sr, si, sg)
% sinit = ss + sr + si + sg
% write(1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
%   (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg

%%%-------------------------- reading rainfall data --------------------%%% 
df = readmatrix(rainfile,'NumHeaderLines',0);
nx_rain = df(1,2);
ny_rain = df(1,3);
tt_max_rain = size(df,1) / (ny_rain + 1);   % fortran + 1
disp([tt_max_rain,nx_rain, ny_rain])        

qp = zeros(nx_rain, ny_rain, tt_max_rain);  % 順番変更
qp_t = zeros(NX,NY);

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

%%%-------------------------- reading evp data -------------------------%%%


% call sub_slo_ij2idx( hs, hs_idx ) 
% call sub_riv_ij2idx( hr, hr_idx ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3 : CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rain_sum = 0; % accumulated rainfall
aevp_sum = 0; % accumulated evaporation
sout = 0;     % accumulated runoff from river mouth

out_dt = maxt / outnum; % output outnum times per run.
out_dt = max(1, out_dt);
out_next = round(out_dt);
TT = 0;

% preparation for RK
load_paramRK
funcr_vr4RK = @(vr_idx, qr_idx, hr_idx) ...
        Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count,...
        domain_riv_idx, zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);
funcr_hs4RK = @(qs_idx, hs_idx, qp_t_idx) ...
        Funcs(qs_idx, hs_idx, qp_t_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx,...
        dis_slo_idx, dis_slo_1d_idx, down_slo_idx, down_slo_1d_idx, len_slo_idx, len_slo_1d_idx,...
        da_idx, dm_idx, beta_idx, dif_slo_idx, soildepth_idx, gammaa_idx, lmax, cellarea);

errmax = 10;
maxt = 36; % practice
for T = 1:maxt
    if mod(T,1)==0; fprintf("%d/%d\n", T, maxt); end

%%%-------------------------- RIVER CALCULATION ------------------------%%%
    
    time = (T - 1) * dt; 
    % time step is initially set to be "dt_riv"
    ddt = dt_riv; % 本来のRRIではRK計算のエラーが大きくなるとddtを小さくすることがある
    ddt_chk_riv = dt_riv;
    
    qr_ave = zeros(NX,NY);
    qr_ave_idx = zeros(riv_count, 1);
    % if dam_switch == 1; dam_vol_temp(:) = 0;end
    
    hr_idx = hr(domAndRiv); % hr_idx = sub_riv_ij2idx(hr, hr_idx, riv_count, riv_idx2i, riv_idx2j);
    
    for K = 1:riv_count
        vr_idx(K) = hr2vr(hr_idx(K), K, cellarea, area_ratio_idx(K));  % Kは不要になりそう
    end
    
    while true
        if rivThresh < 0; break; end        % no river
        % "time + ddt" should be less than "t * dt"
        if time + ddt > T * dt; ddt = T * dt - time; end
        
        %%%% boundary condition
        
        qr_ave_temp_idx = zeros(riv_count, 1);
        % % % Adaptive Runge-Kutta
        [vr_err, vr_temp] = ...
            adaptiveRKvr(ddt, qr_ave_temp_idx, ParamRK, vr_idx, qr_idx, hr_idx, funcr_vr4RK);
    
        hr_err = vr_err ./ (cellarea * area_ratio_idx);
        hr_err(domain_riv_idx==0) = 0;
        [errmax, errmax_loc] =  max(hr_err, [], 'all');

        errmax = errmax / eps;
        % disp([ddt,errmax])
        
        if errmax > 1 && ddt > ddt_min_riv
            ddt = max( safety * ddt * (errmax ^ pshrnk), 0.5 * ddt );
            ddt = max( ddt, ddt_min_riv ); 
            ddt_chk_riv = ddt;
            disp(["shrink (riv): ", ddt, errmax, slo_idx2i(errmax_loc), slo_idx2j(errmax_loc)])
            if ddt==0; error('stepsize underflow'); end  
        else
            if ddt == ddt_min_riv
                kr6 = funcr_vr4RK(vr_idx, qr_idx, hr_idx);
                qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt * 6;
            end
            if time + ddt > T * dt; ddt = T * dt - time; end
            time = time + ddt;
            vr_idx = vr_temp;
            qr_ave_idx = qr_ave_idx + qr_ave_temp_idx;             
        end
        if time >= T * dt; break; end
    end
    qr_ave_idx = qr_ave_idx / dt / 6; 

    for K = 1:riv_count
        hr_idx(K) = vr2hr(vr_idx(K), K, cellarea, area_ratio_idx(K)); 
    end    

    hr(domAndRiv) = hr_idx; 
    qr_ave(domAndRiv) = qr_ave_idx;

%%%-------------------------- SLOPE CALCULATION ------------------------%%%

    time = (T - 1) * dt;
    ddt = dt;
    ddt_chk_slo = dt;

    qs_ave = zeros(i4,NX,NY);
    qs_ave_idx = zeros(i4, slo_count);

    hs_idx = hs(domAndSlo);
    gampt_ff_idx = gampt_ff(domAndSlo);

    while time < T * dt || (errmax > 1 && ddt > ddt_min_slo)
        if time + ddt > T * dt; ddt = T * dt - time; end
    
        % rainfall 
        itemp = -1;
        for jtemp = 1:tt_max_rain
            if t_rain(jtemp) < (time + ddt) && (time + ddt) <= t_rain(jtemp+1); itemp = jtemp; end
        end
        for J = 1:NY
            if itemp<= 0; continue; end  % バグ回避のため追加
            if rain_i(J) < 1 || rain_i(J) > ny_rain; continue; end
            for I = 1:NX
                if rain_j(I) < 1 || rain_j(I) > nx_rain; continue; end
                qp_t(I, J) = qp(rain_j(I), rain_i(J), itemp);
            end
        end
        qp_t_idx = qp_t(domAndSlo); 
        
        %%%% boundary condition
    
        qs_ave_temp_idx = zeros(i4, slo_count);
        % % % Adaptive Runge-Kutta
        [hs_err, hs_temp] = ...
            adaptive_RKhs(ddt, qs_ave_temp_idx, ParamRK, qs_idx, hs_idx, qp_t_idx, funcr_hs4RK);
    
        hs_err(domain_slo_idx==0) = 0;
        [errmax, errmax_loc] =  max(hs_err, [], 'all');
        errmax = errmax / eps;
    
        if errmax > 1 && ddt > ddt_min_slo
            ddt = max( safety * ddt * (errmax ^ pshrnk), 0.5 * ddt );
            ddt = max( ddt, ddt_min_slo ); 
            ddt_chk_slo = ddt;
            disp(["shrink (slo): ", ddt, errmax, slo_idx2i(errmax_loc), slo_idx2j(errmax_loc)])
            if ddt==0; error('stepsize underflow'); end  
        else
            if time + ddt > T * dt; ddt = T * dt - time; end
            time = time + ddt;
            hs_idx = hs_temp;
            qs_ave_idx = qs_ave_idx + qs_ave_temp_idx;
        end

        % cumulative rainfall
        rain_sum = sum(qp_t(domAndSlo)) * cellarea * numOfCell * ddt;
        
        if time >= T * dt; break; end
    end
    qs_ave_idx = qs_ave_idx / dt / 6;

%%%-------------------------- GW CALCULATION ---------------------------%%%


%%%-------------------------- GW Exfiltration  -------------------------%%%


%%%-------------------------- Evapotranspiration  ----------------------%%%

    
    hs(domAndSlo) = hs_idx;


%%%-------------------------- LEVEE BREAK  -----------------------------%%%
% 使用しない

%%%-------------------------- RIVER-SLOPE INTERACTION  -----------------%%%
    
    [hr, hs] = funcrs(hr, hs, NX, NY, domain, ...
        riv, leveeHeight, depth, riv_ij2idx, len_riv_idx, dt, cellarea, area_ratio_idx, width_idx);
    hr_idx = hr(domAndRiv);
    hs_idx = hs(domAndSlo);

%%%-------------------------- INFILTRATION (Green Ampt)  ---------------%%%
    
    hs_idx = infilt(hs_idx, gampt_f_idx, gampt_ff_idx, ksv_idx, faif_idx, gammaa_idx, infilt_limit_idx, dt, slo_count);
    hs(domAndSlo) = hs_idx;
    gampt_f(domAndSlo) = gampt_f_idx;
    gampt_ff(domAndSlo) = gampt_ff_idx;

%%%-------------------------- SET WATER DEPTH 0 AT DOMAIN = 2  ---------%%%
    
    % sout = sout + sum(hs(domain==2)) * cellarea;
    % hs(domain==2) = 0;
    for I = 1:NX
        for J = 1:NY
            if domain(I,J) == 2
                sout = sout + hs(I,J) * cellarea;
                hs(I,J) = 0;
                if riv(I,J) == 1
                    K = riv_ij2idx(I,J);
                    vr_out = hr2vr(hr(I,J), K, cellarea, area_ratio_idx(K));
                    sout = sout + vr_out;
                    hr(I,J) = 0;
                end
            end
        end
    end

    hr_idx = hr(domAndRiv);
    hs_idx = hs(domAndSlo);

    [hr_max, hr_max_loc] = max(hr,[],'all');
    [hs_max, hs_max_loc] = max(hs,[],'all');
    disp(["max hr : ", hr_max, "loc : ", slo_idx2i(hr_max_loc), slo_idx2j(hr_max_loc) ])
    disp(["max hs : ", hs_max, "loc : ", slo_idx2i(hs_max_loc), slo_idx2j(hs_max_loc) ])
%%%-------------------------- OUTPUT -----------------------------------%%%


end

toc