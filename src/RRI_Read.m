%% RRI_Read.f90

lines = readvars('RRI_Input.txt','Range','1:100'); % 修正したい
lines = cellfun(@split, lines, 'UniformOutput', false);

format_version = lines{1}{:};
if format_version ~= "RRI_Input_Format_Ver1_4_2"
    error("This RRI model requires RRI_Input_Format_Ver1_4_2")
end

rainfile = lines{2}{:}; disp(['rainfile : ', rainfile])
demfile  = lines{3}{:}; disp(['demfile : ' , demfile])
accfile  = lines{4}{:}; disp(['accfile : ' , accfile])
dirfile  = lines{5}{:}; disp(['dirfile : ' , dirfile])
disp(' ')

utm             = str_split_cell2mat(lines{6}, 'utm');
eightFlowDir    = str_split_cell2mat(lines{7}, 'eight_dir');
lasth           = str_split_cell2mat(lines{8}, 'lasth');
dt              = str_split_cell2mat(lines{9}, 'dt');
dt_riv          = str_split_cell2mat(lines{10}, 'dt_riv');
outnum          = str_split_cell2mat(lines{11}, 'XLLCORNER_RAIN');
XLLCORNER_RAIN  = str_split_cell2mat(lines{12}, 'XLLCORNER_RAIN');
YLLCORNER_RAIN  = str_split_cell2mat(lines{13}, 'YLLCORNER_RAIN');
CELLSIZE_RAIN_X = str_split_cell2mat(lines{14}(1), 'CELLSIZE_RAIN_X');
CELLSIZE_RAIN_Y = str_split_cell2mat(lines{14}(2), 'CELLSIZE_RAIN_Y');
disp(' ')

ns_river       = str_split_cell2mat(lines{15}, 'ns_river');
num_of_landuse = str_split_cell2mat(lines{16}, 'num_of_landuse');

dif       = str_split_cell2mat(lines{17}, 'dif');
ns_slope  = str_split_cell2mat(lines{18}, 'ns_slope');
soildepth = str_split_cell2mat(lines{19}, 'soildepth');
gammaa    = str_split_cell2mat(lines{20}, 'gammaa');
ksv       = str_split_cell2mat(lines{21}, 'ksv');
faif      = str_split_cell2mat(lines{22}, 'faif');
ka        = str_split_cell2mat(lines{23}, 'ka');
gammam    = str_split_cell2mat(lines{24}, 'gammam');
beta      = str_split_cell2mat(lines{25}, 'beta');
disp(' ')

ksg    = str_split_cell2mat(lines{26}, 'ksg');
gammag = str_split_cell2mat(lines{27}, 'gammag');
kg0    = str_split_cell2mat(lines{28}, 'kg0');
fpg    = str_split_cell2mat(lines{29}, 'fpg');
rgl    = str_split_cell2mat(lines{30}, 'rgl');
disp(' ')



rivfile_switch = str_split_cell2mat(lines{38}, '');
if rivfile_switch == 0
    rivThresh          = str_split_cell2mat(lines{31}, 'rivThresh');
    width_param_c      = str_split_cell2mat(lines{32}, 'width_param_c');
    width_param_s      = str_split_cell2mat(lines{33}, 'width_param_s');
    depth_param_c      = str_split_cell2mat(lines{34}, 'depth_param_c');
    depth_param_s      = str_split_cell2mat(lines{35}, 'depth_param_s');
    height_param       = str_split_cell2mat(lines{36}, 'height_param');
    height_limit_param = str_split_cell2mat(lines{37}, 'height_limit_param');
else
    widthfile  = lines{39}{:}; disp(['widthfile : '  , dirfile])
    depthfile  = lines{40}{:}; disp(['depthfile : '  , dirfile])
    heightfile = lines{41}{:}; disp(['heightfile : ' , dirfile])
end
disp(' ')

% init_slo_switch      = str_split_cell2mat(lines{42}(1), '');
% init_riv_switch      = str_split_cell2mat(lines{42}(2), '');
% init_gw_switch       = str_split_cell2mat(lines{42}(3), '');
% init_gampt_ff_switch = str_split_cell2mat(lines{42}(4), '');
% initfile_slo         = lines{43}{:};
% initfile_riv         = lines{44}{:};
% initfile_gw          = lines{45}{:};
% initfile_gampt_ff    = lines{46}{:};
% if init_slo_switch ~= 0, disp(['initfile_slo : ', initfile_slo]), end
% if init_riv_switch ~= 0, disp(['initfile_riv : ', initfile_riv]), end
% if init_gw_switch  ~= 0, disp(['initfile_gw : ', initfile_gw]), end
% if init_gampt_ff_switch ~= 0, disp(['initfile_gampt_ff : ', initfile_slo]), end
% 
% bound_slo_wlev_switch = str_split_cell2mat(lines{47}(1), '');
% bound_riv_wlev_switch = str_split_cell2mat(lines{47}(2), '');
% boundfile_slo_wlev = lines{48}{:};
% boundfile_riv_wlev = lines{49}{:};
% if bound_slo_wlev_switch ~= 0, disp(['boundfile_slo_wlev : ', boundfile_slo_wlev]), end
% if bound_riv_wlev_switch ~= 0, disp(['boundfile_riv_wlev : ', boundfile_riv_wlev]), end
% 
% bound_slo_disc_switch = str_split_cell2mat(lines{50}(1), '');
% bound_riv_disc_switch = str_split_cell2mat(lines{50}(2), '');
% boundfile_slo_disc = lines{51}{:};
% boundfile_riv_disc = lines{52}{:};
% if bound_slo_disc_switch ~= 0, disp(['boundfile_slo_disc : ', boundfile_slo_disc]), end
% if bound_riv_disc_switch ~= 0, disp(['boundfile_riv_disc : ', boundfile_riv_disc]), end
% disp(' ')

land_switch = str_split_cell2mat(lines{53}, '');
landfile    = lines{54}{:};
if land_switch == 1, disp(['landfile : ', landfile]), end

% dam_switch = str_split_cell2mat(lines{55}, '');
% damfile    = lines{56}{:};
% if dam_switch == 1, disp(['damfile : ', damfile]), end
% 
% div_switch = str_split_cell2mat(lines{57}, '');
% divfile    = lines{58}{:};
% if div_switch == 1, disp(['divfile : ', divfile]), end
% 
% evp_switch = str_split_cell2mat(lines{59}, '');
% if evp_switch ~= 0
%     evpfile = lines{60}{:};
%     xllcorner_evp  = str_split_cell2mat(lines{61}, '');
%     yllcorner_evp  = str_split_cell2mat(lines{62}, '');
%     cellsize_evp_x = str_split_cell2mat(lines{63}(1), '');
%     cellsize_evp_y = str_split_cell2mat(lines{63}(2), '');
% end
% 
% sec_length_switch = str_split_cell2mat(lines{64}, '');
% sec_length_file   = lines{65}{:};
% if sec_length_switch == 1, disp(['sec_length_file : ', sec_length_file]), end
% 
% sec_switch   = str_split_cell2mat(lines{66}, '');
% sec_map_file = lines{67}{:};
% sec_file = cell2mat(lines{68});
% if sec_switch == 1, disp(['sec_map_file : ', sec_map_file]), end
% if sec_switch == 1, disp(['sec_file : ', sec_file]), end
% disp(' ')

% outswitch_hs       = str_split_cell2mat(lines{69}(1), '');
% outswitch_hr       = str_split_cell2mat(lines{69}(2), '');
% outswitch_hg       = str_split_cell2mat(lines{69}(3), '');
% outswitch_qr       = str_split_cell2mat(lines{69}(4), '');
% outswitch_qu       = str_split_cell2mat(lines{69}(5), '');
% outswitch_qv       = str_split_cell2mat(lines{69}(6), '');
% outswitch_gu       = str_split_cell2mat(lines{69}(7), '');
% outswitch_gv       = str_split_cell2mat(lines{69}(8), '');
% outswitch_gampt_ff = str_split_cell2mat(lines{69}(9), '');
% outswitch_storage  = str_split_cell2mat(lines{69}(10), '');
% outfile_hs         = lines{70}{:};
% outfile_hr         = lines{71}{:};
% outfile_hg         = lines{72}{:};
% outfile_qr         = lines{73}{:};
% outfile_qu         = lines{74}{:};
% outfile_qv         = lines{75}{:};
% outfile_gu         = lines{76}{:};
% outfile_gv         = lines{77}{:};
% outfile_gampt_ff   = lines{78}{:};
% outfile_storage    = lines{79}{:};
% if outswitch_hs       ~= 0, disp(['outfile_hs : '     , outfile_hs]), end
% if outswitch_hr       ~= 0, disp(['outfile_hr : '     , outfile_hr]), end
% if outswitch_hg       ~= 0, disp(['outfile_hg : '     , outfile_hg]), end
% if outswitch_qr       ~= 0, disp(['outfile_qr : '     , outfile_qr]), end
% if outswitch_qu       ~= 0, disp(['outfile_qu : '     , outfile_qu]), end
% if outswitch_qv       ~= 0, disp(['outfile_qv : '     , outfile_qv]), end
% if outswitch_gu       ~= 0, disp(['outfile_gu: '      , outfile_gu]), end
% if outswitch_gv       ~= 0, disp(['outfile_gv : '     , outfile_gv]), end
% if outswitch_gampt_ff ~= 0, disp(['outfile_gampt_f : ', outfile_gampt_f]), end
% if outswitch_storage  ~= 0, disp(['outfile_storage : ', outfile_storage]), end
% disp(' ')

hydro_switch  = str_split_cell2mat(lines{80}, '');
location_file = lines{81}{:};
if hydro_switch == 1, disp(['location_file : ', location_file]), end

for I = 1:num_of_landuse
    if (ksv(I) > 0) && (ka(I) > 0), error("Error: both ksv and ka are non-zero."), end
    if gammam(I) > gammaa(I), error("Error: gammam must be smaller than gammaa."), end
end

infilt_limit = zeros(1, num_of_landuse);
da           = zeros(1, num_of_landuse);
dm           = zeros(1, num_of_landuse);
for I = 1:num_of_landuse
    if (soildepth(I) > 0) && (ksv(I) > 0), infilt_limit(I) = soildepth(I) * gammaa(I); end
    if (soildepth(I) > 0) && (ka(I) > 0), da(I) = soildepth(I) * gammaa(I); end
    if (soildepth(I) > 0) && (ksv(I) > 0) && (gammam(I) > 0), dm(I) = soildepth(I) * gammam(I); end
end

%% functions
function out = str_split_cell2mat(lines, s)
    out = split(lines, ' ');
    out = cellfun(@str2double, out);
    out = reshape(out, [], length(out));
    if ~strcmp(s, '')
        disp([s, ' : ', num2str(out)])
    end
end