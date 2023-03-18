%% RRI_Read.f90
tic

lines = readvars('RRI_Input.txt','Range','1:100'); % 修正したい
lines = cellfun(@split, lines, 'UniformOutput', false);

format_version = cell2mat(lines{1});
if format_version ~= "RRI_Input_Format_Ver1_4_2"
    error("This RRI model requires RRI_Input_Format_Ver1_4_2")
end

rainfile = cell2mat(lines{2});
demfile = cell2mat(lines{3});
accfile = cell2mat(lines{4});
dirfile = cell2mat(lines{5});

disp(['rainfile : ', rainfile])
disp(['demfile : ', demfile])
disp(['accfile : ', accfile])
disp(['dirfile : ', dirfile])
disp(' ')

utm = str2double(cell2mat(lines{6}));
eightFlowDir = str2double(cell2mat(lines{7}));
lasth = str2double(cell2mat(lines{8}));
dt = str2double(cell2mat(lines{9}));
dt_riv = str2double(cell2mat(lines{10}));
outnum = str2double(cell2mat(lines{11}));
XLLCORNER_RAIN = str2double(cell2mat(lines{12}));
YLLCORNER_RAIN = str2double(cell2mat(lines{13}));
CELLSIZE_RAIN_X = str2double(cell2mat(lines{14}(1)));
CELLSIZE_RAIN_Y = str2double(cell2mat(lines{14}(2)));

disp(['utm : ' , num2str(utm)]);
disp(['eight_dir : ' , num2str(eightFlowDir)]);
disp(['lasth : ' , num2str(lasth)]);
disp(['dt : ' , num2str(dt)]);
disp(['dt_riv : ' , num2str(dt_riv)]);
disp(['XLLCORNER_RAIN : ' , num2str(XLLCORNER_RAIN)]);
disp(['YLLCORNER_RAIN : ' , num2str(YLLCORNER_RAIN)]);
disp(['CELLSIZE_RAIN_X : ' , num2str(CELLSIZE_RAIN_X), ...
    ',  CELLSIZE_RAIN_Y : ' , num2str(CELLSIZE_RAIN_Y)]);
disp(' ')

ns_river = str2double(cell2mat(lines{15}));
num_of_landuse = str2double(cell2mat(lines{16}));
disp(['ns_river : ' , num2str(ns_river)]);
disp(['num_of_landuse : ' , num2str(num_of_landuse)]);


dif       = str_split_cell2mat(lines{17});
ns_slope  = str_split_cell2mat(lines{18});
soildepth = str_split_cell2mat(lines{19});
gammaa    = str_split_cell2mat(lines{20});
ksv       = str_split_cell2mat(lines{21});
faif      = str_split_cell2mat(lines{22});
ka        = str_split_cell2mat(lines{23});
gammam    = str_split_cell2mat(lines{24});
beta      = str_split_cell2mat(lines{25});


disp(['dif : ', num2str(dif)])
disp(['ns_slope : ', num2str(ns_slope)])
disp(['soildepth : ', num2str(soildepth)])
disp(['gammaa : ', num2str(gammaa)])
disp(['ksv : ', num2str(ksv)])
disp(['faif : ', num2str(faif)])
disp(['ka : ', num2str(ka)])
disp(['gammam : ', num2str(gammam)])
disp(['beta : ', num2str(beta)])
disp(' ')

for I = 1:num_of_landuse   % 土地利用の数だけなら事前割り当て不要？
    ksg(I) = str2double(cell2mat(lines{26}(I)));
    gammag(I) = str2double(cell2mat(lines{27}(I)));
    kg0(I) = str2double(cell2mat(lines{28}(I)));
    fpg(I) = str2double(cell2mat(lines{29}(I)));
    rgl(I) = str2double(cell2mat(lines{30}(I)));
end

disp(['ksg : ', num2str(ksg)])
disp(['gammag : ', num2str(gammag)])
disp(['kg0 : ', num2str(kg0)])
disp(['fpg : ', num2str(fpg)])
disp(['rgl : ', num2str(rgl)])
disp(' ')

rivThresh = str2double(cell2mat(lines{31}));
width_param_c = str2double(cell2mat(lines{32}));
width_param_s = str2double(cell2mat(lines{33}));
depth_param_c = str2double(cell2mat(lines{34}));
depth_param_s = str2double(cell2mat(lines{35}));
height_param = str2double(cell2mat(lines{36}));
height_limit_param = str2double(cell2mat(lines{37}));

rivfile_switch = str2double(cell2mat(lines{38}));
widthfile = cell2mat(lines{39});
depthfile = cell2mat(lines{40});
heightfile = cell2mat(lines{41});

if rivfile_switch == 0
    disp(['riv_thresh : ', num2str(rivThresh)])
    disp(['width_param_c : ', num2str(width_param_c)])
    disp(['width_param_s : ', num2str(width_param_s)])
    disp(['depth_param_c : ', num2str(depth_param_c)])
    disp(['depth_param_s : ', num2str(depth_param_s)])
    disp(['height_param : ', num2str(height_param)])
    disp(['height_limit_param : ', num2str(height_limit_param)])
else
    disp(['widthfile : ', widthfile])
    disp(['depthfile : ', depthfile])
    disp(['heightfile : ', heightfile])
end
disp(' ')

init_slo_switch = str2double(cell2mat(lines{42}(1)));
init_riv_switch = str2double(cell2mat(lines{42}(2)));
init_gw_switch = str2double(cell2mat(lines{42}(3)));
init_gampt_ff_switch = str2double(cell2mat(lines{42}(4)));
initfile_slo = cell2mat(lines{43});
initfile_riv = cell2mat(lines{44});
initfile_gw = cell2mat(lines{45});
initfile_gampt_ff = cell2mat(lines{46});
if init_slo_switch ~= 0, disp(['initfile_slo : ', initfile_slo]), end
if init_riv_switch ~= 0, disp(['initfile_riv : ', initfile_riv]), end
if init_gw_switch ~= 0, disp(['initfile_gw : ', initfile_gw]), end
if init_gampt_ff_switch ~= 0, disp(['initfile_gampt_ff : ', initfile_slo]), end

bound_slo_wlev_switch = str2double(cell2mat(lines{47}(1)));
bound_riv_wlev_switch = str2double(cell2mat(lines{47}(2)));
boundfile_slo_wlev = cell2mat(lines{48});
boundfile_riv_wlev = cell2mat(lines{49});
if bound_slo_wlev_switch ~= 0, disp(['boundfile_slo_wlev : ', boundfile_slo_wlev]), end
if bound_riv_wlev_switch ~= 0, disp(['boundfile_riv_wlev : ', boundfile_riv_wlev]), end

bound_slo_disc_switch = str2double(cell2mat(lines{50}(1)));
bound_riv_disc_switch = str2double(cell2mat(lines{50}(2)));
boundfile_slo_disc = cell2mat(lines{51});
boundfile_riv_disc = cell2mat(lines{52});
if bound_slo_disc_switch ~= 0, disp(['boundfile_slo_disc : ', boundfile_slo_disc]), end
if bound_riv_disc_switch ~= 0, disp(['boundfile_riv_disc : ', boundfile_riv_disc]), end
disp(' ')

land_switch = str2double(cell2mat(lines{53}));
landfile = cell2mat(lines{54});
if land_switch == 1, disp(['landfile : ', landfile]), end

dam_switch = str2double(cell2mat(lines{55}));
damfile = cell2mat(lines{56});
if dam_switch == 1, disp(['damfile : ', damfile]), end

div_switch = str2double(cell2mat(lines{57}));
divfile = cell2mat(lines{58});
if div_switch == 1, disp(['divfile : ', divfile]), end

evp_switch = str2double(cell2mat(lines{59}));
evpfile = cell2mat(lines{60});
xllcorner_evp = str2double(cell2mat(lines{61}));
yllcorner_evp = str2double(cell2mat(lines{62}));
cellsize_evp_x = str2double(cell2mat(lines{63}(1)));
cellsize_evp_y = str2double(cell2mat(lines{63}(2)));
if evp_switch ~= 0
    disp(['evpfile : ', evpfile])
    disp(['xllcorner_evp : ', xllcorner_evp])
    disp(['yllcorner_evp : ', yllcorner_evp])
    disp(['cellsize_evp_x  : ', cellsize_evp_x , 'cellsize_evp_y : ',cellsize_evp_y])
end

sec_length_switch = str2double(cell2mat(lines{64}));
sec_length_file = cell2mat(lines{65});
if sec_length_switch == 1, disp(['sec_length_file : ', sec_length_file]), end

sec_switch = str2double(cell2mat(lines{66}));
sec_map_file = cell2mat(lines{67});
sec_file = cell2mat(lines{68});
if sec_switch == 1, disp(['sec_map_file : ', sec_map_file]), end
if sec_switch == 1, disp(['sec_file : ', sec_file]), end
disp(' ')

outswitch_hs = str2double(cell2mat(lines{69}(1)));
outswitch_hr = str2double(cell2mat(lines{69}(2)));
outswitch_hg = str2double(cell2mat(lines{69}(3)));
outswitch_qr = str2double(cell2mat(lines{69}(4)));
outswitch_qu = str2double(cell2mat(lines{69}(5)));
outswitch_qv = str2double(cell2mat(lines{69}(6)));
outswitch_gu = str2double(cell2mat(lines{69}(7)));
outswitch_gv = str2double(cell2mat(lines{69}(8)));
outswitch_gampt_ff = str2double(cell2mat(lines{69}(9)));
outswitch_storage = str2double(cell2mat(lines{69}(10)));
outfile_hs = cell2mat(lines{70});
outfile_hr = cell2mat(lines{71});
outfile_hg = cell2mat(lines{72});
outfile_qr = cell2mat(lines{73});
outfile_qu = cell2mat(lines{74});
outfile_qv = cell2mat(lines{75});
outfile_gu = cell2mat(lines{76});
outfile_gv = cell2mat(lines{77});
outfile_gampt_ff = cell2mat(lines{78});
outfile_storage = cell2mat(lines{79});
if outswitch_hs ~= 0, disp(['outfile_hs : ', outfile_hs]), end
if outswitch_hr ~= 0, disp(['outfile_hr : ', outfile_hr]), end
if outswitch_hg ~= 0, disp(['outfile_hg : ', outfile_hg]), end
if outswitch_qr ~= 0, disp(['outfile_qr : ', outfile_qr]), end
if outswitch_qu ~= 0, disp(['outfile_qu : ', outfile_qu]), end
if outswitch_qv ~= 0, disp(['outfile_qv : ', outfile_qv]), end
if outswitch_gu ~= 0, disp(['outfile_gu: ', outfile_gu]), end
if outswitch_gv ~= 0, disp(['outfile_gv : ', outfile_gv]), end
if outswitch_gampt_ff ~= 0, disp(['outfile_gampt_f : ', outfile_gampt_f]), end
if outswitch_storage ~= 0, disp(['outfile_storage : ', outfile_storage]), end
disp(' ')

hydro_switch = str2double(cell2mat(lines{80}));
location_file = cell2mat(lines{81});
if hydro_switch == 1, disp(['location_file : ', location_file]), end

for I = 1:num_of_landuse
    if (ksv(I) > 0) && (ka(I) > 0), error("Error: both ksv and ka are non-zero."), end
    if gammam(I) > gammaa(I), error("Error: gammam must be smaller than gammaa."), end
end

infilt_limit = zeros(1,num_of_landuse);
da = zeros(1,num_of_landuse);
dm = zeros(1,num_of_landuse);
for I = 1:num_of_landuse
    if (soildepth(I) > 0) & (ksv(I) > 0), infilt_limit(I) = soildepth(I) * gammaa(I); end
    if (soildepth(I) > 0) & (ka(I) > 0), da(I) = soildepth(I) * gammaa(I); end
    if (soildepth(I) > 0) & (ksv(I) > 0) & (gammam(I) > 0), dm(I) = soildepth(I) * gammam(I), end
end

toc


%% functions
function out = str_split_cell2mat(line)
    out = split(lines, ' ');
    out = cellfun(@str2double, out);
    out = out';
end
