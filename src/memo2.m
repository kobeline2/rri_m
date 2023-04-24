%% ルンゲクッタのパラメータ
a2 = 0.2;
a3 = 0.3;
a4 = 0.6;
a5 = 1.0;
a6 = 0.875;
b21 = 0.2;
b31 = 3/40;
b32 = 9/40;
b41 = 0.3;
b42 = -0.9;
b43 = 1.2;
b51 = -11/54;
b52 = 2.5;
b53 = -70/27;
b54 = 35/27;
b61 = 1631/55296;
b62 = 175/512;
b63 = 575/13824;
b64 = 44275/110592;
b65 = 253/4096;
c1 = 37/378;
c3 = 250/621;
c4 = 125/594;
c6 = 512/1771;
dc1 = c1 - 2825/27648;
dc3 = c3 - 18575/48384;
dc4 = c4 - 13525/55296;
dc5 = -277/14336;
dc6 = c6 - 0.25;

%% その他のパラメータ
safety = 0.9;
pgrow = -0.2;
pshrnk = -0.25;
errcon = 1.89e-4;
eps = 0.01;
ddt_min_riv = 0.1;
ddt_min_slo = 1.0;

% %% calculation
% % Adaptive Runge-Kutta
% % (1)
% % call funcr( vr_idx, fr, qr_idx, t )   
% vr_temp = vr_idx + b21 * ddt * fr
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (2)
% % call funcr( vr_temp, kr2, qr_idx, t )
% vr_temp = vr_idx + ddt * (b31 * fr + b32 * kr2)
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (3)
% % call funcr( vr_temp, kr3, qr_idx, t )
% vr_temp = vr_idx + ddt * (b41 * fr + b42 * kr2 + b43 * kr3)
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (4)
% % call funcr( vr_temp, kr4, qr_idx, t )
% vr_temp = vr_idx + ddt * (b51 * fr + b52 * kr2 + b53 * kr3 + b54 * kr4)
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (5)
% % call funcr( vr_temp, kr5, qr_idx, t )
% vr_temp = vr_idx + ddt * (b61 * fr + b62 * kr2 + b63 * kr3 + b64 * kr4 + b65 * kr5)
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (6)
% % call funcr( vr_temp, kr6, qr_idx, t )
% vr_temp = vr_idx + ddt * (c1 * fr + c3 * kr3 + c4 * kr4 + c6 * kr6)
% where(vr_temp .lt. 0) vr_temp = 0.d0
% qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt
% 
% % (e)
% vr_err = ddt * (dc1 * fr + dc3 * kr3 + dc4 * kr4 + dc5 * kr5 + dc6 * kr6);
% 
% hr_err(:) = vr_err(:) / (area * area_ratio_idx(:));
