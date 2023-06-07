function [hs_err, hs_temp, qs_ave_temp_idx] = ...
    adaptive_RKhs(ddt, qs_ave_temp_idx, ParamRK, qs_idx, hs_idx, qp_t_idx, funcr_hs4RK)
    
% % Adaptive Runge-Kutta
    % % (1)
    [fs, qs_idx] = funcr_hs4RK(qs_idx, hs_idx, qp_t_idx);
    hs_temp = hs_idx + ParamRK.b21 * ddt * fs;
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;

    % % (2)
    [ks2, qs_idx] = funcr_hs4RK(qs_idx, hs_temp, qp_t_idx);
    hs_temp = hs_idx + ddt * (ParamRK.b31 * fs + ParamRK.b32 * ks2);
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;
    % 
    % % (3)
    [ks3, qs_idx] = funcr_hs4RK(qs_idx, hs_temp, qp_t_idx);
    hs_temp = hs_idx + ddt * (ParamRK.b41 * fs + ParamRK.b42 * ks2 + ParamRK.b43 * ks3);
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;
    % 
    % % (4)
    [ks4, qs_idx] = funcr_hs4RK(qs_idx, hs_temp, qp_t_idx);
    hs_temp = hs_idx + ddt * (ParamRK.b51 * fs + ParamRK.b52 * ks2 + ParamRK.b53 * ks3 + ParamRK.b54 * ks4);
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;
    % 
    % % (5)
    [ks5, qs_idx] = funcr_hs4RK(qs_idx, hs_temp, qp_t_idx);
    hs_temp = hs_idx + ddt * (ParamRK.b61 * fs + ParamRK.b62 * ks2 + ParamRK.b63 * ks3 + ParamRK.b64 * ks4 + ParamRK.b65 * ks5);
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;
    % 
    % % (6)
    [ks6, qs_idx] = funcr_hs4RK(qs_idx, hs_temp, qp_t_idx);
    hs_temp = hs_idx + ddt * (ParamRK.c1 * fs + ParamRK.c3 * ks3 + ParamRK.c4 * ks4 + ParamRK.c6 * ks6);
    hs_temp(hs_temp<0) = 0;
    qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt;
    
    hs_err = ddt * (ParamRK.dc1 * fs + ParamRK.dc3 * ks3 + ParamRK.dc4 * ks4 + ParamRK.dc5 * ks5 + ParamRK.dc6 * ks6);
end