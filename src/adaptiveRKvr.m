function [vr_err, vr_temp, qr_ave_temp_idx] = ...
        adaptiveRKvr(ddt, qr_ave_temp_idx, ParamRK, vr_idx, qr_idx, hr_idx, funcr_vr4RK)
    
    % %% Adaptive Runge-Kutta
    % % (1)
    [fr, qr_idx] = funcr_vr4RK(vr_idx, qr_idx, hr_idx);   
    vr_temp = vr_idx + ParamRK.b21 * ddt * fr;
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % disp(max(vr_temp))

    % % (2)
    [kr2, qr_idx] = funcr_vr4RK(vr_temp, qr_idx, hr_idx); 
    vr_temp = vr_idx + ddt * (ParamRK.b31 * fr + ParamRK.b32 * kr2);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (3)
    [kr3, qr_idx] = funcr_vr4RK(vr_temp, qr_idx, hr_idx); 
    vr_temp = vr_idx + ddt * (ParamRK.b41 * fr + ParamRK.b42 * kr2 + ParamRK.b43 * kr3);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (4)
    [kr4, qr_idx] = funcr_vr4RK(vr_temp, qr_idx, hr_idx); 
    vr_temp = vr_idx + ddt * (ParamRK.b51 * fr + ParamRK.b52 * kr2 + ParamRK.b53 * kr3 + ParamRK.b54 * kr4);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (5)
    [kr5, qr_idx] = funcr_vr4RK(vr_temp, qr_idx, hr_idx); 
    vr_temp = vr_idx + ddt * (ParamRK.b61 * fr + ParamRK.b62 * kr2 + ParamRK.b63 * kr3 + ParamRK.b64 * kr4 + ParamRK.b65 * kr5);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (6)
    [kr6, qr_idx] = funcr_vr4RK(vr_temp, qr_idx, hr_idx); 
    vr_temp = vr_idx + ddt * (ParamRK.c1 * fr + ParamRK.c3 * kr3 + ParamRK.c4 * kr4 + ParamRK.c6 * kr6);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    
    vr_err = ddt * (ParamRK.dc1 * fr + ParamRK.dc3 * kr3 + ParamRK.dc4 * kr4 + ParamRK.dc5 * kr5 + ParamRK.dc6 * kr6);
   

end