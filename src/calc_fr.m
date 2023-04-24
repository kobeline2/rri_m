function [fr, vr_temp] = calc_fr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river, ddt)
    
    memo2
    %%%% boundary condition
    
    qr_ave_temp_idx = zeros(riv_count, 1);

    % % Adaptive Runge-Kutta
    % % (1)
    fr = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);   
    vr_temp = vr_idx + b21 * ddt * fr;
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    
    % % (2)
    kr2 = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);
    vr_temp = vr_idx + ddt * (b31 * fr + b32 * kr2);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (3)
    kr3 = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);
    vr_temp = vr_idx + ddt * (b41 * fr + b42 * kr2 + b43 * kr3);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (4)
    kr4 = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);
    vr_temp = vr_idx + ddt * (b51 * fr + b52 * kr2 + b53 * kr3 + b54 * kr4);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (5)
    kr5 = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river);
    vr_temp = vr_idx + ddt * (b61 * fr + b62 * kr2 + b63 * kr3 + b64 * kr4 + b65 * kr5);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    % 
    % % (6)
    kr6 = Funcr(vr_idx, qr_idx, hr_idx, cellarea, area_ratio_idx, riv_count, domain_riv_idx,...
        zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river); 
    vr_temp = vr_idx + ddt * (c1 * fr + c3 * kr3 + c4 * kr4 + c6 * kr6);
    vr_temp(vr_temp<0) = 0;
    qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt;
    
    vr_err = ddt * (dc1 * fr + dc3 * kr3 + dc4 * kr4 + dc5 * kr5 + dc6 * kr6);

    hr_err = vr_err / (cellarea * area_ratio_idx);

end