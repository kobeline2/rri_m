function hs_idx = infilt(hs_idx, gampt_f_idx, gampt_ff_idx, ksv_idx, faif_idx, gammaa_idx, infilt_limit_idx, dt, slo_count)
% infilt  calculating infiltretion by using Green Ampt Model
% subroutine infilt
% hs_idx = infilt(hs_idx, gampt_f_idx, gampt_ff_idx, ksv_idx, faif_idx, gammaa_idx, infilt_limit_idx, dt, slo_count)
%
% 
% [ref]

for K = 1:slo_count

    gampt_f_idx(K) = 0;
    gampt_ff_temp = gampt_ff_idx(K);
    if gampt_ff_temp <= 0.01; gampt_ff_temp = 0.01; end

    % gampt_f_idx(K) : infiltration capacity [m/s]
    % gampt_ff : accumulated infiltration depth [m]
    gampt_f_idx(K) = ksv_idx(K) * (1 + faif_idx(K) * gammaa_idx(K) / gampt_ff_temp);

    % gampt_f_idx(K) : infiltration capacity -> infiltration rate [m/s]
    if gampt_f_idx(K) >= hs_idx(K) / dt; gampt_f_idx(K) = hs_idx(K) / dt; end

    % gampt_ff should not exceeds a certain level
    if infilt_limit_idx(K) >= 0 && gampt_ff_idx(K) >= infilt_limit_idx(K); gampt_f_idx(K) = 0; end

    % update gampt_ff [m]
    gampt_ff_idx(K) = gampt_ff_idx(K) + gampt_f_idx(K) * dt;

    % hs : hs - infiltration rate * dt [m]
    hs_idx(K) = hs_idx(K) - gampt_f_idx(K) * dt;
    if hs_idx(K) <= 0; hs_idx(K) = 0; end

end

end