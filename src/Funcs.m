% subroutine funcs(hs_idx, qp_t_idx, fs_idx, qs_idx )
function fs_idx = Funcs(qs_idx, hs_idx, qp_t_idx, slo_count, zb_slo_idx, ns_slo_idx, ka_idx,...
    dis_slo_idx, dis_slo_1d_idx, down_slo_idx, down_slo_1d_idx, len_slo_idx, len_slo_1d_idx,...
    da_idx, dm_idx, beta_idx, dif_slo_idx, soildepth_idx, gammaa_idx, lmax, area)

% qs_calc(hs_idx, qs_idx)
% lateral discharge (slope)

% $omp parallel do private(kk,zb_p,hs_p,ns_p,ka_p,da_p,dm_p,b_p,dif_p,l,distance,len, &
% $omp                     zb_n,hs_n,ns_n,ka_n,da_n,dm_n,b_n,dif_n,lev_p,lev_n,dh,hw)
for K = 1:slo_count
     zb_p = zb_slo_idx(K);
     hs_p = hs_idx(K);
     ns_p = ns_slo_idx(K);
     ka_p = ka_idx(K);
     da_p = da_idx(K);
     dm_p = dm_idx(K);
     b_p  = beta_idx(K);
     dif_p = dif_slo_idx(K);

    % 8-direction: lmax = 4, 4-direction: lmax = 2
    for L = 1:lmax % (1: right�C2: down, 3: right down, 4: left down)
        if dif_p == 0 && L == 2; break; end % kinematic -> 1-direction
        KK = down_slo_idx(L, K);
        if dif_p == 0; KK = down_slo_1d_idx(K); end
        if KK == -1; continue; end

        distance = dis_slo_idx(L, K);
        len = len_slo_idx(L, K);
        if dif_p == 0; distance = dis_slo_1d_idx(K); end
        if dif_p == 0; len = len_slo_1d_idx(K); end

        zb_n = zb_slo_idx(KK);
        hs_n = hs_idx(KK);
        ns_n = ns_slo_idx(KK);
        ka_n = ka_idx(KK);
        da_n = da_idx(KK);
        dm_n = dm_idx(KK);
        b_n = beta_idx(KK);
        dif_n = dif_slo_idx(KK);

        lev_p = h2lev(hs_p, soildepth_idx(K), gammaa_idx(K));
        lev_n = h2lev(hs_n, soildepth_idx(KK), gammaa_idx(KK));

        % diffusion wave
        dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance;

        % 1-direction : kinematic wave
        if dif_p == 0; dh = max( (zb_p - zb_n) / distance, 0.001 ); end

        % embankment
        % if emb_switch == 1
            % if L == 1; emb = emb_r_idx(K); end
            % if L == 2; emb = emb_b_idx(K); end
            % if L == 3; emb = max( emb_r_idx(k), emb_b_idx(k) ); end
            % if L == 4; emb = max( emb_r_idx(kk), emb_b_idx(k) ); end
        % end

        %  water coming in or going out?
        if dh >= 0
            % going out
            hw = hs_p;
            % if emb > 0; hw = max(hs_p - emb, 0); end
            if zb_p < zb_n; hw = max(0, zb_p + hs_p - zb_n); end
            qs_idx(L,K) = hq(ns_p, ka_p, da_p, dm_p, b_p, hw, dh, len, area);
        else
        % coming in
            hw = hs_n;
            % if emb > 0; hw = max(hs_n - emb, 0); end
            dh = abs(dh);
            if zb_n < zb_p; hw = max(0, zb_n + hs_n - zb_p); end
            qs_idx(L,K) =  - hq(ns_p, ka_p, da_p, dm_p, b_p, hw, dh, len, area);
        end

    end
end
%$omp end parallel do
% end subroutine qs_calc

% boundary condition for slope (discharge boundary)
% if bound_slo_disc_switch >= 1 
%     itemp = -1
%     for jtemp = 1:tt_max_bound_slo_disc
%         if t_bound_slo_disc(jtemp) < (time + ddt) & (time + ddt) <= t_bound_slo_disc(jtemp+1); itemp = jtemp; end
%     end
%     for K = 1:slo_count
%         if bound_slo_disc_idx(itemp, K) <= -100.0; continue; end % not boundary
%         % right
%         if dir(slo_idx2i(K), slo_idx2j(K)) == 1
%             qs_idx(1, K) = bound_slo_disc_idx(itemp, K) / area;
%         % right down
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 2 
%             qs_idx(3, K) = bound_slo_disc_idx(itemp, K) / area;
%         % down
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 4 
%             qs_idx(2, K) = bound_slo_disc_idx(itemp, K) / area;
%         % left down
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 8 
%             qs_idx(4, K) = bound_slo_disc_idx(itemp, K) / area;
%         % left
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 16 
%             qs_idx(1, K) = - bound_slo_disc_idx(itemp, K) / area;
%         % left up
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 32 
%             qs_idx(3, K) = - bound_slo_disc_idx(itemp, K) / area;
%         % up
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 64 
%             qs_idx(2, K) = - bound_slo_disc_idx(itemp, K) / area;
%         % right up
%         elseif dir(slo_idx2i(K), slo_idx2j(K)) == 128 
%             qs_idx(4, K) = - bound_slo_disc_idx(itemp, K) / area;
%         end
%     end
% end

% qs_idx > 0 --> discharge flowing out from a cell

%$omp parallel do
% for K = 1:slo_count
%  fs_idx(K) = qp_t_idx(K) - (qs_idx(1,K) + qs_idx(2,K) + qs_idx(3,K) + qs_idx(4,K));
% end
fs_idx = qp_t_idx - sum(qs_idx,1)';
%$omp end parallel do

for K = 1:slo_count
    for L = 1:lmax
        if dif_slo_idx(K) == 0 && L == 2; break; end % kinematic -> 1-direction
    KK = down_slo_idx(L, K);
    if dif_slo_idx(K) == 0; KK = down_slo_1d_idx(K); end
    if KK == -1; continue; end
    fs_idx(KK) = fs_idx(KK) + qs_idx(L, K);
    end
end

% end subroutine funcs


% water depth (h) to actual water level (lev)
% function lev = h2lev(h, soildepth, gammaa)
% 
% da_temp = soildepth * gammaa;
% 
% if soildepth == 0
%     lev = h;
% elseif h >= da_temp  % including da = 0
%     lev = soildepth + (h - da_temp); % surface water
% else
%     if soildepth > 0; rho = da_temp / soildepth; end
%     lev = h / rho;
% end
% 
% end


end