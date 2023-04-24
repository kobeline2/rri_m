% subroutine funcr( vr_idx, fr_idx, qr_idx )
function fr_idx = Funcr(vr_idx, qr_idx, hr_idx, area, area_ratio_idx, riv_count, domain_riv_idx,...
    zb_riv_idx, dif_riv_idx, dis_riv_idx, down_riv_idx, width_idx, ns_river)

qr_sum_idx = zeros(riv_count, 1);
qr_div_idx = zeros(riv_count, 1);

for K = 1:riv_count
 hr_idx(K) = vr2hr(vr_idx(K),K,area, area_ratio_idx(K));
end

% boundary condition for river (water depth boundary)
% if bound_riv_wlev_switch >= 1 
%  itemp = -1
%  for jtemp = 1: tt_max_bound_riv_wlev
%   if( t_bound_riv_wlev(jtemp-1) < (time + ddt)) & (time + ddt <= t_bound_riv_wlev(jtemp) );0 itemp = jtemp; end
%  end
%  for k = 1,:riv_count
%   if( bound_riv_wlev_idx(itemp, k) <= -100.0 ); continue; end % not boundary
%   hr_idx(k) = bound_riv_wlev_idx(itemp, k)
%   call hr2vr(hr_idx(k), k, vr_idx(k)) % add v1.4
%   end
% end

%%

%%%------------qr_calc(hr_idx, qr_idx)----------------%%%
qr_div_idx = zeros(riv_count, 1);

%$omp parallel do private(kk,zb_p,hr_p,distance,zb_n,hr_n,dh,hw,qr_temp)
for K = 1:riv_count
    if domain_riv_idx(K) == 2; continue; end
    zb_p = zb_riv_idx(K);
    hr_p = hr_idx(K);
    dif_p = dif_riv_idx(K);

    distance = dis_riv_idx(K);

    % information of the destination cell
    KK = down_riv_idx(K);
    zb_n = zb_riv_idx(KK);
    hr_n = hr_idx(KK);
    dif_n = dif_riv_idx(KK);

    % diffusion wave
    dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance; % diffussion

    % kinematic wave
    if dif_p == 0; dh = max( (zb_p - zb_n) / distance, 0.001 ); end

    % the destination cell is outlet (domain = 2)
    if domain_riv_idx(KK) == 2; dh = (zb_p + hr_p - zb_n) / distance; end % kinematic wave (+hr_p)

%     % the destination cell is outlet (domain = 2 ) with water depth boundary water table
%     if domain_riv_idx(kk) == 2 & bound_riv_wlev_switch >= 1 
%         % ver 1.4.2 mod by T.Sayama on June 24, 2015
%         %if( bound_riv_wlev_idx(1, kk) .le. -100.0 ) cycle % not boundary
%         %dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance % diffussion
%         if( bound_riv_wlev_idx(1, kk) > -100.0 ); dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance; end % diffussion
%     end

%     % kinematic wave (for dam and water gate) modified by T.Sayama on Aug 9, 2021
%     if  damflg(k) > 0  
%         dh = max((zb_p - zb_n) / distance, 0.001);
%     end
%     if  damflg(kk) > 0
%         if dam_floodq(damflg(kk)) > 0 
%             dh = max((zb_p - zb_n) / distance, 0.001);
%         end
%     end

    % from a tributary (levee height =< 0) to a main river (levee height > 0) : no reverse flow
    %if( height_idx(k) .le. 0 .and. height_idx(kk) .gt. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )
    
    a = sqrt(abs(dh)) / ns_river;

    if dh >= 0 
        hw = hr_p;
        if zb_p < zb_n; hw = max(0, zb_p + hr_p - zb_n); end
        r = (width_idx(K) * hw) / (width_idx(K) + 2.d0 * hw);
        qr_idx(K) = a * r ^ (2 / 3) * width_idx(K) * hw;
    else
    % reverse flow
        hw = hr_n;
        if zb_n < zb_p; hw = max(0, zb_n + hr_n - zb_p); end
        dh = abs(dh);
        r = (width_idx(K) * hw) / (width_idx(K) + 2.d0 * hw);
        qr_idx(K) = - a * r ^ (2 / 3) * width_idx(K) * hw;
    end

end
%%%---------------------------------------------------%%%

% if( div_switch == 1 ) call RRI_Div(qr_idx, hr_idx, qr_div_idx)

% % boundary condition for river (discharge boundary)
% if bound_riv_disc_switch >= 1 
%  for jtemp = 1:tt_max_bound_riv_disc
%   if( t_bound_riv_disc(jtemp-1) < (time + ddt) & (time + ddt) <= t_bound_riv_disc(jtemp) ); itemp = jtemp; end
%  end
%  for k = 1:riv_count
%   if( bound_riv_disc_idx(itemp, k) <= -100.0 ); continue; end % not boundary
%   %qr_idx(k) = bound_riv_disc_idx(itemp, k) / area  % qr_idx: discharge per unit area
%   qr_idx(k) = bound_riv_disc_idx(itemp, k) % modified v1.4
%   % linear interpolation of the boundary condition
%   %qr_idx(k) = bound_riv_disc_idx(itemp-1, k) * (t_bound_riv_disc(itemp) - (time + ddt)) &
%   %           + bound_riv_disc_idx(itemp, k) * ((time + ddt) - t_bound_riv_disc(itemp-1))
%   %qr_idx(k) = qr_idx(k) / (t_bound_riv_disc(itemp) - t_bound_riv_disc(itemp-1))
%   hr_idx(k) = 0
%   vr_idx(k) = 0 % add v1.4
%  end
% end

% % dam control
% if( dam_switch .eq. 1 ) then
%  call dam_prepare(qr_idx) % calculate inflow to dam
%  do i = 1, dam_num
%   if( dam_volmax(i) .gt. 0 ) then     %
%    % dam
%    call dam_operation( dam_loc(i) )
%    qr_idx( dam_loc(i) ) = dam_qout(i)
%    qr_sum_idx( dam_loc(i) ) = qr_sum_idx( dam_loc(i) ) + dam_qin( dam_loc(i) ) - dam_qout(i)
%   elseif( dam_volmax(i) .eq. 0 ) then
%    % barrage
%    if( hr_idx( dam_loc(i) ) .le. dam_floodq(i) ) then
%     qr_idx( dam_loc(i) ) = 0
%    endif
%   else
%    % water gate (dam_volmax(i) < 0) % added on Aug 7, 2021 by T.Sayama
%    k = dam_loc(i)
%    kk = down_riv_idx(k)
%    call gate_operation( dam_loc(i), hr_idx(k), hr_idx(kk) )
%    qr_idx( dam_loc(i) ) = dam_qout(i)
%    %qr_sum_idx( dam_loc(i) ) = qr_sum_idx( dam_loc(i) ) + dam_qin( dam_loc(i) ) - dam_qout(i)
%   endif
%  enddo
% endif

% qr_sum > 0 --> discharge flowing out from a cell
for K = 1:riv_count
 % outflow from (k)
 qr_sum_idx(K) = qr_sum_idx(K) + qr_idx(K);
 KK = down_riv_idx(K);
 if domain_riv_idx(KK)==0; continue; end
 % qr_sum minus (flowing into) discharge at the destination cell
 qr_sum_idx(KK) = qr_sum_idx(KK) - qr_idx(K);
end

% diversion
% if div_switch == 1 
%  for l = 1:div_id_max
%  % outflow from (k)
%   k = div_org_idx(l)
%   qr_sum_idx(k) = qr_sum_idx(k) + qr_div_idx(k)
%   kk = div_dest_idx(l)
%   if domain_riv_idx(kk) == 0; continue; end
%   % qr_sum minus (flowing into) discharge at the destination cell
%   qr_sum_idx(kk) = qr_sum_idx(kk) - qr_div_idx(k)
%  end
% end

% qr_sum divide by area_ratio
% (area_ratio = ratio of river area against total cell area)
%fr_idx = -qr_sum_idx / area_ratio_idx
fr_idx = - qr_sum_idx; % modified v1.4

% end subroutine funcr
% 
% 
% % lateral discharge (river)
% 
%% subroutine qr_calc(hr_idx, qr_idx)
% 
% qr_idx = zeros(riv_count, 1);
% qr_div_idx = zeros(riv_count, 1);
% 
% %$omp parallel do private(kk,zb_p,hr_p,distance,zb_n,hr_n,dh,hw,qr_temp)
% for K = 1:riv_count
%     if domain_riv_idx(K) == 2; continue; end
%     zb_p = zb_riv_idx(k);
%     hr_p = hr_idx(k);
%     dif_p = dif_riv_idx(k);
% 
%     distance = dis_riv_idx(k);
% 
%     % information of the destination cell
%     kk = down_riv_idx(k);
%     zb_n = zb_riv_idx(kk);
%     hr_n = hr_idx(kk);
%     dif_n = dif_riv_idx(kk);
% 
%     % diffusion wave
%     dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance; % diffussion
% 
%     % kinematic wave
%     if dif_p == 0; dh = max( (zb_p - zb_n) / distance, 0.001 ); end
% 
%     % the destination cell is outlet (domain = 2)
%     if domain_riv_idx(kk) == 2; dh = (zb_p + hr_p - zb_n) / distance; end % kinematic wave (+hr_p)

%     % the destination cell is outlet (domain = 2 ) with water depth boundary water table
%     if domain_riv_idx(kk) == 2 & bound_riv_wlev_switch >= 1 
%         % ver 1.4.2 mod by T.Sayama on June 24, 2015
%         %if( bound_riv_wlev_idx(1, kk) .le. -100.0 ) cycle % not boundary
%         %dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance % diffussion
%         if( bound_riv_wlev_idx(1, kk) > -100.0 ); dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance; end % diffussion
%     end

%     % kinematic wave (for dam and water gate) modified by T.Sayama on Aug 9, 2021
%     if  damflg(k) > 0  
%         dh = max((zb_p - zb_n) / distance, 0.001);
%     end
%     if  damflg(kk) > 0
%         if dam_floodq(damflg(kk)) > 0 
%             dh = max((zb_p - zb_n) / distance, 0.001);
%         end
%     end

    % from a tributary (levee height =< 0) to a main river (levee height > 0) : no reverse flow
    %if( height_idx(k) .le. 0 .and. height_idx(kk) .gt. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )
    
%     a = sqrt(abs(dh)) / ns_river;
% 
%     if dh >= 0 
%         hw = hr_p;
%         if( zb_p .lt. zb_n ); hw = max(0, zb_p + hr_p - zb_n); end
%         r = (width_idx(K) * hw) / (width_idx(K) + 2.d0 * hw);
%         qr_idx(K) = a * r ^ (2 / 3) * width_idx(K) * hw;
%     else
%     % reverse flow
%         hw = hr_n;
%         if zb_n < zb_p; hw = max(0, zb_n + hr_n - zb_p); end
%         dh = abs(dh);
%         r = (width_idx(K) * hw) / (width_idx(K) + 2.d0 * hw);
%         qr_idx(K) = - a * r ^ (2 / 3) * width_idx(K) * hw;
%     end
% 
% end
% %$omp end parallel do
% 
% 
%% water depth and discharge relationship
% %subroutine hq_riv(h, dh, w, q)
% subroutine hq_riv(h, dh, k, w, q)

% a = sqrt(abs(dh)) / ns_river
% r = (w * h) / (w + 2.d0 * h)
% q = a * r ** (2.d0 / 3.d0) * w * h
% 
% if( sec_map_idx(k) .gt. 0 ) then
%  call sec_hq_riv(h, dh, k, q)
% end

end