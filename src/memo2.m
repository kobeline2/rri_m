%% hr2vr（改良）

% for K = 1:riv_count
%     vr_idx(K) = hr2vr_new(hr_idx(K), K, area, area_ratio_idx(K));  % Kは不要になりそう
% end
vr_idx = hr2vr_new(hr_idx, area, area_ratio_idx);

function vr = hr2vr_new(hr, area, area_ratio_idx)  % calculating vr[m^3] by using hr[m]

riv_count = 1484;
sec_map_idx = zeros(riv_count, 1); % 後で削除
id = (sec_map_idx > 0);

vr = hr * area .* area_ratio_idx;

% if id <= 0
%  vr = hr * area * area_ratio_idx; 
% else
% end

end