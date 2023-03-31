%% subroutine hr2vr
function vr = hr2vr(hr,k,area,area_ratio_idx)  % calculating vr[m^3] by using hr[m]

% id = sec_map_idx(k)
id = 0;

if id <= 0
 vr = hr * area * area_ratio_idx; 
else
%  div_max = sec_div(id)
%  if( hr <= sec_hr(id, 1) ) 
%   a = hr * sec_b(id, 1);
%  elseif( hr > sec_hr(id, div_max) )
%   a = sec_area(id, div_max) + (hr - sec_hr(id, div_max)) * sec_b(id, div_max);
%  else
%   for i = 2:div_max
%    if( hr .le. sec_hr(id, i) )
%     a = sec_area(id, i - 1) + (hr - sec_hr(id, i - 1)) * sec_b(id, i);
%     break
%    end
%   end
%  end
%  vr = a * len_riv_idx(k);
end

% return
end