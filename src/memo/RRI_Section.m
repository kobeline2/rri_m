%% subroutine set_section
tic

hr = 3;
hrs = 3;
area = 3;
k = 0;

hr_new = hr_update(hr,hrs,k,area);

b = sec_h2b(hr,k)

toc
%% subroutine sec_hq_riv


%% subroutine hr_update
function hr_new = hr_update(hr,hrs,k,area)

vr_inc = hrs * area;

vr_org = hr2vr(hr,k,area)
vr_new = vr_org + vr_inc
hr_new = vr2hr(vr_new,k,area)

end
%% subroutine hr2vr
function vr = hr2vr(hr,k,area)  % calculating vr[m^3] by using hr[m]

% id = sec_map_idx(k)
id = 0;
area_ratio_idx = 0.5;

if(id <= 0)
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

%% subroutine vr2hr

function hr = vr2hr(vr,k,area)     % calculating hr[m] by using vr[m^3]

% id = sec_map_idx(k);
id = 0;
area_ratio_idx = 0.5;

if( id <= 0 )
 hr = vr / ( area * area_ratio_idx ) 
% else
%  div_max = sec_div(id);
%  a = vr / len_riv_idx(k);
%  if( a <= sec_area(id, 1) ) 
%   hr = a / sec_b(id, 1);
%  elseif( a > sec_area(id, div_max) ) 
%   hr = (a - sec_area(id, div_max)) / sec_b(id, div_max) + sec_hr(id, div_max);
%  else
%   for i = 2:div_max
%    if( a <= sec_area(id, i) ) 
%     hr = (a - sec_area(id, i-1)) / sec_b(id, i) + sec_hr(id, i - 1);
%     break
%    end
%   end
%  end
end

% return
end
%% subroutine sec_h2b

function b = sec_h2b(h,k)

width_idx = 1;
% id = sec_map_idx(k)
id = 0;
if id <= 0 
    b = width_idx
else
%     div_max = sec_div(id);
%     for I = 1:div_max;
%         if h <= sec_hr
%             b = sec_b;
%             break
%         end
%         b = sec_b;
%     end
end

end