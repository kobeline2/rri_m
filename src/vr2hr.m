function hr = vr2hr(vr,k,area,area_ratio_idx)   
% vr2hr  calculating hr[m] by using vr[m^3]
% subroutine vr2hr
% hr = vr2hr(vr,k,area,area_ratio_idx)
%
% 
% [ref]

% id = sec_map_idx(k);
id = 0;

if id <= 0 
    hr = vr / ( area * area_ratio_idx );
else
%     div_max = sec_div(id);
%     a = vr / len_riv_idx(k);
%     if( a <= sec_area(id, 1) ) 
%         hr = a / sec_b(id, 1);
%     elseif( a > sec_area(id, div_max) ) 
%         hr = (a - sec_area(id, div_max)) / sec_b(id, div_max) + sec_hr(id, div_max);
%     else
%     for i = 2:div_max
%         if( a <= sec_area(id, i) ) 
%             hr = (a - sec_area(id, i-1)) / sec_b(id, i) + sec_hr(id, i - 1);
%             break
%         end
%     end
end

    
end