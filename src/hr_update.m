function hr_new = hr_update(hr,hrs,k,area, area_ratio_idx)
% hr_update  updating hr[m]
% subroutine hr_update
% hr_new = hr_update(hr,hrs,k,area, area_ratio_idx)
%
% 
% [ref]

vr_inc = hrs * area;

vr_org = hr2vr(hr,k,area,area_ratio_idx);
vr_new = vr_org + vr_inc;
hr_new = vr2hr(vr_new,k,area,area_ratio_idx);

end