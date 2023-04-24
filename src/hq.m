function q = hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, area)
% water depth and discharge relationship
% subroutine hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, q)
if b_p > 0
    km = ka_p / b_p;
else
    km = 0;
end
vm = km * dh;

if da_p > 0
    va = ka_p * dh;
else
    va = 0;
end

if dh < 0; dh = 0; end
al = sqrt(dh) / ns_p;
m = 5 / 3;

if h < dm_p 
    q = vm * dm_p * (h / dm_p) ^ b_p;
elseif h < da_p 
    q = vm * dm_p + va * (h - dm_p);
else
    q = vm * dm_p + va * (h - dm_p) + al * (h - da_p) ^ m;
end

% discharge per unit area
% (q multiply by width and divide by area)
q = q * len / area;

% water depth limitter (1 mm)
% note: it can be set to zero
% if( h.le.0.001 ) q = 0.d0

end 