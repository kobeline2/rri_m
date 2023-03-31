function lev = h2lev(h, soildepth, gammaa)
% h2lev
% subroutine h2lev
% lev = h2lev(h, soildepth, gammaa)
% 
% [ref]

da_temp = soildepth * gammaa;

if soildepth == 0
    lev = h;
elseif h >= da_temp  % including da = 0
    lev = soildepth + (h - da_temp); % surface water
else
    if soildepth > 0; rho = da_temp / soildepth; end
    lev = h / rho;
end

end