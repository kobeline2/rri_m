%%
% syms theta(t) theta_t(t) omega_0
% eqs = [diff(theta)   == theta_t;
%        diff(theta_t) == -omega_0^2*sin(theta)];
% 
% eqs  = subs(eqs,omega_0,omega_0Value);
% vars = [theta, theta_t];
% 
% [M,F] = massMatrixForm(eqs,vars);
% 
% f = M\F;
% 
% f = odeFunction(f, vars);
% 
% x0 = [0; 1.99*omega_0Value];
% tInit  = 0;
% tFinal = 10;
% 
% sols = ode45(f,[tInit tFinal],x0);


%% river

% KK = down_riv_idx(K);
% zb_p = zb_riv_idx(K);
% % hr_p = hr_idx(K);
% hr_p = 3;
% zb_n = zb_riv_idx(KK);
% % hr_n = hr_idx(KK);
% hr_n = 3.1;
% distance = dis_riv_idx(K);
% 
%  dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance; 
%  h = hr_p;
% 
% opts = odeset('AbsTol',eps*2);
% [t,h]=ode45(@(t,h) odefun_r(t,h,dh,ns_river,width_idx(K),area,area_ratio_idx(K)),[time time+ddt],h,opts);
% plot(t,h,'-o')
% disp(h(end))
% 
% function dhdt = odefun_r(t,h,dh,n,w,A1,A2) 
% A = sqrt(abs(dh)) / n;
% R = (w * h) / (w + 2 * h);
% dhdt = A * R^(2/3) * w * h  / (A1 * A2);
% end

%% slope

K = 100;

zb_p = zb_slo_idx(K);
hs_p = hs_idx(K);
ns_p = ns_slo_idx(K);
ka_p = ka_idx(K);
da_p = da_idx(K);
dm_p = dm_idx(K);
b_p  = beta_idx(K);
dif_p = dif_slo_idx(K);

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

    %  water coming in or going out?
    if dh >= 0
        % going out
        h = hs_p;
        if zb_p < zb_n; h = max(0, zb_p + hs_p - zb_n); end
%         qs_idx(L,K) = hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, area);
        [t,h]=ode45(@(t,h) odefun_s(t,h,dh,ns_river,width_idx(K),area,area_ratio_idx(K)),[time time+ddt],h);
    else
    % coming in
%         h = hs_n;
%         dh = abs(dh);
%         if zb_n < zb_p; h = max(0, zb_n + hs_n - zb_p); end
%         qs_idx(L,K) =  - hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, area);
    end

end

function q = hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, area)

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

end

