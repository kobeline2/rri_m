%% subroutine riv_idx_setting
tic

riv_count = sum(domain > 0 & riv == 1,'all'); % number of river cell
  
down_riv_idx = zeros(riv_count,1);
domain_riv_idx = zeros(riv_count,1);
dis_riv_idx = zeros(riv_count,1); % rivセルの距離


domAndRiv = (domain>0 & riv==1);                   % 1:rivセル， 0:範囲外 or sloのみ                        
[riv_idx2i,riv_idx2j] = find(domAndRiv');          % rivセルのx,y座標
riv_ij2idx = zeros(size(domAndRiv));              % rivセルの座標
riv_ij2idx(domAndRiv) = 1:sum(domAndRiv, 'all');  % rivセルの番号を振る  
domain_riv_idx         = domain(domAndRiv);  % domain
width_idx              = width(domAndRiv);   % 川幅
depth_idx              = depth(domAndRiv);   % 河道深さ
% height_idx = height(domAndRiv);   % 堤防高
% area_ratio_idx = area_ratio(domAndRiv);   % セルにおける河川の割合
% zb_riv_idx = zb_riv(domAndRiv);   % 不透水層の標高
% dif_riv_idx = dif(domAndRiv);   %%%%%　修正　%%%%%%%
% sec_map_idx = sec_map(domAndRiv);   % 
% len_riv_idx = len_riv(domAndRiv);   % 

%%
riv_count = 0;
for J = 1:NY
    for I = 1:NX
        if domain(I,J) == 0 | riv(I,J) ~= 1; continue; end
        riv_count = riv_count + 1;
        % east (down) 
        if flowDir(I,J) == 1 
            II = I + 1;
            JJ = J;
            distance = dx;
        % south east (right down)
        elseif flowDir(I,J) == 2 
            II = I + 1;
            JJ = J + 1;
            distance = sqrt(dx*dx + dy*dy);
        % south (right)
        elseif flowDir(I,J) == 4 
            II = I;
            JJ = J + 1;
            distance = dy;
        % south west (right up)
        elseif flowDir(I,J) == 8 
            II = I - 1;
            JJ = J + 1;
            distance = sqrt(dx*dx + dy*dy);
        % west (up)
        elseif flowDir(I,J) == 16 
            II = I - 1;
            JJ = J;
            distance = dx;
        % north west (left up)
        elseif flowDir(I,J) == 32 
            II = I - 1;
            JJ = J - 1;
            distance = sqrt(dx*dx + dy*dy);
        % north (left)
        elseif flowDir(I,J) == 64 
            II = I;
            JJ = J - 1;
            distance = dx;
        % north east (left down)
        elseif flowDir(I,J) == 128
            II = I + 1;
            JJ = J - 1;
            distance = sqrt(dx*dx + dy*dy);
        end

        
        %  If the downstream cell is outside the domain, set domain(i, j) = 2
        if II < 1 || II > NY || JJ < 1 || JJ > NX || domain(II,JJ) == 0
            domain(I,J) = 2;
            flowDir(I,J) = 0;
            II = I;
            JJ = J;
        end
        
        if riv(II,JJ) == 0
            %error("riv(II, JJ) should be 1 (", I, J, ")", "(", II, JJ, ")")
            disp('error')
        end
        
        dis_riv_idx(riv_count) = distance;
        down_riv_idx(riv_count) = riv_ij2idx(II,JJ);
    end
end

disp(size(dis_riv_idx))
disp(size(down_riv_idx))

toc

% if sec_length_switch == 1 
%     for k = 1:riv_count
%         kk = down_riv_idx(k);
%         dis_riv_idx(k) = ( len_riv_idx(k) + len_riv_idx(kk) ) / 2.d0;
%     end
% end

%% subroutine slo_idx_setting
tic
i4 = 1;

slo_count = sum(domain ~= 0,'all'); % number of slope cell
        
down_slo_idx = zeros(i4,slo_count);   %
zb_slo_idx = zeros(slo_count,1);      % 不透水層の標高
dis_slo_idx = zeros(slo_count,1);     % sloセルの距離
len_slo_idx = zeros(i4,slo_count);    %
down_slo_1d_idx = zeros(slo_count,1); %
dis_slo_1d_idx = zeros(slo_count,1);  %
len_slo_1d_idx = zeros(slo_count,1);  %
land_idx = zeros(slo_count,1);        %

dif_slo_idx = zeros(slo_count,1);   %
ns_slo_idx = zeros(slo_count,1);    %
soildepth_idx = zeros(slo_count,1); %
gammaa_idx = zeros(slo_count,1);    %

% ksv_idx = zeros(slo_count,1);          %
% faif_idx = zeros(slo_count,1);         %
% infilt_limit_idx = zeros(slo_count,1); %
% ka_idx = zeros(slo_count,1);           %
% gammam_idx = zeros(slo_count,1);       %
% beta_idx = zeros(slo_count,1);         %
% da_idx = zeros(slo_count,1);           %
% dm_idx = zeros(slo_count,1);           %
% ksg_idx = zeros(slo_count,1);          %
% gammag_idx = zeros(slo_count,1);       %
% kg0_idx = zeros(slo_count,1);          %
% fpg_idx = zeros(slo_count,1);          %
% rgl_idx = zeros(slo_count,1);          %

%%

slo_count = 0;
for J = 1:NY
    for I = 1:NX
        if domain(I,J)==0; continue; end
        slo_count = slo_count + 1;     
%         zb_slo_idx(slo_count) = zb(I, J);                                 
%         land_idx(slo_count) = land(I, J);            
% 
%         dif_slo_idx(slo_count) = dif(land(I, J));
%         ns_slo_idx(slo_count) = ns_slope(land(I,J));
%         soildepth_idx(slo_count) = soildepth(land(I, J));
%         gammaa_idx(slo_count) = gammaa(land(I, J));
% 
%         ksv_idx(slo_count) = ksv(land(I,J));
%         faif_idx(slo_count) = faif(land(I,J));
%         infilt_limit_idx(slo_count) = infilt_limit(land(I,J));
% 
%         ka_idx(slo_count) = ka(land(I,J));
%         gammam_idx(slo_count) = gammam(land(I,J));
%         beta_idx(slo_count) = beta(land(I,J));
%         da_idx(slo_count) = da(land(I,J));
%         dm_idx(slo_count) = dm(land(I,J));
% 
%         ksg_idx(slo_count) = ksg(land(I, J));
%         gammag_idx(slo_count) = gammag(land(I, J));
%         kg0_idx(slo_count) = kg0(land(I, J));
%         fpg_idx(slo_count) = fpg(land(I, J));
%         rgl_idx(slo_count) = rgl(land(I, J));

    end
end
 
domAndSlo = (domain>0);                            % 1:範囲内， 0:範囲外                        
[slo_idx2i,slo_idx2j] = find(domAndSlo');          % sloセルのx,y座標
slo_ij2idx = zeros(size(domAndSlo));              % sloセルの座標
slo_ij2idx(domAndSlo) = 1:sum(domAndSlo, 'all');  % sloセルの番号を振る  
domain_slo_idx         = domain(domAndSlo);  % domain
acc_slo_idx = flowAcc(domAndSlo);  

%%
eightFlowDir = 1;
if eightFlowDir == 1
    lmax = 4;
    l1 = dy / 2;
    l2 = dx / 2;
    l3 = sqrt(dx^2 + dy^2)/4;
elseif eightFlowDir == 0
    lmax =2;
    l1 = dy;
    l2 = dx;
    l3 = 0;
else
    error("error: eight_dir should be 0 or 1")
end



