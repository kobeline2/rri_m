%% subroutine storage calc
% storageCalc : srtorage 計算
% subroutine storage calc
% 
% 
%
% [ref]
tic

ss = 0;    %
sr = 0;    %
si = 0;    %
sg = 0;    %

vr_temp = 0;

for J = 1:NY
    for I = 1:NX
        if domain(I,J) == 0;continue;end

        if rivThresh >= 0 & riv(I,J) == 1
            vr_temp = hr2vr(hr(I,J), riv_ij2idx(I,J), area, area_ratio_idx);
            sr = sr + vr_temp;
        end
        si = si + gampt_ff(I,J) * area;
        sg = sg - hg(I,J) * gammag_idx(slo_ij2idx(I,J)) * area;
    end
end

toc
%%
function vr = hr2vr(hr,k,area,area_ratio_idx)  % calculating vr[m^3] by using hr[m]
% 簡略版

id = 0;

if(id <= 0)
 vr = hr * area * area_ratio_idx(k); 
end

end