%% RRI_RivSlo.f90
tic

% 変数の準備
% hs = [1 2;1 2];
% hr = [2 4;5 4.6];
% height = [0 6;1 1];
% depth = [3 3;3 3];
% domain = [1 1;1 1];
% riv = [1 1;1 1];
% k = 0;
% nx = 2;
% ny = 2;
% qrs = zeros(2,2);
% 
% dt = 1;
% len = 1;
% area = 10;

%% subroutine funcrs

Mu1 = (2/3) ^ (3/2);
Mu2 = 0.35;
Mu3 = 0.91;

for I = 1:NX
    for J = 1:NY
        if domain(I,J) == 0 | riv(I,J) == 0;continue;end
        
        hs_top = hs(I,J);
        hr_top = hr(I,J) - depth(I,J);
        
        % k = riv_ij2idx(i, j)
        % len = len_riv_idx(k)
        
        %%
        if (height(I,J) == 0 & hr_top < 0)|(height(I,J) > 0 & hr_top < 0 & hs_top <= height(I,J))
            % From slope to river (hrs > 0)
            disp('Case A')
            hrs = Mu1 * hs_top * sqrt(9.81 * hs_top) * dt * len * 2 /area;
            if hrs > hs(I,J)
                hrs = hs(I,J);
            end
            hs(I,J) = hs(I,J) - hrs;
        
            hr_new = hr_update(hr(I,J),hrs,k,area);
            hr(I,J) = hr_new;
            qrs(I,J) = hrs;
        
            % avoid the situation of hr_top > hs_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            %%
            if hr_top > -0.00001 & hr_top > hs_top;
                for count = 1:10;
                    b = sec_h2b(hr(I,J),k);
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(2/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,k,area);
                    hr(I,J) = hr_new;
                    qrs(I,J) = qrs(I,J) + hrs;
                    if abs(hs(I,J) - (hr(I,J) - depth(I,J))) < 0.00001
                        break
                    end
                    hs_top = hs(I,J);
                    hr_top = hr(I,J) - depth(I,J);
                end
                hr(I,J) = hs(I,J) + depth(I,J);
            end
        
        
        elseif height(I,J) > 0 & hs_top <= height(I,J) & hr_top <= height(I,J) & hr_top >= 0
            % No exchange
            disp('case B')
            qrs(I,J) = 0;
        
        elseif hs_top <= hr_top & hr_top >= height(I,J)
            % From river to slope (hrs < 0)
            disp('case C')
            h1 = hr_top - height(I,J);
            h2 = hs_top - height(I,J);
            if h2/h1 <= 2/3
                hrs = -Mu2 * h1 * sqrt(2 * 9.81 * h1) * dt * len * 2 / area;
                disp('C-1')
            else
                hrs = -Mu3 * h2 * sqrt(2 * 9.81 * (h1-h2) ) * dt * len * 2 / area;
                disp('C-2')
            end
        
            b = sec_h2b(hr(I,J),k);
            ar = len * b / area;
            if abs(hrs/ar) > (hr_top - height(I,J))
                hrs = -(hr_top - height(I,J)) * ar;
            end
            qrs(I,J) = hrs;
        
            hs(I,J) = hs(I,J) - hrs;
            hr_new = hr_update(hr(I,J),hrs,k,area);
            hr(I,J) = hr_new;
            % avoid the situation of hs_top > hr_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            %%
            if hs_top > hr_top
                for count = 1:10
                    b = sec_h2b(hr(I,J),k);
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(2/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,k,area);
                    hr(I,J) = hr_new;
                    qrs(I,J) = qrs(I,J) + hrs;
                    if abs(hs(I,J) - (hr(I,J) - depth(I,J))) < 0.00001
                        break
                    end
                    hs_top = hs(I,J);
                    hr_top = hr(I,J) - depth(I,J);
                end
                hr(I,J) = hs(I,J) + depth(I,J);
            end
        
        elseif hs_top>= hr_top & hs_top >= height(I,J)
            % From slope to river (hrs > 0)
            disp('CaseD')
            h1 = hs_top - height(I,J);
            h2 = hr_top - height(I,J);
            if h2/h1 <= 2/3
                hrs = Mu2 * h1 * sqrt(2 * 9.81 * h1) * dt * len * 2 / area;
                disp('D-1');
            else
                hrs = Mu3 * h2 * sqrt(2 * 9.81 * (h1-h2) ) * dt * len * 2 / area;
                disp('D-2');
            end    
        
            if hrs > (hs_top - height(I,J))
                hrs = hs - height(I,J);
            end
            qrs(I,J) = hrs;
        
            hs(I,J) = hs(I,J) - hrs;
            hr_new = hr_update(hr(I,J),hrs,k,area);
            hr(I,J) = hr_new;
        
             % avoid the situation of hr_top > hs_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            if hr_top > -0.00001 & hr_top > hs_top
                for count = 1:10
                    b = sec_h2b(hr(I,J),k);
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(2/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,k,area);
                    hr(I,J) = hr_new;
                    qrs(I,J) = qrs(I,J) + hrs;
                    if abs(hs(I,J) - (hr(I,J) - depth(I,J))) < 0.00001
                        break
                    end
                    hs_top = hs(I,J);
                    hr_top = hr(I,J) - depth(I,J);
                end
                hr(I,J) = hs(I,J) + depth(I,J);
            end   
        
        else
            disp("Error : Rivslo");
            quit
        
        end
        
        qrs(I,J) = qrs(I,J) / dt;
    end
end  

 toc

%% RRI_Section.f90
function hr_new = hr_update(hr,hrs,k,area)

vr_inc = hrs * area;

vr_org = hr2vr(hr,k,area);
vr_new = vr_org + vr_inc;
hr_new = vr2hr(vr_new,k,area);

end

%%
function vr = hr2vr(hr,k,area)  % calculating vr[m^3] by using hr[m]
% 簡略版

id = 0;
area_ratio_idx = 0.5;

if(id <= 0)
 vr = hr * area * area_ratio_idx; 
end

end

%%
function hr = vr2hr(vr,k,area)     % calculating hr[m] by using vr[m^3]
% 簡略版

id = 0;
area_ratio_idx = 0.5;

if( id <= 0 )
 hr = vr / ( area * area_ratio_idx ) ;
end

end

%%
function b = sec_h2b(h,k)
% 簡略版

width_idx = 1;
id = 0;
if id <= 0 
    b = width_idx;
end

end
