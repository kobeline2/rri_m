function [hr, hs] = funcrs(hr, hs, NX, NY, domain, riv, height, depth, ...
    riv_ij2idx, len_riv_idx, dt, area, area_ratio_idx, width_idx)
% funcrs  calculating river-slope interactions
% subroutine funcrs
% [hr, hs] = funcrs(hr, hs, NX, NY, domain, riv, height, depth, ...
%    riv_ij2idx, len_riv_idx, dt, area, area_ratio_idx)
%
% 
% [ref]

Mu1 = (2/3) ^ (3/2);
Mu2 = 0.35;
Mu3 = 0.91;
qrs = zeros(NX,NY);

for I = 1:NX
    for J = 1:NY
        if domain(I,J) == 0 || riv(I,J) == 0;continue;end
        
        hs_top = hs(I,J);
        hr_top = hr(I,J) - depth(I,J);
        
        K = riv_ij2idx(I, J);
        len = len_riv_idx(K);
        
        %%
        if (height(I,J) == 0 && hr_top < 0)||(height(I,J) > 0 && hr_top < 0 && hs_top <= height(I,J))
            % From slope to river (hrs > 0)
            hrs = Mu1 * hs_top * sqrt(9.81 * hs_top) * dt * len * 2 /area;
            if hrs > hs(I,J)
                hrs = hs(I,J);
            end
            hs(I,J) = hs(I,J) - hrs;
        
            hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
            hr(I,J) = hr_new;
            qrs(I,J) = hrs;
        
            % avoid the situation of hr_top > hs_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            %%
            if hr_top >= -0.00001 && hr_top > hs_top
                for count = 1:10
                    b = sec_h2b(hr(I,J),width_idx(K));
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(1 + 1/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
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
        
        
        elseif height(I,J) > 0 && hs_top <= height(I,J) && hr_top <= height(I,J) && hr_top >= 0
            % No exchange
            qrs(I,J) = 0;
        
        elseif hs_top <= hr_top && hr_top >= height(I,J)
            % From river to slope (hrs < 0)
            h1 = hr_top - height(I,J);
            h2 = hs_top - height(I,J);
            if h2/h1 <= 2/3
                hrs = -Mu2 * h1 * sqrt(2 * 9.81 * h1) * dt * len * 2 / area;
            else
                hrs = -Mu3 * h2 * sqrt(2 * 9.81 * (h1-h2) ) * dt * len * 2 / area;
            end
        
            b = sec_h2b(hr(I,J),width_idx(K));
            ar = len * b / area;
            if abs(hrs/ar) > (hr_top - height(I,J))
                hrs = -(hr_top - height(I,J)) * ar;
            end
            qrs(I,J) = hrs;
        
            hs(I,J) = hs(I,J) - hrs;
            hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
            hr(I,J) = hr_new;
            % avoid the situation of hs_top > hr_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            %%
            if hs_top > hr_top
                for count = 1:10
                    b = sec_h2b(hr(I,J),width_idx(K));
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(1 + 1/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
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
        
        elseif hs_top>= hr_top && hs_top >= height(I,J)
            % From slope to river (hrs > 0)
            h1 = hs_top - height(I,J);
            h2 = hr_top - height(I,J);
            if h2/h1 <= 2/3
                hrs = Mu2 * h1 * sqrt(2 * 9.81 * h1) * dt * len * 2 / area;
            else
                hrs = Mu3 * h2 * sqrt(2 * 9.81 * (h1-h2) ) * dt * len * 2 / area;
            end    
        
            if hrs > (hs_top - height(I,J))
                hrs = hs(I,J) - height(I,J);
            end
            qrs(I,J) = hrs;
        
            hs(I,J) = hs(I,J) - hrs;
            hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
            hr(I,J) = hr_new;
        
             % avoid the situation of hr_top > hs_top
            hs_top = hs(I,J);
            hr_top = hr(I,J) - depth(I,J);
            if hr_top > -0.00001 && hr_top > hs_top
                for count = 1:10
                    b = sec_h2b(hr(I,J),width_idx(K));
                    ar = len * b / area;
                    hrs = (hs_top - hr_top)/(1 + 1/ar);
                    hs(I,J) = hs(I,J) - hrs;
        
                    hr_new = hr_update(hr(I,J),hrs,K,area, area_ratio_idx(K));
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
            error("Error : Rivslo")
        end
        
        qrs(I,J) = qrs(I,J) / dt;
    end
end  

end



