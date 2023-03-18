function [dis_riv_idx, down_riv_idx] = rivIdxSetting(riv_count,NX,NY,domain,riv,flowDir,dx,dy,riv_ij2idx)
% readIdxSetting
% subroutine riv_idx_seting
% [dis_riv_idx, down_riv_idx] = rivIdxSetting(riv_count,NX,NY,domain,riv,flowDir,dx,dy,riv_ij2idx)
% 
% [ref]

down_riv_idx = zeros(riv_count,1); % 行先
dis_riv_idx  = zeros(riv_count,1); % rivセルの距離

distance  = 0;
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
            distance = dy;
        % north east (left down)
        elseif flowDir(I,J) == 128
            II = I + 1;
            JJ = J - 1;
            distance = sqrt(dx*dx + dy*dy);
        elseif flowDir(I,J) == 0 | flowDir(I,J) == -1
            II = I;
            JJ = J;
        else
            error(['dir(i, j) is error (', num2str(I),',', num2str(J), ')"',num2str(flowDir(I,J))])
        end
        

        %  If the downstream cell is outside the domain, set domain(i, j) = 2
        if II < 1 || II > NX || JJ < 1 || JJ > NY || domain(II,JJ) == 0
            domain(I,J) = 2;
            flowDir(I,J) = 0;
            II = I;
            JJ = J;
        end
        
        if riv(II,JJ) == 0
            formatSpec = "riv(II, JJ) should be 1 (%d, %d)  (%d, %d)";
            error(formatSpec, I,J,II,JJ);
        end
        
        dis_riv_idx(riv_count) = distance;
        down_riv_idx(riv_count) = riv_ij2idx(II,JJ);
    end
end



% if sec_length_switch == 1 
%     for k = 1:riv_count
%         kk = down_riv_idx(k);
%         dis_riv_idx(k) = ( len_riv_idx(k) + len_riv_idx(kk) ) / 2.d0;
%     end
% end