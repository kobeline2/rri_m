function [down_slo_idx,dis_slo_idx,len_slo_idx,down_slo_1d_idx,dis_slo_1d_idx,len_slo_1d_idx] ...
    = sloIdxSetting(slo_count,NX,NY,domain,flowDir,dx,dy,slo_ij2idx,i4,eightFlowDir)
% readIdxSetting
% subroutine slo_idx_seting
% [down_slo_idx,dis_slo_idx,len_slo_idx,down_slo_1d_idx,dis_slo_1d_idx,len_slo_1d_idx] ...
%    = sloIdxSetting(slo_count,NX,NY,domain,flowDir,dx,dy,slo_ij2idx,i4,eightFlowDir)
% 
% [ref]

down_slo_idx = zeros(i4,slo_count);   %
dis_slo_idx = zeros(i4,slo_count);    % sloセルの距離
len_slo_idx = zeros(i4,slo_count);    % sloセルの距離
down_slo_1d_idx = zeros(slo_count,1); %
dis_slo_1d_idx = zeros(slo_count,1);  %
len_slo_1d_idx = zeros(slo_count,1);  %

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

%  search for downstream gridcell (down_slo_idx)
slo_count = 0;
down_slo_idx(:,:) = -1;
for J = 1:NY
    for I = 1:NX
        if domain(I,J) == 0 ; continue; end
        slo_count = slo_count + 1;   
        
        for L = 1:lmax
            % east (down) 
            if L == 1
               II = I + 1;
               JJ = J;
               distance = dx;
               len = l1;
            elseif L == 2
                % south (right)
                II = I;
                JJ = J + 1;
                dustance = dy;
                len = l2;
            elseif L == 3
                % south east (right down)
                II = I + 1;
                JJ = J + 1;
                distance = sqrt(dx^2 + dy^2);
                len = l3;
            else % L == 4
                % south west (right up)
                II = I - 1;
                JJ = J + 1;
                distance = sqrt(dx^2 + dy^2);
                len = l3;                            
            end
        
        if II > NX; continue; end
        if JJ > NY; continue; end
        if II < 1; continue; end
        if JJ < 1; continue; end
        if domain(II,JJ) == 0; continue; end
        
        down_slo_idx(L, slo_count) = slo_ij2idx(II,JJ);
        dis_slo_idx(L, slo_count) = distance;
        len_slo_idx(L, slo_count) = len;
            
        end
    end
end

% search for downstream gridcell (down_slo_1d_idx) (used only for kinematic with 1-direction)
slo_count = 0;
down_slo_1d_idx(:) = -1;
dis_slo_1d_idx(:) = l1;
len_slo_1d_idx(:) = dx;

l1_kin = dy;
l2_kin = dx;
l3_kin = dx * dy / sqrt( dx ^ 2 + dy ^ 2 );

for J = 1: NY
    for I = 1: NX

        if domain(I,J) == 0; continue; end
        % domain(i, j) = 1 or 2
        slo_count = slo_count + 1;

        % east (down) 
        if flowDir(I,J) == 1 
            II = I + 1;
            JJ = J;
            distance = dx;
            len = l1_kin;
        % south east (right down)
        elseif flowDir(I,J) == 2 
            II = I + 1;
            JJ = J + 1;
            distance = sqrt(dx*dx + dy*dy);
            len = l3_kin;
        % south (right)
        elseif flowDir(I,J) == 4 
            II = I;
            JJ = J + 1;
            distance = dy;
            len = l2_kin;
        % south west (right up)
        elseif flowDir(I,J) == 8 
            II = I - 1;
            JJ = J + 1;
            distance = sqrt(dx*dx + dy*dy);
            len = l3_kin;
        % west (up)
        elseif flowDir(I,J) == 16 
            II = I - 1;
            JJ = J;
            distance = dx;
            len = l1_kin;
        % north west (left up)
        elseif flowDir(I,J) == 32 
            II = I - 1;
            JJ = J - 1;
            distance = sqrt(dx*dx + dy*dy);
            len = l3_kin;
        % north (left)
        elseif flowDir(I,J) == 64 
            II = I;
            JJ = J - 1;
            distance = dy;
            len = l2_kin;
        % north east (left down)
        elseif flowDir(I,J) == 128
            II = I + 1;
            JJ = J - 1;
            distance = sqrt(dx*dx + dy*dy);
            len = l3_kin;
        elseif flowDir(I,J) == 0 || flowDir(I,J) == -1
            II = I;
            JJ = J;
            distance = dx;
            len = l1_kin;            
        else
            error(['dir(i, j) is error (', num2str(I),',', num2str(J), ')"',num2str(flowDir(I,J))])
        end

        if II > NX; continue; end
        if JJ > NY; continue; end
        if II < 1; continue; end
        if JJ < 1; continue; end
        if domain(II,JJ) == 0; continue; end

        down_slo_1d_idx(slo_count) = slo_ij2idx(II,JJ);
        dis_slo_idx(L, slo_count) = distance;
        len_slo_idx(L, slo_count) = len;

    end
end
