function a = sub_riv_idxij(a_idx, a, riv_count, riv_idx2i, riv_idx2j)
% sub_riv_idx2ij
% subroutine sub_riv_idx2ij
% a = sub_riv_idx2ij(a_idx, a, riv_count, riv_idx2i, riv_idx2j)
% 
% [ref]


for k = 1:riv_count
    a(riv_idx2j(k), riv_idx2i(k)) = a_idx(k); % 転置してあるのでFortranと逆
end

end