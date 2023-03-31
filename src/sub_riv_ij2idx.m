function a_idx = sub_riv_ij2idx(a, a_idx, riv_count, riv_idx2i, riv_idx2j)
% sub_riv_ij2idx
% subroutine sub_riv_ij2idx
% a_idx = sub_riv_ij2idx(a, a_idx, riv_count, riv_idx2i, riv_idx2j)
% 
% [ref]

for k = 1:riv_count
    a_idx(k) = a(riv_idx2j(k), riv_idx2i(k)); % 転置してあるのでFortranと逆
end

end