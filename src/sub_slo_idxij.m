function a = sub_slo_idxij(a_idx, a, slo_count, slo_idx2i, slo_idx2j)
% sub_slo_idx2ij
% subroutine sub_slo_idx2ij
% a = sub_slo_idx2ij(a_idx, a, slo_count, slo_idx2i, slo_idx2j)
% 
% [ref]


for k = 1:slo_count
    a(slo_idx2j(k), slo_idx2i(k)) = a_idx(k); % 転置してあるのでFortranと逆
end

end