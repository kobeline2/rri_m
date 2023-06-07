function a_idx = sub_slo_ij2idx(a, a_idx, slo_count, slo_idx2i, slo_idx2j)
% sub_slo_ij2idx
% subroutine sub_slo_ij2idx
% a_idx = sub_slo_ij2idx(a, a_idx, slo_count, slo_idx2i, slo_idx2j)
% 
% [ref]

for k = 1:slo_count
    a_idx(k) = a(slo_idx2j(k), slo_idx2i(k)); % 転置してあるのでFortranと逆
end

end