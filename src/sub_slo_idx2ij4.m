function a = sub_slo_idx2ij4(a_idx, a, slo_count, slo_idx2i, slo_idx2j, i4)
% sub_slo_idx2ij4
% subroutine sub_slo_idx2ij4
% a = sub_slo_idx2ij4(a_idx, a, riv_count, riv_idx2i, riv_idx2j, i4)
% 
% [ref]

for I = 1:i4
    for K = 1:slo_count
        a(slo_idx2j(K), slo_idx2i(K), I) = a_idx(I, K); 
    end
end

end