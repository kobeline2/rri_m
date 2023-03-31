%% rainfall

itemp = -1;
for jtemp = 1:tt_max_rain
    if t_rain(jtemp) < (time + ddt) && (time + ddt) <= t_rain(jtemp+1); itemp = jtemp; end
end
for J = 1:NY
    if itemp<= 0; continue; end  % バグ回避のため追加
    if rain_i(J) < 1 || rain_i(J) > ny_rain; continue; end
    for I = 1:NX
        if rain_j(I) < 1 || rain_j(I) > nx_rain; continue; end
        qp_t(I, J) = qp(rain_j(I), rain_i(J), itemp);
    end
end