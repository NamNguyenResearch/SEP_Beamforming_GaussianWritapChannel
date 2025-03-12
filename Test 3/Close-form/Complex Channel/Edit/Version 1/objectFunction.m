function sum = objectFunction(M, H, s, N_0, W)
    sum = 0;

    for i = 1:M
        for j = 1:M
            if j ~= i
                diff_s = s(:, i) - s(:, j);
                argument = ((diff_s.') * (W.') * (H.') * H * W * diff_s) / (2 * N_0);
                sum = sum + qfunc(sqrt(argument)) / M;
            end
        end
    end
end