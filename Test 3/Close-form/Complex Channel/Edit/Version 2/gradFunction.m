function sum = gradFunction(M, H, s, N_0, W)
    sum = 0;

    for i = 1:M
        for j = 1:M
            if j ~= i
                diff_s = s(:, i) - s(:, j);
                term1 = exp(-((diff_s') * (W') * (H') * H * W * diff_s) / (4 * N_0)) / sqrt(2 * pi);
                term2_denom = 2 * sqrt(((diff_s') * (W') * (H') * H * W * diff_s) / (2 * N_0));
                term2 = 1 / term2_denom;
                term3_num = (H') * H * W * diff_s * (diff_s');
                term3_denom = N_0;

                sum = sum + 1/M * term1 * term2 * term3_num / term3_denom;
            end
        end
    end
end