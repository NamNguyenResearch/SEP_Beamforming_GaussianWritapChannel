function sum = objectFunction(M, H, s, N_0, W)
    sum = 0;

    for i = 1:M
        for j = 1:M
            if j ~= i
                diff_s = s(:, i) - s(:, j);
                argument = norm(H * W * diff_s)/(2*sqrt(N_0/2));
                sum = sum + qfunc(argument) / M;
            end
        end
    end
end