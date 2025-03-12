function [ErrorProbability_Bob, ErrorProbability_Eve] = calculate(N, P, N_B, N_E, H_B, H_E, a, constraint)
    w_optimal = optimizationFunc(H_B, H_E, P, N, constraint);

    norm_term_Bob = norm(H_B * w_optimal * abs(a));
    norm_term_Eve = norm(H_E * w_optimal * abs(a));

    ErrorProbability_Bob = qfunc(norm_term_Bob / sqrt(N_B / 2));
    ErrorProbability_Eve = qfunc(norm_term_Eve / sqrt(N_E / 2));
end
