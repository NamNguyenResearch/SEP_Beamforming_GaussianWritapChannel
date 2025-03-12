function matrix = inverse_vectorizationCal(N,vector)
    matrix = zeros(N);

    for i=1:N
        for j=1:N
            matrix(i,j) = vector((i-1)*N + j);
        end
    end
end