function matrix = inverse_vectorizationCal(N,L,vector)
    matrix = zeros(N,L);

    for i=1:N
        for j=1:L
            matrix(i,j) = vector((i-1)*L + j);
        end
    end
end