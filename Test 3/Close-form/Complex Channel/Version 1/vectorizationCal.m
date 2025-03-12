function vector = vectorizationCal(N,matrix)
    vector = zeros(N*N,1);

    for i=1:N
        for j=1:N
            vector((i-1)*N + j) = matrix(i,j);
        end
    end   
end