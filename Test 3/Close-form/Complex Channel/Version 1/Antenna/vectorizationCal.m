function vector = vectorizationCal(N,L,matrix)
    vector = zeros(N*L,1);

    for i=1:N
        for j=1:L
            vector((i-1)*L + j) = matrix(i,j);
        end
    end   
end