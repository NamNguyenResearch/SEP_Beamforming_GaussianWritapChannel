function [sum] = objectFunction(M,H,s,N_0,W) 

sum = 0;

    for i = 1:M 
        for j = 1:M 
            if j ~= i 
                sum = sum + 1/M * qfunc(norm(H*W*(s(:,i) - s(:,j)))/(2*sqrt(N_0/2)));
            end
        end
    end
end

