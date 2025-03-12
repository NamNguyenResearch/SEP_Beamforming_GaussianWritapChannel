function [sum] = objectFunction(M,H,s,N_0,W) 

sum = 0;

    for i = 1:M 
        for j = 1:M 
            if j ~= i 
                sum = sum + 1/M * qfunc(sqrt((((s(:,i) - s(:,j)).')*(W.') ...
                      *(H.')*H*W*(s(:,i) - s(:,j)))/(2*N_0)));
            end
        end
    end
end

