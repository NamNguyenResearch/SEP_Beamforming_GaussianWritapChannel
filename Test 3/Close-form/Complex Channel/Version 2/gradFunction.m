function [sum] = gradFunction(M,H,s,N_0,W) 
    
    sum = 0;

    for i = 1:M 
        for j = 1:M 
            if j ~= i 
                sum = sum + 1/M * 1/sqrt(2*pi)*exp(-(((s(:,i) - s(:,j)).')*(W.')*(H.')*H*W*(s(:,i) - s(:,j)))/(4*N_0))...
                          * 1/(2*sqrt((((s(:,i) - s(:,j)).')*(W.')*(H.')*H*W*(s(:,i) - s(:,j)))/(2*N_0)))...
                          *((H.')*H*W*(s(:,i) - s(:,j))*((s(:,i) - s(:,j)).'))/(N_0);
            end
        end
    end
end