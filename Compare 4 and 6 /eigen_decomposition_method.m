function [w_optimal_final] = eigen_decomposition_method(H_B, H_E, N, N_B, N_E, P, D, a)
    %==================================================================
        % Search through case 4, 3, and 2 to find the optimal beamforming vector
        
        % Case 4: For lambda1 > 0 and lambda2 > 0
        epsilon4 = 0.00000000001;
        
        lambda1 = 1:0.001:50; 
        
        for j = 1:length(lambda1) 
            [eigenVectors4,eigenValues4] = eig((H_B')*H_B - lambda1(j)*(H_E')*H_E);
        
            index_set4 = [];
        
            for i = 1:N
                condition4 = (norm(H_E*eigenVectors4(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;
        
                if abs(condition4) <= epsilon4 
                    index_set4 = [index_set4, i];
                end
            end
      
            if length(index_set4) >= 2
                lambda1_optimal4 = lambda1(j);
        
                for k = 1:length(index_set4) 
                    objectiveValue4(k) = qfunc(norm(H_B*eigenVectors4(:,index_set4(k))*a)*sqrt(2*P/N_B));
                end
        
                [minOjectiveValue4,minIndex4] = min(objectiveValue4);
        
                lambda2_optimal4 = eigenValues4(minIndex4,minIndex4);
        
                if (lambda2_optimal4 > 0) 
                    w_optimal4 = eigenVectors4(:,index_set4(minIndex4));
        
                    ErrorProbability_Bob4 = minOjectiveValue4;
                    ErrorProbability_Eve4 = qfunc(norm(H_E*w_optimal4*a)*sqrt(2*P/N_E));

                    w_optimal_Array(:,1) = w_optimal4;

                    % flag4 = 1
                    break
                else
                    ErrorProbability_Bob4 = NaN;
                    ErrorProbability_Eve4 = NaN;

                    w_optimal_Array(:,1) = zeros(N,1);
                end
            else
                ErrorProbability_Bob4 = NaN;
                ErrorProbability_Eve4 = NaN;

                w_optimal_Array(:,1) = zeros(N,1);
            end
        end
        
        
        %==================================================================
        % Case 3: For lambda1 > 0 and lambda2 = 0
        
        epsilon3 = 0.000001;
        constant3 = 0.9999;
       
        %==================================================================
        [eigenVectors3,eigenValues3] = eig(((H_B')*H_B),((H_E')*H_E));
        
        index_set3 = [];
        
        for i = 1:N
            objectiveValue3 = qfunc(norm(H_B*eigenVectors3(:,i)*a)*sqrt(2*P/N_B));
        
            % Normalize the eigenvectors
            eigenVectors3(:,i) = eigenVectors3(:,i)/norm(eigenVectors3(:,i))*constant3;
        
            condition3 = (norm(H_E*eigenVectors3(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;
        
            if (abs(condition3) <= epsilon3) 
                index_set3 = [index_set3, i];
            end
        end
        
        if length(index_set3) >= 2
            for k = 1:length(index_set4) 
                objectiveValue3(k) = qfunc(norm(H_B*eigenVectors3(:,index_set3(k))*a)*sqrt(2*P/N_B));
            end
        
            [minOjectiveValue3,minIndex3] = min(objectiveValue3);
        
            lambda1_optimal3 = eigenValues3(minIndex3,minIndex3);
        
            if (lambda1_optimal3 > 0) 
                w_optimal3 = eigenVectors3(:,minIndex3);
                ErrorProbability_Bob3 = minOjectiveValue3;
                ErrorProbability_Eve3 = qfunc(norm(H_E*w_optimal3*a)*sqrt(2*P/N_E));

                w_optimal_Array(:,2) = w_optimal3;

                % flag3 = 1
            else
                ErrorProbability_Bob3 = NaN;
                ErrorProbability_Eve3 = NaN;

                w_optimal_Array(:,2) = zeros(N,1);
            end       
        else
            ErrorProbability_Bob3 = NaN;
            ErrorProbability_Eve3 = NaN;

            w_optimal_Array(:,2) = zeros(N,1);
        end
        
        
        %==================================================================
        % Case 2: For lambda1 = 0 and lambda2 > 0
        [eigenVectors2,eigenValues2] = eig((H_B')*H_B);
        
        index_set2 = [];
            
        for i = 1:N
            objectiveValue2 = qfunc(norm(H_B*eigenVectors2(:,i)*a)*sqrt(2*P/N_B));

            condition2 = (norm(H_E*eigenVectors2(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;
            
            if (condition2 < 0) 
                index_set2 = [index_set2, i];
            end
        end
         
        if length(index_set2) >= 2
            for k = 1:length(index_set2) 
                objectiveValue2(k) = qfunc(norm(H_B*eigenVectors2(:,index_set2(k))*a)*sqrt(2*P/N_B));
            end
            
            [minOjectiveValue2,minIndex2] = min(objectiveValue2);
                
            lambda2_optimal2 = eigenValues2(minIndex2,minIndex2);
                
            if (lambda2_optimal2 > 0) 
                w_optimal2 = eigenVectors2(:,minIndex2);
                ErrorProbability_Bob2 = minOjectiveValue2;
                ErrorProbability_Eve2 = qfunc(norm(H_E*w_optimal2*a)*sqrt(2*P/N_E));

                w_optimal_Array(:,3) = w_optimal2;

                % flag2 = 1
            else
                ErrorProbability_Bob2 = NaN;
                ErrorProbability_Eve2 = NaN;

                w_optimal_Array(:,3) = zeros(N,1);
            end       
        else
            ErrorProbability_Bob2 = NaN;
            ErrorProbability_Eve2 = NaN;

            w_optimal_Array(:,3) = zeros(N,1);
        end
        
        
        %==================================================================
        % Final result
       
        ErrorProbability_Bob_Array = [ErrorProbability_Bob4, ErrorProbability_Bob3, ErrorProbability_Bob2];
        
        [ErrorProbability_Bob,minIndex] = min(ErrorProbability_Bob_Array); 

        w_optimal_final = w_optimal_Array(:,minIndex);
end