clear;
clc;

rng('shuffle'); %Initiate the random number generators with a random seed

%==========================================================================
N = 2;
K = 2;
N_B = 1;
N_E = 1;
a = 1; % BPSK modulation constellations {-1, 1}
channel_var_B = 1;
channel_var_E = 1;

num_MonteCarlo = 100; 

P_dBm = 10;
P = 10.^(P_dBm/10);

D = 0:0.1:0.5;

for num = 1:num_MonteCarlo
    % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
    H_B = sqrt(channel_var_B/2)*(randn(K,N)+1i*randn(K,N));
    H_E = sqrt(channel_var_E/2)*(randn(K,N)+1i*randn(K,N));

    for iter = 1:length(D) 
        %======================================================================
        % Search through case 4, 3, and 2 to find the optimal beamforming vector
        
        % Case 4: For lambda1 > 0 and lambda2 > 0
        epsilon4 = 0.00000000001;
        
        lambda1 = 1:0.001:50; 
        
        for j = 1:length(lambda1) 
            [eigenVectors4,eigenValues4] = eig((H_B')*H_B - lambda1(j)*(H_E')*H_E);
        
            index_set4 = [];
        
            for i = 1:N
                condition4 = (norm(H_E*eigenVectors4(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D(iter))/abs(a))^2;
        
                if abs(condition4) <= epsilon4 
                    index_set4 = [index_set4, i];
                end
            end
      
            if sum(index_set4) > 0 
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
    
                    break
                else
                    ErrorProbability_Bob4 = NaN;
                    ErrorProbability_Eve4 = NaN;
                end
            else
                ErrorProbability_Bob4 = NaN;
                ErrorProbability_Eve4 = NaN;
            end
        end
        
        
        %======================================================================
        % Case 3: For lambda1 > 0 and lambda2 = 0
        
        epsilon3 = 0.000001;
        constant3 = 0.9999;
       
        %======================================================================
        [eigenVectors3,eigenValues3] = eig(((H_B')*H_B),((H_E')*H_E));
        
        index_set3 = [];
        
        for i = 1:N
            objectiveValue3 = qfunc(norm(H_B*eigenVectors3(:,i)*a)*sqrt(2*P/N_B));
        
            % Normalize the eigenvectors
            eigenVectors3(:,i) = eigenVectors3(:,i)/norm(eigenVectors3(:,i))*constant3;
        
            condition3 = (norm(H_E*eigenVectors3(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D(iter))/abs(a))^2;
        
            if (abs(condition3) <= epsilon3) 
                index_set3 = [index_set3, i];
            end
        end
        
        if sum(index_set3) > 0
            for k = 1:length(index_set4) 
                objectiveValue3(k) = qfunc(norm(H_B*eigenVectors3(:,index_set3(k))*a)*sqrt(2*P/N_B));
            end
        
            [minOjectiveValue3,minIndex3] = min(objectiveValue3);
        
            lambda1_optimal3 = eigenValues3(minIndex3,minIndex3);
        
            if (lambda1_optimal3 > 0) 
                w_optimal3 = eigenVectors3(:,minIndex3);
                ErrorProbability_Bob3 = minOjectiveValue3;
                ErrorProbability_Eve3 = qfunc(norm(H_E*w_optimal3*a)*sqrt(2*P/N_E));
            else
                ErrorProbability_Bob3 = NaN;
                ErrorProbability_Eve3 = NaN;
            end       
        else
            ErrorProbability_Bob3 = NaN;
            ErrorProbability_Eve3 = NaN;
        end
        
        
        %======================================================================
        % Case 2: For lambda1 = 0 and lambda2 > 0
        [eigenVectors2,eigenValues2] = eig((H_B')*H_B);
        
        index_set2 = [];
            
        for i = 1:N
            objectiveValue2 = qfunc(norm(H_B*eigenVectors2(:,i)*a)*sqrt(2*P/N_B));
            condition2 = (norm(H_E*eigenVectors2(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D(iter))/abs(a))^2;
            
            if (condition2 < 0) 
                index_set2 = [index_set2, i];
            end
        end
            
        if sum(index_set2) > 0
            for k = 1:length(index_set2) 
                objectiveValue2(k) = qfunc(norm(H_B*eigenVectors2(:,index_set2(k))*a)*sqrt(2*P/N_B));
            end
            
            [minOjectiveValue2,minIndex2] = min(objectiveValue2);
            
            lambda2_optimal2 = eigenValues2(minIndex2,minIndex2);
            
            if (lambda2_optimal2 > 0) 
                w_optimal2 = eigenVectors2(:,minIndex2);
                ErrorProbability_Bob2 = minOjectiveValue2;
                ErrorProbability_Eve2 = qfunc(norm(H_E*w_optimal2*a)*sqrt(2*P/N_E));
            else
                ErrorProbability_Bob2 = NaN;
                ErrorProbability_Eve2 = NaN;
            end       
        else
            ErrorProbability_Bob2 = NaN;
            ErrorProbability_Eve2 = NaN;
        end
        
        
        %======================================================================
        % Final result
        ErrorProbability_Bob_Array = [ErrorProbability_Bob4, ErrorProbability_Bob3, ErrorProbability_Bob2];
        ErrorProbability_Eve_Array = [ErrorProbability_Eve4, ErrorProbability_Eve3, ErrorProbability_Eve2];
        
        [ErrorProbability_Bob(iter,num),minIndex] = min(ErrorProbability_Bob_Array); 
        % minIndex
        ErrorProbability_Eve(iter,num) = ErrorProbability_Eve_Array(minIndex);
    end
end

for iter = 1:length(D) 
    ErrorProbability_Bob_average(iter) = sum(ErrorProbability_Bob(iter,:))/num_MonteCarlo;
    ErrorProbability_Eve_average(iter) = sum(ErrorProbability_Eve(iter,:))/num_MonteCarlo;
end


%==========================================================================
% Figure
figure(1)
plot(D,ErrorProbability_Bob_average,'b-o','MarkerSize',7,'LineWidth',1.5);
hold on
plot(D,ErrorProbability_Eve_average,'r--square','MarkerSize',7,'LineWidth',1.5);
hold off
grid on
xlabel('D');
ylabel('Error Probability');
legend('Bob','Eve','Location','SouthEast');