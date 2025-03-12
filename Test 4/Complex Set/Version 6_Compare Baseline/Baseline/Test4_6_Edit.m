clear;
clc;

rng('shuffle'); %Initiate the random number generators with a random seed

%==========================================================================
N = 2;
K = 2;
N_B = 0.01;
N_E = 0.01;
a = 1; % BPSK modulation constellations {-1, 1}

num_MonteCarlo = 1; % Deterministic!
D = 0.3;

SNR_dB = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 11.5, 12];
SNR = 10.^(SNR_dB/10);

P = SNR*N_B;

% % Othorgonal direction deterministic channel
% H_B = [0.21 0.011; 0.09 0.3]; 
% H_E = [-0.21 0.011; -0.09 0.3]; 

% H_B = sqrt(channel_var_B)*randn(K,N)
% H_E = sqrt(channel_var_E)*randn(K,N)

% % Test 1
% N_B = N_E = 0.01, var_B = var_E = 0.01, D = 0.3
% H_B = [0.1458 0.0633; 0.0974 0.0366]; 
% H_E = [-0.0803 -0.0827; 0.0058 0.0530]; 

% Test 2
H_B = [0.0262 0.0049; -0.1598  -0.2414]; 
H_E = [0.0498 0.0194; -0.0446 -0.0758]; 

% % Test 4
% H_B = [0.0262 0.0049 0.0262; -0.1598 -0.1598  -0.2414; -0.1598  -0.2414 0.0049]; 
% H_E = [0.0498 0.0498 0.0194; 0.0498 0.0194 -0.0446; -0.0446 0.0498 0.0194];

for num = 1:num_MonteCarlo 

    for iter = 1:length(P) 
        %==================================================================
        % Search through case 4, 3, and 2 to find the optimal beamforming vector
        
        % Case 4: For lambda1 > 0 and lambda2 > 0
        epsilon4 = 0.00000000001;
        
        lambda1 = 1:0.001:50; 
        
        for j = 1:length(lambda1) 
            [eigenVectors4,eigenValues4] = eig((H_B')*H_B - lambda1(j)*(H_E')*H_E);
        
            index_set4 = [];
        
            for i = 1:N
                condition4 = (norm(H_E*eigenVectors4(:,i)))^2 - (sqrt(N_E/(2*P(iter)))*qfuncinv(D)/abs(a))^2;
        
                if abs(condition4) <= epsilon4 
                    index_set4 = [index_set4, i];
                end
            end
      
            if length(index_set4) >= 2
                lambda1_optimal4 = lambda1(j);
        
                for k = 1:length(index_set4) 
                    objectiveValue4(k) = qfunc(norm(H_B*eigenVectors4(:,index_set4(k))*a)*sqrt(2*P(iter)/N_B));
                end
        
                [minOjectiveValue4,minIndex4] = min(objectiveValue4);
        
                lambda2_optimal4 = eigenValues4(minIndex4,minIndex4);
        
                if (lambda2_optimal4 > 0) 
                    w_optimal4 = eigenVectors4(:,index_set4(minIndex4));
        
                    ErrorProbability_Bob4 = minOjectiveValue4;
                    ErrorProbability_Eve4 = qfunc(norm(H_E*w_optimal4*a)*sqrt(2*P(iter)/N_E));

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
            objectiveValue3 = qfunc(norm(H_B*eigenVectors3(:,i)*a)*sqrt(2*P(iter)/N_B));
        
            % Normalize the eigenvectors
            eigenVectors3(:,i) = eigenVectors3(:,i)/norm(eigenVectors3(:,i))*constant3;
        
            condition3 = (norm(H_E*eigenVectors3(:,i)))^2 - (sqrt(N_E/(2*P(iter)))*qfuncinv(D)/abs(a))^2;
        
            if (abs(condition3) <= epsilon3) 
                index_set3 = [index_set3, i];
            end
        end
        
        if length(index_set3) >= 2
            for k = 1:length(index_set4) 
                objectiveValue3(k) = qfunc(norm(H_B*eigenVectors3(:,index_set3(k))*a)*sqrt(2*P(iter)/N_B));
            end
        
            [minOjectiveValue3,minIndex3] = min(objectiveValue3);
        
            lambda1_optimal3 = eigenValues3(minIndex3,minIndex3);
        
            if (lambda1_optimal3 > 0) 
                w_optimal3 = eigenVectors3(:,minIndex3);
                ErrorProbability_Bob3 = minOjectiveValue3;
                ErrorProbability_Eve3 = qfunc(norm(H_E*w_optimal3*a)*sqrt(2*P(iter)/N_E));

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
            objectiveValue2 = qfunc(norm(H_B*eigenVectors2(:,i)*a)*sqrt(2*P(iter)/N_B));

            condition2 = (norm(H_E*eigenVectors2(:,i)))^2 - (sqrt(N_E/(2*P(iter)))*qfuncinv(D)/abs(a))^2;
            
            if (condition2 < 0) 
                index_set2 = [index_set2, i];
            end
        end
         
        if length(index_set2) >= 2
            for k = 1:length(index_set2) 
                objectiveValue2(k) = qfunc(norm(H_B*eigenVectors2(:,index_set2(k))*a)*sqrt(2*P(iter)/N_B));
            end
            
            [minOjectiveValue2,minIndex2] = min(objectiveValue2);
                
            lambda2_optimal2 = eigenValues2(minIndex2,minIndex2);
                
            if (lambda2_optimal2 > 0) 
                w_optimal2 = eigenVectors2(:,minIndex2);
                ErrorProbability_Bob2 = minOjectiveValue2;
                ErrorProbability_Eve2 = qfunc(norm(H_E*w_optimal2*a)*sqrt(2*P(iter)/N_E));

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
        ErrorProbability_Eve_Array = [ErrorProbability_Eve4, ErrorProbability_Eve3, ErrorProbability_Eve2];
        
        [ErrorProbability_Bob(iter,num),minIndex] = min(ErrorProbability_Bob_Array); 
        ErrorProbability_Eve(iter,num) = ErrorProbability_Eve_Array(minIndex);

        w_optimal_final = w_optimal_Array(:,minIndex);
        
        % =================================================================
        % Comparison: A. Mukherjee and A. L. Swindlehurst, "Robust Beamforming for Security in MIMO Wiretap Channels With Imperfect CSI," in IEEE Transactions on Signal Processing, vol. 59, no. 1, pp. 351-361, Jan. 2011.
        
        % For case: K_E >= N
        scaleFactor = 0.067;
        
        [eigenVectors_add,eigenValues_add] = eig(((H_B')*H_B),((H_E')*H_E));
        
        [max_value_add, max_index_add] = max(max(eigenValues_add));
        
        w_optimal_add = eigenVectors_add(:,max_index_add);

        if (norm(w_optimal_add) > 1) 
            w_optimal_add = w_optimal_add*scaleFactor;
        end

        ErrorProbability_Bob_add(iter,num) = qfunc(norm(H_B*w_optimal_add*a)*sqrt(2*P(iter)/N_B));
        ErrorProbability_Eve_add(iter,num) = qfunc(norm(H_E*w_optimal_add*a)*sqrt(2*P(iter)/N_E));
    end
end

for iter = 1:length(P) 
    ErrorProbability_Bob_average(iter) = sum(ErrorProbability_Bob(iter,:))/num_MonteCarlo;
    ErrorProbability_Eve_average(iter) = sum(ErrorProbability_Eve(iter,:))/num_MonteCarlo;

    ErrorProbability_Bob_add_average(iter) = sum(ErrorProbability_Bob_add(iter,:))/num_MonteCarlo;
    ErrorProbability_Eve_add_average(iter) = sum(ErrorProbability_Eve_add(iter,:))/num_MonteCarlo;
end


%==========================================================================
% Figure
figure(1)
semilogy(SNR_dB,ErrorProbability_Bob_average,'b-o','MarkerSize',7,'LineWidth',1.5);
hold on
semilogy(SNR_dB,ErrorProbability_Bob_add_average,'b-','MarkerSize',7,'LineWidth',1.5);
semilogy(SNR_dB,ErrorProbability_Eve_average,'r--square','MarkerSize',7,'LineWidth',1.5);
semilogy(SNR_dB,ErrorProbability_Eve_add_average,'r--','MarkerSize',7,'LineWidth',1.5);
hold off
grid on
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Bob - our proposed','Bob - [33]','Eve - our proposed ','Eve - [33]','Location','SouthWest');
axis([-10,11.5,0,0.5]);