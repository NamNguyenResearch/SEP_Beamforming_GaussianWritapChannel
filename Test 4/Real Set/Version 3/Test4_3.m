clear;
clc;

%========================================================================
N = 2;
K = 2;
N_B = 0.02;
N_E = 0.01;
a = 1; % BPSK modulation constellations {-1, 1}
channel_var_B = 1;
channel_var_E = 1;

% Setup 2 + 5
H_B = [0.21 0.011; 0.09 0.3]; % N_B = 0.02; N_E = 0.01
H_E = [0.01 0.02; 0.017 0.01];

% % Setup 3
% H_B = [0.21 0.011; 0.09 0.3]; 
% H_E = [-0.01 0.02; 0.01 0.01];

% % Setup 4 -> Investigate!
% H_B = [0.21 0.015; 0.1 0.12];
% H_E = [0.01 0.071; 0.01 0.01];

% Generate the channel matrices and the modulated symbol vectors
% Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices

% H_B = sqrt(channel_var_B/2)*(randn(K,N)+1i*randn(K,N));
% H_E = sqrt(channel_var_E/2)*(randn(K,N)+1i*randn(K,N));

D = 0.2;
SNR = -10:1:25;
P = 10.^(SNR/10)*N_E;

for iter = 1:length(P) 
    %====================================================================
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
  
        if sum(index_set4) > 0 
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

                break
            else
                lambda1_optimal4 = 0.5;
                lambda2_optimal4 = 0.5;
                w_optimal4 = 10^-1;
                ErrorProbability_Bob4 = 0.5;
                ErrorProbability_Eve4 = 0.5;
            end
        else
            lambda1_optimal4 = 0.5;
            lambda2_optimal4 = 0.5;
            w_optimal4 = 10^0;
            ErrorProbability_Bob4 = 0.5;
            ErrorProbability_Eve4 = 0.5;
        end
    end
    
    % Result for case 4

    % lambda1_optimal4
    % lambda2_optimal4
    % w_optimal4
    % ErrorProbability_Bob4
    % ErrorProbability_Eve4
    
    
    %====================================================================
    % Case 3: For lambda1 > 0 and lambda2 = 0
    
    epsilon3 = 0.000001;
    constant3 = 0.9999;
   
    %====================================================================
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
    
    if sum(index_set3) > 0
        for k = 1:length(index_set4) 
            objectiveValue3(k) = qfunc(norm(H_B*eigenVectors3(:,index_set3(k))*a)*sqrt(2*P(iter)/N_B));
        end
    
        [minOjectiveValue3,minIndex3] = min(objectiveValue3);
    
        lambda1_optimal3 = eigenValues3(minIndex3,minIndex3);
    
        if (lambda1_optimal3 > 0) 
            w_optimal3 = eigenVectors3(:,minIndex3);
            ErrorProbability_Bob3 = minOjectiveValue3;
            ErrorProbability_Eve3 = qfunc(norm(H_E*w_optimal3*a)*sqrt(2*P(iter)/N_E));
        else
            lambda1_optimal3 = 0.5;
            w_optimal3 = 0.5;
            ErrorProbability_Bob3 = 0.5;
            ErrorProbability_Eve3 = 0.5;
        end       
    else
        lambda1_optimal3 = 0.5;
        w_optimal3 = 0.5;
        ErrorProbability_Bob3 = 0.5;
        ErrorProbability_Eve3 = 0.5;
    end
    
    % Result for case 3

    % lambda1_optimal3 
    % w_optimal3 
    % ErrorProbability_Bob3 
    % ErrorProbability_Eve3 
    
    
    %====================================================================
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
        
    if sum(index_set2) > 0
        for k = 1:length(index_set2) 
            objectiveValue2(k) = qfunc(norm(H_B*eigenVectors2(:,index_set2(k))*a)*sqrt(2*P(iter)/N_B));
        end
        
    [minOjectiveValue2,minIndex2] = min(objectiveValue2);
        
    lambda2_optimal2 = eigenValues2(minIndex2,minIndex2);
        
    if (lambda2_optimal2 > 0) 
        w_optimal2 = eigenVectors2(:,minIndex2);
        ErrorProbability_Bob2 = minOjectiveValue2;
        ErrorProbability_Eve2 = qfunc(norm(H_E*w_optimal2*a)*sqrt(2*P(iter)/N_E));
    else
        lambda2_optimal2 = 0.5;
        w_optimal2 = 0.5;
        ErrorProbability_Bob2 = 0.5;
        ErrorProbability_Eve2 = 0.5;
    end       
    else
        lambda2_optimal2 = 0.5;
        w_optimal2 = 0.5;
        ErrorProbability_Bob2 = 0.5;
        ErrorProbability_Eve2 = 0.5;
    end
    
    % Result for case 2

    % lambda2_optimal2 
    % w_optimal2 
    % ErrorProbability_Bob2 
    % ErrorProbability_Eve2 
    
    
    %======================================================================
    % Final result
    ErrorProbability_Bob_Array = [ErrorProbability_Bob4, ErrorProbability_Bob3, ErrorProbability_Bob2];
    ErrorProbability_Eve_Array = [ErrorProbability_Eve4, ErrorProbability_Eve3, ErrorProbability_Eve2];
    
    [ErrorProbability_Bob(iter),minIndex] = min(ErrorProbability_Bob_Array); 
    % minIndex
    ErrorProbability_Eve(iter) = ErrorProbability_Eve_Array(minIndex);
end


%==========================================================================
% Figure

figure(1)
plot(SNR,ErrorProbability_Bob,'b--o','LineWidth',1);
hold on
plot(SNR,ErrorProbability_Eve,'r--square','LineWidth',1);
hold off
grid on
xlabel('Average SNR (dB)');
ylabel('Error Probability');
legend('Bob','Eve','Location','SouthWest');