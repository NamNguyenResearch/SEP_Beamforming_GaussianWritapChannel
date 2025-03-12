clear;
clc;

%==========================================================================
N = 2;
K = 2;
P = 1;
N_B = 0.01;
N_E = 0.01;
a = 1; % BPSK modulation constellations {-1, 1}
channel_var_B = 1;
channel_var_E = 1;

% Generate the channel matrices and the modulated symbol vectors

% % Real value determined channel matrices

% % Setup 1
% H_B = [0 0.1];
% H_E = [0.1 0];

% Setup 2 + 5
H_B = [0.21 0.011; 0.09 0.3]; 
H_E = [0.01 0.02; 0.017 0.01];

% % Setup 3
% H_B = [0.21 0.011; 0.09 0.3]; 
% H_E = [-0.01 0.02; 0.01 0.01];

% % Setup 4
% H_B = [0.21 0.015; 0.1 0.12];
% H_E = [0.01 0.071; 0.01 0.01];

% %--------------------------------------------------------------------------
% Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
% H_B = sqrt(channel_var_B/2)*(randn(K,N)+1i*randn(K,N));
% H_E = sqrt(channel_var_E/2)*(randn(K,N)+1i*randn(K,N));

% % Setup 1
% H_B = [0.1176 + 1.4314i -0.1605 - 0.3606i; 0.2661 - 1.6684i -0.8124 - 0.9345i]; 
% H_E = [-0.4498 + 0.5494i 0.0976 + 0.4578i; 0.2248 + 0.4401i -0.5026 - 0.3010i];

% % Setup 2
% H_B = [-0.5212 - 0.0565i 0.6438 + 0.1299i; -1.2374 - 0.6353i 0.6131 + 0.2056i]; 
% H_E = [0.0799 + 0.8250i 0.0719 + 0.8066i; 0.3111 - 1.3112i 1.9709 - 0.7731i];

D = 0.3555;

% %==========================================================================
% % Case 4: For \lambda1 > 0 and \lambda2 > 0
% epsilon4 = 0.0000001;
% 
% % % Find D corresponding to lambda1_setup4
% % lambda1_setup4 = 1.5;
% % 
% % [eigenVectors4_inverse,eigenValues4_inverse] = eig((H_B')*H_B - lambda1_setup4*(H_E')*H_E);
% % 
% % for i = 1:N
% %     objectiveValue4_inverse(i) = qfunc(norm(H_B*eigenVectors4_inverse(:,i)*a)*sqrt(2*P/N_B));
% % end
% % 
% % [minOjectiveValue4_inverse,minIndex4_inverse] = min(objectiveValue4_inverse);
% % w_optimal4_inverse = eigenVectors4_inverse(:,minIndex4_inverse);
% % 
% % ErrorProbability_Eve4_inverse = qfunc(norm(H_E*w_optimal4_inverse*a)*sqrt(2*P/N_E));
% % 
% % % We choose D sastify the first constraint in case of \lambda1 > 0 and \lambda2 > 0
% % D =  ErrorProbability_Eve4_inverse;
% 
% %==========================================================================
% lambda1 = 1:0.001:50; 
% 
% for j = 1:length(lambda1) 
%     [eigenVectors4,eigenValues4] = eig((H_B')*H_B - lambda1(j)*(H_E')*H_E);
% 
%     index_set = [];
% 
%     for i = 1:N
%         condition4 = (norm(H_E*eigenVectors4(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;
% 
%         if abs(condition4) <= epsilon4 
%             index_set = [index_set, i];
%         end
%     end
% 
%     if size(index_set) > 0 
%         lambda1_optimal4 = lambda1(j);
% 
%         for k = 1:length(index_set) 
%             objectiveValue4(k) = qfunc(norm(H_B*eigenVectors4(:,index_set(k))*a)*sqrt(2*P/N_B));
%         end
% 
%         [minOjectiveValue4,minIndex4] = min(objectiveValue4);
% 
%         lambda2_optimal4 = eigenValues4(minIndex4,minIndex4);
% 
%         w_optimal4 = eigenVectors4(:,index_set(minIndex4));
% 
%         ErrorProbability_Bob4 = minOjectiveValue4;
%         ErrorProbability_Eve4 = qfunc(norm(H_E*w_optimal4*a)*sqrt(2*P/N_E));
% 
%         break
%     end
% end
% 
% lambda1_optimal4
% lambda2_optimal4
% w_optimal4
% ErrorProbability_Bob4
% ErrorProbability_Eve4


%==========================================================================
% Case 3: For \lambda1 > 0 and \lambda2 = 0

epsilon3 = 0.00001;
minOjectiveValue3 = 10^10;
constant3 = 0.9999;

% % Find D corresponding to lambda1_setup3
% 
% [eigenVectors3_inverse,eigenValues3_inverse] = eig(((H_B')*H_B),((H_E')*H_E));
% 
% for i = 1:N
%     % Normalize the eigenvectors
%     eigenVectors3_inverse(:,i) = eigenVectors3_inverse(:,i)/norm(eigenVectors3_inverse(:,i))*constant3;
% 
%     objectiveValue3_inverse(i) = qfunc(norm(H_B*eigenVectors3_inverse(:,i)*a)*sqrt(2*P/N_B));
% end
% 
% [minOjectiveValue3_inverse,minIndex3_inverse] = min(objectiveValue3_inverse);
% w_optimal3_inverse = eigenVectors3_inverse(:,minIndex3_inverse);
% 
% ErrorProbability_Eve3_inverse = qfunc(norm(H_E*w_optimal3_inverse*a)*sqrt(2*P/N_E));
% 
% % We choose D sastify the first constraint in case of \lambda1 > 0 and \lambda2 = 0
% D =  ErrorProbability_Eve3_inverse


%==========================================================================
[eigenVectors3,eigenValues3] = eig(((H_B')*H_B),((H_E')*H_E));

for i = 1:N
    objectiveValue3 = qfunc(norm(H_B*eigenVectors3(:,i)*a)*sqrt(2*P/N_B));

    % Normalize the eigenvectors
    eigenVectors3(:,i) = eigenVectors3(:,i)/norm(eigenVectors3(:,i))*constant3;

    condition3 = (norm(H_E*eigenVectors3(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;

    if (abs(condition3) <= epsilon3) && (objectiveValue3 < minOjectiveValue3)
        minOjectiveValue3 = objectiveValue3;
        minIndex3 = i;
    end
end

lambda1_optimal3 = eigenValues3(minIndex3,minIndex3);
w_optimal3 = eigenVectors3(:,minIndex3);
ErrorProbability_Bob3 = minOjectiveValue3;
ErrorProbability_Eve3 = qfunc(norm(H_E*w_optimal3*a)*sqrt(2*P/N_E));

lambda1_optimal3
w_optimal3
ErrorProbability_Bob3
ErrorProbability_Eve3


%==========================================================================
% Case 2: For \lambda1 = 0 and \lambda2 > 0

% D = 2e-04;

minOjectiveValue2 = 10^10;

[eigenVectors2,eigenValues2] = eig((H_B')*H_B);

for i = 1:N
    objectiveValue2 = qfunc(norm(H_B*eigenVectors2(:,i)*a)*sqrt(2*P/N_B));
    condition2 = (norm(H_E*eigenVectors2(:,i)))^2 - (sqrt(N_E/(2*P))*qfuncinv(D)/abs(a))^2;

    if (condition2 < 0) && (objectiveValue2 < minOjectiveValue2)
        minOjectiveValue2 = objectiveValue2;
        minIndex2 = i;
    end
end

lambda2_optimal2 = eigenValues2(minIndex2,minIndex2);

lambda1_optimal2 = 0;
w_optimal2 = eigenVectors2(:,minIndex2);
ErrorProbability_Bob2 = minOjectiveValue2;
ErrorProbability_Eve2 = qfunc(norm(H_E*w_optimal2*a)*sqrt(2*P/N_E));

lambda1_optimal2
lambda2_optimal2
w_optimal2
ErrorProbability_Bob2
ErrorProbability_Eve2