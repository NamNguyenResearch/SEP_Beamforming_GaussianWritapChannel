clear;
clc;

rng('shuffle'); %Initiate the random number generators with a random seed

% N_B = 0.01, N_E = 0.03, D = 0.2, var_B = var_E = 0.01 -> Fig.1
% N_B = 0.01, N_E = 0.03, D = 0.3, var_B = var_E = 0.01 -> Fig.2
% N_B = 0.01, N_E = 0.03, D = 0.3, var_B = var_E = 0.01 -> Fig.3
% N_B = 0.01, N_E = 0.01, D = 0.3, var_B = var_E = 0.01 -> Fig.4

%==========================================================================
K = 2;
N_B = 0.01;
N_E = 0.01;
a = 1; % BPSK modulation constellations {-1, 1}
channel_var_B = 0.01;
channel_var_E = 0.01;

num_MonteCarlo = 1; 
D = 0.3;

SNR_dB = -10:2.5:5;
SNR = 10.^(SNR_dB/10);

P = SNR*N_B;

N = [2, 5];

% % Test 1
% H_B = [-0.0980 + 0.0517i  -0.1005 - 0.0514i; -0.0189 - 0.0835i  -0.0824 - 0.0608i]
% H_E = [-0.0280 + 0.0596i   0.0170 - 0.0553i; 0.1377 + 0.1336i   0.0156 - 0.0072i]
% 
% H_B = [0.0767 + 0.1229i  -0.0562 - 0.0293i   0.0145 + 0.0939i  -0.0588 - 0.0956i  -0.0358 - 0.0170i; -0.0631 + 0.0972i  -0.0335 + 0.0308i  -0.0334 + 0.0672i  -0.0284 + 0.0380i  -0.0088 + 0.0061i]
% H_E = [0.1977 - 0.0740i  -0.0566 + 0.0685i  -0.0369 - 0.0236i   0.0120 + 0.0755i  -0.0578 - 0.0671i; -0.0266 - 0.0467i   0.0870 - 0.1119i  -0.0094 + 0.0034i   0.0206 - 0.0661i   0.0133 + 0.0990i]

% % Test 2
% H_B = [0.0578 - 0.0746i  -0.0272 + 0.0709i; 0.0107 - 0.0635i   0.0038 + 0.0275i]
% H_E = [0.0487 + 0.0785i  -0.0057 - 0.0455i; 0.0108 + 0.0409i   0.0208 + 0.0335i]
% 
% H_B = [-0.0947 + 0.0119i  -0.0606 - 0.1411i  -0.1166 + 0.0315i  -0.1289 + 0.0850i  -0.0343 - 0.0430i; 0.0464 - 0.0002i   0.0273 - 0.0979i  -0.0404 - 0.0734i   0.0310 + 0.0485i  -0.0991 + 0.0528i]
% H_E = [-0.0265 - 0.0013i  -0.0830 + 0.0723i   0.0953 + 0.0792i   0.0393 - 0.0754i   0.0289 + 0.0506i; -0.0719 + 0.0281i   0.1470 - 0.0686i  -0.0016 + 0.0042i  -0.0673 + 0.0552i   0.0036 + 0.0353i]


for num = 1:num_MonteCarlo 

    for numAten = 1:length(N)

        % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
        H_B = sqrt(channel_var_B/2)*(randn(K,N(numAten))+1i*randn(K,N(numAten)))
        H_E = sqrt(channel_var_E/2)*(randn(K,N(numAten))+1i*randn(K,N(numAten)))

        for iter = 1:length(P) 
            %==============================================================
            % Search through case 4, 3, and 2 to find the optimal beamforming vector
            
            % Case 4: For lambda1 > 0 and lambda2 > 0
            epsilon4 = 0.00000000001;
            
            lambda1 = 1:0.001:100; 
            
            for j = 1:length(lambda1) 
                [eigenVectors4,eigenValues4] = eig((H_B')*H_B - lambda1(j)*(H_E')*H_E);
            
                index_set4 = [];
            
                for i = 1:N(numAten)
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
                        ErrorProbability_Bob4 = NaN;
                        ErrorProbability_Eve4 = NaN;
                    end
                else
                    ErrorProbability_Bob4 = NaN;
                    ErrorProbability_Eve4 = NaN;
                end
            end
            
            
            %==============================================================
            % Case 3: For lambda1 > 0 and lambda2 = 0
            
            epsilon3 = 0.000001;
            constant3 = 0.9999;
           
            %==============================================================
            [eigenVectors3,eigenValues3] = eig(((H_B')*H_B),((H_E')*H_E));
            
            index_set3 = [];
            
            for i = 1:N(numAten)
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
                    ErrorProbability_Bob3 = NaN;
                    ErrorProbability_Eve3 = NaN;
                end       
            else
                ErrorProbability_Bob3 = NaN;
                ErrorProbability_Eve3 = NaN;
            end
            
            
            %==============================================================
            % Case 2: For lambda1 = 0 and lambda2 > 0
            [eigenVectors2,eigenValues2] = eig((H_B')*H_B);
            
            index_set2 = [];
                
            for i = 1:N(numAten)
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
                    ErrorProbability_Bob2 = NaN;
                    ErrorProbability_Eve2 = NaN;
                end       
            else
                ErrorProbability_Bob2 = NaN;
                ErrorProbability_Eve2 = NaN;
            end
            
            
            %==============================================================
            % Final result
            ErrorProbability_Bob_Array = [ErrorProbability_Bob4, ErrorProbability_Bob3, ErrorProbability_Bob2];
            ErrorProbability_Eve_Array = [ErrorProbability_Eve4, ErrorProbability_Eve3, ErrorProbability_Eve2];
            
            [ErrorProbability_Bob(iter,numAten,num),minIndex] = min(ErrorProbability_Bob_Array); 
            ErrorProbability_Eve(iter,numAten,num) = ErrorProbability_Eve_Array(minIndex);
        end
    end
end

for numAten = 1:length(N) 
    for iter = 1:length(P) 
        ErrorProbability_Bob_average(iter,numAten) = sum(ErrorProbability_Bob(iter,numAten,:))/num_MonteCarlo;
        ErrorProbability_Eve_average(iter,numAten) = sum(ErrorProbability_Eve(iter,numAten,:))/num_MonteCarlo;
    end
end


%==========================================================================
% Figure
figure(1)
semilogy(SNR_dB,ErrorProbability_Bob_average(:,1),'b--',SNR_dB,ErrorProbability_Eve_average(:,1),'r--','MarkerSize',7,'LineWidth',1.5);
hold on
semilogy(SNR_dB,ErrorProbability_Bob_average(:,2),'b-o',SNR_dB,ErrorProbability_Eve_average(:,2),'r-square','MarkerSize',7,'LineWidth',1.5);
hold off
grid on
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Bob, N = 2','Eve, N = 2','Bob, N = 5','Eve, N = 5','Location','SouthWest');
% axis([-10,5,0.2,0.5]);
% saveas(gcf,'Figure2.fig');