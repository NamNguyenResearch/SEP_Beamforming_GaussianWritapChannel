clear;
clc;
rng('shuffle');

%==========================================================================
% Parameters
N_Comp = 2; % Number of the transmitted antennas
K_Comp = 2; % Number of the received antennas => Not depend on for this case! 

N = 2*N_Comp;
K = 2*K_Comp;
M = 4; % M-ary detection
lambda = 1; % Tradeoff Parameter
P = 1.5; % Given power
N_B = 0.005; % Noise power
N_E = 0.005; % Noise power
channel_var = 1;

num_iter = 4500; 
max_iteration = 35;

EbN0dB_vector = [-20 -15 -10 -5 0 5 10 15 20 25 30 35]; 
Eb = 1; % Because we set the symbol by s=sign(rand(Nt,1)-0.5) below

%==========================================================================
% Initial value 

s_Comp = [1+1i, -1-1i, -1+1i, -1-1i;...
          1-1i,  1-1i,  1-1i, -1+1i]; % => QPSK modulation constellations {+- 1 +- 1i}

s = [real(s_Comp(:,1)) real(s_Comp(:,2)) real(s_Comp(:,3)) real(s_Comp(:,4));...
     imag(s_Comp(:,1)) imag(s_Comp(:,2)) imag(s_Comp(:,3)) imag(s_Comp(:,4))];

% % Define QPSK constellation points
% s1 =  1 + 1j;
% s2 =  1 - 1j;
% s3 = -1 + 1j;
% s4 = -1 - 1j;
% 
% % Create an array of QPSK constellation points
% QPSK_constellation = [s1, s2, s3, s4];
% 
% % Number of QPSK points
% num_points = length(QPSK_constellation);
% 
% % Preallocate a matrix to hold all 2-dimensional vectors
% s_Comp = zeros(num_points^2, 2);
% 
% % Generate all possible 2-dimensional vectors
% index = 1;
% for i = 1:num_points
%     for j = 1:num_points
%         s_Comp(index, :) = [QPSK_constellation(i), QPSK_constellation(j)];
%         index = index + 1;
%     end
% end
% 
% s_Comp = s_Comp.';
% 
% % 16 cases
% s = [real(s_Comp(:,1)) real(s_Comp(:,2)) real(s_Comp(:,3)) real(s_Comp(:,4)) real(s_Comp(:,5)) real(s_Comp(:,6)) real(s_Comp(:,7)) real(s_Comp(:,8)) real(s_Comp(:,9)) real(s_Comp(:,10)) real(s_Comp(:,11)) real(s_Comp(:,12)) real(s_Comp(:,13)) real(s_Comp(:,14)) real(s_Comp(:,15)) real(s_Comp(:,16));...
%      imag(s_Comp(:,1)) imag(s_Comp(:,2)) imag(s_Comp(:,3)) imag(s_Comp(:,4)) imag(s_Comp(:,5)) imag(s_Comp(:,6)) imag(s_Comp(:,7)) imag(s_Comp(:,8)) imag(s_Comp(:,9)) imag(s_Comp(:,10)) imag(s_Comp(:,11)) imag(s_Comp(:,12)) imag(s_Comp(:,13)) imag(s_Comp(:,14)) imag(s_Comp(:,15)) imag(s_Comp(:,16))];

% Preallocate matrices
objective_function_Bob = zeros(max_iteration, num_iter);
objective_function_Eve = zeros(max_iteration, num_iter);

objective_function_Bob_SNR = zeros(length(EbN0dB_vector), num_iter);
objective_function_Eve_SNR = zeros(length(EbN0dB_vector), num_iter);

for iter = 1:num_iter 

    W_0_Comp = complex(rand(N_Comp,N_Comp),rand(N_Comp,N_Comp)); % Generate randomly

    W_0 = [real(W_0_Comp), -imag(W_0_Comp); imag(W_0_Comp), real(W_0_Comp)];

    % Generate the channel matrices and the modulated symbol vectors

    % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
    H_B_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp));
    H_E_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp));
    
    % Real channel matrices
    H_B = [real(H_B_Comp), -imag(H_B_Comp); imag(H_B_Comp), real(H_B_Comp)];
    H_E = [real(H_E_Comp), -imag(H_E_Comp); imag(H_E_Comp), real(H_E_Comp)];
    
    % Optimization process
    W = W_0;
    
    for i = 1:max_iteration 

        W = update(W,lambda,M,N,P,H_B,H_E,s,N_B,N_E); 
        
        objective_function_Bob(i,iter) = objectFunction(M,H_B,s,N_B,W); % Error Probability of Bob
    
        objective_function_Eve(i,iter) = objectFunction(M,H_E,s,N_E,W); % Error Probability of Eve
    end
    
    % Beamforming matrix normalization
    W = W/norm(W,"fro");

    for snr_i = 1:length(EbN0dB_vector) 
        EbN0dB = EbN0dB_vector(snr_i);
        EbN0 = 10^(EbN0dB / 10);
        N0 = Eb / EbN0;

        objective_function_Bob_SNR(snr_i,iter) = objectFunction(M,H_B,s,N0,W); % Error Probability of Bob
    
        objective_function_Eve_SNR(snr_i,iter) = objectFunction(M,H_E,s,N0,W); % Error Probability of Eve
    end
end

% Find mean of probability values
objective_function_Bob_average = mean(objective_function_Bob, 2);
objective_function_Eve_average = mean(objective_function_Eve, 2);

objective_function_Bob_SNR_average = mean(objective_function_Bob_SNR, 2);
objective_function_Eve_SNR_average = mean(objective_function_Eve_SNR, 2);

%==========================================================================
% Figure of the objective value versus iteration 
figure(1);
iterations = 1:max_iteration;
plot(iterations,objective_function_Bob_average,'-b','LineWidth',1.5);
hold on 
plot(iterations,objective_function_Eve_average,'--r','LineWidth',1.5);
hold off
grid on
xlabel('Iterations');
ylabel('Bound of symbol error probability');
legend('Bob','Eve','Location','SouthEast');

figure(2);
semilogy(EbN0dB_vector,objective_function_Bob_SNR_average,'-ob','LineWidth',1.5);
hold on
semilogy(EbN0dB_vector,objective_function_Eve_SNR_average,'-*r','LineWidth',1.5);
hold off
grid on
xlabel('SNR');
ylabel('Symbol Error Probability');
legend('Bob','Eve','Location','SouthWest');
% saveas(gcf,'Figure2_2.fig');