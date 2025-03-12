clear;
clc;

%==========================================================================
% Parameters
N = 2; % Number of the transmitted antennas
K = 2; % Number of the received antennas => Not depend on for this case!   
N_B = 0.01;
N_E = 0.01;
channel_var_B = 1;
channel_var_E = 1;

s_0 = [1; 1]; % BPSK modulation constellations {-1, 1}]
s_1 = [-1; 1];

%==========================================================================

num_MonteCarlo = 10; 
D = 0.00001;

SNR_dB = -10:2.5:10;
SNR = 10.^(SNR_dB/10);

P = SNR*N_B;

constraint = 2*sqrt(N_B/2)*qfuncinv(D);

for iter = 1:num_MonteCarlo
    % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
    H_B1 = sqrt(channel_var_B/2)*(randn(1,N)+1i*randn(1,N));
    H_E1 = sqrt(channel_var_E/2)*(randn(K,N)+1i*randn(K,N));

    for i = 1:length(P) 
        [W_optimal1] = optimizationFunc(H_B1,H_E1,s_0,s_1,P(i),N,constraint);

        ErrorProbability_Bob_1(i,iter) = qfunc(norm(H_B1*W_optimal1*(s_0-s_1))/(2*sqrt(N_B/2)));
        ErrorProbability_Eve_1(i,iter) = qfunc(norm(H_E1*W_optimal1*(s_0-s_1))/(2*sqrt(N_E/2))); 
    end
end

for i = 1:length(D) 
    ErrorProbability_Bob_average1(i) = sum(ErrorProbability_Bob_1(i,:))/num_MonteCarlo;
    ErrorProbability_Eve_average1(i) = sum(ErrorProbability_Eve_1(i,:))/num_MonteCarlo;
end

%==========================================================================
figure(1)
plot(SNR_dB,ErrorProbability_Bob_average1);
grid on
hold on
plot(SNR_dB,ErrorProbability_Eve_average1);
hold off
xlabel('SNR (dB)');
ylabel('Error Probability');
legend('Bob','Eve');