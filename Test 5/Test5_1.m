clear;
clc;

%==========================================================================
% Parameters
N = 2; % Number of the transmitted antennas
K = 2; % Number of the received antennas => Not depend on for this case! 

P = 1;   
N_B = 0.01;
N_E = 0.1;
channel_var = 1;

H_B2 = [0  1i]; 
H_E2 = [1i -1i; 0 0];

s_0 = [1; 1]; % BPSK modulation constellations {-1, 1}]
s_1 = [-1; 1];

%==========================================================================
D = 0.00000000001:0.01:(0.5+0.00000000001); % Take care the condtion for D since we use Q-function
num_iter = 100; 

constraint = 2*sqrt(N_B/2)*qfuncinv(D);

for iter = 1:num_iter
    % Generate the channel matrices and the modulated symbol vectors

    % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
    H_B1 = sqrt(channel_var/2)*(randn(1,N)+1i*randn(1,N));
    H_E1 = sqrt(channel_var/2)*(randn(K,N)+1i*randn(K,N));

    for i = 1:length(D) 
        [W_optimal1] = optimizationFunc(H_B1,H_E1,s_0,s_1,P,N,constraint(i));

        ErrorProbability_Bob_1(i,iter) = qfunc(norm(H_B1*W_optimal1*(s_0-s_1))/(2*sqrt(N_B/2)));
        ErrorProbability_Eve_1(i,iter) = qfunc(norm(H_E1*W_optimal1*(s_0-s_1))/(2*sqrt(N_E/2))); 
    end
end

for i = 1:length(D) 
    ErrorProbability_Bob_average1(i) = sum(ErrorProbability_Bob_1(i,:))/num_iter;
    ErrorProbability_Eve_average1(i) = sum(ErrorProbability_Eve_1(i,:))/num_iter;

    [W_optimal2] = optimizationFunc(H_B2,H_E2,s_0,s_1,P,N,constraint(i));
    
    ErrorProbability_Bob_2(i) = qfunc(norm(H_B2*W_optimal2*(s_0-s_1))/(2*sqrt(N_B/2)));
    ErrorProbability_Eve_2(i) = qfunc(norm(H_E2*W_optimal2*(s_0-s_1))/(2*sqrt(N_E/2)));
end

%==========================================================================
figure(1)
subplot(2,1,1)
plot(D,ErrorProbability_Bob_average1);
grid on
hold on
plot(D,ErrorProbability_Eve_average1);
hold off
xlabel('D');
ylabel('Error Probability');
legend('Bob','Eve');
title('(a)');
axis([0,0.5,0,0.5]);

subplot(2,1,2)
plot(D,ErrorProbability_Bob_2);
grid on
hold on
plot(D,ErrorProbability_Eve_2);
hold off
xlabel('D');
ylabel('Error Probability');
legend('Bob','Eve');
title('(b)');
axis([0,0.5,0,0.5]);