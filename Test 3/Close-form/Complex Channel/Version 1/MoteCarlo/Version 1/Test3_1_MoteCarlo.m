clear;
clc;
%rng('shuffle');

%==========================================================================
% Parameters
N_Comp = 2; % Number of the transmitted antennas
K_Comp = 2; % Number of the received antennas => Not depend on for this case! 

N = 2*N_Comp;
K = 2*K_Comp;
M = 4; % M-ary detection, QPSK modulation
lambda = 1; % Tradeoff Parameter
P = 1.5; % Given power
N_B = 0.005; % Noise power
N_E = 0.005; % Noise power
channel_var = 1;

num_iter = 4500; 
max_iteration = 30;

% global W_optimal;

%==========================================================================
% Initial value 

% First approach
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

        [W_update] = update(W,lambda,M,N,P,H_B,H_E,s,N_B,N_E); 

        objective_function_Bob(i,iter) = objectFunction(M,H_B,s,N_B,W_update); % Error Probability of Bob

        objective_function_Eve(i,iter) = objectFunction(M,H_E,s,N_E,W_update); % Error Probability of Eve

        objective_function(i,iter) = objective_function_Bob(i,iter) - lambda * objective_function_Eve(i,iter); % Objective Value

        W = W_update;
    end 

    % W_optimal(:,:,iter) = W;
end

% Find mean of probability values
for i=1:max_iteration

    objective_function_Bob_average(i) = sum(objective_function_Bob(i,:))/num_iter;
    objective_function_Eve_average(i) = sum(objective_function_Eve(i,:))/num_iter;
    objective_function_average(i) = sum(objective_function(i,:))/num_iter;
end

% W_optimal_average = zeros(N,N);
% 
% for i=1:num_iter 
%     W_optimal_average = W_optimal(:,:,i) + W_optimal_average;
% end
% 
% W_optimal_average = W_optimal_average/num_iter;


%==========================================================================
% Monte Carlo simulation

% EbN0dB_vector = [-10 -8 -6 -4 -2 0 2 4 6 8 12 14 16 18 20 22 25]; 
EbN0dB_vector = [0 2 4 5 8 10 14 15 18 20];  

Eb = 1; % Because we set the symbol by s=sign(rand(Nt,1)-0.5) below
Nerrs_stop = 10; 
BER = zeros(1,length(EbN0dB_vector)); % Preallocate BER array
SER = zeros(1,length(EbN0dB_vector)); % Preallocate SER array

% Beamforming matrix 
% Extract the real and imaginary parts from the block matrix
realPart_optimal = W_optimal_average(1:N_Comp, 1:N_Comp);
imaginaryPart_optimal = W_optimal_average(N_Comp+1:end, 1:N_Comp);

% Reconstruct the complex matrix
W_optimal_complex = realPart_optimal + 1i * imaginaryPart_optimal;

for snr_i = 1:length(EbN0dB_vector)
    EbN0dB = EbN0dB_vector(snr_i);
    EbN0 = 10^(EbN0dB / 10);
    N0 = Eb / EbN0;

    Nerrs_ml = 0;
    Nbits = 0;

    Nerrs_ml_symb = 0;
    Nsymbols = 0;

    % Precompute constants
    noise_scale = sqrt(N0/2);
    H_scale = sqrt(channel_var/K_Comp);

    while Nerrs_ml < Nerrs_stop

        % Generate random symbols 

        % random_variable = randi([1, 4]);
        % syms_tr = s_Comp(:,random_variable);

        syms_tr = sign(rand(N_Comp, 1) - 0.5) + 1j * sign(rand(N_Comp, 1) - 0.5);

        %Alcie - Bob
        % Generate the complex Gaussian noise for Bob
        realPart_B = randn(K_Comp, 1);
        imaginaryPart_B = randn(K_Comp, 1);
        n_B = noise_scale*(realPart_B + 1i * imaginaryPart_B);

        % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel matrices
        H_B_Comp_MC = H_scale*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp));

        y_B = H_B_Comp_MC*syms_tr + n_B;

        % Compute estimated symbols
        syms_tr_hat_ml = ML_Detection_QPSK(y_B, H_B_Comp_MC);

        % Calculate bit errors
        Nerrs_ml = Nerrs_ml + sum(real(syms_tr_hat_ml) ~= real(syms_tr))...
                 + sum(imag(syms_tr_hat_ml) ~= imag(syms_tr));
        Nbits = Nbits + 2 * N_Comp;

        % Calculate symbol errors
        Nerrs_ml_symb = Nerrs_ml_symb + sum(syms_tr_hat_ml ~= syms_tr);
        Nsymbols = Nsymbols + 1;

    end

    % Calculate Bit Error Rate
    BER(snr_i) = Nerrs_ml / Nbits; 

    % Calculate Symbol Error Rate
    SER(snr_i) = Nerrs_ml_symb / Nsymbols; 
end

%==========================================================================
%Figure of the objective value versus iteration 
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
semilogy(EbN0dB_vector,BER,'-*b','LineWidth',1.5)
semilogy(EbN0dB_vector,SER,'-ok','LineWidth',1.5)
grid on
xlabel('SNR');
ylabel('BER');