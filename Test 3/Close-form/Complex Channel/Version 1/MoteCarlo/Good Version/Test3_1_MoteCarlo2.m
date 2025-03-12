clear;
clc;
rng('shuffle');

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

% Monte-Carlo Parameters  
EbN0dB_vector = [0 5 10 15 20 25 30]; 
Eb = 1; % Because we set the symbol by s=sign(rand(Nt,1)-0.5) below
channel_var_Sim = 0.5;
Nsyms = 100; % Number of symbols 

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

% Preallocate arrays for speed
objective_function_Bob = zeros(max_iteration, num_iter);
objective_function_Eve = zeros(max_iteration, num_iter);

% objective_function_Bob_noBeam = zeros(num_iter);
% objective_function_Eve_noBeam = zeros(num_iter);

BER_B = zeros(length(EbN0dB_vector), num_iter);
BER_E = zeros(length(EbN0dB_vector), num_iter);

% BER_B_noBeam = zeros(length(EbN0dB_vector), num_iter);
% BER_E_noBeam = zeros(length(EbN0dB_vector), num_iter);

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

        W = W_update;
    end

    %==========================================================================
    % Monte Carlo simulation
   
    % Extract the real and imaginary parts from the block matrix
    realPart_optimal = W(1:N_Comp, 1:N_Comp);
    imaginaryPart_optimal = W(N_Comp+1:end, 1:N_Comp);

    % Reconstruct the complex matrix
    W_optimal_complex = realPart_optimal + 1i * imaginaryPart_optimal;

    % % Beamforming matrix normalization
    % W_optimal_complex = W_optimal_complex/norm(W_optimal_complex,"fro");
    
    for snr_i = 1:length(EbN0dB_vector)
        EbN0dB = EbN0dB_vector(snr_i);
        EbN0 = 10^(EbN0dB / 10);
        N0 = Eb / EbN0;
    
        Nerrs_ml_B = 0;
        Nerrs_ml_E = 0;

%         Nerrs_ml_B_noBeam = 0;
%         Nerrs_ml_E_noBeam = 0;

        Nbits = 0;
    
        % Precompute constants
        noise_scale = sqrt(N0/2);
    
        while Nbits < Nsyms
            % Generate random symbols 
            
            random_pick = randi([1, 4]);
            syms_tr = s_Comp(:,random_pick);

            % % Generate random symbols
            % syms_tr = sign(rand(N_Comp, 1) - 0.5) + 1j * sign(rand(N_Comp, 1) - 0.5);
    
            % Alice - Bob channel
            % Generate the complex Gaussian noise 
            n_B = noise_scale*(randn(K_Comp,1) + 1i * randn(K_Comp,1));

            H_B_Comp_Sim = H_B_Comp/(sqrt(channel_var/2))*(sqrt(channel_var_Sim/2));

            y_B = H_B_Comp_Sim*W_optimal_complex*syms_tr + n_B;
        
%             y_B_noBeam = H_B_Comp*syms_tr + n_B;
    
            % Compute estimated symbols
            syms_tr_hat_ml_B = ML_Detection_QPSK(y_B, H_B_Comp_Sim, W_optimal_complex);

%             syms_tr_hat_ml_B_noBeam = ML_Detection_QPSK_noBeam(y_B_noBeam, H_B_Comp);
    
            % Calculate errors
            Nerrs_ml_B = Nerrs_ml_B + sum(real(syms_tr_hat_ml_B) ~= real(syms_tr))...
                       + sum(imag(syms_tr_hat_ml_B) ~= imag(syms_tr));

%             Nerrs_ml_B_noBeam = Nerrs_ml_B_noBeam + sum(real(syms_tr_hat_ml_B_noBeam) ~= real(syms_tr))...
%                               + sum(imag(syms_tr_hat_ml_B_noBeam) ~= imag(syms_tr));

            % Alice - Eve channel
            % Generate the complex Gaussian noise 
            n_E = noise_scale*(randn(K_Comp,1) + 1i * randn(K_Comp,1));
            
            H_E_Comp_Sim = H_E_Comp/(sqrt(channel_var/2))*(sqrt(channel_var_Sim/2));

            y_E = H_E_Comp_Sim*W_optimal_complex*syms_tr + n_E;

%             y_E_noBeam = H_E_Comp*syms_tr + n_E;
    
            % Compute estimated symbols
            syms_tr_hat_ml_E = ML_Detection_QPSK(y_E, H_E_Comp_Sim, W_optimal_complex);

%             syms_tr_hat_ml_E_noBeam = ML_Detection_QPSK_noBeam(y_E_noBeam, H_E_Comp);
    
            % Calculate errors
            Nerrs_ml_E = Nerrs_ml_E + sum(real(syms_tr_hat_ml_E) ~= real(syms_tr))...
                       + sum(imag(syms_tr_hat_ml_E) ~= imag(syms_tr));

%             Nerrs_ml_E_noBeam = Nerrs_ml_E_noBeam + sum(real(syms_tr_hat_ml_E_noBeam) ~= real(syms_tr))...
%                               + sum(imag(syms_tr_hat_ml_E_noBeam) ~= imag(syms_tr));

            Nbits = Nbits + 2 * N_Comp;
        end
    
        % Calculate Bit Error Rate
        BER_B(snr_i,iter) = Nerrs_ml_B / Nbits; 
        BER_E(snr_i,iter) = Nerrs_ml_E / Nbits; 

%         BER_B_noBeam(snr_i,iter) = Nerrs_ml_B_noBeam / Nbits; 
%         BER_E_noBeam(snr_i,iter) = Nerrs_ml_E_noBeam / Nbits; 
    end 
end

% Find mean of probability values
objective_function_Bob_average = mean(objective_function_Bob, 2);
objective_function_Eve_average = mean(objective_function_Eve, 2);


BER_average_B = mean(BER_B, 2);
BER_average_E = mean(BER_E, 2);

% BER_average_B_noBeam = mean(BER_B_noBeam, 2);
% BER_average_E_noBeam = mean(BER_E_noBeam, 2);

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
legend('Bob','Eve','Bob, No Beam','Eve, No Beam','Location','SouthEast');
saveas(gcf,'Figure1_2.fig');

figure(2);
semilogy(EbN0dB_vector,BER_average_B,'-ob','LineWidth',1.5);
hold on
semilogy(EbN0dB_vector,BER_average_E,'-*r','LineWidth',1.5);
% semilogy(EbN0dB_vector,BER_average_B_noBeam,'--b','LineWidth',1.5);
% semilogy(EbN0dB_vector,BER_average_E_noBeam,'--r','LineWidth',1.5);
hold off
grid on
xlabel('SNR');
ylabel('BER');
% legend('Bob, M-ary Beam','Eve, M-ary Beam','Bob, No Beam','Eve, No Beam','Location','SouthWest');
legend('Bob','Eve','Location','SouthWest');
saveas(gcf,'Figure2_2.fig');