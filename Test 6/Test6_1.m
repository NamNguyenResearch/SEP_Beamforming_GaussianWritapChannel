clear;
clc;

% Parameters
N = 2;
K = 2;
P = 1;
N_B = 0.01;
N_E = 0.1;
channel_var = 1;
a = 1;

D = 0.000001:0.01:(0.5+0.000001);
constraint = (sqrt(N_B/2)*qfuncinv(D)/abs(a)).^2;
num_iter = 100;

H_B2 = [0   1i; 0 -1i];
H_E2 = [1i -1i; 0   0];

% Preallocate arrays
ErrorProbability_Bob1 = zeros(length(D), num_iter);
ErrorProbability_Eve1 = zeros(length(D), num_iter);
ErrorProbability_Bob_average1 = zeros(size(D));
ErrorProbability_Eve_average1 = zeros(size(D));
ErrorProbability_Bob2 = zeros(size(D));
ErrorProbability_Eve2 = zeros(size(D));

for iter = 1:num_iter
    % Generate the channel matrices and the modulated symbol vectors
    H_B1 = sqrt(channel_var/2)*(randn(K,N)+1i*randn(K,N));
    H_E1 = sqrt(channel_var/2)*(randn(K,N)+1i*randn(K,N));

    for i = 1:length(D)
        [ErrorProbability_Bob1(i,iter), ErrorProbability_Eve1(i,iter)] = calculate(N,P,N_B,N_E,H_B1,H_E1,a,constraint(i));
    end
end

for i = 1:length(D)
    ErrorProbability_Bob_average1(i) = mean(ErrorProbability_Bob1(i,:));
    ErrorProbability_Eve_average1(i) = mean(ErrorProbability_Eve1(i,:));

    [ErrorProbability_Bob2(i), ErrorProbability_Eve2(i)] = calculate(N,P,N_B,N_E,H_B2,H_E2,a,constraint(i));
end

% Plotting
figure(1)
subplot(2,1,1)
plot(D, ErrorProbability_Bob_average1, D, ErrorProbability_Eve_average1);
grid on
xlabel('D');
ylabel('Error Probability');
title('(a)');
legend('Bob','Eve','Location','SouthEast');
axis([0.000001, 0.5, 0, 0.5]);

subplot(2,1,2)
plot(D, ErrorProbability_Bob2, D, ErrorProbability_Eve2);
grid on
xlabel('D');
ylabel('Error Probability');
title('(b)');
legend('Bob','Eve','Location','SouthEast');
axis([0.000001, 0.5, 0, 0.5]);