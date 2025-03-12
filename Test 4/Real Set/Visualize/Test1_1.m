clear;
clc;

P = 1;
N_0 = 0.01;
a = 1;

D = 0.346;

constraint = (sqrt(N_0/(2*P))*qfuncinv(D)/sqrt(a))^2;

% % % Setup 2 + 3 + 5
% H_B = [0.21 0.011; 0.09 0.3]; 
% H_E = [0.01 0.02; 0.017 0.01];

% H_B = [0.21 0.011; 0.09 0.3]; 
% H_E = [-0.21 0.011; -0.09 0.3]; 

H_B = sqrt(0.01)*randn(2,2)
H_E = [-H_B(1,1) H_B(1,2); -H_B(2,1) H_B(2,2)]


%==========================================================================
objec_func1 = @(w_1,w_2) -(H_B(1,1)*w_1+H_B(1,2)*w_2).^2 - (H_B(2,1)*w_1+H_B(2,2)*w_2).^2;
[X1,Y1] = meshgrid(-1:0.01:1);
Z1 = objec_func1(X1,Y1); 

objec_func2 = @(w_1,w_2) w_1.^2 + w_2.^2;
[X2,Y2] = meshgrid(-1:0.01:1);
Z2 = objec_func2(X2,Y2); 

objec_func3 = @(w_1,w_2) (H_E(1,1)*w_1+H_E(1,2)*w_2).^2 + (H_E(2,1)*w_1+H_E(2,2)*w_2).^2;
[X3,Y3] = meshgrid(-1:0.01:1);
Z3 = objec_func3(X3,Y3); 

%==========================================================================
figure(1)

% contour(X1,Y1,Z1,'b','ShowText','on');
contour(X1,Y1,Z1,'b');

hold on 

% contour(X2,Y2,Z2,1,'r','ShowText','on');
contour(X2,Y2,Z2,1,'r');

% contour(X3,Y3,Z3,'k','ShowText','on');
contour(X3,Y3,Z3,'k');

hold off
grid on
xlabel('w_1');
ylabel('w_2');
% legend('-||H_Bw||^2','||w||^2 = 1','||H_Ew||^2');