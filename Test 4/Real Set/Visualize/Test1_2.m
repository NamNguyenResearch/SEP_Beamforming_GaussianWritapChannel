clear;
clc;

P = 1;
N_0 = 0.01;
a = 1;

% Setup 2 + 5
H_B_2 = [0.21 0.011; 0.09 0.3]; 
H_E_2 = [0.01 0.02; 0.017 0.01];

% Setup 3
H_B_3 = [0.21 0.011; 0.09 0.3]; 
H_E_3 = [-0.01 0.02; 0.01 0.01];

% Setup 4
H_B_4 = [0.21 0.015; 0.1 0.12];
H_E_4 = [0.01 0.071; 0.01 0.01];

% Setup 5
H_B_5 = [0.21 0.011; 0.09 0.3]; 
H_E_5 = [0.01 0.02; 0.017 0.01];


%==========================================================================
% Setup 2
D2 = 0.346;
constraint2 = (sqrt(N_0/(2*P))*qfuncinv(D2)/sqrt(a))^2;
range2 = [-3 3 -3 3];
w_optimal2 = [0.475;  0.88];
optimalValue2 = -norm(H_B_2*w_optimal2)^2;

figure(1)
hold on
Objective2 = ezplot(@(w_1,w_2) -(H_B_2(1,1)*w_1+H_B_2(1,2)*w_2).^2 - (H_B_2(2,1)*w_1+H_B_2(2,2)*w_2).^2 - optimalValue2,range2);
set(Objective2,'Color','#77AC30','LineWidth', 1);

Constraint2 = ezplot(@(w_1,w_2) (H_E_2(1,1)*w_1+H_E_2(1,2)*w_2)^2 + (H_E_2(2,1)*w_1+H_E_2(2,2)*w_2)^2 - constraint2,range2);
set(Constraint2,'Color','black','LineWidth',1);

PowerConstraint2 = ezplot(@(w_1,w_2) w_1^2 + w_2^2 - 1,range2);
set(PowerConstraint2,'Color','#A2142F','LineWidth',1);

plot(w_optimal2(1),w_optimal2(2),'ro','MarkerSize',6,'Linewidth',1);
hold off
grid on
axis equal
xlabel('w_1');
ylabel('w_2');
legend(['-||H_Bw||^2_2 =' num2str(round(optimalValue2,3))], ['||H_Ew||^2_2 = ' num2str(round(constraint2,3))],'||w||^2_2 = 1','Optimal');


%==========================================================================
% Setup 3
D3 = 0.2;
constraint3 = (sqrt(N_0/(2*P))*qfuncinv(D3)/sqrt(a))^2;
range3 = [-6 6 -6 6];
w_optimal3 = [0.4779;  0.8784];
optimalValue3 = -norm(H_B_3*w_optimal3)^2;

figure(2)
hold on
Objective3 = ezplot(@(w_1,w_2) -(H_B_3(1,1)*w_1+H_B_3(1,2)*w_2).^2 - (H_B_3(2,1)*w_1+H_B_3(2,2)*w_2).^2 - optimalValue3,range3);
set(Objective3,'Color','#77AC30','LineWidth',1);

Constraint3 = ezplot(@(w_1,w_2) (H_E_3(1,1)*w_1+H_E_3(1,2)*w_2)^2 + (H_E_3(2,1)*w_1+H_E_3(2,2)*w_2)^2 - constraint3,range3);
set(Constraint3,'Color','black','LineWidth',1);

PowerConstraint3 = ezplot(@(w_1,w_2) w_1^2 + w_2^2 - 1,range3);
set(PowerConstraint3,'Color','#A2142F','LineWidth',1);

plot(w_optimal3(1),w_optimal3(2),'ro','MarkerSize',5,'Linewidth',1);
hold off
grid on
axis equal
xlabel('w_1');
ylabel('w_2');
legend(['-||H_Bw||^2_2 =' num2str(round(optimalValue3,3))], ['||H_Ew||^2_2 = ' num2str(round(constraint3,3))],'||w||^2_2 = 1','Optimal');


%==========================================================================
% Setup 4
D4 = 0.3246;
constraint4 = (sqrt(N_0/(2*P))*qfuncinv(D4)/sqrt(a))^2;
range4 = [-4 4 -4 4];
w_optimal4 = [-0.9592;  -0.2828];
optimalValue4 = -norm(H_B_4*w_optimal4)^2;

figure(3)
hold on
Objective4 = ezplot(@(w_1,w_2) -(H_B_4(1,1)*w_1+H_B_4(1,2)*w_2).^2 - (H_B_4(2,1)*w_1+H_B_4(2,2)*w_2).^2 - optimalValue4,range4);
set(Objective4,'EdgeColor','#77AC30','LineWidth',1);

Constraint4 = ezplot(@(w_1,w_2) (H_E_4(1,1)*w_1+H_E_4(1,2)*w_2).^2 + (H_E_4(2,1)*w_1+H_E_4(2,2)*w_2).^2 - constraint4,range4);
set(Constraint4,'EdgeColor','black','LineWidth',1);

PowerConstraint4 = ezplot(@(w_1,w_2) w_1^2 + w_2^2 - 1,range4);
set(PowerConstraint4,'EdgeColor','#A2142F','LineWidth',1);

plot(w_optimal4(1),w_optimal4(2),'ro','MarkerSize',5,'Linewidth',1);
hold off
grid on
axis equal
xlabel('w_1');
ylabel('w_2');
legend(['-||H_Bw||^2_2 =' num2str(round(optimalValue4,3))], ['||H_Ew||^2_2 = ' num2str(round(constraint4,3))],'||w||^2_2 = 1','Optimal');

%==========================================================================
% Setup 5
D5 = 0.3555;
constraint5 = (sqrt(N_0/(2*P))*qfuncinv(D5)/sqrt(a))^2;
range5 = [-3 3 -3 3];
w_optimal5 = [-0.9106;  -0.4132];
optimalValue5 = -norm(H_B_5*w_optimal5)^2;

figure(4)
hold on
Objective5 = ezplot(@(w_1,w_2) -(H_B_5(1,1)*w_1+H_B_5(1,2)*w_2).^2 - (H_B_5(2,1)*w_1+H_B_5(2,2)*w_2).^2 - optimalValue5,range5);
set(Objective5,'Color','#77AC30','LineWidth',1);

Constraint5 = ezplot(@(w_1,w_2) (H_E_5(1,1)*w_1+H_E_5(1,2)*w_2)^2 + (H_E_5(2,1)*w_1+H_E_5(2,2)*w_2)^2 - constraint5,range5);
set(Constraint5,'Color','black','LineWidth',1);

PowerConstraint5 = ezplot(@(w_1,w_2) w_1^2 + w_2^2 - 1,range5);
set(PowerConstraint5,'Color','#A2142F','LineWidth',1);

plot(w_optimal5(1),w_optimal5(2),'ro','MarkerSize',5,'Linewidth',1);
hold off
grid on
axis equal
xlabel('w_1');
ylabel('w_2');
legend(['-||H_Bw||^2_2 =' num2str(round(optimalValue5,3))], ['||H_Ew||^2_2 = ' num2str(round(constraint5,3))],'||w||^2_2 = 1','Optimal');