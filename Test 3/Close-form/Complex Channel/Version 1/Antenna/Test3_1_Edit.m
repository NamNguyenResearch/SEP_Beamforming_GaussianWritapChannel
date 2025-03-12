clear;
clc;

N_Comp1 = 2; % Number of the transmitted antennas
K_Comp1 = 2; % Number of the received antennas => Not depend on for this case! 
L_Comp1 = 2;
N_Comp2 = 5; % Number of the transmitted antennas
K_Comp2 = 5; % Number of the received antennas => Not depend on for this case!
L_Comp2 = 2;

max_iteration = 45;
    
[objective_function_Bob_average1, objective_function_Eve_average1] = optimization_case(N_Comp1,K_Comp1,L_Comp1);
[objective_function_Bob_average2, objective_function_Eve_average2] = optimization_case(N_Comp2,K_Comp2,L_Comp2);
    
%==========================================================================
%Figure of the objective value versus iteration 
figure(1);
iterations = 1:max_iteration;
plot(iterations,objective_function_Bob_average1,'b-','LineWidth',1.5);
hold on 
plot(iterations,objective_function_Eve_average1,'r--','LineWidth',1.5);
plot(iterations,objective_function_Bob_average2,'b-','LineWidth',1.5);
plot(iterations,objective_function_Eve_average2,'r--','LineWidth',1.5);
hold off
grid on
xlabel('Iterations');
ylabel('Bound of symbol error probability');
legend('Bob','Eve','Location','NorthWest');