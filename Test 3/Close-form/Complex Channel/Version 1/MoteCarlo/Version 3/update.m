function [W_update] = update(W,lambda,M,N,P,H_B,H_E,s,N_B,N_E)

    [grad_W] = gradient(W,lambda,M,H_B,H_E,s,N_B,N_E);
    
    % % Backtracking line search
    % alpha_linesearch = 0.4; % (0,0.5)
    % beta_linesearch = 0.8; % (0,1)
    % 
    % alpha = 1;
    % 
    % objective_function_Bob1 = objectFunction(M,H_B,s,N_B,W); % Error Probability of Bob
    % 
    % objective_function_Eve1 = objectFunction(M,H_E,s,N_E,W); % Error Probability of Eve
    % 
    % objective_function1 = objective_function_Bob1 - lambda * objective_function_Eve1; % Objective Value
    % 
    % % objective_function1_tmp = objective_function1 + alpha_linesearch*alpha*(grad_W')*(-grad_W);
    % 
    % objective_function1_tmp = objective_function1 + alpha_linesearch*alpha*(-norm(grad_W,"fro"));
    % 
    % tmp = W - alpha * grad_W;
    % 
    % objective_function_Bob2 = objectFunction(M,H_B,s,N_B,tmp); % Error Probability of Bob
    % 
    % objective_function_Eve2 = objectFunction(M,H_E,s,N_E,tmp); % Error Probability of Eve
    % 
    % objective_function2 = objective_function_Bob2 - lambda * objective_function_Eve2; % Objective Value
    % 
    % 
    % while objective_function2 > objective_function1_tmp
    % 
    %     alpha = beta_linesearch*alpha;
    % 
    %     % objective_function1_tmp = objective_function1 + alpha_linesearch*alpha*(grad_W')*(-grad_W);
    % 
    %     objective_function1_tmp = objective_function1 + alpha_linesearch*alpha*(-norm(grad_W,"fro"));
    % 
    %     tmp = W - alpha * grad_W;
    % 
    %     objective_function_Bob2 = objectFunction(M,H_B,s,N_B,tmp); % Error Probability of Bob
    % 
    %     objective_function_Eve2 = objectFunction(M,H_E,s,N_E,tmp); % Error Probability of Eve
    % 
    %     objective_function2 = objective_function_Bob2 - lambda * objective_function_Eve2; % Objective Value
    % end

    alpha = 0.005;

    W = W - alpha * grad_W;
    [W_update] = projection(W,N,P);
end