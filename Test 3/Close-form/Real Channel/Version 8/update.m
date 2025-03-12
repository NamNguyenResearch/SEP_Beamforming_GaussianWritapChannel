function [W_update] = update(W_value,P)
    global alpha % Step size
    
    [grad_W] = gradient(W_value);  

    W_value = W_value - alpha * grad_W;
    [W_update] = projection(W_value,P);
end