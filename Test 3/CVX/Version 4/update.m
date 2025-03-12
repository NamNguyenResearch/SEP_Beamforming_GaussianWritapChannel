function [W_update] = update(W)
    global alpha % Step size
    
    [grad_W] = gradient(W);  

    W = W - alpha * grad_W;
    [W_update] = projection(W);
end