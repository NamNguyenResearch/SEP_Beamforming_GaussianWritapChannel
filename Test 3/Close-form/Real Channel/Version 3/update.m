function [W_update] = update(W,lambda)
    global alpha % Step size
    
    [grad_W] = gradient(W,lambda);  

    W = W - alpha * grad_W;
    [W_update] = projection(W);
end