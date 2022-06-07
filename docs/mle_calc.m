function [L, CS] = mle_calc(omega, x)

    x = x(:);
    T = length(x);

    nu = (T-1)/2;
    t = [0:T-1] - nu;
    c = cos(omega*t)';
    s = sin(omega*t)';
    
    CS = [c'*c c'*s; s'*c s'*s];
    
    L = [c'*x s'*x]*pinv(CS)*[c'*x; s'*x];

return


