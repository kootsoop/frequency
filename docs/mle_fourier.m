function Ls = mle_fourier(x, T)

    if (~exist('T','var'))
        T = length(x);
    end
    Ls = zeros(T,1);
    indx = 1;
    for omega = [0:(T-1)]/T*2*pi 
        Ls(indx) = mle_calc(omega, x); 
        indx = indx + 1;
    end
return
