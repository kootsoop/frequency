%
function [x] = quin(pk3vect)

% Adapted by Eric Jacobsen, 1994.
%
% pk3vect is a three-element complex vector with the
% DFT output magnitude maximizer as the center element.
% 
% Returns -0.5 < x < 0.5, which is the fraction of the sample
% spacing (i.e., bin width) about the center element where the
% peak is estimated to be.
%
% Ref: Quinn, BG, "Estimating frequency by interpolation using
% Fourier coefficients," IEEE Trans. Sig. Proc. Vol 42 No 5,
% May 1994, pp1264-1268.

% $Id: Quinn.m,v 1.1 1999/02/21 12:27:45 PeterK Exp PeterK $

alpha1=real(pk3vect(1)/pk3vect(2));
alpha2=real(pk3vect(3)/pk3vect(2));

delta1= alpha1/(1-alpha1);
delta2=-alpha2/(1-alpha2);

if ((delta1>0) & (delta2>0))
	x=delta2;
else
        x=delta1;

end
% 
