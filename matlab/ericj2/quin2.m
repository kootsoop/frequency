%
function [x] = quin2(Y)

% Adapted by Eric Jacobsen, June, 2002.
%
% Y is a three-element complex vector with the
% DFT output magnitude maximizer as the center element.
% 
% Returns -0.5 < x < 0.5, which is the fraction of the sample
% spacing (i.e., bin width) about the center element where the
% peak is estimated to be.
%
% Ref: Quinn, BG, "Estimation of frequency, amplitude and phase
% from the DFT of a time series," IEEE Trans. Sig. Proc. Vol 45, No 3,
% Mar 1997, pp814-817.


betam=real(Y(1)/Y(2));
betap=real(Y(3)/Y(2));

dm=-betam/(betam-1);
dp= betap/(betap-1);

kappap=(1/4)*log(3*(dp^4)+6*(dp^2)+1) - (sqrt(6)/24)*log(((dp^2)+1-sqrt(2/3))/((dp^2)+1+sqrt(2/3)));
kappam=(1/4)*log(3*(dm^4)+6*(dm^2)+1) - (sqrt(6)/24)*log(((dm^2)+1-sqrt(2/3))/((dm^2)+1+sqrt(2/3)));

delta=(dm+dp)/2 +kappap-kappam;

x=delta;

% 
