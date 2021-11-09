%
function [x] = macleod(Y)

% Adapted by Eric Jacobsen, June 2002.
%
% Y is a three-element complex vector with the
% DFT output magnitude maximizer as the center element.
% 
% Returns -0.5 < x < 0.5, which is the fraction of the sample
% spacing (i.e., bin width) about the center element where the
% peak is estimated to be.
%
% Ref: Macleod, M.D., "Fast Nearly ML Estimation of the Parameters
% of Real or Complex Single Tones or Resolved Multiple Tones,"
% IEEE Trans. Sig. Proc. Vol 46 No 1,
% January 1998, pp141-148.

ref = Y(2);                   % Isolate phase reference.

R = real(Y.*conj(ref));       % Generate phase corrected coefficient vector.

gamma = (R(1)-R(3))/((2*R(2))+R(1)+R(3));       % Calculate offset.

delta = (sqrt(1 + 8*gamma*gamma)-1)/(4*gamma);  % Final estimate.

x = delta;
% 
