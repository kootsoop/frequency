
%
function x=quadterp(y)
%
%   Eric Jacobsen 1994
%
% Performs a quadratic-fit peak location interpolation on a three-element input vector y.
% Returns -0.5 < x < 0.5, which is the fraction of the sample spacing about the center
% element where the peak is estimated to be.
%
%   Modified:
%
%   June 2002 - eaj - Integrated real() operation for consistency with other estimators.
%

x=real((y(1)-y(3))/((2*y(2))-y(1)-y(3)));

%
