% <PRE>
function x=quadterp(y)
%
% Eric Jacobsen, EF Data Corp., 1994
%
% Performs a quadratic-fit peak location interpolation on a three-element input vector y.
% Returns -0.5 < x < 0.5, which is the fraction of the sample spacing about the center
% element where the peak is estimated to be.
%
% Obtained from Robert Bristow-Johnson, robert@audioheads.win.net
%

% $Id: Quadterp.m,v 1.1 1999/02/21 12:27:45 PeterK Exp PeterK $

x=(y(1)-y(3))/((2*y(2))-y(1)-y(3));
% </PRE>
