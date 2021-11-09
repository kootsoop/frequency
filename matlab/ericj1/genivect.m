% <PRE>
function M=genivect(k,NS,NL)
%
% M = genivect( k, NS, NL )     12/1/94 Eric Jacobsen
%
% Generates a length NS interpolating vector M for NS:NL interpolation.
% k is the index of the output in the interpolated length NL vector.
% See the MathCad file intrpltr.mcd for the derivation.
%

% 
% For MATLAB, sign of exp_j should be (-), for MATHCAD, (+).

% $Id: Genivect.m,v 1.1 1999/02/21 12:27:45 PeterK Exp PeterK $

M(1:NS)=zeros(size(1:NS));
for n = 0:(NS-1)
 for i = 0:(NS-1)
  M(n+1)=M(n+1)+exp(-j*2*pi*i*((k/NL)-(n/NS))); % Sign of exponent changes between MATLAB/MathCAD.
  end
 M(n+1)=M(n+1)/NS;
 end

% </PRE>

