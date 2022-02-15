%
function omegahat = discperiod(signal)

%DISCPERIOD Find the discrete-frequency periodogram maximizer 
%           frequency estimates of the signals in each column 
%           of signal.
%
%          omegahat = discperiod(signal)
%
% $Id: discperiod.m,v 1.2 1999/01/09 11:10:35 PeterK Exp PeterK $

% File: discperiod.m
%
% Copyright (C) 1999 Peter J. Kootsookos
% 
% This software provided under the GNU General Public Licence, which
% is available from the website from which this file originated. For
% a copy, write to the Free Software Foundation, Inc., 675 Mass Ave, 
% Cambridge, MA 02139, USA.

%Type cast it to double.
signal = double(signal);

%--------------------------------------------------------
% Error condition checks

[T,N] = size(signal);

if any([T N]==0)
  error('discperiod: zero size data not allowed.');
end

if T==1
   signal = signal(:);
   [T,N] = size(signal);
end

SS = abs(fft(signal,T));

[mx,ix] = max(SS);

omegahat = 2*pi*(ix-1)/T;

% Author: Peter J. Kootsookos (kootsoop@ieee.org)
%
% Based on: P.J. Kootsookos, S.J. Searle and B.G. Quinn, 
% "Frequency Estimation Algorithms," CRC for Robust and 
% Adaptive Systems Internal Report, June 1993.
%
