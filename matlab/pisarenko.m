%
function omegahat = pisarenko(signal)

%PISARENKO Find the Pisarenko frequency estimates of the signals 
%          in each column of signal.
%
%          omegahat = pisarenko(signal)
%
% [1] V. F.  Pisarenko, "On the Estimation of Spectra by 
%     Means of Non-linear Functions of the Covariance Matrix,"
%     Geophys. J. Roy. astr. Soc., Vol. 28, pp. 511-531, 1972. 
%
% [2] V. F. Pisarenko, ``The Retrieval of Harmonics from a 
%     Covariance Function," Geophys. J. Roy. astr. Soc., 
%     Vol.   33, pp. 347-366, 1973. 
%
% $Id: pisarenko.m,v 1.2 1999/01/09 11:10:35 PeterK Exp PeterK $

% File: pisarenko.m
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
  error('pisarenko: zero size data not allowed.');
end;

if T==1
   signal = signal(:);
   [T,N] = size(signal);
end;

Rss1 = sum(signal(2:T,:).*signal(1:T-1,:));
Rss2 = sum(signal(3:T,:).*signal(1:T-2,:));

alpha = ( Rss2 + sqrt(Rss2.^2 + 8*Rss1.^2) ) ./ ( Rss1 + eps ) / 2;

omegahat = acos(alpha/2);

% Author: Peter J. Kootsookos (kootsoop@ieee.org)
%
% Based on: P.J. Kootsookos, S.J. Searle and B.G. Quinn, 
% "Frequency Estimation Algorithms," CRC for Robust and 
% Adaptive Systems Internal Report, June 1993.
%
