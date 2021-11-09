%
function omegahat = wlp(signal,window, SNR)

%WLP  Find the Weighted Linear Predictor frequency estimates of 
%     the signals in each column of signal.
%
%          omegahat = wlp( signal [,window[, SNR]] )
%
%     where window is one of:
%       'ckq'    the Clarkson-Kootsookos-Quinn estimator [1]
%       'kay'    Kay's estimator [2] <------- (DEFAULT)
%       'lovell' the Lovell-Williamson estimator [3]
%       'lpr'    Lank-Reed-Pollon estimator [4]
%     and SNR is the ratio of the sinusoid amplitude squared
%     to the noise variance (i.e. NOT the usual SNR 
%     definition). The SNR is only required for the 'ckq' 
%     option.
%    
% [1]  V. Clarkson, P. J.  Kootsookos and B. G.  Quinn,
%      ``Variance Analysis of Kay's Weighted Linear Predictor 
%      Frequency Estimator,'' IEEE Transactions on Signal 
%      Processing, Vol. 42, pp. 2370-2379, 1994. 
%
% [2] S. M.  Kay, `` A Fast and Accurate Single Frequency 
%     Estimator,'' IEEE Transactions on Acoustics, Speech 
%     and Signal Processing, Vol. 37(12), pp. 1987-1989,
%     1989. 
%
% [3]  B.C. Lovell and R.C. Williamson,
%      "The Statistical Performance of Some Instantaneous
%      Frequency Estimators," IEEE Trans. on Acoustics, 
%      Speech and Signal Processing, Vol. 40, pp. 1708-1723,
%      1992.
%
% [4] G. W. Lank, I. S.  Reed and G. E.  Pollon, 
%     `` A Semicoherent Detection Statistic and Doppler 
%     Estimation Statistic,'' IEEE Transactions on 
%     Aerospace and Electronic Systems, Vol. AES-9(2), 
%     pp. 151-165, 1973. 
%
% $Id: wlp.m,v 1.4 1999/01/09 11:10:35 PeterK Exp PeterK $

% File: wlp.m
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
  error('wlp: zero size data not allowed.');
end;

if T==1
   signal = signal(:);
   [T,N] = size(signal);
end;

if ~exist('window')
  window = 'kay';
end;

% Need the analytic signal
if isreal(signal)
  signal = conj(hilbert(signal));
end;

T = T-1;
t = [1:T-1]';
switch (window)
  case 'ckq'
    if ~exist('SNR')
      SNR = 1;
    end;
    Theta = log(1 + SNR + sqrt(SNR^2 + SNR));
    WN = ( sinh(T*Theta) - sinh(t*Theta) - sinh((T-t)*Theta) ) ...
	/ ( (T-1)*sinh(T*Theta) - 2*sinh(T*Theta/2)*sinh((T-1)/2*Theta) / ...
	sinh(Theta/2));
  case 'kay'
    WN = (6*t.*(T-t));
  case 'lovell'
    WN = (6*t.*(T-t));
    signal = signal./abs(signal);
  case 'lrp'
    WN = ones(T-1,1);
end;

WN = WN*ones(1,N);

omegahat = angle(sum([ signal(1:T-1,:).*conj(signal(3:T+1,:)) ].*WN)) / 2;

% Author: Peter J. Kootsookos (kootsoop@ieee.org)
%
% Based on: P.J. Kootsookos, S.J. Searle and B.G. Quinn, 
% "Frequency Estimation Algorithms," CRC for Robust and 
% Adaptive Systems Internal Report, June 1993.

%
