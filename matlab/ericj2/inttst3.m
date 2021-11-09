%
%
%       1/5/94  Eric Jacobsen
%
%	Interpolator comparison test simulator.
%
%	This version compares the quadratic/parabolic and Quinn's
%	interpolator across a bin with fractionally spaced samples.
%
%	Variables NZ and NOISE must be intialized.
%	NZ = 0/1 turns the noise generator off/on. NOISE, set
%	to 1.0, 0.7071, and 0.5, provides -3, 0, and +3dB SNR.
%
%	
%   Last modified:
%
%   6/14/02 - eaj  Added tests for Macleod's and Quinn's 2nd method.
%   1/28/99 - eaj  Cleaned up a little for readability.
%



tstlen=64;              % Effective length of test vector.

NZ=1;			% Switch to enable (1) or disable (0) noise.
NOISE=0.5;		% Set noise level.

ACW=1;			% Amplitude of tone.
N=10000;			% Number of trials in test.

%hw=Hanning(tstlen);	% Hanning window for Grandke's method.
			% Requires Signal Processing toolbox.

%bin=9.25;		% Desired bin number of tone relative to long test length.

%fid = fopen('intdata.doc','w');		% Open file for results.

quinerr=zeros(size(1:N));	% Allocate vector for results.
quaderr=zeros(size(1:N));	% Allocate vector for results.
quin2err=zeros(size(1:N));	% Allocate vector for results.
maclderr=zeros(size(1:N));	% Allocate vector for results.
%jainerr=zeros(size(1:N));	% Allocate vector for results.
%granerr=zeros(size(1:N));	% Allocate vector for results.
quinest=zeros(size(1:N));	% Allocate vector for results.
quadest=zeros(size(1:N));	% Allocate vector for results.
quin2est=zeros(size(1:N));	% Allocate vector for results.
macldest=zeros(size(1:N));	% Allocate vector for results.
%jainest=zeros(size(1:N));	% Allocate vector for results.
%granest=zeros(size(1:N));	% Allocate vector for results.
SNR=zeros(size(1:N));		% Allocate vector for results.

K=10;				% Number of bins to test.
quinr=zeros(size(1:K));		% Allocate vector for results.
quadr=zeros(size(1:K));		% Allocate vector for results.
quin2r=zeros(size(1:K));		% Allocate vector for results.
macldr=zeros(size(1:K));		% Allocate vector for results.
%jainr=zeros(size(1:K));		% Allocate vector for results.
%granr=zeros(size(1:K));		% Allocate vector for results.
quinvar=zeros(size(1:K));	% Allocate vector for results.
quadvar=zeros(size(1:K));	% Allocate vector for results.
quin2var=zeros(size(1:K));	% Allocate vector for results.
macldvar=zeros(size(1:K));	% Allocate vector for results.
%jainvar=zeros(size(1:K));	% Allocate vector for results.
%granvar=zeros(size(1:K));	% Allocate vector for results.
quinbias=zeros(size(1:K));	% Allocate vector for results.
quadbias=zeros(size(1:K));	% Allocate vector for results.
quin2bias=zeros(size(1:K));	% Allocate vector for results.
macldbias=zeros(size(1:K));	% Allocate vector for results.
%jainbias=zeros(size(1:K));	% Allocate vector for results.
%granbias=zeros(size(1:K));	% Allocate vector for results.

targ=zeros(size(1:K));		% Allocate vector for results.

M=0;				% Current test number.

binstrt = 9.0;          % Starting bin number.
binstep = 0.1;          % Bin delta step size.
binend = 9.9;           % Ending bin number.

for bin = binstrt: binstep: binend,

M=M+1;				% Current test number.

target=bin;			% Calculate desired target result for comparison.
targ(M)=bin;

%fprintf(fid,'Target is %f.\n',target);
%fprintf(fid,'N = %f.\n',N);

fprintf('Peak is at %f.\n',bin);


%for C = 1.0: -0.25: 0.25,	% Run these tests for NZ.

%NOISE=C;

for I = 1:N,		% Run trials.

 phz=2*pi*rand(1);

 cw=ACW*exp(j*((2*pi*(1:tstlen)*bin/tstlen)+phz));		% Generate signal.
 sigp=sum(abs(cw(1:tstlen).^2));			% Calculate signal power.

 cwn=cw;

 if(NZ==1)			% Generate channel noise.

  % Set signal level and add noise (noise variance=1 on I and on Q)                        
  nzi=NOISE*randn(1,tstlen);
  nzq=NOISE*randn(1,tstlen);
  nzv=nzi+j*nzq;

  nzp=sum(abs(nzv(1:tstlen).^2));
  
  if(nzp>0)
   SNR(I)=10*log10(sigp/nzp);	% Calculate SNR=CNR.
  else
   SNR(I)=100.0;		% Use this for no noise.
  end

  cwn=cw+nzv;  

 end				% End channel noise.


 dftshrt(1:tstlen)=fft(cwn(1:tstlen));      	% DFT.
 magshrt(1:tstlen)=abs(dftshrt);      		% DFT magnitude.
 [rawmag,rawind]=max(magshrt);			% Find raw peak magnitude and location.


 pk3vect(1:3)=dftshrt(rawind-1:rawind+1);	% Isolated 3 samples around peak.
 
 quinest(I)=rawind-1+quin(pk3vect);		% Do Quinn's first estimation.
 quinerr(I)=target-quinest(I);			% Calculate and save error.

 quadest(I)=rawind-1+quadterp(pk3vect);	% Do modified quadratic estimation.
 quaderr(I)=target-quadest(I);			% Calculate and save error.

 quin2est(I)=rawind-1+quin2(pk3vect);	% Do Quinn's second estimation.
 quin2err(I)=target-quin2est(I);		% Calculate and save error.

 macldest(I)=rawind-1+macleod(pk3vect);	% Do Macleod's estimation.
 maclderr(I)=target-macldest(I);		% Calculate and save error.
 
 
 % The following is Jain's method.

% if(pk3vect(1)>pk3vect(3))			% Find which adjacent bin is greatest.
%  alpha=abs(pk3vect(2))/abs(pk3vect(1));  
%  delta=alpha/(1+alpha);
%  jainest(I)=rawind-2+delta;			% Jain's method.
% else
%  alpha=abs(pk3vect(3))/abs(pk3vect(2));  
%  jainest(I)=rawind-1+delta;
% end

% jainerr(I)=target-jainest(I);			% Calculate and save error.


 % The following is Grandke's method.

% hcwn=cwn.*hw';					% Apply Hanning window to DFT input.
% gdftshrt(1:tstlen)=fft(hcwn(1:tstlen));      	% Weighted DFT.
% gmagshrt(1:tstlen)=abs(gdftshrt);      	% Weighted DFT magnitude.
% [grawmag,grawind]=max(gmagshrt);		% Find raw peak magnitude and location.

% gpk3vect(1:3)=gdftshrt(grawind-1:grawind+1);	% Isolated 3 samples around peak.

% if(gpk3vect(1)>gpk3vect(3))			% Find which adjacent bin is greatest.
%  alpha=abs(gpk3vect(2))/abs(gpk3vect(1));  
%  delta=(2*alpha-1)/(alpha+1);
%  grandkest(I)=grawind-2+delta;			% Grandke's method.
% else
%  alpha=abs(gpk3vect(3))/abs(gpk3vect(2));  
%  delta=(2*alpha-1)/(alpha+1);
%  grandkest(I)=grawind-1+delta;
% end

% granerr(I)=target-grandkest(I);		% Calculate and save error.

end						% End inner for loop.

quinvar(M)=sqrt(mean(quinerr.^2));	% Calculate result rms error.
quadvar(M)=sqrt(mean(quaderr.^2));	% Calculate result rms error.
quin2var(M)=sqrt(mean(quin2err.^2));	% Calculate result rms error.
macldvar(M)=sqrt(mean(maclderr.^2));	% Calculate result rms error.
%jainvar(M)=sqrt(mean(jainerr.^2));	% Calculate result rms error.
%granvar(M)=sqrt(mean(granerr.^2));	% Calculate result rms error.
SNRmn=mean(SNR);		% Calculate result mean.

quinbias(M)=mean(quinerr);		% Calculate bias.
quadbias(M)=mean(quaderr);		% Calculate bias.
quin2bias(M)=mean(quin2err);		% Calculate bias.
macldbias(M)=mean(maclderr);		% Calculate bias.
%jainbias(M)=mean(jainerr);		% Calculate bias.
%granbias(M)=mean(granerr);		% Calculate bias.

quinr(M)=mean(quinest);		% Calculate average result.
quadr(M)=mean(quadest);		% Calculate average result.
quin2r(M)=mean(quin2est);		% Calculate average result.
macldr(M)=mean(macldest);		% Calculate average result.
%jainr(M)=mean(jainest);		% Calculate average result.
%granr(M)=mean(grandkest);	% Calculate average result.


%end;						% End outer for loop.

end;						% End bin loop.

%fclose(fid);

xs(1:M)=binstrt+(((1:M)-1)*binstep);    % Generate scale for plot horizontal axis.

figure(1);
plot(xs(1:M),targ(1:M),'r-',xs(1:M),quadr(1:M),'mo',xs(1:M),quinr(1:M),'gs',xs(1:M),quin2r(1:M),'k+',xs(1:M),macldr(1:M),'bx');
title('Average Peak Location Estimates (bin) vs Bin Number');

figure(2);
plot(xs(1:M),quadvar(1:M),'mo:',xs(1:M),quinvar(1:M),'gs:',xs(1:M),quin2var(1:M),'k+:',xs(1:M),macldvar(1:M),'bx:');
title('Estimator Variance vs Bin Number');

figure(3);
plot(xs(1:M),quadbias(1:M),'mo:',xs(1:M),quinbias(1:M),'gs:',xs(1:M),quin2bias(1:M),'k+:',xs(1:M),macldbias(1:M),'bx:');
title('Estimator Bias vs Bin Number');

%end;
% 
