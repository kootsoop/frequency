% <PRE>
%
%       1/5/94  Eric Jacobsen
%
% Interpolator comparison test simulator.
%
% This version compares the quadratic/parabolic and Quinn's
% interpolator across a bin with fractionally spaced samples.
%
% Variables NZ and NOISE must be intialized.
% NZ = 0/1 turns the noise generator off/on. NOISE, set
% to 1.0, 0.7071, and 0.5, provides -3, 0, and +3dB SNR.
%
% Last modified 1/28/99 - eaj  Cleaned up a little for readability.
% 

% $Id: Inttst2.m,v 1.1 1999/02/21 12:27:45 PeterK Exp PeterK $

tstlen=64;              % Effective length of test vector.

NZ=1; % Switch to enable (1) or disable (0) noise.
NOISE=0.7071; % Set noise level.

ACW=1; % Amplitude of tone.
N=1000; % Number of trials in test.

%hw=Hanning(tstlen); % Hanning window for Grandke's method.
 % Requires Signal Processing toolbox.

%bin=9.25; % Desired bin number of tone relative to long test length.

%fid = fopen('intdata.doc','w'); % Open file for results.

quinerr=zeros(size(1:N)); % Allocate vector for results.
quaderr=zeros(size(1:N)); % Allocate vector for results.
%jainerr=zeros(size(1:N)); % Allocate vector for results.
%granerr=zeros(size(1:N)); % Allocate vector for results.
quinest=zeros(size(1:N)); % Allocate vector for results.
quadest=zeros(size(1:N)); % Allocate vector for results.
%jainest=zeros(size(1:N)); % Allocate vector for results.
%granest=zeros(size(1:N)); % Allocate vector for results.
SNR=zeros(size(1:N)); % Allocate vector for results.

K=10; % Number of bins to test.
quinr=zeros(size(1:K)); % Allocate vector for results.
quadr=zeros(size(1:K)); % Allocate vector for results.
%jainr=zeros(size(1:K)); % Allocate vector for results.
%granr=zeros(size(1:K)); % Allocate vector for results.
quinvar=zeros(size(1:K)); % Allocate vector for results.
quadvar=zeros(size(1:K)); % Allocate vector for results.
%jainvar=zeros(size(1:K)); % Allocate vector for results.
%granvar=zeros(size(1:K)); % Allocate vector for results.
quinbias=zeros(size(1:K)); % Allocate vector for results.
quadbias=zeros(size(1:K)); % Allocate vector for results.
%jainbias=zeros(size(1:K)); % Allocate vector for results.
%granbias=zeros(size(1:K)); % Allocate vector for results.

targ=zeros(size(1:K)); % Allocate vector for results.

M=0; % Current test number.

for bin = 9.0: 0.1: 9.9,

M=M+1; % Current test number.

target=bin; % Calculate desired target result for comparison.
targ(M)=bin;

%fprintf(fid,'Target is %f.\n',target);
%fprintf(fid,'N = %f.\n',N);

fprintf('Peak is at %f.\n',bin);


%for C = 1.0: -0.25: 0.25, % Run these tests for NZ.

%NOISE=C;

for I = 1:N, % Run trials.

 phz=2*pi*rand(1);

 cw=ACW*exp(j*((2*pi*(1:tstlen)*bin/tstlen)+phz)); % Generate signal.
 sigp=sum(abs(cw(1:tstlen).^2)); % Calculate signal power.

 cwn=cw;

 if(NZ==1) % Generate channel noise.

  % Set signal level and add noise (noise variance=1 on I and on Q)                        
  nzi=NOISE*randn(1,tstlen);
  nzq=NOISE*randn(1,tstlen);
  nzv=nzi+j*nzq;

  nzp=sum(abs(nzv(1:tstlen).^2));
  
  if(nzp>0)
   SNR(I)=10*log10(sigp/nzp); % Calculate SNR=CNR.
  else
   SNR(I)=100.0; % Use this for no noise.
  end

  cwn=cw+nzv;  

 end % End channel noise.


 dftshrt(1:tstlen)=fft(cwn(1:tstlen));       % DFT.
 magshrt(1:tstlen)=abs(dftshrt);       % DFT magnitude.
 [rawmag,rawind]=max(magshrt); % Find raw peak magnitude and location.


 pk3vect(1:3)=dftshrt(rawind-1:rawind+1); % Isolated 3 samples around peak.
 
 quinest(I)=rawind-1+quin(pk3vect); % Do Quinn's estimation.
 quinerr(I)=target-quinest(I); % Calculate and save error.

 quadest(I)=rawind-1+real(quadterp(pk3vect)); % Do quadratic estimation.
 quaderr(I)=target-quadest(I); % Calculate and save error.

 % The following is Jain's method.

% if(pk3vect(1)>pk3vect(3)) % Find which adjacent bin is greatest.
%  alpha=abs(pk3vect(2))/abs(pk3vect(1));  
%  delta=alpha/(1+alpha);
%  jainest(I)=rawind-2+delta; % Jain's method.
% else
%  alpha=abs(pk3vect(3))/abs(pk3vect(2));  
%  jainest(I)=rawind-1+delta;
% end

% jainerr(I)=target-jainest(I); % Calculate and save error.


 % The following is Grandke's method.

% hcwn=cwn.*hw'; % Apply Hanning window to DFT input.
% gdftshrt(1:tstlen)=fft(hcwn(1:tstlen));       % Weighted DFT.
% gmagshrt(1:tstlen)=abs(gdftshrt);       % Weighted DFT magnitude.
% [grawmag,grawind]=max(gmagshrt); % Find raw peak magnitude and location.

% gpk3vect(1:3)=gdftshrt(grawind-1:grawind+1); % Isolated 3 samples around peak.

% if(gpk3vect(1)>gpk3vect(3)) % Find which adjacent bin is greatest.
%  alpha=abs(gpk3vect(2))/abs(gpk3vect(1));  
%  delta=(2*alpha-1)/(alpha+1);
%  grandkest(I)=grawind-2+delta; % Grandke's method.
% else
%  alpha=abs(gpk3vect(3))/abs(gpk3vect(2));  
%  delta=(2*alpha-1)/(alpha+1);
%  grandkest(I)=grawind-1+delta;
% end

% granerr(I)=target-grandkest(I); % Calculate and save error.

end % End inner for loop.

quinvar(M)=sqrt(mean(quinerr.^2)); % Calculate result rms error.
quadvar(M)=sqrt(mean(quaderr.^2)); % Calculate result rms error.
%jainvar(M)=sqrt(mean(jainerr.^2)); % Calculate result rms error.
%granvar(M)=sqrt(mean(granerr.^2)); % Calculate result rms error.
SNRmn=mean(SNR); % Calculate result mean.

quinbias(M)=mean(quinerr); % Calculate bias.
quadbias(M)=mean(quaderr); % Calculate bias.
%jainbias(M)=mean(jainerr); % Calculate bias.
%granbias(M)=mean(granerr); % Calculate bias.

quinr(M)=mean(quinest); % Calculate average result.
quadr(M)=mean(quadest); % Calculate average result.
%jainr(M)=mean(jainest); % Calculate average result.
%granr(M)=mean(grandkest); % Calculate average result.


%end; % End outer for loop.

end; % End bin loop.

%fclose(fid);

figure(1);
plot((1:M),targ(1:M),'r+',(1:M),quadr(1:M),'yo',(1:M),quinr(1:M),'g-');
title('Average Peak Location Estimates (bin) vs Trial Number');

figure(2);
plot((1:M),quadvar(1:M),'yx',(1:M),quinvar(1:M),'g-');
title('Estimator Variance vs Trial Number');

figure(3);
plot((1:M),quadbias(1:M),'yx',(1:M),quinbias(1:M),'g-');
title('Estimator Bias vs Trial Number');

end;
% </PRE>

