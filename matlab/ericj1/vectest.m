%
%
%       8/15/94  Eric Jacobsen
%
%	Interpolator comparison test simulator.
%	
%	This one takes statistics on the EF Data vector interpolator.
%

% $Id: Vectst.m,v 1.1 1999/02/21 12:27:45 PeterK Exp PeterK $


NZ=1;			% Noise switch, 1/0 = on/off.

M=genivect(3,4,8);	% Generate EF Data interpolation vector.

tstlen=128;             % Effective length of test vector.
shtlen=tstlen/2;	% Length of short DFT to be interpolated.

ACW=1;			% Amplitude of tone.
N=1000;			% Number of trials in test.


%bin=9.25;		% Desired bin number of tone relative to long test length.

fid = fopen('vecdata.doc','w');		% Open file for results.

LDFTerr=zeros(size(1:N));	% Allocate vector for results.
SDFTerr=zeros(size(1:N));	% Allocate vector for results.
ZPerr=zeros(size(1:N));		% Allocate vector for results.
EFerr=zeros(size(1:N));		% Allocate vector for results.

LDFTmag=zeros(size(1:N));	% Allocate vector for results.
SDFTmag=zeros(size(1:N));	% Allocate vector for results.
ZPmag=zeros(size(1:N));		% Allocate vector for results.
EFmag=zeros(size(1:N));		% Allocate vector for results.

SNR=zeros(size(1:N));		% Allocate vector for results.


for bin = 9.0: 0.125: 9.625,

target=bin;		% Calculate desired target result for comparison.

fprintf(fid,'Target is %f.\n',target);
fprintf(fid,'N = %f.\n',N);

fprintf('\nPeak is at %f.\n',bin);


for C = 1.125: -0.125: 0.25,	% Run these tests for NZ.

NOISE=C;

for I = 1:N,		% Run trials.

 phz=2*pi*rand(1);					% Random phase.

 cw=ACW*exp(j*((4*pi*(1:tstlen)*bin/tstlen)+phz));	% Generate signal.
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


 
 dftlong(1:tstlen)=fft(cwn(1:tstlen));      	% Long (reference) DFT.
 maglong(1:tstlen)=abs(dftlong);	      	% Long DFT magnitude.

 [longmag,longind]=max(maglong);
 maxlong=dftlong(longind);

 LDFTerr(I)=target-((longind-1)/2);
 LDFTmag(I)=longmag;

 %fprintf('Long Reference:  Mag = %f, Phase = %f radians at %f.\n', longmag, angle(maxlong), (longind-1)/2);

 dftshrt(1:shtlen)=fft(cwn(1:shtlen));      	% Short DFT.
 magshrt(1:shtlen)=abs(dftshrt);      		% Short DFT magnitude.
 [rawmag,rawind]=max(magshrt);			% Find raw peak magnitude and location.
 
 if(rawind>62)					% Correct for bogus case.
  magshrt(rawind)=0;
  [rawmag,rawind]=max(magshrt);			% Find raw peak magnitude and location.
  end
 if(rawind<3)					% Correct for bogus case.
  magshrt(rawind)=0;
  [rawmag,rawind]=max(magshrt);			% Find raw peak magnitude and location.
  end

 maxshrt=dftshrt(rawind);

 SDFTerr(I)=target-(rawind-1);
 SDFTmag(I)=rawmag;

 %fprintf('Raw peak is %f at bin number %d.\n', rawmag, rawind-1);
 %fprintf('Short Reference: Mag = %f, Phase = %f radians at %f.\n', rawmag, angle(maxshrt), rawind-1);
  
 invshrt=ifft(dftshrt);				% Inverse transform short DFT.
 zpdftshrt=zeros(size(dftlong));		% Create zero padding vector for dftshort.
 zpdftshrt(1:shtlen)=invshrt(1:shtlen);

 zidftshrt=fft(zpdftshrt);			% Forward transform zero padded vector.
 magzi(1:tstlen)=abs(zidftshrt);      		% Zero padded DFT magnitude.
 [zpmag,zpind]=max(magzi);

 maxzi=zidftshrt(zpind);

 ZPerr(I)=target-((zpind-1)/2);
 ZPmag(I)=zpmag;
 
 %fprintf('Zero padded: Mag = %f, Phase = %f radians.\n', abs(maxzi), angle(maxzi));
 %fprintf('Zero Padded: Mag = %f, Phase = %f radians at %f.\n', zpmag, angle(maxzi), (zpind-1)/2);

 if(magshrt(rawind+1)>magshrt(rawind-1))	% Which side of raw peak to interpolate.
  pkvect(1:4)=dftshrt(rawind-1:rawind+2);	% Isolated 4 samples around peak.
%  fprintf('Interpolating high side of bin %d.\n', rawind-1);
  efint=sum(pkvect.*M);				% Calculate our interpolation method. 
  efimag=1.2*abs(efint);
  if(efimag>rawmag)				% Make freq estimate.
   efind=rawind+0.5;
  else
   efind=rawind;
   efimag=rawmag;
  end
 else
  pkvect(1:4)=dftshrt(rawind-2:rawind+1);	% Isolated 4 samples around peak.
%  fprintf('Interpolating low side of bin %d.\n', rawind-1);
  efint=sum(pkvect.*M);				% Calculate our interpolation method. 
  efimag=1.2*abs(efint);
  if(efimag>rawmag)				% Make freq estimate.
   efind=rawind-0.5;
  else
   efind=rawind;
   efimag=rawmag;
  end
 end

 EFerr(I)=target-(efind-1);
 EFmag(I)=efimag; 

% fprintf('EF Interpolated: Mag = %f, Phase = %f radians.\n', abs(efint), angle(efint));




end						% End inner for loop.


LDFTmn=sqrt(mean(LDFTerr.^2));	% Calculate result rms error.
SDFTmn=sqrt(mean(SDFTerr.^2));	% Calculate result rms error.
ZPmn=sqrt(mean(ZPerr.^2));	% Calculate result rms error.
EFmn=sqrt(mean(EFerr.^2));	% Calculate result rms error.
LDMmn=mean(LDFTmag);		% Calculate result mean magnitude.
SDMmn=mean(SDFTmag);		% Calculate result mean magnitude.
ZPMmn=mean(ZPmag);		% Calculate result mean magnitude.
EFMmn=mean(EFmag);		% Calculate result mean magnitude.
SNRmn=mean(SNR);		% Calculate test mean SNR.

fprintf('\nTest case noise variance = %d.\n\n',C);
fprintf('Bin number rms error from target.\n');
fprintf('LDFT rms=%f\n',LDFTmn);
fprintf('SDFT rms=%f\n',SDFTmn);
fprintf('ZP   rms=%f\n',ZPmn);
fprintf('EF   rms=%f\n',EFmn);
fprintf('Result mean peak magnitude estimates.\n');
%fprintf('LDM rms=%f\n',LDMmn);
fprintf('SDM =%f\n',SDMmn);
fprintf('ZPM =%f\n',ZPMmn);
fprintf('EFM =%f\n',EFMmn);
fprintf('SNR     X=%f\n',SNRmn);

fprintf(fid,'\nTest case %d.\n',C);
fprintf(fid,'Bin number rms error from target.\n');
fprintf(fid,'LDFT rms=%f\n',LDFTmn);
fprintf(fid,'SDFT rms=%f\n',SDFTmn);
fprintf(fid,'ZP   rms=%f\n',ZPmn);
fprintf(fid,'EF   rms=%f\n',EFmn);
fprintf(fid,'Result mean peak magnitude estimates.\n');
%fprintf(fid,'LDM rms=%f\n',LDMmn);
fprintf(fid,'SDM rms=%f\n',SDMmn);
fprintf(fid,'ZPM rms=%f\n',ZPMmn);
fprintf(fid,'EFM rms=%f\n',EFMmn);
fprintf(fid,'SNR     X=%f\n',SNRmn);

end;						% End outer for loop.

end;						% End bin loop.

fclose(fid);

end;
% 
