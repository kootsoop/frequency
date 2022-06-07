LLs = [];
FFs = [];
for phase_1 = [0:.5:2*pi]
    %phase_1 = 2;
    disp(phase_1);
    %for phase_2 = [0:.5:2*pi]
    phase_2 = 0;

        x = cos(2*pi*0.123456789*[0:9] + phase_1) + cos(2*pi*(0.123456789+1/15)*[0:9] + phase_2);
        LLs = [LLs mle_fourier(x)]; 
        FFs = [FFs; abs(fft(x, 10001))];
    end
%end

