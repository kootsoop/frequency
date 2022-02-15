%
% Simulate over different frequencies, different initial phases, and different noise levels.
%
noise_realizations = 100;
signal_length = 100;
SNRdb = [-20:.1:10];
frequencies = 2*pi*0.25 + 2*pi/signal_length/2; % + linspace(-2*pi/signal_length, 2*pi/signal_length, 100);
phase = pi/2;


signal_amplitude = 1;
signal_power = signal_amplitude*signal_amplitude/2;
noise = randn(signal_length, noise_realizations);

results = [];

for snr = SNRdb
    snr
    noise_variance = signal_power / 10^(snr/10); % calculated noise power
    scaled_noise = noise_variance*noise;
    for frequency = frequencies
        frequency
        noisy_signal = cos(frequency*kron([0:signal_length-1]',ones(1,noise_realizations)) + phase) + scaled_noise;
        discperiod_err = sum(abs(discperiod(noisy_signal)-frequency).^2);
        pisarenko_err = sum(abs(pisarenko(noisy_signal)-frequency).^2);
        qnf_err = sum(abs(qnf(noisy_signal)-frequency).^2);
        
        results = [ results; snr discperiod_err pisarenko_err qnf_err];
    end
end

figure(1);
clf
semilogy(results(:,1),results(:,2))
hold on;
semilogy(results(:,1),results(:,3),'m')
semilogy(results(:,1),results(:,4),'g')