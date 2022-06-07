rng(123456789);
phase_1 = pi/4;
T = 32;
NT = 8192;
x = cos(2*pi*0.0123456789*[0:T-1] + phase_1) + cos(2*pi*(0.0123456789 + 1/200)*[0:T-1] + pi/4) + 0.5*randn(1,T);
clf
periodogram = abs(fft(x, NT)).^2;
mle = mle_fourier(x, NT);

figure(1);
clf;
semilogy([0:(NT-1)]/NT*2*pi, periodogram/max(periodogram), ':', 'LineWidth',10)
hold on
semilogy([0:(NT-1)]/NT*2*pi, mle/max(mle),'r', 'LineWidth',5)
xlabel('Normalized frequency', 'FontSize',32)
ax = gca; 
ax.FontSize = 24; 

figure(2);
clf;
semilogy([0:(NT-1)]/NT*2*pi, periodogram/max(periodogram), ':', 'LineWidth',10)
hold on
semilogy([0:(NT-1)]/NT*2*pi, mle/max(mle),'r', 'LineWidth',5)
xlabel('Normalized frequency', 'FontSize',32)
axis([0 0.5  0.0001 1])
ax2 = gca; 
ax2.FontSize = 24; 

figure(3);
clf;
semilogy([0:(NT-1)]/NT*2*pi, periodogram/max(periodogram), ':', 'LineWidth',10)
hold on
semilogy([0:(NT-1)]/NT*2*pi, mle/max(mle),'r', 'LineWidth',5)
xlabel('Normalized frequency', 'FontSize',32)
axis([pi-0.1 pi+0.1  0.000001 1])
ax3 = gca; 
ax3.FontSize = 24; 