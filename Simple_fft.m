DeltaT = 0.0001;
Fs = 1/DeltaT;
L = length(x);
t = (0:L-1)*DeltaT;

figure;
plot(x,y);

Y = fft(y);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
plot(f,P1); 
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');