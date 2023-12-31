%% Problem 2:

% Use the function generator to generate a sinusoidal wave having 
% 500 Hz frequency, 2 Vpp amplitude and no offset. 
% Use the measure function to verify all properties. 
% Take a hard copy in time domain.

F = 500;
T = 1/F;
Fs = 1E6;
Ts = 1/Fs;
w = 2*pi*F;

Vpp = 2;
Vamp = Vpp / 2;
t = 0:Ts:T - Ts;
signal = Vamp * sin(w*t);

subplot(2, 1, 1);
plot(t, signal);
title("Signal");
ylabel("Amplitude(in V)");
xlabel("Time(in s)")

N = length(signal);

y = fft(signal);
y = 2 * (abs(y) / N);
y = y(1:floor(N/2));

f = linspace(0, Fs/2, length(y));
subplot(2, 1, 2);
plot(f, mag2db(y / sqrt(2)));
title("Frequency Spectrum");

xlim([-0, 1e4]);
ylabel("Amplitude(in VdB)");
xlabel("Frequency(in Hz)")

%% Problem 3
close all
clear

% Generate a sinusoidal wave having 0 dB spectrum peak, 2 KHz frequency, 
% without a dc offset. 
% What is the amplitude value? 
% Use the measure function and the cursors. 
% Take hard copies of time and frequency domain.

F = 2000;
T = 1/F;
Fs = 1E6;
Ts = 1/Fs;
w = 2*pi*F;

t = 0:Ts:T - Ts;
V_amp = sqrt(2);
signal = V_amp*sin(w*t);

subplot(2, 1, 1);
plot(t, signal);
title("Signal");
ylabel("Amplitude(in V)");
xlabel("Time(in s)")

N = length(signal);

y = fft(signal);
y = 2 * (abs(y) / N);
y = y(1:floor(N/2));

f = linspace(0, Fs/2, length(y));
subplot(2, 1, 2);
plot(f, floor(mag2db(y / sqrt(2))));
title("Frequency Spectrum");

xlim([-0, 1e5]);
ylabel("Amplitude(in VdB)");
xlabel("Frequency(in Hz)")
