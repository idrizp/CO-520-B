%% Prelab:
close all
clear all

% Prelab Problem 2:
Fs = 1e6;
t = 0:1/Fs:0.01 - 1/Fs; % 5 seconds
f = 20E3; % 20 Khz

A_c = 5;

f_m = 5E2;
x = sin(2*pi*f_m*t);

f_c = 20E3;
y = A_c*(1+0.5*x).*cos(2*pi*f_c*t);

subplot(3, 2, 1);
plot(t, x);
title('Modulating Signal');
xlabel('Time(in s)');
ylabel('Voltage(in V)');

subplot(3, 2, 2);
plot(t, y);
title('Amplitude Modulated Signal');
xlabel('Time (in s)');
ylabel('Voltage (in V)');

% Frequency spectrum
N = length(y);
Y = fft(y);

spectrum = abs(Y/N);
spectrum_single = spectrum(1:N/2+1);
spectrum_single(1:end-1) = 2*spectrum_single(1:end-1);

F = Fs * (0:(N/2)) / N;
subplot(3, 2, 3);
plot(F, spectrum_single);
xlim([1.9E4, 2.1E4]);
title('Amplitude Modulated Signal Frequency Spectrum');
xlabel("Frequency(in Hz)")
ylabel("Amplitude(in V)")

frequencies = logspace(3, 5, 100); % 100 Hz to 100 KHz
Wn = 1000/(Fs/2); % Normalized cutoff frequency

[b1, a1] = butter(1, Wn); % Butterworth filter of first order
[b3, a3] = butter(3, Wn); % Butterworth filter of third order

% Rectify the signal

% We can use an envelope detector for this, the simplest analogue to a
% diode would be the absolute value of the signal, so that's what will be
% used

% (Ac + abs(Ac)) / 2;
% Or -0.5 for the diode effect
rectified = (A_c + abs(y)) / 2 - 0.5;
filtered = filter(b1, a1, rectified);

% figure;
subplot(3, 2, 4);
plot(t, filtered);
title('First-Order Demodulated Signal');
xlabel('Time(in s)');
ylabel('Voltage(in V)');

filtered3 = filter(b3, a3, rectified);
subplot(3, 2, 5);
plot(t, filtered3);
title('Third-Order Demodulated Signal');
xlabel('Time(in s)');
ylabel('Voltage(in V)');

% Bode plot for the filters
figure;
freqs(b1, a1);
title("First Order Butterworth Filter Bode Plot");

figure;
freqs(b3, a3);
title("Third Order Butterworth Filter Bode Plot");


% Lab copies

%% Part 1:

% First Hardcopy - 50% Modulation Index Amplitude
% Second Hardcopy - 50% Modulation Index Frequency
% Third Hardcopy - 70% Modulation Index Frequency
% Fourth Hardcopy - 70% Modulation Index Amplitude
% Fifth Hardcopy - 120% Modulation Index Amplitude
% Sixth Hardcopy - 120% Modulation Index Frequency
% Seventh Hardcopy - 70% Modulation Index Frequency and Amplitude of FFT
% Eighth Hardcopy - Circuit 1st order demodulation hardcopy - just a view
% Ninth Hardcopy - Circuit 3rd order demodulation hardcopy - amplitude
% Tenth Hardcopy - Circuit 3rd order demodulation hardcopy - FFT amplitude
% Eleventh Hardcopy - Circuit 3rd order demodulation hardcopy - 20kHz
% component

%% Evaluation: Part 2
close all
clc

w_c = 20E3 * 2 * pi;
w_m = 5E2 * 2 * pi;
c = 5*sin(w_c*t);
y = (1 + 0.70*cos(w_m * t)) .* c;

N = length(y);
Y = fft(y);

spectrum = abs(Y/N);
spectrum_single = spectrum(1:N/2+1);
spectrum_single(1:end-1) = 2*spectrum_single(1:end-1);

F = Fs * (0:(N/2)) / N;
plot(F, 20*log10(spectrum_single / sqrt(2)));
xlim([1.9E4, 2.1E4]);
title('Amplitude Modulated Signal Frequency Spectrum');
xlabel("Frequency(in Hz)")
ylabel("Amplitude(in V)")