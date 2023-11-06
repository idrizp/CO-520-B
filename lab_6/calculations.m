%% Prelab 6:

% Part 2:
% Problem 2: FM signal in the frequency domain
% Plot a frequency modulated signal in the frequency domain. 
% The signal exhibits 2.5 V peak carrier amplitude, 40 kHz carrier frequency, 5 kHz modulation frequency. 
% Vary βf between 0.2, 1, and 2. Display the magnitudes in dBrms! Calculate the bandwidth using Carlsons rule. 
% Tabulate the peak magnitudes inside the bandwidth from the three plots.

Fs = 1E6;
t = 0:1/Fs:0.1-(1/Fs);
w_c = 40E3 *2*pi;
w_m = 5E3*2*pi;

% For m = 0.2
m = 0.2;
y = 2.5.*cos(w_c.*t + m.*sin(w_m.*t));

% subplot(3, 1, 1);
% plot(t, y);

N = length(y);
Y = fft(y);

spectrum = abs(Y/N);
spectrum_single = spectrum(1:N/2+1);
spectrum_single(1:end-1) = 2*spectrum_single(1:end-1);

F = Fs * (0:(N/2)) / N;
subplot(3, 1, 1);

BT = w_m/pi*(0.2+1);
disp("Bandwidth for m=0.2");
disp(BT);

m_02_dB = 20*log10(spectrum_single / sqrt(2));
disp("For m=0.2, peaks are");
display_peaks(m_02_dB);

plot(F, m_02_dB);
xlabel("Frequencies(in Hz)");
ylabel("Amplitudes(in VdBrms)");
xlim([0, 75000]);
ylim([-40, 40]);

% For m = 1
m = 1;
y = 2.5.*cos(w_c.*t + m.*sin(w_m.*t));

N = length(y);
Y = fft(y);

spectrum = abs(Y/N);
spectrum_single = spectrum(1:N/2+1);
spectrum_single(1:end-1) = 2*spectrum_single(1:end-1);

F = Fs * (0:(N/2)) / N;
subplot(3, 1, 2);

m_1_dB = 20*log10(spectrum_single / sqrt(2));

BT = w_m/pi*(1+1);
disp("Bandwidth for m=1");
disp(BT);
disp("For m=1, peaks are");
display_peaks(m_1_dB);

plot(F, m_1_dB);
xlabel("Frequencies(in Hz)");
ylabel("Amplitudes(in VdBrms)");
xlim([0, 75000]);
ylim([-40, 40]);

% For m = 2
m = 2;
y = 2.5.*cos(w_c.*t + m.*sin(w_m.*t));

N = length(y);
Y = fft(y);

% BT ∼=2fm(βf +1)
BT = w_m/pi*(2+1);
disp("Bandwidth for m=2");
disp(BT);

spectrum = abs(Y/N);
spectrum_single = spectrum(1:N/2+1);
spectrum_single(1:end-1) = 2*spectrum_single(1:end-1);

F = Fs * (0:(N/2)) / N;
m_2_dB = 20*log10(spectrum_single / sqrt(2));

% Peaks:
disp("For m=2, peaks are:");
display_peaks(m_2_dB);

subplot(3, 1, 3);
plot(F, m_2_dB);
xlabel("Frequencies(in Hz)");
ylabel("Amplitudes(in VdBrms)");
xlim([0, 75000]);
ylim([-40, 40]);

function [] = display_peaks(m_dB)
    frequencies = find(m_dB > -40); % more than -40dB
    decibels = m_dB(frequencies);
    for i = 1:length(frequencies)
        fprintf("Frequency: %d Hz, Decibel: %f dB\n", frequencies(i), round(decibels(i)));
    end
end

