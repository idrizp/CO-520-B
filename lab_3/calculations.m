% Prelab:

%% Problem 2:
t = -0.5:0.01:0.5;
f = 4*t.^2;
f_fs = 1/3;

% Fourier series representation of the function
for i=1:5
    % w_0 = 2pi
    % T = 1
    a_i = (4*(-1)^i)/(pi^2*i^2);
    f_fs = f_fs + (a_i .* cos(i.*2.*pi.*t));
end
hold on

plot(t, f, "blue");
plot(t, f_fs, "black");
xlim([-0.5, 0.5]);

legend({"Original Function", ...
    "Fourier Series" ...
    }, "Location", "southwest")

%% Problem 3:

%% Time plots

[t_fifty, signal_fifty, f_fifty, y_fifty] = get_square(50);
[t_thirty, signal_thirty, f_thirty, y_thirty] = get_square(30);
[t_twenty, signal_twenty, f_twenty, y_twenty] = get_square(20);

subplot(3, 1, 1);
plot(t_fifty, signal_fifty, "blue", "LineWidth", 2);
title("Duty Cycle: 50%");

subplot(3, 1, 2);
plot(t_thirty, signal_thirty, "blue", "LineWidth", 2);
title("Duty Cycle: 30%");

subplot(3, 1, 3);
plot(t_twenty, signal_twenty, "blue", "LineWidth", 2);
title("Duty Cycle: 20%");

%% Frequency plots
subplot(3, 1, 1);
plot(f_fifty, y_fifty, "blue", "LineWidth", 1);
title("FT: Duty Cycle: 50%");


subplot(3, 1, 2);
plot(f_thirty, y_thirty, "blue", "LineWidth", 1);
title("FT: Duty Cycle: 30%");


subplot(3, 1, 3);
plot(f_twenty, y_twenty, "blue", "LineWidth", 1);
title("FT: Duty Cycle: 20%");


%% Problem 4:

[y, Fs] = audioread("s_samp.wav");

N_samples = Fs * 10E-3;
y_first = y(1:N_samples);
t = 0:1/Fs:((10E-3)-1/Fs);
plot(t, y_first);
%% 

N_samples = length(y);
Y = fft(y);
rms = sqrt(mean(abs(Y.^2)));

Fs_nyquist = Fs / 2; 

Y_single = 2 * abs(Y) / N_samples;
Y_single = Y_single(1:floor(N_samples/2));
f = linspace(0,Fs_nyquist,length(Y_single));

Y_db = 20*log10(Y_single ./ rms);
Y_db(Y_db == -Inf) = 0;

% plot(f, Y_single);

plot(f, Y_db);
ylabel("Amplitude(in dBVrms)")
xlabel("Frequency (in Hz)");

%% Lab start
clear
clc

% Problem 1:

% Oscilloscope found to be incalibrated perhaps, 5Hz error, 1% error when
% it should be 0.1% error.
% Vpp used for 0dB: 2.7Vpp, delta dB = 46.8dB, Cursor dB: -148mdB

% Problem 2:

% First harmonic: -989mdB, 1050Hz
% Second harmonic: -10.1dB, 2950Hz
% Third harmonic: -14.9dB, 5000Hz
% Fourth harmonic: -17.7dB, 7000Hz

% 250kS/s was used

% First non-DC harmonic: -5.39dB, 1010Hz
% Second non-DC harmonic: -7.39dB, 2010Hz
% Third non-DC harmonic: -10.9dB, 2970Hz
% Fourth non-DC harmonic: -17.3dB, 3980Hz

%50kS/s

% Problem 3:

% First harmonic: -26.9dB, 1000Hz
% Second harmonic: -23.7dB, 9900Hz
%% Functions

function [t, signal, f, y_single_db] = get_square(duty_cycle)
    period = 1e-3;
    Fs = 200e3;
    frequency = 1/period;
    Vpp = 2;
    duration = period * 6;
    
    % The signal
    t = 0:1/Fs:duration;
    signal = Vpp*square((2*pi*frequency)*t, duty_cycle);
    plot(t, signal, "blue", "LineWidth", 2);
    
    % The fourier transform of the signal
    rms_value = sqrt(mean(signal.^2));
    N = length(signal); % The length of the signal
    
    y = fft(signal, N);
    y_mag = 2*(abs(y)/N); % Magnitudes of y
    
    y_single = y_mag(1:floor(N/2)) * 2;
    f_nyquist = Fs / 2;
    
    y_single_db = 20*log10(y_single / rms_value);
    f = linspace(0, f_nyquist, length(y_single));
end