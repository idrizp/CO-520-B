% Prelab:

% Given is a series RLC resonant circuit with R = 390Ω, C = 270 nF 
% and L = 10 mH.
% 1. Name the filter characteristic measured over the different components, com-
% ponent combinations.

% 2. Show the Bode magnitude plot across the resistor, 
% the capacitor, the inductor and across the capacitor 
% and the inductor together. Use a 5 V amplitude and vary the
% frequency starting at 100 Hz up to 100 KHz.

% Develop a Matlab script to plot the four characteristic in one graph.
% Attach the script to the prelab!
% 3. Taking the magnitude across the resistance represents a band-pass filter. 
% Calculate the bandwidth and the Q factor of the circuit. 
% Extract the bandwidth from the Matlab plot and compare.

R = 390; % in ohm
C = 270E-9; % in F
L = 10E-3; % in H

w = 100:100:100E3;
jw = 1i*w;

% Voltage over the resistor
H_r = (jw .* (R .* C)) ./ ((jw.^2 .* L .* C) + (jw .* R .* C) + 1);

semilogx(w, 20 * log10(abs(H_r)), "blue");

hold on

% Voltage over the inductor
H_l = (jw.^2 .* L .* C) ./ ((jw.^2 .* L .* C) + (jw .* R .* C) + 1);
semilogx(w, 20 * log10(abs(H_l)), "red");

% Voltage over the capacitor
H_c = 1./((jw.^2 .* L .* C) + (jw .* R .* C) + 1);
semilogx(w, 20 * log10(abs(H_c)), "green");

% Voltage over the inductor and the capacitor
H_lc = ((jw.^2 .* L .* C) + 1) ./ ((jw.^2 .* L .* C) + jw .* R .*C + 1);
semilogx(w, 20 * log10(abs(H_lc)), "black");
ylim([-50, 0]);

legend("H_R", "H_L", "H_C", "H_LC", "Location", "southwest")

% Calculated bandwidth and quality factor:
B_calculated = R/L;
w_0 = 1/sqrt(L*C);
X_0 = sqrt(L/C);
Q_s = X_0/R;

% Obtained from the plot, respectively at their 11dB cutoff
% points: -20log10(1/sqrt(2) * 5) and 20log10(1/sqrt(2)*5)
w_2 = 4.69E4;
w_1 = 7.9E3;

B_plot = w_2 - w_1;

disp("Plot Bandwidth: " + (w_2 - w_1));

disp("Calculated Quality Factor: " + Q_s);
disp("Calculated Bandwidth: " + B_calculated);

%% Lab data results

% V_C hardcopy first - Low pass filter
% V_L hardcopy second - High pass filter 
% V_R hardcopy third - Band pass filter
% V_LC hardcopy fourth - Band-stop/Notch filter
% At resonance, the Lissajou is linear(a line), at 3100 Hz we see the
% Lissajou turn linear, so that must be the resonant frequency

% Phase shift: -500mDeg to 0 Deg(more conceretely, 0 deg)

% Corner frequencies: 1. 1290Hz  2. 7640Hz
