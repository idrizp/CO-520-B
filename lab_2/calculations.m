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

H_R = tf([R*C, 0], [L*C, R*C, 1]);
% The transfer function of the voltage taken across the resistor shows a
% bandpass filter.

H_L = tf([L*C, 0, 0], [L*C, R*C, 1]);
% The transfer function of the voltage across the inductor shows a
% high-pass filter.

H_C = tf(1,  [L*C, R*C, 1]);
% The transfer function of the voltage across the capacitor shows a
% low-pass filter.

H_LC = tf([L*C, 0, 1], [L*C, R*C, 1]);
% The transfer function of the voltage taken across the inductor and the
% capacitor shows a band-stop filter.

% Calculated bandwidth and quality factor:
B_calculated = R/L;
w_0 = 1/sqrt(L*C);
X_0 = sqrt(L/C);
Q_s = X_0/R;

bodemag(H_R,H_L,H_C,H_LC,{1E2, 1E5});
ylim([-100, 50])

% Obtained from the plot, respectively at their 11dB cutoff
% points: -20log10(1/sqrt(2) * 5) and 20log10(1/sqrt(2)*5)
w_2 = 4.73E4;
w_1 = 7.81E3;

B_plot = w_2 - w_1;

legend("H_R", "H_L", "H_C", "H_LC", 'Location', 'southwest');

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
