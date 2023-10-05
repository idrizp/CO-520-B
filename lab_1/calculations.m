%% Theory

% w_n = 1/sqrt(L*C); % The undamped natural frequency (in rad/s)
% zeta = (R / 2) * sqrt(C/L); % The damping ratio
% 2*zeta = R*sqrt(C/L);
% R = 2*zeta/sqrt(C/L);
% k = 1; % The gain
% overshoot:
% percentage_overshoot = (v_max - v_steady_state) / (v_steady_state) * 100;

% CASES: 
% zeta < 1: Underdamped
% zeta = 1: Critically damped
% zeta > 1: Overdamped

%% Lab Calculations

% % Obtained from hard-copy
t_peak = 0; % The time at which the sinusoid peaks
t_peak_2 = 1.5; % The time at which the sinusoid peaks again.

% Calculated
period = t_peak_2 - t_peak;

f = 1/period; % The frequency
R_i = 50; % R-internal
R_decade = 1e2; % R-decade
L = 10E-3; % Inductor
C = 6.8E-9; % Capacitor

R = R_decade + R_i; % Account for both the Resistor Decade and the internal Resistance

w_n = 1/sqrt(L*C); % The undamped natural frequency(in rad/s)
zeta = (R / 2) * sqrt(C/L); % Damping ratio

% Oscilloscope measurements:

period = 52E-6; % The period between the oscillations of the sinusoid caused by the damping.
f_d = 1/period; % The frequency of the oscillations of the sinusoid caused by the damping.

w_d_calculated = 1 / sqrt(L*C); % The damped angular frequency calculated using LC
w_d_measured = (2*pi) * f_d; % The damped 

w_d_diff = (w_d_calculated - w_d_measured) / (w_d_calculated) * 100; % The difference in the damped angular frequency.

% Taken from R = (2*zeta) / sqrt(C/L) where our optimal "zeta" is one.
optimal_r = 2*(1 / sqrt(C/L)) - 50; % The optimal resistance to achieve critically damped behaviour.
% Subtract 50 to get rid of internal resistance

% Other data:

% Resistance by playing around: 1905
% Varied resistance: 2375
% Varied resistance: 2425 + 10K

% 10% tolerance for the capacitor
% 20% tolerance for the inductor

 
%% Evaluation
close all

% y(t) = exp (−ζωnt)(C1 cos (ωdt) + C2 sin (ωdt))

% Underdamped
% C1: y(0) = C1 + 0, 1 = C1 + 0, C1 = 1.
% C2: y(infty) = 0: C2 = 0

% For the underdamped case, C1 = -1, C2 = 0
t = 0:1E-9:1E-3;
y = -exp(-zeta.*w_n.*t).*cos(w_d_measured.*t);
plot(t, y, 'r', 'DisplayName', 'Underdamped');

hold on

% For the critically damped case, C1 = -1, C2 = 0
zeta_optimal = (optimal_r / 2) * sqrt(C/L); % Damping ratio with optimal resistance
y = -exp(-zeta_optimal.*w_n.*t);
plot(t, y, 'b', 'DisplayName', 'Critically Damped')

legend('Location', 'best')
