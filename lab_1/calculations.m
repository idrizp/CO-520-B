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
clear 

% Calculated
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

w_d_calculated = w_n * sqrt(1 - zeta^2); % The damped angular frequency calculated using LC
f_d_nominal = w_d_calculated / (2*pi);

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

%% Part 1:
clear

% For the case where R = 100Ohm, C = 6.8nF and L=10mH
R = 100;
C = 6.8E-9;
L = 10E-3;

zeta = R/2 * sqrt(C/L); % approximately 0.0412, so underdamped
w_n = 1/sqrt(L*C);
w_d = w_n * sqrt(1 - zeta^2);

% We know that the circuit is underdamped because zeta < 1

% We obtained C1 = -1, C2 = -(w_n/w_d) * zeta
C1 = -1;
C2 = -(w_n/w_d) * zeta;

t = 0:1E-6:1E-3;
y = exp(-zeta * w_n .* t) .* (C1*cos(w_d.*t) + C2 * sin(w_d.*t)) + 1;
plot(t, y, 'red');

hold on

%% Critically damped

% For the critically damped case:
R = 2 * (1 / sqrt(C/L));
zeta = R/2 * sqrt(C/L);
w_n = 1/sqrt(L*C);

% We obtained that C1 = -1, and C2 = -w_n.
C2 = -w_n;
C1 = -1;
y = (C1*exp(-zeta * w_n .* t) + C2.*t.*exp(-zeta * w_n .* t)) + 1;
plot(t, y, 'blue');

legend({'Underdamped','Critically Damped'},'Location','southwest')

%% Part 3:
clear
close all

R1 = 25;
R2 = 56;
L = 20E-3;
C = 2E-6;
Vin = 16.2;

a2 = R1*L*C;
a1 = R1*R2*C + L;
a0 = R1+R2;

w_n = sqrt(a0 / a2);
zeta = a1 / (2*sqrt(a0 * a2));
K = 1/a0;

iL_p = Vin / (R1 + R2);
C2 = iL_p * (zeta - sqrt(zeta^2 - 1)) / (2 * sqrt(zeta^2 - 1));
C1 = -C2 - iL_p;

% The current over the inductor.
t = 0:1E-6:1E-2;
% It is overdamped.
y = C1 .* exp((-zeta + sqrt(zeta^2 - 1)) .* t .* w_n) + ...
        C2 .* exp((-zeta - sqrt(zeta^2 - 1)) .* t .* w_n) + iL_p;

plot(t, y, "blue", "LineWidth", 2);
ylim([0, 0.3])
legend("Overdamped Behavior")

