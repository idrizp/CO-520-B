%% Prelab:

% first part:

f = 50; % frequency, in hertz
f_n = 2*f; % nyquist frequency, in hertz

t = 0:1E-5:1; % smaller iterations lead to closer to continuous behavior

V_pp = 5;
w = 2*pi*f;
x = V_pp * sin(w .* t);

% Original signal:
subplot(3, 2, 1);
plot(t, x);
ylabel("Voltage")
xlabel("Time")
title("Original signal");
ylim([-10, 10]);
xlim([0, 1/f]);

% Dirac comb
T = 1/f_n;
subplot(3, 2, 2);
p = zeros(size(t));
p(mod(t, T) == 0) = 1;
stem(t, p);
ylabel("Voltage")
xlabel("Time")
title("Dirac comb");
ylim([0, 1]);
xlim([0, 10*T]);

% Undersampled signal:
fs = 48; % In hertz
T = 1/fs;

p = zeros(size(t));
p(mod(t, T) == 0) = 1;

subplot(3, 2, 3);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Undersampled signal");
ylim([-10, 10]);
xlim([0, 1/2]);

% Nyquist-Frequency sampled signal
fs = 100; % In hertz
T = 1/fs;

p = zeros(size(t));
p(mod(t, T) == 0) = 1;

subplot(3, 2, 4);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Nyquist-Frequency signal");

ylim([-10, 10]);
xlim([0, 1/50]);


% Oversampled signal
fs = 1000; % In hertz
T = 1/fs;

p = zeros(size(t));
p(mod(t, T) == 0) = 1;

subplot(3, 2, 5);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Oversampled signal");

ylim([-10, 10]);
xlim([0, 1/50]);

%%

clear
close all

f = 50;
t = 0:1E-5:1;
V_pp = 5;
w = 2*pi*f;
x = V_pp * sin(w .* t);

% Underdamped signal:
fs = 48; % In hertz
T = 1/fs;
T0 = 0.5 * 1/fs;
p = gen_pulse_train(t, T0, T);

subplot(3, 2, 1);
plot(t, p);
title("Rectangular Pulse Train(Undersampled)");
ylabel("Voltage")
xlabel("Time")
xlim([0, 1/2]);
ylim([-0.5, 2]);

subplot(3, 2, 2);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Undersampled Signal")
xlim([0, 1/2]);
ylim([-6, 6]);


% Nyquist Damped signal:
fs = 100; % In hertz
T = 1/fs;
T0 = 0.5 * 1/fs;
p = gen_pulse_train(t, T0, T);

subplot(3, 2, 3);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Nyquist Sampled Signal")
xlim([0, 1/50]);
ylim([-6, 6]);

% Nyquist Sampled signal:
fs = 1000; % In hertz
T = 1/fs;
T0 = 0.5 * 1/fs;
p = gen_pulse_train(t, T0, T);

subplot(3, 2, 4);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Oversampled Signal")
xlim([0, 1/50]);
ylim([-6, 6]);

subplot(3, 2, [5 6]);
plot(t, x);
ylabel("Voltage")
xlabel("Time")
title("Original Signal")
xlim([0, 1/50]);
ylim([-6, 6]);

% Function to generate pulse train
function [p] = gen_pulse_train(t, T0, T)
    p = zeros(size(t));
    for k=1:length(t)
        if mod(t(k), T) == 0
           left = t(k) - T0;
           right = t(k) + T0;
           p((t > left) & (t < right)) = 1;
        else
            p(k) = 0;
        end
    end
end

