%% Prelab:

% first part:

f = 50; % frequency, in hertz
f_n = 2*f; % nyquist frequency, in hertz

samples=1E5;
t = 0:1/samples:1; % 100kS

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
xlim([0, 0.5]);

% Dirac comb
subplot(3, 2, 2);
p = zeros(size(t));
p(1:10750:length(t)) = 1;
stem(t, p);
ylabel("Voltage")
xlabel("Time")
title("Dirac comb");
ylim([0, 1]);
xlim([0, 1]);

% Undersampled signal:
fs = 48; % In hertz

p = zeros(size(t));
p(1:samples/fs:length(t)) = 1;

subplot(3, 2, 3);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Undersampled signal");
ylim([-10, 10]);
xlim([0, 0.5]);

% Nyquist-Frequency sampled signal
fs = 100; % In hertz

p = zeros(size(t));
p(1:samples/fs:length(t)) = 1;

subplot(3, 2, 4);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Nyquist-Frequency signal");

ylim([-10, 10]);
xlim([0, 0.5]);

% Oversampled signal
fs = 1000; % In hertz

p = zeros(size(t));
p(1:samples/fs:length(t)) = 1;

subplot(3, 2, 5);
plot(t, p .* x);
ylabel("Voltage")
xlabel("Time")
title("Oversampled signal");

ylim([-10, 10]);
xlim([0, 0.5]);

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


% Nyquist Sampled signal:
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

% Overly Sampled signal:
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


%% 
% Used: 2.5VPP at function generator, saw 5Vpp on oscilloscope

% First: 24900Hz Display

% Second: Dots display
% Third: Vector display
% Fourth: 25000Hz - If oscilloscope and function generator were
% synchronized perfectly, you would only see a straight line.
% Fifth: 25020Hz
% Sixth: 25500Hz

% Seventh: 150Hz - Synced with the sync output and source set to Ext.
% Eighth: 200Hz 

% At 24900Hz, 100Hz Frequency Displayed on Oscilloscope

% At 25000Hz, only a line displayed on the oscilloscope, but we also
% observe that the line has a very very small frequency to it(this is an
% error) - an aliasing difference between the oscilloscope and function
% generator frequency

% Pre-lab: The reason Nyquist can be both zero and 2.5 according to the formula depending on the rounding,
% is it depends on the
% phase shift, but for others it does not matter

% We observed a square wave for the DC input that followed the
% variation(for 2V - approximately 216mVpp, for 3V - approximately 312mVpp)

% The VPP was set to 0.75Vpp on the function generator.
