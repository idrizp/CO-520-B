R = 390; % in ohm
C = 270E-9; % in F
L = 10E-3; % in H

w = 100:100:100E3;
jw = 1i*w;

% Voltage over the resistor
H_r = (jw .* (R .* C)) ./ ((jw.^2 .* L .* C) + (jw .* R .* C) + 1);

semilogx(w, 20 * log10(abs(H_r)), "blue");

hold on

H_l = (jw.^2 .* L .* C) ./ ((jw.^2 .* L .* C) + (jw .* R .* C) + 1);
semilogx(w, 20 * log10(abs(H_l)), "red");

H_c = 1./((jw.^2 .* L .* C) + (jw .* R .* C) + 1);
semilogx(w, 20 * log10(abs(H_c)), "green");

H_lc = ((jw.^2 .* L .* C) + 1) ./ ((jw.^2 .* L .* C) + jw .* R .*C + 1);
semilogx(w, 20 * log10(abs(H_lc)), "black");
ylim([-50, 0]);

legend("H_R", "H_L", "H_C", "H_LC", "Location", "southwest")