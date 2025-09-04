% Small rotor thrust calculation with hub cutoff and 12 blades
clear; close all; clc;

% ---------------- Parameters ----------------
N = 12;            % number of blades
R = 0.034;         % blade radius (m)
Rhub = 0.015;      % hub radius (m)
RPM = 34000;       % rotational speed
rho = 1.225;       % air density (kg/m^3)

% Nondimensional span range (from hub to tip)
r = linspace(Rhub/R, 1.0, 2000);

% ---------------- Alpha distribution ----------------
r_alpha_points = [0.3, 0.51, 1.0];       % given nondimensional stations
alpha_points_deg = [45, 20, 10];         % corresponding alpha (deg)
alpha_r_deg = interp1(r_alpha_points, alpha_points_deg, r, 'linear', 'extrap');
alpha_r = deg2rad(alpha_r_deg);          % convert to radians

% ---------------- Chord distribution (fit hyperbola) ----------------
r_chord_points = [0.3, 0.51, 1.0];
C_points = [0.02, 0.025, 0.02];   % chord values (m)
A_fit = [1./r_chord_points(:), ones(size(r_chord_points(:)))];
x_fit = A_fit \ C_points(:);      % least-squares [a; b]
a = x_fit(1); b = x_fit(2);
C_r = a ./ r + b;
C_r(C_r < 0) = eps;  % avoid negative chord

% ---------------- Lift coefficient (piecewise) ----------------
Cl_r = zeros(size(r));
Cl_r(r <= 0.51) = 0.39 + 0.687  .* alpha_r(r <= 0.51);
Cl_r(r > 0.51)  = 0.31 + 0.0916 .* alpha_r(r > 0.51);

% ---------------- Differential thrust coefficient integrand ----------------
prefactor = 0.5 * (N / (pi * R));
dCT = prefactor .* C_r .* Cl_r .* (r.^2);

% ---------------- Integrate to get CT ----------------
CT = trapz(r, dCT);

% ---------------- Rotor disc area ----------------
A_rotor = pi * R^2;

% ---------------- Tip speed ----------------
omega = 2*pi*RPM/60;   % rad/s
Vtip = omega * R;

% ---------------- Total thrust ----------------
T = rho * A_rotor * Vtip^2 * CT;

% ---------------- Display results ----------------
fprintf('Fitted hyperbola: C(r) = %.6f./r + %.6f\n', a, b);
fprintf('Thrust coefficient C_T = %.6f\n', CT);
fprintf('Total thrust T = %.4f N (%.2f grams)\n', T, T*1000/9.81);

% ---------------- Plots ----------------
figure('Units','normalized','Position',[0.12 0.12 0.65 0.7]);

subplot(2,2,1);
plot(r, rad2deg(alpha_r), 'LineWidth', 1.4);
hold on; plot(r_alpha_points, alpha_points_deg, 'ro','MarkerFaceColor','r');
xlabel('r'); ylabel('\alpha (deg)'); grid on; title('Alpha vs r');

subplot(2,2,2);
plot(r, C_r, 'LineWidth', 1.4);
hold on; plot(r_chord_points, C_points, 'ro','MarkerFaceColor','r');
xlabel('r'); ylabel('Chord C(r) [m]'); grid on; title('Hyperbolic Chord Distribution');
legend('C_{fit}','Given data','Location','best');

subplot(2,2,3);
plot(r, Cl_r, 'LineWidth', 1.4);
xlabel('r'); ylabel('C_l'); grid on; title('Piecewise C_l(r)');

subplot(2,2,4);
plot(r, dCT, 'LineWidth', 1.4);
xlabel('r'); ylabel('dC_T/dr'); grid on; title('Differential integrand dC_T');

sgtitle(sprintf('C_T = %.6f   |   Thrust = %.4f N (%.2f grams)', ...
    CT, T, T*1000/9.81),'FontWeight','Bold');
