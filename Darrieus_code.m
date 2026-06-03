% Darrieus VAWT simulation - Double Multiple Streamtube BEM
% FIRST-PRINCIPLES LEARNING EXERCISE. This hand-written solver does NOT
% produce physically valid results for this rotor (see project writeup).
% It is kept as a documented attempt; the validated analysis was done in QBlade,
% which diverged in the same way - establishing that the rotor's solidity
% (Nc/R ~ 0.45) is beyond the reliable range of streamtube theory.
%
% Geometry from Darrieus_Motor.FCStd (confirmed visually + bounding-box derived):
%   N = 2 (confirmed: curved-blade troposkein "egg-beater" rotor)
%   c = 0.067 m (measured), H = 1.0 m and R = 0.30 m 
% Note: R treated as constant (straight-blade approximation); real blades curve.
% NACA 0018 polar, base MATLAB, no toolboxes.

clear; clc; close all;

% Turbine parameters
R   = 0.300;
H   = 1.000;
N   = 2;
c   = 0.067;
rho = 1.225;
sigma = N * c / (2 * pi * R);

fprintf('Darrieus VAWT: R=%.3f m, H=%.3f m, N=%d, c=%.3f m, solidity=%.4f\n', ...
        R, H, N, c, sigma);
fprintf('Swept area: %.4f m^2\n\n', 2*R*H);

% NACA 0018 polar (representative, Sheldahl & Klimas-style, full +/-180 deg)
alpha_data = [-180,-170,-160,-150,-140,-130,-120,-110,-100,-90,...
               -80, -70, -60, -50, -40, -30, -20, -15, -10,  -8,...
                -6,  -4,  -2,   0,   2,   4,   6,   8,  10,  12,...
                14,  16,  18,  20,  30,  40,  50,  60,  70,  80,...
                90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

Cl_data = [0, 0.40, 0.70, 0.85, 0.80, 0.65, 0.45, 0.20,-0.10,-0.35,...
          -0.55,-0.70,-0.75,-0.72,-0.58,-0.40,-0.65,-0.95,-0.65,-0.52,...
          -0.38,-0.24,-0.12, 0.00, 0.12, 0.24, 0.38, 0.52, 0.65, 0.78,...
           0.85, 0.88, 0.80, 0.65, 0.40, 0.58, 0.72, 0.75, 0.70, 0.55,...
           0.35, 0.10,-0.20,-0.45,-0.65,-0.80,-0.85,-0.70,-0.40, 0.00];

Cd_data = [0.020,0.040,0.080,0.140,0.220,0.300,0.370,0.420,0.440,0.440,...
           0.420,0.390,0.340,0.280,0.220,0.160,0.060,0.035,0.025,0.021,...
           0.016,0.013,0.011,0.010,0.011,0.013,0.016,0.021,0.025,0.032,...
           0.040,0.055,0.070,0.090,0.160,0.220,0.280,0.340,0.390,0.420,...
           0.440,0.440,0.420,0.370,0.300,0.220,0.140,0.080,0.040,0.020];

% Cp vs TSR at several wind speeds
TSR_vec = linspace(0.5, 5.0, 60);
V_cases = [5, 7, 9, 11];
colors  = lines(numel(V_cases));

figure('Name','Cp vs TSR','Position',[100 100 700 450]);
hold on; grid on;
Cp_matrix = zeros(numel(V_cases), numel(TSR_vec));

for vi = 1:numel(V_cases)
    Vinf = V_cases(vi);
    for ti = 1:numel(TSR_vec)
        Cp = dmst_solve(TSR_vec(ti), Vinf, R, H, N, c, rho, alpha_data, Cl_data, Cd_data);
        Cp_matrix(vi, ti) = max(Cp, 0);
    end
    plot(TSR_vec, Cp_matrix(vi,:), 'LineWidth', 2, 'Color', colors(vi,:), ...
         'DisplayName', sprintf('V = %d m/s', Vinf));
end

yline(0.593, '--k', 'Betz Limit', 'LabelHorizontalAlignment','right', ...
      'HandleVisibility','off');
xlabel('Tip Speed Ratio'); ylabel('Power Coefficient C_p');
title('Darrieus VAWT: C_p vs TSR');
legend('Location','northeast'); ylim([0 0.65]);

% Power curve and optimal TSR vs wind speed
V_sweep = linspace(2, 15, 40);
P_out   = zeros(1, numel(V_sweep));
TSR_opt = zeros(1, numel(V_sweep));

for vi = 1:numel(V_sweep)
    Vinf = V_sweep(vi);
    Cp_sweep = zeros(1, numel(TSR_vec));
    for ti = 1:numel(TSR_vec)
        Cp_sweep(ti) = max(dmst_solve(TSR_vec(ti), Vinf, R, H, N, c, rho, ...
                            alpha_data, Cl_data, Cd_data), 0);
    end
    [Cp_max, idx] = max(Cp_sweep);
    TSR_opt(vi) = TSR_vec(idx);
    P_out(vi)   = Cp_max * 0.5 * rho * Vinf^3 * (2*R*H);
end

figure('Name','Power Curve','Position',[820 100 700 450]);
subplot(2,1,1);
plot(V_sweep, P_out, 'b-o', 'LineWidth', 2, 'MarkerSize', 4); grid on;
xlabel('Wind Speed [m/s]'); ylabel('Power [W]'); title('Power Curve');
subplot(2,1,2);
plot(V_sweep, TSR_opt, 'r-s', 'LineWidth', 2, 'MarkerSize', 4); grid on;
xlabel('Wind Speed [m/s]'); ylabel('Optimal TSR'); title('Optimal TSR vs Wind Speed');

% 3D Cp surface over TSR and wind speed
TSR_surf = linspace(0.5, 5.0, 30);
V_surf   = linspace(3, 13, 25);
[TSR_g, V_g] = meshgrid(TSR_surf, V_surf);
Cp_surf = zeros(size(TSR_g));

for vi = 1:numel(V_surf)
    for ti = 1:numel(TSR_surf)
        Cp_surf(vi,ti) = max(dmst_solve(TSR_g(vi,ti), V_g(vi,ti), R, H, N, c, rho, ...
                              alpha_data, Cl_data, Cd_data), 0);
    end
end

figure('Name','Cp Surface','Position',[100 580 700 450]);
surf(TSR_g, V_g, Cp_surf, 'EdgeColor','none'); colorbar; colormap(jet);
xlabel('TSR'); ylabel('Wind Speed [m/s]'); zlabel('C_p');
title('C_p Surface'); view(45, 30);

% Torque ripple: single blade vs combined N blades
TSR_ripple  = 2.5;
Vinf_ripple = 9;
[~, torque_1] = dmst_solve(TSR_ripple, Vinf_ripple, R, H, N, c, rho, ...
                           alpha_data, Cl_data, Cd_data);

torque_total = zeros(size(torque_1));
for b = 1:N
    offset = round((b-1) * numel(torque_1) / N);
    torque_total = torque_total + circshift(torque_1, offset);
end

theta_deg = rad2deg(linspace(0, 2*pi, numel(torque_1)));
ripple = (max(torque_total) - min(torque_total)) / (mean(abs(torque_total)) + eps);

figure('Name','Torque Ripple','Position',[820 580 700 450]);
subplot(2,1,1);
plot(theta_deg, torque_1, 'b-', 'LineWidth', 1.5); grid on;
xlabel('Azimuth [deg]'); ylabel('Torque [N.m]');
title(sprintf('Single Blade (TSR=%.1f, V=%d m/s)', TSR_ripple, Vinf_ripple));
subplot(2,1,2);
plot(theta_deg, torque_total, 'r-', 'LineWidth', 1.5); grid on;
xlabel('Azimuth [deg]'); ylabel('Total Torque [N.m]');
title(sprintf('%d Blades Combined (ripple = %.2f)', N, ripple));

% Animated rotor
figure('Name','Rotor Animation','Position',[100 100 500 500]);
axis equal; axis([-0.55 0.55 -0.55 0.55]); grid on; hold on;
title('Darrieus Rotor'); xlabel('x [m]'); ylabel('y [m]');

th_circ = linspace(0, 2*pi, 200);
plot(R*cos(th_circ), R*sin(th_circ), '--', 'Color', [0.7 0.7 0.7]);

blade_len = c * 8;
blade_h   = c * 0.3;
h_blades  = gobjects(N,1);
for b = 1:N
    h_blades(b) = fill([0 0 0 0], [0 0 0 0], [0.2 0.4 0.8]);
end
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

omega = TSR_ripple * Vinf_ripple / R;
dt = 0.05;
for frame = 1:120
    psi = omega * frame * dt;
    for b = 1:N
        ang = psi + (b-1)*2*pi/N;
        cx = R * cos(ang);
        cy = R * sin(ang);
        corners_x = blade_len/2 * [-1 1 1 -1];
        corners_y = blade_h/2  * [-1 -1 1 1];
        rot = [cos(ang+pi/2) -sin(ang+pi/2); sin(ang+pi/2) cos(ang+pi/2)];
        pts = rot * [corners_x; corners_y];
        set(h_blades(b), 'XData', cx + pts(1,:), 'YData', cy + pts(2,:));
    end
    drawnow limitrate;
    pause(dt);
end

% Summary
[Cp_best, ~] = max(Cp_matrix(:));
fprintf('Peak Cp: %.4f\n', Cp_best);
fprintf('Power at 9 m/s: %.2f W\n', P_out(find(V_sweep >= 9, 1)));


function [Cp, torque_blade] = dmst_solve(TSR, Vinf, R, H, N, c, rho, alpha_data, Cl_data, Cd_data)
    % Double multiple streamtube solver (first-principles attempt).
    % NOTE: this implementation does not produce physically valid Cp for the
    % project rotor; see writeup. Kept as a learning artifact.
    n_theta = 72;
    theta   = linspace(0, 2*pi, n_theta+1);
    theta   = theta(1:end-1);
    dtheta  = 2*pi / n_theta;
    sigma_local = N * c / R;

    % Upwind half
    a_up = 0;
    for iter = 1:50
        CT = 0;
        for i = 1:n_theta/2
            th = theta(i);
            W_u = TSR + sin(th);
            W_v = cos(th);
            W   = hypot(W_u, W_v);
            alpha = atan2(W_v, W_u);
            Cl = interp1(alpha_data, Cl_data, rad2deg(alpha), 'linear', 0);
            Cd = interp1(alpha_data, Cd_data, rad2deg(alpha), 'linear', 0.02);
            Cn = Cl*cos(alpha) + Cd*sin(alpha);
            CT = CT + (sigma_local/(4*pi)) * Cn * W^2 * dtheta;
        end
        if CT > 0.96
            a_new = (1/3)*(2 - CT);
        else
            a_new = 0.5*(1 - sqrt(max(0, 1 - CT)));
        end
        if abs(a_new - a_up) < 1e-5, break; end
        a_up = a_up + 0.3*(a_new - a_up);
    end

    % Downwind half sees the wake of the upwind half
    Vdown = Vinf * (1 - 2*a_up);
    a_dn = 0;
    for iter = 1:50
        CT = 0;
        for i = n_theta/2+1:n_theta
            th = theta(i);
            W_u = TSR + sin(th);
            W_v = cos(th);
            W   = hypot(W_u, W_v);
            alpha = atan2(W_v, W_u);
            Cl = interp1(alpha_data, Cl_data, rad2deg(alpha), 'linear', 0);
            Cd = interp1(alpha_data, Cd_data, rad2deg(alpha), 'linear', 0.02);
            Cn = Cl*cos(alpha) + Cd*sin(alpha);
            CT = CT + (sigma_local/(4*pi)) * Cn * W^2 * dtheta;
        end
        if CT > 0.96
            a_new = (1/3)*(2 - CT);
        else
            a_new = 0.5*(1 - sqrt(max(0, 1 - CT)));
        end
        if abs(a_new - a_dn) < 1e-5, break; end
        a_dn = a_dn + 0.3*(a_new - a_dn);
    end

    % Torque around the full rotation
    torque_blade = zeros(1, n_theta);
    for i = 1:n_theta
        th = theta(i);
        if th <= pi
            V_loc = Vinf * (1 - a_up);
        else
            V_loc = Vdown * (1 - a_dn);
        end
        W_u = TSR + sin(th);
        W_v = cos(th);
        alpha = atan2(W_v, W_u);
        Cl = interp1(alpha_data, Cl_data, rad2deg(alpha), 'linear', 0);
        Cd = interp1(alpha_data, Cd_data, rad2deg(alpha), 'linear', 0.02);
        Ct_local = Cl*sin(alpha) - Cd*cos(alpha);
        dF = 0.5 * rho * V_loc^2 * c * H * Ct_local;
        torque_blade(i) = dF * R;
    end

    Cp = trapz(theta, torque_blade) / (0.5 * rho * Vinf^3 * 2*R*H) * N / (2*pi) * TSR;
end