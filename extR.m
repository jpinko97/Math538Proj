%% Total-population return map: y(t) vs y(t+1), and fit r_eff
clear; clc;

% --- Parameters (edit) ---
r1 = 3.7;
r2 = 3.9;
D  = 0.25;

T_total = 8000;
T_burn  = 2000;

% Initial conditions
x1 = zeros(1, T_total);
x2 = zeros(1, T_total);
x1(1) = 0.3;
x2(1) = 0.7;

% --- Simulate coupled maps ---
for t = 1:T_total-1
    u1 = r1 * x1(t) * (1 - x1(t));
    u2 = r2 * x2(t) * (1 - x2(t));
    x1(t+1) = (1 - D) * u1 + D * u2;
    x2(t+1) = (1 - D) * u2 + D * u1;
end

% --- Build totals y(t) and y(t+1) after transient ---
t_idx = (T_burn):(T_total-1);
y_t   = x1(t_idx)   + x2(t_idx);
y_tp1 = x1(t_idx+1) + x2(t_idx+1);

%% 1) Plot y(t) vs y(t+1)
figure;
plot(y_t, y_tp1, '.', 'MarkerSize', 3);
xlabel('y(t) = x_1(t) + x_2(t)');
ylabel('y(t+1) = x_1(t+1) + x_2(t+1)');
title(sprintf('Total return map, r_1=%.2f, r_2=%.2f, D=%.2f', r1, r2, D));
xlim([0 2]); ylim([0 2]);
grid on; box on;

%% 2) Fit a single effective r for the total map
% Model: y(t+1) ≈ r_eff * y(t) * (1 - y(t)/2)
phi = y_t .* (1 - y_t/2);

% (Optional) ignore points where phi is tiny (near 0 or 2) to reduce noise
mask = abs(phi) > 1e-4;

% Least-squares best fit for scalar r_eff
r_eff = (phi(mask) * y_tp1(mask)') / (phi(mask) * phi(mask)');

fprintf('Estimated r_eff (LS fit to y(t+1)=r_eff*y(t)*(1-y(t)/2)): %.6f\n', r_eff);

% Overlay fitted curve on the scatter (for visual check)
hold on;
ygrid = linspace(0, 2, 400);
yfit  = r_eff * ygrid .* (1 - ygrid/2);
plot(ygrid, yfit, 'k-', 'LineWidth', 1.5);
legend('data', sprintf('fit: r_{eff}=%.4f', r_eff), 'Location', 'best');
