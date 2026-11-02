%% Bifurcation diagram for IN-PHASE solution (logistic map)
clear; clc;

% r range
r_min = 3.0;
r_max = 4.0;
n_r   = 800;              % resolution in r

% iteration settings
n_transient = 500;        % throw away
n_keep      = 200;        % plot points

r_values = linspace(r_min, r_max, n_r);

R_plot = zeros(1, n_r * n_keep);
X_plot = zeros(1, n_r * n_keep);
idx = 1;

for ir = 1:n_r
    r = r_values(ir);

    % in-phase populations => single logistic variable x
    x = 0.5;              % initial condition in (0,1)

    % transient
    for k = 1:n_transient
        x = r * x * (1 - x);
    end

    % collect points for the diagram
    for k = 1:n_keep
        x = r * x * (1 - x);
        R_plot(idx) = r;
        X_plot(idx) = 2 * x;   % total population x1 + x2, with x1 = x2 = x
        idx = idx + 1;
    end
end

% trim
R_plot = R_plot(1:idx-1);
X_plot = X_plot(1:idx-1);

% plot
figure;
plot(R_plot, X_plot, '.', 'MarkerSize', 1);
xlabel('r');
ylabel('Total Population');
title('In-Phase Bifurcation Diagram, D = 0.1');
xlim([r_min r_max]);
ylim([0 2]);   % since 0 <= x <= 1 and x1 + x2 = 2x
grid on;
