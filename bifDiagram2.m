%% Bifurcation diagram for OUT-OF-PHASE solutions (signed difference)
clear; clc;

D = 0.15;                % diffusion parameter
r_min = 3.0;
r_max = 4.0;
n_r   = 800;            % number of r values

n_transient = 500;      % steps to discard (transient)
n_keep      = 100;      % steps to sample after transient
eps_diff    = 1e-3;     % threshold to define out-of-phase behavior

r_values = linspace(r_min, r_max, n_r);

% Preallocate (upper bound)
R_plot     = zeros(1, n_r * n_keep);
XDIFF_plot = zeros(1, n_r * n_keep);
idx = 1;

for ir = 1:n_r
    r = r_values(ir);

    % random initial conditions
    x1 = rand;
    x2 = rand;

    % transient
    for k = 1:n_transient
        x1_tilde = r * x1 * (1 - x1);
        x2_tilde = r * x2 * (1 - x2);

        x1 = x1_tilde + D * (x2_tilde - x1_tilde);
        x2 = x2_tilde + D * (x1_tilde - x2_tilde);
    end

    % after transient, record
    for k = 1:n_keep
        x1_tilde = r * x1 * (1 - x1);
        x2_tilde = r * x2 * (1 - x2);

        x1 = x1_tilde + D * (x2_tilde - x1_tilde);
        x2 = x2_tilde + D * (x1_tilde - x2_tilde);

        xdiff = x1 - x2;   % signed difference

        % keep only out-of-phase states
        if abs(xdiff) > eps_diff
            R_plot(idx)     = r;
            XDIFF_plot(idx) = xdiff;
            idx = idx + 1;
        end
    end
end

% trim
R_plot     = R_plot(1:idx-1);
XDIFF_plot = XDIFF_plot(1:idx-1);

% Plot bifurcation diagram
figure;
plot(R_plot, XDIFF_plot, '.', 'MarkerSize', 8);
xlabel('r_1 = r_2');
ylabel('Difference in Populations');
title('Out-of-phase bifurcation diagram, D = 0.15');

ylim([-1, 1]);      % <-- requested y-limits
xlim([3.77 3.83]);
grid on;
