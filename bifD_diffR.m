%% Bifurcation diagram for OUT-OF-PHASE solutions (signed difference) FOR D

%bifurcation diagram for b with different r

clear; clc;

r1 = 3.7;                % diffusion parameter
r2 = 3.9;
d_min = 0;
d_max = 0.5;
n_d   = 800;            % number of r values

n_transient = 500;      % steps to discard (transient)
n_keep      = 100;      % steps to sample after transient
eps_diff    = 1e-3;     % threshold to define out-of-phase behavior

d_values = linspace(d_min, d_max, n_d);

% Preallocate (upper bound)
D_plot     = zeros(1, n_d * n_keep);
XDIFF_plot = zeros(1, n_d * n_keep);
idx = 1;

for id = 1:n_d
    d = d_values(id);

    % random initial conditions
    x1 = rand;
    x2 = rand;

    % transient
    for k = 1:n_transient
        x1_tilde = r1 * x1 * (1 - x1);
        x2_tilde = r2 * x2 * (1 - x2);

        x1 = x1_tilde + d * (x2_tilde - x1_tilde);
        x2 = x2_tilde + d * (x1_tilde - x2_tilde);
    end

    % after transient, record
    for k = 1:n_keep
        x1_tilde = r1 * x1 * (1 - x1);
        x2_tilde = r2 * x2 * (1 - x2);

        x1 = x1_tilde + d * (x2_tilde - x1_tilde);
        x2 = x2_tilde + d * (x1_tilde - x2_tilde);

        xdiff = x1 - x2;   % signed difference

        % keep only out-of-phase states
        if abs(xdiff) > eps_diff
            D_plot(idx)     = d;
            XDIFF_plot(idx) = xdiff;
            idx = idx + 1;
        end
    end
end

% trim
D_plot     = D_plot(1:idx-1);
XDIFF_plot = XDIFF_plot(1:idx-1);

% Plot bifurcation diagram
figure;
plot(D_plot, XDIFF_plot, '.', 'MarkerSize', 3);
xlabel('D');
ylabel('Difference in Populations');
title('Out-of-phase bifurcation diagram, r_1=3.9 ,r_2=3.7');

ylim([-1, 1]);      % <-- requested y-limits
xlim([d_min d_max]);
grid on;
