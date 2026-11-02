% Parameters
r1 = 3.6;
r2 = 3.6;
D  = 0.1;       % 0 ≤ D ≤ 0.5
T  = 50;       % number of iterations
x1_0 = 0.9;
x2_0 = 0.5;

%% for figure 2, r1=r2=3.8 and D=0.15. 
    % for A, x10 = 0.5 and x20 = 0.9
    % for B, x10 = 0.3 and x20 = 0.5

% Run simulation
[t, x1, x2] = coupled_logistic(x1_0, x2_0, r1, r2, D, T);

% Plot
figure;
plot(t, x1, 'LineWidth', 1.5); hold on;
plot(t, x2, 'LineWidth', 1.5); hold on;
%plot(t, x1+x2, 'LineWidth', 1.5)
xlabel('Generation');
ylabel('Number of individuals');
legend('Population 1', 'Population 2', 'Total Population');
title('Out of Phase Coupled Logistic Trajectories,','x_1(0) = 0.9, x_2(0) = 0.5, r_1 = 3.7, r_2 = 3.9, D = 0.05');
grid on;


function [t, x1, x2] = coupled_logistic(x1_0, x2_0, r1, r2, D, T)
% COUPLED_LOGISTIC  Simulates the system:
%
%   x̃1(t) = r1 * x1(t) * (1 - x1(t))
%   x̃2(t) = r2 * x2(t) * (1 - x2(t))
%
%   x1(t+1) = x̃1(t) + D * (x̃2(t) - x̃1(t))
%   x2(t+1) = x̃2(t) + D * (x̃1(t) - x̃2(t))
%
% INPUTS:
%   x1_0, x2_0 : initial values in [0,1]
%   r1, r2     : logistic map parameters
%   D          : diffusion parameter (0 ≤ D ≤ 0.5)
%   T          : number of time steps
%
% OUTPUTS:
%   t          : time vector
%   x1, x2     : trajectories of x1(t) and x2(t)

    % Preallocate
    t  = 0:T;
    x1 = zeros(1, T+1);
    x2 = zeros(1, T+1);

    % Initial conditions
    x1(1) = x1_0;
    x2(1) = x2_0;

    for k = 1:T
        % Compute intermediate values x̃i(t)
        x1_tilde = r1 * x1(k) * (1 - x1(k));
        x2_tilde = r2 * x2(k) * (1 - x2(k));

        % Update coupled dynamics
        x1(k+1) = x1_tilde + D * (x2_tilde - x1_tilde);
        x2(k+1) = x2_tilde + D * (x1_tilde - x2_tilde);
    end
end

