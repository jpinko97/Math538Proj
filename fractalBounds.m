%% Fractal basin boundaries for coupled logistic map
clear; clc;

r1 = 3.8;
r2 = 3.8;
D = 0.155;

Ngrid = 100;          % 100 x 100 grid
n_transient = 2000;   % iterations to approach attractor
n_post      = 64;     % iterations used for classification

% tolerances used in classification
tol_const = 1e-5;     % constancy within each phase for period-2
tol_swap  = 1e-3;     % swap x1<->x2 condition for out-of-phase
tol_total = 1e-3;     % constancy of total population for period-2
tol_per4  = 1e-4;     % error allowed for period-4 orbit


%% -------------------- Figure A: full [0,1]x[0,1] -----------------------

x1_vals_A = linspace(0, 1, Ngrid);
x2_vals_A = linspace(0, 1, Ngrid);

basinA = zeros(Ngrid, Ngrid);   % 1 = P2, 2 = P4, 3 = other

for i = 1:Ngrid
    for j = 1:Ngrid

        x1 = x1_vals_A(i);
        x2 = x2_vals_A(j);

        % transient
        for k = 1:n_transient
            [x1, x2] = coupled_step(x1, x2, r1, r2, D);
        end

        % collect series
        x1_series = zeros(1, n_post);
        x2_series = zeros(1, n_post);
        for k = 1:n_post
            [x1, x2] = coupled_step(x1, x2, r1, r2, D);
            x1_series(k) = x1;
            x2_series(k) = x2;
        end

        % classify
        if is_period2_out_of_phase(x1_series, x2_series, ...
                                   tol_const, tol_swap, tol_total)
            basinA(j,i) = 1;            % white
        elseif is_period4_orbit(x1_series, x2_series, tol_per4)
            basinA(j,i) = 2;            % black
        else
            basinA(j,i) = 3;            % red (other / chaos)
        end
    end
end

% plot Figure A
figure;
imagesc(x1_vals_A, x2_vals_A, basinA);
set(gca,'YDir','normal');
axis square;
xlabel('Population 1 initial condition');
ylabel('Population 2 initial condition');
title('Basin of attraction, r=3.8, D = 0.155');

% colormap: 1=white (P2), 2=black (P4), 3=red (other)
cmapA = [1 1 1;    % white
         0 0 0;    % black
         1 0 0];   % red
colormap(cmapA);
colorbar off;

% ---- Legend for colors ----
hold on;

% Create invisible (NaN) scatter points with your colors
hP2 = scatter(NaN, NaN, 80, [1 1 1], 'filled', 'MarkerEdgeColor','k');  % white
hP4 = scatter(NaN, NaN, 80, [0 0 0], 'filled');                         % red
hCh = scatter(NaN, NaN, 80, [1 0 0], 'filled');                         % black

legend([hP2, hP4, hCh], ...
       {'Period 2 stable', 'Period 4 stable', 'Chaotic/other'}, ...
       'Location', 'northwest');


%% -------------------- Figure B: zoom [0,1e-3]x[0,1e-3] -----------------

x1_vals_B = linspace(0, 1e-3, Ngrid);
x2_vals_B = linspace(0, 1e-3, Ngrid);

basinB = zeros(Ngrid, Ngrid);

for i = 1:Ngrid
    for j = 1:Ngrid

        x1 = x1_vals_B(i);
        x2 = x2_vals_B(j);

        % transient
        for k = 1:n_transient
            [x1, x2] = coupled_step(x1, x2, r1, r2, D);
        end

        % collect series
        x1_series = zeros(1, n_post);
        x2_series = zeros(1, n_post);
        for k = 1:n_post
            [x1, x2] = coupled_step(x1, x2, r1, r2, D);
            x1_series(k) = x1;
            x2_series(k) = x2;
        end

        % classify
        if is_period2_out_of_phase(x1_series, x2_series, ...
                                   tol_const, tol_swap, tol_total)
            basinB(j,i) = 1;
        elseif is_period4_orbit(x1_series, x2_series, tol_per4)
            basinB(j,i) = 2;
        else
            basinB(j,i) = 3;
        end
    end
end

% plot Figure B
figure;
imagesc(x1_vals_B, x2_vals_B, basinB);
set(gca,'YDir','normal');
axis square;
xlabel('Population 1 initial condition');
ylabel('Population 2 initial condition');
title('Basin of attraction zoomed, r_1=3.9, r_2=3.7, D = 0.13');

% ---- Legend for colors ----
hold on;

% Create invisible (NaN) scatter points with your colors
hP2 = scatter(NaN, NaN, 80, [1 1 1], 'filled', 'MarkerEdgeColor','k');  % white
hP4 = scatter(NaN, NaN, 80, [0 0 0], 'filled');                         % red
hCh = scatter(NaN, NaN, 80, [1 0 0], 'filled');                         % black

legend([hP2, hP4, hCh], ...
       {'Period 2 stable', 'Period 4 stable', 'Chaotic/other'}, ...
       'Location', 'northwest');

colormap(cmapA);
colorbar off;



%% =============== Helper functions ======================================

function [x1_next, x2_next] = coupled_step(x1, x2, r1, r2, D)
    % One step of the coupled logistic maps
    u1 = r1 * x1 * (1 - x1);
    u2 = r2 * x2 * (1 - x2);
    x1_next = (1 - D) * u1 + D * u2;
    x2_next = (1 - D) * u2 + D * u1;
end


function isP2 = is_period2_out_of_phase(x1s, x2s, tol_const, tol_swap, tol_total)
    % Detect out-of-phase period-2 solution with constant total population

    n = numel(x1s);
    odd  = 1:2:n;
    even = 2:2:n;
    if isempty(even)
        isP2 = false;
        return;
    end

    x1o = x1s(odd);  x1e = x1s(even);
    x2o = x2s(odd);  x2e = x2s(even);

    % each phase should be nearly constant
    if std(x1o) > tol_const || std(x1e) > tol_const || ...
       std(x2o) > tol_const || std(x2e) > tol_const
        isP2 = false; return;
    end

    m1o = mean(x1o); m1e = mean(x1e);
    m2o = mean(x2o); m2e = mean(x2e);

    % out-of-phase: x1_even ≈ x2_odd and x1_odd ≈ x2_even
    if abs(m1e - m2o) > tol_swap || abs(m1o - m2e) > tol_swap
        isP2 = false; return;
    end

    % exclude in-phase behavior
    if abs(m1e - m2e) < tol_swap && abs(m1o - m2o) < tol_swap
        isP2 = false; return;
    end

    % total population almost constant between the two phases
    if abs((m1e + m2e) - (m1o + m2o)) > tol_total
        isP2 = false; return;
    end

    isP2 = true;
end


function isP4 = is_period4_orbit(x1s, x2s, tol_per4)
    % Detect (approximate) period-4 orbit (not period-2)
    n = numel(x1s);
    if n < 8
        isP4 = false; return;
    end

    % Compare state at k and k+4
    dx1 = x1s(1:(n-4)) - x1s(5:n);
    dx2 = x2s(1:(n-4)) - x2s(5:n);
    err4 = mean(sqrt(dx1.^2 + dx2.^2));

    isP4 = (err4 < tol_per4);
end
