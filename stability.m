%% Stability diagram: period-2 vs period-4 vs other
clear; clc;

%% Parameter ranges (edit if you like)
r_min = 3.0; r_max = 4.0;
D_min = 0.0; D_max = 0.25;

Nr = 201;                    % # r values
ND = 201;                    % # D values
r_vec = linspace(r_min, r_max, Nr);
D_vec = linspace(D_min, D_max, ND);

%% Simulation settings
n_transient = 2000;          % burn-in
n_post      = 80;            % samples for classification

% tolerances for period-2 (out-of-phase, ~constant total)
tol_const = 1e-5;
tol_swap  = 1e-3;
tol_total = 1e-3;

% tolerance for period-4 detection
tol_p4 = 1e-3;               % try 1e-3–2e-3 if needed

%% Codes: 1 = P2, 2 = P4, 0 = other/chaotic
region = zeros(ND, Nr);      % rows = D, cols = r

%% Sweep over (r,D)
for ir = 1:Nr
    r = r_vec(ir);

    for id = 1:ND
        D = D_vec(id);

        % --- simulate from an out-of-phase IC ---
        x1 = 0.5;
        x2 = 0.9;

        % burn transient
        for k = 1:n_transient
            [x1, x2] = coupled_step(x1, x2, r, D);
        end

        % collect series
        x1s = zeros(1, n_post);
        x2s = zeros(1, n_post);
        for k = 1:n_post
            [x1, x2] = coupled_step(x1, x2, r, D);
            x1s(k) = x1;
            x2s(k) = x2;
        end

        % --- classify dynamics ---
        if is_period2_out_of_phase(x1s, x2s, ...
                                   tol_const, tol_swap, tol_total)
            region(id, ir) = 1;          % period-2 (white)
        elseif is_period4(x1s, x2s, tol_p4)
            region(id, ir) = 2;          % period-4 (red)
        else
            region(id, ir) = 0;          % other / chaotic (black)
        end
    end
end

%% Plot stability diagram
figure;
imagesc(r_vec, D_vec, region);
set(gca,'YDir','normal');
axis tight;
xlabel('r_1 = r_2 = r');
ylabel('D');
title('Stability of Period 2 solutions x_1(t)=0.5, x_2(t)=0.9,', 'found at D=0.15, r=3.8');

% colormap rows correspond to values 0,1,2 (region) via +1 offset
cmap = [1 0 0;    % region==0  -> red (other/chaotic)
        1 1 1;    % region==1  -> white (P2)
        0 0 0];   % region==2  -> black   (P4)
colormap(cmap);
caxis([0 2]);
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


box on;


%% ================= Helper functions =================

function [x1_next, x2_next] = coupled_step(x1, x2, r, D)
    % One time step of the coupled logistic maps with r1=r2=r.
    u1 = r * x1 * (1 - x1);
    u2 = r * x2 * (1 - x2);
    x1_next = (1 - D) * u1 + D * u2;
    x2_next = (1 - D) * u2 + D * u1;
end

function isP2 = is_period2_out_of_phase(x1s, x2s, ...
                                        tol_const, tol_swap, tol_total)
    % Detect out-of-phase period-2 orbit with (approximately) constant total.

    n = numel(x1s);
    odd  = 1:2:n;
    even = 2:2:n;
    if isempty(even)
        isP2 = false; return;
    end

    x1o = x1s(odd);  x1e = x1s(even);
    x2o = x2s(odd);  x2e = x2s(even);

    % each phase nearly constant
    if std(x1o) > tol_const || std(x1e) > tol_const || ...
       std(x2o) > tol_const || std(x2e) > tol_const
        isP2 = false; return;
    end

    m1o = mean(x1o); m1e = mean(x1e);
    m2o = mean(x2o); m2e = mean(x2e);

    % out-of-phase: x1_even ≈ x2_odd, x1_odd ≈ x2_even
    if abs(m1e - m2o) > tol_swap || abs(m1o - m2e) > tol_swap
        isP2 = false; return;
    end

    % exclude in-phase period-2 (x1≈x2 each step)
    if abs(m1e - m2e) < tol_swap && abs(m1o - m2o) < tol_swap
        isP2 = false; return;
    end

    % total population nearly constant between phases
    tot_even = m1e + m2e;
    tot_odd  = m1o + m2o;
    if abs(tot_even - tot_odd) > tol_total
        isP2 = false; return;
    end

    isP2 = true;
end

function isP4 = is_period4(x1s, x2s, tol)
    % Detect period-4 orbit (not period-2), more generally (in- or out-of-phase).
    % Uses 4-step and 2-step recurrences.

    n = numel(x1s);
    if n < 12
        isP4 = false; return;
    end

    X = [x1s(:), x2s(:)];

    % 4-step match: X(t+4) ≈ X(t)
    d4 = vecnorm(X(1:n-4,:) - X(5:n,:), 2, 2);
    err4 = max(d4);

    % 2-step mismatch: X(t+2) should NOT line up too well (to rule out P2)
    d2 = vecnorm(X(1:n-2,:) - X(3:n,:), 2, 2);
    err2 = max(d2);

    tol2 = 5*tol;
    isP4 = (err4 < tol) && (err2 > tol2);
end
