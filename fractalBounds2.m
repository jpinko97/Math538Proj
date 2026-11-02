
% Classic Hastings-like parameters (same r):
fractal_basin_two_r(3.8, 3.8, 0.1);

% Different r in each patch:
%fractal_basin_two_r(3.4, 3.9, 0.15);


function fractal_basin_two_r(r1, r2, D)
% FRACTAL_BASIN_TWO_R  Basins of attraction for coupled logistic maps
% with possibly different r values in the two patches.
%
%   fractal_basin_two_r(r1, r2, D)
%
%   r1, r2 : logistic parameters for patch 1 and 2
%   D      : diffusion parameter

    % ----------------- Parameters -----------------
    Ngrid       = 100;   % grid resolution (Ngrid x Ngrid)
    n_transient = 2000;  % iterations to approach attractor
    n_post      = 64;    % iterations used for classification

    % tolerances
    tol_const = 1e-5;    % constancy within each phase for period-2
    tol_swap  = 1e-3;    % swap condition for out-of-phase P2
    tol_total = 1e-3;    % constancy of total population for P2
    tol_per4  = 1e-3;    % tolerance for P4 detection

    % colormap: 1=white (P2), 2=black (P4), 3=red (other/chaotic)
    cmap = [1 1 1;    % white
            0 0 0;    % black
            1 0 0];   % red

    %% ---------------- Figure A: [0,1] x [0,1] ----------------
    x1_vals_A = linspace(0, 1, Ngrid);
    x2_vals_A = linspace(0, 1, Ngrid);

    basinA = zeros(Ngrid, Ngrid);   % 1=P2, 2=P4, 3=other

    for i = 1:Ngrid
        for j = 1:Ngrid

            x1 = x1_vals_A(i);
            x2 = x2_vals_A(j);

            % transient
            for k = 1:n_transient
                [x1, x2] = coupled_step_two_r(x1, x2, r1, r2, D);
            end

            % collect series
            x1_series = zeros(1, n_post);
            x2_series = zeros(1, n_post);
            for k = 1:n_post
                [x1, x2] = coupled_step_two_r(x1, x2, r1, r2, D);
                x1_series(k) = x1;
                x2_series(k) = x2;
            end

            % classify
            if is_period2_out_of_phase(x1_series, x2_series, ...
                                       tol_const, tol_swap, tol_total)
                basinA(j,i) = 1;   % white: P2
            elseif is_period4_orbit(x1_series, x2_series, tol_per4)
                basinA(j,i) = 2;   % black: P4
            else
                basinA(j,i) = 3;   % red: other/chaotic
            end
        end
    end

    % plot Figure A
    figure;
    imagesc(x1_vals_A, x2_vals_A, basinA);
    set(gca,'YDir','normal');
    axis square;
    colormap(cmap);
    colorbar off;
    xlabel('x_1(0)');
    ylabel('x_2(0)');
    title(sprintf('Basins, r_1 = %.3f, r_2 = %.3f, D = %.3f (full)', r1, r2, D));


    %% ------------- Figure B: [0,1e-3] x [0,1e-3] --------------
    x1_vals_B = linspace(0, 1e-3, Ngrid);
    x2_vals_B = linspace(0, 1e-3, Ngrid);

    basinB = zeros(Ngrid, Ngrid);   % 1=P2, 2=P4, 3=other

    for i = 1:Ngrid
        for j = 1:Ngrid

            x1 = x1_vals_B(i);
            x2 = x2_vals_B(j);

            % transient
            for k = 1:n_transient
                [x1, x2] = coupled_step_two_r(x1, x2, r1, r2, D);
            end

            % collect series
            x1_series = zeros(1, n_post);
            x2_series = zeros(1, n_post);
            for k = 1:n_post
                [x1, x2] = coupled_step_two_r(x1, x2, r1, r2, D);
                x1_series(k) = x1;
                x2_series(k) = x2;
            end

            % classify
            if is_period2_out_of_phase(x1_series, x2_series, ...
                                       tol_const, tol_swap, tol_total)
                basinB(j,i) = 1;   % white: P2
            elseif is_period4_orbit(x1_series, x2_series, tol_per4)
                basinB(j,i) = 2;   % black: P4
            else
                basinB(j,i) = 3;   % red: other/chaotic
            end
        end
    end

    % plot Figure B
    figure;
    imagesc(x1_vals_B, x2_vals_B, basinB);
    set(gca,'YDir','normal');
    axis square;
    colormap(cmap);
    colorbar off;
    xlabel('x_1(0)');
    ylabel('x_2(0)');
    title(sprintf('Basins, r_1 = %.3f, r_2 = %.3f, D = %.3f (zoom)', r1, r2, D));
end


%% ==================== Helper functions ====================

function [x1_next, x2_next] = coupled_step_two_r(x1, x2, r1, r2, D)
    % One step of coupled logistic maps with different r values
    u1 = r1 * x1 * (1 - x1);
    u2 = r2 * x2 * (1 - x2);
    x1_next = (1 - D) * u1 + D * u2;
    x2_next = (1 - D) * u2 + D * u1;
end


function isP2 = is_period2_out_of_phase(x1s, x2s, ...
                                        tol_const, tol_swap, tol_total)
    % Detect out-of-phase period-2 orbit with constant total population.

    n = numel(x1s);
    odd  = 1:2:n;
    even = 2:2:n;

    if isempty(even)
        isP2 = false; return;
    end

    x1o = x1s(odd);  x1e = x1s(even);
    x2o = x2s(odd);  x2e = x2s(even);

    % each phase nearly constant ⇒ period-2
    if std(x1o) > tol_const || std(x1e) > tol_const || ...
       std(x2o) > tol_const || std(x2e) > tol_const
        isP2 = false; return;
    end

    m1o = mean(x1o); m1e = mean(x1e);
    m2o = mean(x2o); m2e = mean(x2e);

    % out-of-phase: x1_even ≈ x2_odd AND x1_odd ≈ x2_even
    if abs(m1e - m2o) > tol_swap || abs(m1o - m2e) > tol_swap
        isP2 = false; return;
    end

    % exclude in-phase period-2 (x1 ≈ x2)
    if abs(m1e - m2e) < tol_swap && abs(m1o - m2o) < tol_swap
        isP2 = false; return;
    end

    % total population nearly constant
    tot_even = m1e + m2e;
    tot_odd  = m1o + m2o;
    if abs(tot_even - tot_odd) > tol_total
        isP2 = false; return;
    end

    isP2 = true;
end


function isP4 = is_period4_orbit(x1s, x2s, tol4)
    % Robust detection of (approximate) period-4 orbit, not period-2.

    n = numel(x1s);
    if n < 10
        isP4 = false; 
        return;
    end

    X = [x1s(:), x2s(:)];

    % 4-step repetition: X(t+4) ≈ X(t)
    d4 = vecnorm(X(1:n-4,:) - X(5:n,:), 2, 2);
    err4 = max(d4);

    % 2-step mismatch must be larger (to rule out period-2)
    d2 = vecnorm(X(1:n-2,:) - X(3:n,:), 2, 2);
    err2 = max(d2);

    tol2 = 5*tol4;    % threshold separation

    isP4 = (err4 < tol4) && (err2 > tol2);
end
