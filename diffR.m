%% Coupled logistic model – phase plots like Hastings Fig.
clear; clc;

% Parameters
r1 = 3.7;
r2 = 3.9;
D  = 0.13;

T_total   = 6000;   % total iterations
T_burn    = 1000;   % discard as transient

% Preallocate
x1 = zeros(1, T_total);
x2 = zeros(1, T_total);

% Initial conditions (can be changed)
x1(1) = 0.3;
x2(1) = 0.7;

% Iterate the coupled logistic system
for t = 1:T_total-1
    % local logistic growth
    u1 = r1 * x1(t) * (1 - x1(t));
    u2 = r2 * x2(t) * (1 - x2(t));

    % diffusive coupling
    x1(t+1) = (1 - D) * u1 + D * u2;
    x2(t+1) = (1 - D) * u2 + D * u1;
end

% Indices we actually plot (after transient and with t-1 defined)
idx = (T_burn+2):T_total;

x1_t   = x1(idx);
x1_tm1 = x1(idx-1);
x2_t   = x2(idx);
x2_tm1 = x2(idx-1);

y_t    = x1_t + x2_t;
y_tm1  = x1_tm1 + x2_tm1;

%% ---- Make 4-panel figure ----
figure;

% 1) x2(t) vs x1(t)
subplot(2,2,1);
plot(x1_t, x2_t, '.', 'MarkerSize', 3);
xlabel('x_1(t)');
ylabel('x_2(t)');
title('Population 1 vs Population 2');
xlim([0 1]); ylim([0 1]);

% 2) total population y(t) vs y(t-1)
subplot(2,2,2);
plot(y_tm1, y_t, '.', 'MarkerSize', 3);
xlabel('x_{tot}(t-1) = x_1(t-1) + x_2(t-1)');
ylabel('x_{tot}(t) = x_1(t) + x_2(t)');
title('Total population');
xlim([0 2]); ylim([0 2]);

% 3) population 1: x1(t) vs x1(t-1)
subplot(2,2,3);
plot(x1_tm1, x1_t, '.', 'MarkerSize', 3);
xlabel('x_1(t-1)');
ylabel('x_1(t)');
title('Population 1');
xlim([0 1]); ylim([0 1]);

% 4) population 2: x2(t) vs x2(t-1)
subplot(2,2,4);
plot(x2_tm1, x2_t, '.', 'MarkerSize', 3);
xlabel('x_2(t-1)');
ylabel('x_2(t)');
title('Population 2');
xlim([0 1]); ylim([0 1]);

sgtitle('Coupled logistic maps, r_1=3.7, r_2=3.9, D=0.13');

figure;
hold on;
plot(x1_tm1, x1_t, '.', 'MarkerSize', 3, 'Color', [0 0.4470 0.7410]); % blue
plot(x2_tm1, x2_t, '.', 'MarkerSize', 3, 'Color', [0.8500 0.3250 0.0980]); % red
xlabel('x(t-1)');
ylabel('x(t)');
title('Overlay of return maps: Pop 1 (blue), Pop 2 (red)');
xlim([0 1]); ylim([0 1]);
legend('Population 1','Population 2','Location','best');
grid on;
