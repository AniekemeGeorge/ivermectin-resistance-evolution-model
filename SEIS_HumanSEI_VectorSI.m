clear; close all; clc;

% Simulation time (days)
simTime = 365; % 1 year

% Example parameter values (Can be adjusted as needed)
mu_h     = 1/(70*365); % human natural mortality (~70-year life)
sigma_h  = 1/10;     % incubation ~10 days
gamma_h  = 1/20;     % recovery untreated with ivermectin and without ACT (natural recovery) ~20 days
gamma_T  = 1/20;     % recovery rate of those treated with ivermectin ~20 days without ACT
m_rate   = 0.0;     % treatment rate per day
w_rate   = 1/3;    % waning treatment ~3 days

mu_v     = 1/10;     % mosquito baseline mortality (days^-1)
a        = 0.4;      % biting rate per mosquito per day
a_u      = 0.4;      %biting rate susceptible unexposed to Ivermectin
a_e      = 0.4;
beta_vh  = 0.2;      % prob human infected per infectious bite
phi      = 1;      % prob bite on treated human exposes susceptible genotype
phi_r    = 1;      % prob bite on treated human exposes resistant genotype

rho_v    = 0;      % rate leaving 'exposed' back to unexposed
d_s      = 2.0;     % extra mortality exposed susceptible
d_r      = 0.01;     % extra mortality exposed resistant

rs       = 0.2;      % per-capita birth rate susceptible genotype
rr       = rs*0.8;     % per-capita birth rate resistant genotype
K        = 100000;   % mosquito carrying capacity - 10x human population

beta_hv_u = 0.2;     % prob unexposed mosquito infected per bite on infectious human
beta_hv_e = 0.2;    % prob exposed mosquito infected per bite on infectious human
beta_hv_r = 0.2;    % prob resistant mosquito infected per bite on infectious human

% Call the model
[T,Y] = HumanSEI_VectorSI (mu_h, sigma_h, gamma_h, gamma_T, m_rate, w_rate, ...
     mu_v, a, a_u,a_e, beta_vh, phi, phi_r, ...
    rho_v, d_s, d_r, rs, rr, K, beta_hv_u, beta_hv_e, simTime);

% Y columns correspond to:
% 1 S_h_U  2 E_h_U  3 I_h_U  4 S_h_T  5 E_h_T  6 I_h_T
% 7 S_s_u  8 I_s_u  9 S_s_e 10 I_s_e 11 S_r_u 12 I_r_u 13 S_r_e 14 I_r_e

% Compute totals for plotting
I_h_total = Y(:,3) + Y(:,6);
S_h_total = Y(:,1) + Y(:,4);
E_h_total = Y(:,2) + Y(:,5);
N_h = Y(:,1) + Y(:,2) + Y(:,3) + Y(:,4)+ Y(:,5)+ Y(:,6);
Nv_total  = sum(Y(:,7:14), 2);
Iv_total = Y(:,8) + Y(:,10) + Y(:,12) + Y(:,14); % infected mosquitoes (all categories)

% Plot human infectious over time
figure;
plot(T, I_h_total, 'LineWidth', 1.6);
xlabel('Time (days)');
ylabel('Number of infectious humans');
title('Human Infectious (I_h_U + I_h_T)');
grid on;

% Stacked human compartments
figure;
plot(T, S_h_total, '-', T, E_h_total, '--', T, I_h_total, '-.', 'LineWidth', 1.4);
legend('Susceptible (humans)','Exposed (humans)','Infectious (humans)');
xlabel('Time (days)');
ylabel('Counts');
title('Human compartments (total treated+untreated)');
grid on;

% Quick printout of final-day values
fprintf('--- Final day (t=%.1f days) ---\n', T(end));
fprintf('Total humans = %.1f\n', sum(Y(end,1:6)));
fprintf('Infectious humans = %.2f\n', I_h_total(end));
fprintf('Total vectors = %.1f\n', Nv_total(end));
fprintf('Infectious vectors = %.2f\n', Iv_total(end));
