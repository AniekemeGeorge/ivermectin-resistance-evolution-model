clear; close all; clc;

simTime = 365;    % Total simulation time

% -------- FIXED PARAMETERS --------
mu_h      = 1/(70*365);
sigma_h   = 1/10;
gamma_T   = 1/20;
gamma_h   = 1/20;
w_rate    = 1/3;
mu_v      = 1/10;
a_u       = 0.4;
a_e       = 0.4;
beta_vh   = 0.2;
K         = 100000;
beta_hv_u = 0.2;
beta_hv_e = 0.2;
phi_r     = 1;
phi       = 1;
rs        = 0.2;
d_r       = 0.01;
rho_v     = 0;
m_rate    = 0;

% -------- SAMPLING SETTINGS --------
nSamples = 3*24;
rng(1);   % Reproducibility

% Results:
% columns = [a, d_s, rr, prop_resistant]
results = zeros(nSamples, 7);  %added extra column to keep track of malaria Ro and mosquito numbers before ivermectin

d_s = 2;      % susceptibles die  0.5 days days after exposure


k=1;

for rr_adjust=0.75:0.05:0.95

    rr  = rr_adjust* rs;       % high rr (low cost of resistance)


    for a = 0.08:0.04:0.9
  
        %bringing in code from pulsed model

        init_pop(1)  = 10000-10;   % S_h_U
        init_pop(2)  = 9;      % E_h_U
        init_pop(3)  = 1;      % I_h_U
        init_pop(4)  = 0;   % S_h_T
        init_pop(5)  = 0;      % E_h_T
        init_pop(6)  = 0;      % I_h_T
        % mosquitoes: simple allocation across compartments
        init_pop(7)  = 100000;  % S_s_u Susceptible at birth
        init_pop(8)  = 0;                    % I_s_u
        init_pop(9)  = 0;                    % S_s_e
        init_pop(10) = 0;                    % I_s_e
        init_pop(11) = 0;     % S_r_u
        init_pop(12) = 0;                    % I_r_u
        init_pop(13) = 0;                    % S_r_e
        init_pop(14) = 0;

        [Tdeteq,Ydeteq] = HumanSEI_VectorSI(mu_h, sigma_h, gamma_h, gamma_T, ...
            m_rate, w_rate, mu_v, a, a_u, a_e, beta_vh, ...
            phi, phi_r, rho_v, d_s, d_r, ...
            rs, rr, K, beta_hv_u, beta_hv_e, simTime,init_pop);

        TotalMosquitoes=sum(Ydeteq(end,7:14));
        TotalHumans=sum(Ydeteq(end,1:6));

        ratioMosquitoesToHumans=TotalMosquitoes/TotalHumans;

        R0=ratioMosquitoesToHumans*(((a^2)*beta_vh*beta_hv_u)/((mu_h+gamma_h)*mu_v))*(sigma_h/(sigma_h+mu_h));

        init_pop=Ydeteq(end,:);

        prop_start_resistant=0.0001;

        init_pop(11)=prop_start_resistant*init_pop(7);
        init_pop(7)=(1-prop_start_resistant)*init_pop(7);

        %end of brought in code from pulsed model


        % ---- Run ODE model for actual simulation----
        [~, Y] = HumanSEI_VectorSI(mu_h, sigma_h, gamma_h, gamma_T, ...
            m_rate, w_rate, mu_v, a, a_u, a_e, beta_vh, ...
            phi, phi_r, rho_v, d_s, d_r, ...
            rs, rr, K, beta_hv_u, beta_hv_e, simTime);

        % ---- Resistant proportion at final time ----
        Nv = sum(Y(end, 7:14));
        Nr = sum(Y(end, 11:14));

        if Nv <= 0
            prop_resistant = NaN;
        else
            prop_resistant = Nr / Nv;
        end

        % ---- Store results ----
        results(k,:) = [a, d_s, rr_adjust, prop_resistant,prop_resistant-0.0001,R0,TotalMosquitoes];

        k=k+1;
    end

end


% Remove rows with NaN outcome, if any
valid = ~isnan(results(:,4));
results = results(valid,:);

ResultsTable=array2table(results,'VariableNames',{'a';'d_s';'rr_adjust';'prop_resistant';'changepropresistant';'R0';'EquMosqBeforeIvermectin'});

scatter(ResultsTable.R0,ResultsTable.changepropresistant,[],(1-ResultsTable.rr_adjust),'filled')

cb = colorbar;
cb.Label.String = '% fecundity cost';

xlabel ('Malaria R_0')
ylabel ({'Change in proportion resistant over year post';'1 ivermectin MDA'})