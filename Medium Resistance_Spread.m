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
nSamples = 200;
rng(1);   % Reproducibility

% Results:
% columns = [a, d_s, rr, prop_resistant]
results = zeros(nSamples, 7);  %added extra column to keep track of malaria Ro and mosquito numbers before ivermectin

for k = 1:nSamples %Medium resistance spread Varied parameters
    a   = 0.05 + 0.45*rand;              % biting rate always varying between 0.05 and 0.5, means R0 varies between less than 1 and up to 10
    rr  = (0.8 + 0.1*rand) * rs;       % moderate rr
    d_s = d_r + (0.50 + 2*rand);      % susceptibles die between 0.5 days and 2 days after exposure


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
    results(k,:) = [a, d_s, rr, prop_resistant,prop_resistant-0.0001,R0,TotalMosquitoes];
end

% Remove rows with NaN outcome, if any
valid = ~isnan(results(:,4));
results = results(valid,:);

ResultsTable=array2table(results,'VariableNames',{'a';'d_s';'rr';'prop_resistant';'changepropresistant';'R0';'EquMosqBeforeIvermectin'});

% -------------------------------------------------------
% SCATTER / HISTOGRAM MATRIX
% -------------------------------------------------------
X = results(:,1:3);
severity = results(:,4);
changeinprop=results(:,5);
nVars = size(X,2);
varNames = {'a', 'd_s', 'rr'};

% Spearman correlations
rho_vals = zeros(1, nVars);
p_vals   = zeros(1, nVars);
for k = 1:nVars
    [rho_vals(k), p_vals(k)] = corr(X(:,k), severity, 'Type', 'Spearman', 'Rows', 'complete');
end

figure('Color','w');
for i = 1:nVars
    for j = 1:nVars
        subplot(nVars, nVars, (i-1)*nVars + j)

        if i == j

mask=changeinprop>0;  %filtering so that only the values of the relevant parameter where resistance actually spread in the population are displayed

            histogram(X(mask,i), 20);
            xlabel(varNames{i}, 'Interpreter', 'none');
            ylabel('Count');

%annotate with Spearman correlation for the x-variable in that panel
             text(0.05, 0.90, sprintf('\\rho(%s, y)=%.2f', varNames{j}, rho_vals(j)), ...
                 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
        else
            scatter(X(:,j), X(:,i), 10, severity, 'filled');
            clim([min(severity) max(severity)]);
            box on;
                          
        end

        if i == nVars
            xlabel(varNames{j}, 'Interpreter', 'none');
        end

        if j == 1
            ylabel(varNames{i}, 'Interpreter', 'none');
        end

colormap("turbo");

        if  ((i-1)*nVars + j)==3

cb = colorbar;
cb.Label.String = 'Proportion resistant';
        end

    end
end



% Print correlation ranking
[~, idx] = sort(abs(rho_vals), 'descend');
fprintf('\nSpearman correlation ranking (absolute value):\n');
for k = 1:nVars
    i = idx(k);
    fprintf('%s : rho = %.3f, p = %.3g\n', varNames{i}, rho_vals(i), p_vals(i));
end

% -------------------------------------------------------
% BINNED MEAN TRENDS
% -------------------------------------------------------
nbins = 6;
figure('Color','w');

for k = 1:size(X,2)
    subplot(1, 3, k)

    edges = linspace(min(X(:,k)), max(X(:,k)), nbins+1);
    [~,~,bin] = histcounts(X(:,k), edges);

    meanY = accumarray(bin(bin > 0), severity(bin > 0), [nbins,1], @mean, NaN);
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    plot(centers, meanY, '-o', 'LineWidth', 1.5);
    xlabel(varNames{k}, 'Interpreter', 'none');
    ylabel('Mean prop\_resistant');
    title(sprintf('%s: mean per bin', varNames{k}));
    grid on;
    box on;
end

sgtitle('Binned Mean Resistance (Numeric Trend by Parameter)');

% -------------------------------------------------------
% SUMMARY STATISTICS
% -------------------------------------------------------
fprintf('\nResistance outcome summary:\n');
fprintf('Mean = %.3f\n', mean(severity, 'omitnan'));
fprintf('Std  = %.3f\n', std(severity, 'omitnan'));
fprintf('Min  = %.3f\n', min(severity));
fprintf('Max  = %.3f\n', max(severity));

% % -------------------------------------------------------
% % SAVE RESULTS
% % -------------------------------------------------------
% ResultsTable = array2table(results, ...
%     'VariableNames', {'a', 'd_s', 'rr', 'prop_resistant'});
% 
% writetable(ResultsTable, 'ResistanceSweep_RandomUniform.csv');