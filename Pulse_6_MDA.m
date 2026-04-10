clear; close all; clc;

% -------- FIXED PARAMETERS (unchanged) --------
%human parameters
mu_h     = 1/(70*365);
sigma_h  = 1/10;
gamma_h  = 1/20;
gamma_T  = gamma_h;
mu_v     = 1/10;
a_u    = 0.4;
a_e    = 0.4;
beta_vh  = 0.2;
d_r    = 0.01;

%mosquito parameters
K        = 100000;
beta_hv_u = 0.2;
beta_hv_e = beta_hv_u;


rs      = 0.2;  % The fecundity number of 0.2 here means that essentially every
%  5 days, one of the offspring of a single mosquito survives to
%  adulthood...since the mosquitoes live for ~10 days in this model, each
%  mosquito has 2 surviving offspring that make it to adulthood (which is
%  plenty to sustain the population!  

rho_v = 0;  %setting to zero for simplicity for now
w_rate =0;  %also set to zero for simplicity here, changed again later when we introduce the pulsed MDA
m_rate = 0; %setting to zero for simplicity
phi_r  = 1; %setting to 1 for simplicity
phi    = 1; %setting to 1 for simplicity

% -------------------------------------------------------
nSamples = 300;      % Number of model runs (can be adjust as needed)
rng(1);               % Fix random seed for reproducibility

a      = 0.4; %to be varied in the "one hit" analysis
rr     = rs*0.8;   %making the maximum value of the fecundity of the resistant mosquito 80% that of the susceptible mosquito- but
%the value here is something to vary in the one hit analysis
 %making sure that d_r is always smaller than d_s (because the cost of ivermectin is lower for d_r),  vary by how much in the one hit analysis
d_s   = d_r + 2; %setting to 1 because we assume susceptible mosquitoes die within half a day
% ---- Run ODE model ----

init_pop = zeros(14,1);
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

simTime=1000;

[Tdeteq,Ydeteq] = HumanSEI_VectorSI(mu_h, sigma_h, gamma_h, gamma_T, ...
    m_rate, w_rate, mu_v, a, a_u, a_e, beta_vh, ...
    phi, phi_r, rho_v, d_s, d_r, ...
    rs, rr, K, beta_hv_u, beta_hv_e, simTime,init_pop);

TotalMosquitoes=sum(Ydeteq(end,7:14));

TotalHumans=sum(Ydeteq(end,1:6));

ratioMosquitoesToHumans=TotalMosquitoes/TotalHumans;

R0=ratioMosquitoesToHumans*(((a^2)*beta_vh*beta_hv_u)/((mu_h+gamma_h)*mu_v))*(sigma_h/(sigma_h+mu_h));


init_pop=Ydeteq(end,:);  %equilibrium behaviour of infection

PulseInterval=30;  %note change of name
time_before_drop=2; %number of days before waning starts
maxpulse=6;
prop_start_resistant=0.0001;

init_pop(11)=prop_start_resistant*init_pop(7);
init_pop(7)=(1-prop_start_resistant)*init_pop(7);

proportion_treated=0.8;

OverallT=[];
OverallY=[];

for pulse=1:1:maxpulse

    simTime=time_before_drop;
    w_rate=0;
    init_pop(4)=init_pop(4)+proportion_treated*init_pop(1);
    init_pop(5)= init_pop(5)+proportion_treated*init_pop(2);
    init_pop(6)=init_pop(6)+proportion_treated*init_pop(3);

    init_pop(1)=(1-proportion_treated)*init_pop(1);
    init_pop(2)=(1-proportion_treated)*init_pop(2);
    init_pop(3)=(1-proportion_treated)*init_pop(3);

    [T,Y] = HumanSEI_VectorSI(mu_h, sigma_h, gamma_h, gamma_T, ...
        m_rate, w_rate, mu_v, a, a_u, a_e, beta_vh, ...
        phi, phi_r, rho_v, d_s, d_r, ...
        rs, rr, K, beta_hv_u, beta_hv_e, simTime,init_pop);

    init_pop=Y(end,:);

    OverallT=[OverallT;T(2:end,:)+PulseInterval*(pulse-1)];
    OverallY=[OverallY;Y(2:end,:)];

    simTime=PulseInterval-time_before_drop;
    w_rate=1/3;

    [T,Y] = HumanSEI_VectorSI(mu_h, sigma_h, gamma_h, gamma_T, ...
        m_rate, w_rate, mu_v, a, a_u, a_e, beta_vh, ...
        phi, phi_r, rho_v, d_s, d_r, ...
        rs, rr, K, beta_hv_u, beta_hv_e, simTime,init_pop);


    init_pop=Y(end,:);

    OverallT=[OverallT;T(2:end,:)+PulseInterval*(pulse-1)+time_before_drop];
    OverallY=[OverallY;Y(2:end,:)];


end


% ---- Resistant proportion Computation----
Nv = sum(OverallY(end,7:14));
Nr = sum(OverallY(end,11:14));
prop_resistant = Nr /Nv;

TotalHumansInfected=OverallY(:,2)+OverallY(:,3)+OverallY(:,5)+OverallY(:,6);

TotalHumansTreated=OverallY(:,4)+OverallY(:,5)+OverallY(:,6);

TotalMosquitoPopulation=sum(OverallY(:,7:end),2);

TotalResistantMosquitoPopulation=sum(OverallY(:,11:end),2);

ProportionResistantMosquitoes=TotalResistantMosquitoPopulation./TotalMosquitoPopulation;

figure('Name','MDA simulation results','Color','w','Units','normalized','Position',[0.2 0.2 0.45 0.6]);

% top: PLOT humans infected & treated
subplot(3,1,1);
h1 = plot(OverallT, TotalHumansInfected, '-','LineWidth',1.6); hold on;
h2 = plot(OverallT, TotalHumansTreated, '--','LineWidth',1.6);
grid on;
xlabel('Time (days)');
ylabel('Total humans');
title('Humans: infected vs treated');
legend([h1 h2], {'Total infected (E+I)','Total treated (S_T+E_T+I_T)'}, 'Location','best');
set(gca,'FontSize',11);

% middle: total mosquito population
subplot(3,1,2);
h3 = plot(OverallT, TotalMosquitoPopulation, '-','LineWidth',1.6);
grid on;
xlabel('Time (days)');
ylabel('Number of mosquitoes');
title('Total mosquito population');
set(gca,'FontSize',11);

subplot(3,1,3)

% bottom: proportion resistant mosquitoes
subplot(3,1,3);
h4 = plot(OverallT, ProportionResistantMosquitoes, '-','LineWidth',1.6);               
grid on;
xlabel('Time (days)');
ylabel('Prop Resistant');
title('Proportion of resistant mosquitoes');
set(gca,'FontSize',11);

% tight layout
drawnow;

% Save figure (optional)
saveas(gcf,'MDA_simulation_plots.png');   % PNG in current folder
% or saveas(gcf,'MDA_simulation_plots.pdf'); % save as PDF

%% --- Plotting: 3 stacked subplots ---------------------------------------------------
 

 