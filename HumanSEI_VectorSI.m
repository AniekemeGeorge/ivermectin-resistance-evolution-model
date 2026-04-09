function [T,Y] = HumanSEI_VectorSI (mu_h, sigma_h, gamma_h, gamma_T, m_rate, w_rate, ...
     mu_v, a, a_u,a_e, beta_vh, phi, phi_r, ...
    rho_v, d_s, d_r, rs, rr, K, beta_hv_u, beta_hv_e, time, init_pop_in)

%  This wrapper:
%   - builds pars = [ ... ] (same order as the equations function expects)
%   - sets sensible example init_pop (14 states)
%   - calls ode45(@equations, 0:0.01:time, init_pop, odeset(...), pars)
%
%  Model state ordering (y):
%   1 S_h_U  2 E_h_U  3 I_h_U  4 S_h_T  5 E_h_T  6 I_h_T
%   7 S_s_u  8 I_s_u  9 S_s_e 10 I_s_e 11 S_r_u 12 I_r_u 13 S_r_e 14 I_r_e

% Parameter vector is built in the order the equations expect
pars = [mu_h,sigma_h, gamma_h, gamma_T, m_rate, w_rate, ...
         mu_v, a, a_u, a_e, beta_vh, phi, phi_r, ...
         rho_v, d_s, d_r, rs, rr, K, beta_hv_u, beta_hv_e];

% Initial conditions (to be edited as needed)
if nargin >= 23 && ~isempty(init_pop_in)
    % use the provided initial state vector (must be 14x1)
    init_pop = init_pop_in(:);           % ensure column
else

% Initial conditions (to be edited as needed)
init_pop = zeros(14,1);
init_pop(1)  = 4800;   % S_h_U
init_pop(2)  = 50;      % E_h_U
init_pop(3)  = 30;      % I_h_U
init_pop(4)  = 5000;   % S_h_T                            
init_pop(5)  = 0;      % E_h_T
init_pop(6)  = 0;      % I_h_T

% mosquitoes: simple allocation across compartments
totalMosq = 100000;
init_pop(7)  = totalMosq*0.999;  % S_s_u Susceptible at birth  
init_pop(8)  = 0;                    % I_s_u
init_pop(9)  = 0;                    % S_s_e
init_pop(10) = 0;                    % I_s_e
init_pop(11) = totalMosq*0.001;     % S_r_u  
init_pop(12) = 0;                    % I_r_u
init_pop(13) = 0;                    % S_r_e
init_pop(14) = 0;                    % I_r_e
end

% Solve with ode45  
tspan = 0:0.01:time;
options = odeset('NonNegative', 1:14, 'RelTol', 1e-6, 'AbsTol', 1e-8);
[T,Y] = ode45(@equations, tspan, init_pop, options, pars);

end


 
% Equations function (same file) - formatted like BasicSIR's "equations"
% -------------------------------------------------------------------------
function dydt = equations(~, y, pars)
% dydt = equations(t,y,pars)
% y = state vector (14) inputs
% pars = parameter vector (19)

% Unpack states for clarity
S_h_U = y(1); E_h_U = y(2); I_h_U = y(3);
S_h_T = y(4); E_h_T = y(5); I_h_T = y(6);

S_s_u = y(7);  I_s_u = y(8);
S_s_e = y(9);  I_s_e = y(10);
S_r_u = y(11); I_r_u = y(12);
S_r_e = y(13); I_r_e = y(14);

% Unpack parameters in the same order as the wrapper
 
mu_h    = pars(1); %(human natural mortality, day^-1)
sigma_h = pars(2); %(human progression rate E->I, day^-1)
gamma_h = pars(3); %(human recovery rate untreated, day^-1)
gamma_T = pars(4); %(human recovery rate treated, day^-1)
m_rate  = pars(5); %(rate untreated->treated, day^-1)    % m(t) in eqn doc
w_rate  = pars(6); % (rate waning treated->untreated, day^-1) % w(t) in eqn doc
mu_v    = pars(7); %(baseline mosquito mortality, day^-1)
a       = pars(8); %(mosquito biting rate, bites per mosquito per day)
a_u     = pars(9); % biting rate unexposed mosquitoes
a_e     = pars(10); % biting rate exposed-to-ivermectin mosquitoes      
beta_vh = pars(11); % (prob human infected per bite from infectious mosquito)
phi     = pars(12);%(prob bite on treated human exposes susceptible vector genotype)
phi_r   = pars(13); %(prob bite on treated human exposes resistant vector genotype)
rho_v   = pars(14); %(rate leaving 'exposed' back to unexposed due to adaptation, detoxification, recovery; can be kept at 0)
d_s      = pars(15); %(extra mortality for exposed susceptible mosquitoes (day^-1))
d_r      = pars(16); %(extra mortality for exposed resistant mosquitoes (day^-1))
rs       = pars(17); %per-capita birth rate of susceptible genotype
rr       = pars(18); %per-capita birth rate of resistant genotype
K        = pars(19); %mosquito carrying capacity
beta_hv_u = pars(20) ;%(prob unexposed mosquito infected per bite on infectious human)
beta_hv_e = pars(21) ;%(prob exposed mosquito infected per bite on infectious human)
 


% Totals and forces
%Human force of infection
%I_h = I_h_U + I_h_T;
N_h = S_h_U + E_h_U + I_h_U + S_h_T + E_h_T + I_h_T; % Total Human population
Ns = S_s_u + I_s_u + S_s_e + I_s_e;   % total susceptible-genotype adults
Nr = S_r_u + I_r_u + S_r_e + I_r_e;   % total resistant-genotype adults
Nv = Ns + Nr;  %total Vector Population - Susceptible + Resistant populations
B_h = mu_h * N_h; %Constant human population


if Nv > 0 && N_h > 0 
    lambda_h = a * beta_vh * (I_s_u + I_s_e + I_r_u + I_r_e) / N_h; % Force of Infection_Vector
    c_frac = (S_h_T + E_h_T + I_h_T) / N_h; % fraction of humans treated with Ivermectin
else
     lambda_h = 0;
     c_frac = 0;
end
%Vector general force of Infection
%lambda_v = a * beta_hv * (I_h_U + I_h_T) / N_h;

%specific force of infections
if N_h > 0 
    lambda_v1 = a_u * beta_hv_u * (I_h_U) / N_h; %Force of infection_Human_Untreated with Ivermectin
    lambda_v2 = a_e * beta_hv_e * (I_h_T) / N_h; %Force of infection_Human_Treated with Ivermectin
    lambda_v3 =  lambda_v1 + lambda_v2;
else
    lambda_v1 = 0;
    lambda_v2 = 0;
    lambda_v3 = 0;

end

if Nv > 0
   B_s = (rs*Ns + rr*Nr) * (1 - Nv/K) * ((rs*Ns)/(rs*Ns + rr*Nr));% Susceptible Mosquito birth/Carrying capacity in the population
   B_r = (rs*Ns + rr*Nr) * (1 - Nv/K) * ((rr*Nr)/(rs*Ns + rr*Nr));% Resistant Mosquito birth/Carrying capacity in the population
else
    B_s = 0;
    B_r = 0;
end   

% Human equations: SEIS Model in Human Compartment

% 1. Untreated susceptible
dS_h_U = B_h - lambda_h * S_h_U - mu_h * S_h_U - m_rate * S_h_U + w_rate * S_h_T + gamma_h * I_h_U; 
% 2. Untreated exposed
dE_h_U = lambda_h * S_h_U - sigma_h * E_h_U - mu_h * E_h_U - m_rate * E_h_U + w_rate * E_h_T;
% 3. Untreated infected
dI_h_U = sigma_h * E_h_U - gamma_h * I_h_U - mu_h * I_h_U - m_rate * I_h_U + w_rate * I_h_T;

% 4. Treated susceptible
dS_h_T = m_rate * S_h_U - w_rate * S_h_T - lambda_h * S_h_T - mu_h * S_h_T + gamma_T * I_h_T;
% 5. Treated exposed
dE_h_T = m_rate * E_h_U + lambda_h * S_h_T - sigma_h * E_h_T - w_rate * E_h_T - mu_h * E_h_T;
% 6. Treated infected
dI_h_T = m_rate * I_h_U - w_rate * I_h_T + sigma_h * E_h_T - gamma_T * I_h_T - mu_h * I_h_T;



% Mosquito equations: SI Model Used for the Moaquito compartment

% 7. S_s^u  (susceptible genotype, unexposed)
dS_s_u = B_s - lambda_v1 * S_s_u - mu_v * S_s_u - a * c_frac * phi * S_s_u;
% 8. I_s^u  (susceptible genotype, unexposed, infectious)
dI_s_u = lambda_v1 * S_s_u - mu_v * I_s_u - a * c_frac * phi * I_s_u;


% 9. S_s^e (susceptible genotype, exposed to ivermectin, susceptible)
dS_s_e = a * c_frac * phi * S_s_u - rho_v * S_s_e - lambda_v3 * S_s_e - (mu_v + d_s) * S_s_e;
% 10. I_s^e (susceptible genotype, exposed, infectious)
dI_s_e = a * c_frac * phi * I_s_u - rho_v * I_s_e + lambda_v3 *S_s_e +  lambda_v2 * S_s_u - (mu_v + d_s) * I_s_e;


% 11. S_r^u (resistant genotype, unexposed)
dS_r_u = B_r - a * c_frac * phi_r * S_r_u - lambda_v3 * S_r_u - mu_v * S_r_u;
% 12. I_r^u (resistant genotype, unexposed, infectious)
dI_r_u =  lambda_v1 * S_r_u - mu_v * I_r_u - a * c_frac * phi_r * I_r_u;


% 13. S_r^e (resistant genotype, exposed to ivermectin, susceptible)
dS_r_e = a * c_frac * phi_r * S_r_u - lambda_v3 * S_r_e - (mu_v + d_r) * S_r_e;
% 14. I_r^e (resistant genotype, exposed, infectious)
dI_r_e = a * c_frac * phi_r * I_r_u + lambda_v3 *S_r_u +  lambda_v2 *S_r_e - (mu_v + d_r) * I_r_e;

% Pack into column vector (derivative)  
dydt = [
    dS_h_U;
    dE_h_U;
    dI_h_U;
    dS_h_T;
    dE_h_T;
    dI_h_T;
    dS_s_u;
    dI_s_u;
    dS_s_e;
    dI_s_e;
    dS_r_u;
    dI_r_u;
    dS_r_e;
    dI_r_e
    ];

end