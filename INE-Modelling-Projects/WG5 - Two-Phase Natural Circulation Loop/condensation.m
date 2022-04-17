function [T_6, output_data] = condensation(p_in, h_in, m_dot)

data = initialize_data();
g = 9.81;

N = 20;

T_sat = XSteam('Tsat_p', p_in);
h_LS = XSteam('hL_p', p_in);
h_VS = XSteam('hV_p', p_in);
h_evap = h_VS-h_LS; 
h_evap_mod = h_evap + (3/8)*XSteam('CpL_p', p_in)*(T_sat-data.T_pool);
rho_LS = XSteam('rhoL_p', p_in);
rho_VS = XSteam('rhoV_p', p_in);
k_LS = XSteam('tcV_p', p_in);
my_LS = XSteam('my_ph', p_in, h_LS);



k_steel = 16.3;
k_LS = XSteam('tcL_p', p_in);

k_liquid = zeros(1,N);
k(1) = k_LS;

W = zeros(1,N);
T = zeros(1,N);
T(1) = T_sat;

h = zeros(1,N);
h(1) = h_in;

rho = zeros(1,N);

L_bi = 0;
W_cond = 0;

%% compute the length of the tube in which the fluid is biphase and condensation power
if h_in > h_LS % check if fluid is biphase
    % compute power to exchange to condensate the mixture
    W_cond = m_dot*(h_in - XSteam('hL_p', p_in))*1000; % [W], defined positive 
    
    % compute alpha convection 
    Re_vapor = XSteam('x_ph', p_in, h_in)*m_dot*data.D_condenser/(data.S_condenser*XSteam('my_ph', p_in, h_in));

    alpha_condensation = 0.555*(g*(rho_LS-rho_VS)*rho_LS*(k_LS^3)*(h_evap_mod*1000)/(data.D_condenser*my_LS*(T_sat-data.T_pool)))^0.25; % Chato formula [W/m2*°K]

    % compute length of condensation
    L_bi = (W_cond*(1/(alpha_condensation*data.D_condenser) + log(0.073025/0.059)/(2*k_steel)))/(pi*(T_sat - data.T_pool)); 
    
    if L_bi > data.L_condenser
        T_6 = T_sat;
        output_data.exit_condenser_biphase = 1;
        output_data.mass = XSteam('rho_ph', p_in, h_in)*data.V_condenser;
        return
    end
    h(1) = XSteam('hL_p', p_in);
else
    %display('Entered condenser condensated')
end

rho_in = XSteam('rho_ph', p_in, h_in);
rho_mid = XSteam('rho_ph', p_in, h(1));
rho_mean = 0.5*(rho_in + rho_mid);
mass_bi = rho_mean*data.S_condenser*L_bi;

rho(1) = rho_mid;
%% heat transfer for liquid fluid

L_mono = data.L_condenser - L_bi;
delta_L = L_mono/N;

R_th_cond = log(0.073025/0.059)/(2*pi*delta_L*16.3);


for i = 1:N-1
    %compute alpha convection and total thermal resistance
    Re = (m_dot/(pi*(data.D_condenser/2)^2))*data.D_condenser/XSteam('my_ph', p_in, h(i)); 

    alpha = 0.023*(Re^0.8)*(7.56^0.3)*(k(i)/data.D_condenser); % Dittus-Boelter (Pr_water = 7.56)

    R_tot = 1/(pi*data.D_condenser*delta_L*alpha) + R_th_cond;

    % compute W and new enthalpy
    W(i) = (T(i) - data.T_pool)/R_tot;
    h(i+1) = h(i) - W(i)/(1000*m_dot);
    rho(i+1) = XSteam('rho_ph', p_in, h(i));

    % update thermodynamic variables
    k(i+1) = XSteam('tc_ph', p_in, h(i+1));
    T(i+1) = XSteam('T_ph', p_in, h(i+1));
end

%% compute return quantities
T_6 = T(N);
W_sub = sum(W); % [W] defined positive
W_tot = W_cond + W_sub;

output_data.W_tot = W_tot;
output_data.L_bi = L_bi;
output_data.h_in = h_in;
output_data.enthalpies = h;
output_data.temperatures = T;
output_data.mass = mass_bi + sum(rho*data.S_condenser*delta_L);

end
