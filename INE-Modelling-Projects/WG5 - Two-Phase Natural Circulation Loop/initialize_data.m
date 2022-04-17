function data = initialize_data()

%% CONDENSER T/H      
data.T_pool = 100;        % [°C]

%% GEOMETRY DATA [m, m^2, m^3]
% Steam Generator
data.D_sg = 0.01253;
data.S_sg = 0.25*data.D_sg^2*pi;
data.Lh_sg= 24; 
data.Luh_sg = 8;
data.L_sg = data.Lh_sg + data.Luh_sg;
data.V_sg= data.S_sg*data.L_sg;

% Riser
data.D_riser=0.02093;
data.S_riser= 0.25*data.D_riser^2*pi;
data.L_rup= 10.7;
data.L_ror= 9.45;
data.L_rtot= data.L_rup + data.L_ror;
data.V_riser= data.S_riser*data.L_rtot;

% Condenser
data.D_condenser= 0.059;
data.L_condenser= 1;
data.S_condenser= 0.25*pi*data.D_condenser^2;
data.V_condenser= data.S_condenser*data.L_condenser; 

% Downcomer
data.L_dcup= 9.45 + 9.32;
data.L_dcor= 8;
data.L_dctot= data.L_dcup + data.L_dcor;
data.D_dc= 0.02093;
data.S_dc=0.25*data.D_dc^2*pi;
data.V_dc= data.S_dc*data.L_dctot;

% Filling Ratio
data.V_tot= data.V_dc + data.V_condenser + data.V_riser + data.V_sg;


data.C = 20;
data.K = 2;
end