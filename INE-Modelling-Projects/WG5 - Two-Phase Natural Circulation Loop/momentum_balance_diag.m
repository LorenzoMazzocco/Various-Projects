function [] = momentum_balance_diag(p_1, T_1, Q)

data = initialize_data();

K = data.K;
C = data.C;
g = 9.81;
N = 400;



pds = zeros(1,N);
pgs = zeros(1,N);

pds_sg = zeros(1,N);

pds_riser_or = zeros(1,N);
pds_riser_up = zeros(1,N);

m_dot = 0.005;
difference = 1;
i = 1;
for i=1:N
    i
    m_dot = m_dot + 0.0005; %step minimo deve essere 0.0005

% Pressure drops in STEAM GENERATOR
    h_1 = XSteam('h_pT', p_1, T_1);
    rho_1 = XSteam('rhoL_p', p_1);
    my_1 = XSteam('my_pT', p_1,T_1);
    
    pd_conc_in_sg = (- (m_dot^2)*(K/(2*(rho_1)*data.S_sg^2)))*1E-5;
    
    %heated
    [p_3_, heated_sg_data] = biphase_pressure_drops(p_1, XSteam('h_pT', p_1, T_1), m_dot ,data.Lh_sg, data.Lh_sg*sin(14.3/180*pi), data.S_sg, data.D_sg, Q/(data.Lh_sg*pi*data.D_sg), false);
    pd_distr_h_sg = -p_1 + p_3_;
    %unheated
    [p_3, unheated_sg_data] = biphase_pressure_drops(p_3_, XSteam('h_pT', p_1, T_1) + Q/(m_dot*1000), m_dot ,data.Luh_sg, data.Luh_sg*sin(14.3/180*pi), data.S_sg, data.D_sg, 0, false);
    pd_distr_uh_sg = -p_3_ + p_3;

    pd_distr_sg = pd_distr_h_sg + pd_distr_uh_sg;
    pds_sg(i) = pd_distr_sg; 

    p_3 = p_3 + pd_conc_in_sg;
    
    h_3 = h_1 + Q/(m_dot*1000);
    rho_3 = XSteam('rho_ph', p_3, h_3);
    
    pd_conc_out_sg = (- (m_dot^2)*K/(2*rho_3*data.S_sg^2))*1E-5;
    
    p_3 = p_3 + pd_conc_out_sg;
    
    % Pressure drops in RISER
    [p_4_, horizontal_riser_data] = biphase_pressure_drops(p_3, h_3, m_dot, data.L_ror, 0, data.S_riser, data.D_riser, 0, false);
    pd_distr_riser_or = -p_3 + p_4_;
    [p_4, vertical_riser_data] = biphase_pressure_drops(p_4_, h_3, m_dot, data.L_rup, data.L_rup, data.S_riser, data.D_riser, 0, false);
    pd_distr_riser_up = -p_4_ + p_4; 
    
    pd_distr_riser = pd_distr_riser_or + pd_distr_riser_up;

    h_4 = h_3;
    rho_4 = XSteam('rho_ph', p_4, h_4);
    
    pd_conc_riser = (m_dot^2)*(- K/(2*rho_4*data.S_riser^2))*1E-5;
    p_4 = p_4 + pd_conc_riser;
    
    % Pressure drops/gains in DOWNCOMER
    my_out_dc = XSteam('my_pT', p_1, T_1);
    my_in_dc = XSteam('my_pT', p_4,T_1); 
    my_mean_dc= (my_out_dc + my_in_dc)*0.5;
    pd_distr_dc = (-((m_dot^2)*(0.184*((m_dot*data.D_dc)/(my_mean_dc*data.S_dc))^(-0.2))*data.L_dctot)/(2*rho_1*data.D_dc*data.S_dc^2))*1E-5;
    p_6 = p_4;
    rho_6 = XSteam('rho_pT', p_6, T_1);
    
    pd_conc_in_dc = (- (m_dot^2)*(K/(2*(rho_6)*data.S_dc^2)))*1E-5;
    
    pg_grav_dc = ( +((rho_1+rho_6)/2)*g*data.L_dcup)*1E-5;
    
    
    % ALL LOOP PRESSURE DROPS
    pd_tot = pd_conc_in_sg + pd_distr_sg + pd_conc_out_sg + pd_distr_riser + pd_conc_riser + pd_distr_dc + pd_conc_in_dc;
    pg_tot = pg_grav_dc;
    
    pds(i) = pd_tot;
    pgs(i) = pg_tot;
    pdd(i) = pd_distr_dc;
    difference = pd_tot + pg_tot;
    
    if isreal(pd_tot) && isreal(pg_tot) && ~isnan(pd_tot) && ~isnan(pg_tot)
        difference = abs(pd_tot+pg_tot);
    else
        pd_tot;
        pg_tot;
        difference = 1;
    end
    
    pd_tot_sg = p_3 - p_1;
    pd_tot_riser = p_4 - p_3;

end

m_dots = linspace(0,m_dot, N);
hold on
plot(m_dots, abs(pds))
plot(m_dots, pgs)
hold off

% hold on
% plot(m_dots,abs(pds_sg))
% hold off

%hold on
%plot(m_dots, abs(pds_riser_or))
%plot(m_dots, abs(pds_riser_up))
%hold off

end