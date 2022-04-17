function [p_final, output_data] = biphase_pressure_drops(p_in, h_in, m_dot, L, height, S, D, surface_heat_flux, output)

data = initialize_data();

N=100;
delta_L = L/N;
delta_z = height/N;

C= data.C;
g=9.81;

% INITIALIZE ARRAYS
h_LS = zeros(1,N);
h_VS = zeros(1,N);
h = zeros(1,N);
my_LS = zeros(1,N);
my_VS = zeros(1,N);

rho_LS = zeros(1,N);
rho_VS = zeros(1,N);
rho = zeros(1,N);

pd_l = zeros(1,N);
martinelli_factor = zeros(1,N);

p = zeros(1,N);

p(1) = p_in;
h(1) = h_in;

%Compute initial values

for k = 1:N
    h_LS(k) = XSteam('hL_p', p(k));
    h_VS(k) = XSteam('hV_p', p(k));
    my_LS(k) = XSteam('my_ph', p(k), h_LS(k));
    my_VS(k) = XSteam('my_ph', p(k), h_VS(k));
    rho_LS(k)= XSteam('rhoL_p', p(k));
    rho_VS(k) = XSteam('rhoV_p', p(k));

    if h(k) < h_LS(k) || h(k) > h_VS(k) % if the fluid is liquid or gas
        h(k+1) = h(k) + (surface_heat_flux*pi*D*delta_L)/(m_dot*1000);
    
        my_l = XSteam('my_ph', p(k), h(k));
        rho_l = XSteam('rho_ph', p(k), h(k));
        pd_tot(k) = (-((m_dot^2)*(0.184*((m_dot*D)/(my_l*S))^(-0.2))*delta_L)/(2*rho_l*D*S^2))*1E-5 - (rho_l*g*delta_z)*1E-5;

        rho(k) = rho_l;
        x(k) = (h(k)-h_LS(k))/(h_VS(k)-h_LS(k));
    
    else % if the fluid is biphase
        h(k+1) = h(k) + (surface_heat_flux*pi*D*delta_L)/(m_dot*1000);

        x(k) = (h(k) - h_LS(k))/(h_VS(k)-h_LS(k));

        rho(k) = XSteam('rho_ph', p(k), h(k));
        
        X = ((1-x(k))/x(k))^0.9 * (rho_VS(k)/rho_LS(k))^0.5 * (my_LS/my_VS)^0.1;

        martinelli_factor(k) = 1 + C/X + 1/(X^2); % Separeted flows

        %martinelli_factor(k) = 1 + x(k)*(rho_LS(k)/rho_VS(k) - 1); % HEM

        %pd_l(k) = (-(((m_dot)^2)*(0.184*((m_dot*D)/(my_LS(k)*S))^(-0.2))*delta_L*2)/(rho_LS(k)*D*S^2))*1E-5; %da usare con HEM
        
        pd_l(k) =(-(((m_dot*(1-x(k)))^2)*(0.184*((m_dot*(1-x(k))*D)/(my_LS(k)*S))^(-0.2))*delta_L)/(2*rho_LS(k)*D*S^2))*1E-5;
        %da usare con Martinelli
        
        pd_tot(k) = pd_l(k)*martinelli_factor(k) - (rho(k)*g*delta_z)*1E-5;
    end

    p(k+1) = p(k) + pd_tot(k);
end

p_final = p(N);

%outputs all info of interest
if output
    output_data.pressure = p;
    output_data.quality = x;
    output_data.enthalpy = h;
    output_data.MF = martinelli_factor;
    output_data.density = rho;
    output_data.discretized_volume = S*delta_L;
else
    output_data = 0;
end

end