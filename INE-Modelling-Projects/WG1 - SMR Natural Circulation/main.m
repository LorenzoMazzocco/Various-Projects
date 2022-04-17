
%% MODEL INPUT DATA

% CORE
input_data.W= 300*10^6;       % Wth
input_data.Tin= 292;          % °C
input_data.hin= 1300*10^3;    % J/kg
input_data.rho_in= 727.6;     % kg/m^3
input_data.Tout= 329;         % °C
input_data.hout= 1520*10^3;   % J/kg
input_data.rho_out= 643.5;    % kg/m^3
input_data.Pc= 155*10^5;      % Pa
input_data.Nfa= 121;          % number of fuel assemblies
input_data.Nfr= 17;           % number of fuel rods per fuel assembly
input_data.Dfr= 0.0095;       % [m] fuel rod diameter
input_data.Pfr= 12.67*10^(-3);% [m] fuel rod pitch
input_data.Lfr= 3;            % [m] fuel rod length

% STEAM GENERATOR
input_data.OD= 0.010;         % [m] sg outer diameter
input_data.ID= 0.0085;        % [m] sg inner diameter
input_data.Psg= 0.015;       % [m] sg tube pitch
input_data.Ps=65*10^5;        % steam generator pressure
input_data.Tsg=280.82;        % degrees saturation temperature for 65 bars 
input_data.Ksg= 30;           % W/mK tube thermal conductivity
input_data.RPV= 3.75;         % reactor pressure vessel
input_data.B = 2.75;          % [m] barrel outer diameter
input_data.alpha= 5100;       % W/m^2K global heat transfer coefficient

% PRESSURE DROP COEFFICIENTS
input_data.K1= 4;             % core support plate
input_data.K2= 4;             % core upper plate
input_data.K3= 3.5;           % steam generator area inlet
input_data.K4= 3.5;           % steam generator area outlet
input_data.e= 4*10^-6;        % [m] sg tube roughness
input_data.mu= 8.284*10^(-5); % [kg/ m s] fluid viscosity



%% PARAMETRIC INPUTS

core_powers = [300, 400, 500, 600]*10^6;
secondary_pressures = input_data.Ps;
outer_diameters = input_data.OD;



%% FIND L
L = compute_L(input_data, core_powers, outer_diameters)
figure(1)
plot(core_powers, L)
xlabel('core powers [MWth]')
ylabel('SG lenght [m]')
title('SG-Core linear relation')

%% FIND H
H = compute_H(input_data, core_powers, outer_diameters)
figure(2)
plot(core_powers, H)
xlabel('core powers [Mwth]')
ylabel('Clearance value [m]')
title('Core-SG clearance as a function of core powers')
