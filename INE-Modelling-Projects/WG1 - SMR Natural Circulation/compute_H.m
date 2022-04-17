function[H]=compute_H(input_data, core_powers,outer_diameters)
%preprocessing - geometric parameters
    %CORE
    sigma_core= (input_data.Pfr)^2 - 0.25*pi*(input_data.Dfr)^2                %calculation of the power channel area
    D_eq_core= (4*sigma_core)/(pi*input_data.Dfr)                          %calculation of the Deq for the core
    sigma_totcore= input_data.Nfa*(input_data.Nfr)^2*sigma_core            %calculation of the totale area of the core
    L_core= input_data.Lfr;
    
    %STEAM GENERATOR
    sigma_sg= (input_data.Psg)^2 - 0.25*pi*(input_data.OD)^2             %calculation of the sg channel area
    D_eq_sg= (4*sigma_sg)/(pi*input_data.OD)                           %calculation of the Deq for the SG
    A_anu = pi*(input_data.RPV/2)^2 - pi*(input_data.B/2)^2
    p = input_data.Psg
    Nt= 0.8* (A_anu)/(p^2)
    sigma_totsg= Nt*sigma_sg; %calculation of the total area on the secondary side for SG
    L= compute_L(input_data, core_powers, outer_diameters);
    
    %RISER
    riser_area= 0.25*pi*input_data.B^2;
    
    %DOWNCAMER
    D_eq_dc= 4*A_anu/(pi*input_data.RPV + pi*input_data.B)
    
    g=9.81;
    
%preprocessing - thermal parameters   
    h_out=1520
    h_in=1300
    W= core_powers;
    Q= (W/(h_out-h_in))*10^(-3);
    rho_mean= (input_data.rho_in + input_data.rho_out)*0.5;
    
%preprocessing - velocities
    v_core= Q/(rho_mean*sigma_totcore)
    v_sg= Q/(rho_mean*sigma_totsg)
    v_riser= Q/(input_data.rho_out*riser_area)
    v_downcamer= Q/(input_data.rho_in*A_anu)

%preprocessing - calculation of the reynolds number
    Re_core= (rho_mean*v_core*D_eq_core)/input_data.mu
    Re_sg=(rho_mean*v_sg*D_eq_sg)/input_data.mu
    Re_riser=(input_data.rho_out*v_riser*input_data.B)/input_data.mu
    Re_downcamer= (input_data.rho_in*v_downcamer*D_eq_dc)/input_data.mu
    
%calculation of the friction factor
    f_core= (3.8*log10(10./Re_core + 0.2*input_data.e/D_eq_core)).^(-2)
    f_sg= (3.8*log10(10./Re_sg + 0.2*input_data.e/D_eq_sg)).^(-2)
    f_riser= (3.8*log10(10./Re_riser + 0.2*input_data.e/input_data.B)).^(-2)
    f_downcamer= (3.8*log10(10./Re_downcamer + 0.2*input_data.e/D_eq_dc)).^(-2)

%calculation of the pressure drops in the core
    deltac_pg= -rho_mean*g*L_core                                                                              %gravitational pressure drops
    deltac_pf= -L_core*2*f_core.*(Q.^2)/(rho_mean*D_eq_core*(sigma_totcore).^2)                                       %frictional pressure drops
    deltac_pc= -input_data.K1*input_data.rho_in*v_core.^2*0.5 - input_data.K2*input_data.rho_out*v_core.^2*0.5              %concentrated pressure drops


%calculation of the pressure drops in the riser
    deltar_pg= -input_data.rho_out*g*L
    deltar_pf= -L.*2.*f_riser.*(Q.^2)./(input_data.rho_out*input_data.B*riser_area.^2)
    deltar_pc= 0 
    deltar_pfl= deltar_pf/L
    deltar_pgl= deltar_pg/L
    
%calculation of the pressure drops in the steam generator
    deltas_pg= +rho_mean*g*L
    deltas_pf= -L.*2.*f_sg.*(Q.^2)./(rho_mean*D_eq_sg*sigma_totsg.^2)
    deltas_pc= -input_data.K3*input_data.rho_out*v_sg.^2*0.5 - input_data.K4*input_data.rho_in*v_sg.^2*0.5
    
%calculation of the pressure drops in the downcamer 
    deltad_pg= +input_data.rho_in*g*L_core
    deltad_pf= -L_core*2.*f_downcamer.*(Q.^2)/(input_data.rho_in*D_eq_dc*A_anu^2)
    deltad_pc=0;
    deltad_pgl=deltad_pg/L_core
    deltad_pfl=deltad_pf/L_core
    
%total pressure drops
    deltatot_pg= deltac_pg + deltar_pg + deltas_pg + deltad_pg;
    deltatot_pf= deltac_pf + deltar_pf + deltas_pf + deltad_pf;
    deltatot_pc= deltac_pc + deltar_pc + deltas_pc + deltad_pc;
    deltatot_diff= deltatot_pg + deltatot_pf + deltatot_pc
   
    H= -deltatot_diff./(deltar_pfl + deltar_pgl + deltad_pfl + deltad_pgl)
end