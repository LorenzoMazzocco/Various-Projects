clear all
clc

continuation = true;


load('inputs_to_model.mat');

n = size(input_to_model, 2);

if continuation
    load('output_model.mat');
    k = size(output, 2) + 1; % continues from the last result
else
    k=1;
end

for i=k:n

        p_1 = input_to_model(1,i)
        Q = input_to_model(2,i)*1000

        T_1 = find_temperature(p_1,Q);
        [m_dot, momentum_converged] = momentum_loop(p_1,T_1,Q);
        
        data = initialize_data();
        
        K = data.K;
        g=9.81;
        
        
        if momentum_converged
        
            %% STEAM GENERATOR
            
            h_1 = XSteam('h_pT', p_1, T_1);    
            rho_1 = XSteam('rhoL_p', p_1);
            my_1 = XSteam('my_pT', p_1,T_1);
            
            pd_conc_in_sg = (- (m_dot^2)*(K/(2*(rho_1)*data.S_sg^2)))*1E-5; 
            
            %heated
            [p_2, heated_sg_data] = biphase_pressure_drops(p_1, h_1, m_dot ,data.Lh_sg, data.Lh_sg*sin(14.3/180*pi), data.S_sg, data.D_sg, Q/(data.Lh_sg*pi*data.D_sg), true);
            pd_distr_h_sg = -p_1 + p_2;
            h_2 = h_1 + Q/(m_dot*1000);
            %unheated
            [p_3, unheated_sg_data] = biphase_pressure_drops(p_2, h_2, m_dot ,data.Luh_sg, data.Luh_sg*sin(14.3/180*pi), data.S_sg, data.D_sg, 0, true);
            pd_distr_uh_sg = -p_2 + p_3;
            
            pd_distr_sg = pd_distr_h_sg + pd_distr_uh_sg;
            
            p_3 = p_3 + pd_conc_in_sg;
            h_3 = h_1 + Q/(m_dot*1000);
            T_3 = XSteam('T_ph', p_3, h_3);
            rho_3 = XSteam('rho_ph', p_3, h_3);
            
            pd_conc_out_sg = (- (m_dot^2)*K/(2*rho_3*data.S_sg^2))*1E-5;
            
            p_3 = p_3 + pd_conc_out_sg;
            
            
            %% RISER
            
            [p_4_, horizontal_riser_data] = biphase_pressure_drops(p_3, h_3, m_dot, data.L_ror, 0, data.S_riser, data.D_riser, 0, true);
            pd_distr_riser_or = -p_3 + p_4_;
            [p_4, vertical_riser_data] = biphase_pressure_drops(p_4_, h_3, m_dot, data.L_rup, data.L_rup, data.S_riser, data.D_riser, 0, true);
            pd_distr_riser_up = -p_4_ + p_4; 
            
            pd_distr_riser = pd_distr_riser_or + pd_distr_riser_up;
            
            h_4 = h_3;
            rho_4 = XSteam('rho_ph', p_4, h_4);
            
            pd_conc_riser = (m_dot^2)*(- K/(2*rho_4*data.S_riser^2))*1E-5;
            p_4 = p_4 + pd_conc_riser;
            
            %% CONDENSER
            [T_6, condenser_data] = condensation(p_4, h_4, m_dot);
            x_4 = XSteam('x_ph', p_4, h_4);
            
            
            %% DOWNCOMER
            my_out_dc = XSteam('my_pT', p_1, T_1);
            my_in_dc = XSteam('my_pT', p_4,T_1); 
            my_mean_dc= (my_out_dc + my_in_dc)*0.5; 
            pd_distr_dc = (-((m_dot^2)*(0.184*((m_dot*data.D_dc)/(my_mean_dc*data.S_dc))^(-0.2))*data.L_dctot)/(2*rho_1*data.D_dc*data.S_dc^2))*1E-5;                       
            p_6 = p_4;
            rho_6 = XSteam('rho_pT', p_6, T_1); 
            
            pd_conc_in_dc = (- (m_dot^2)*(K/(2*(rho_6)*data.S_dc^2)))*1E-5;
            
            pg_grav_dc = ( +((rho_1+rho_6)/2)*g*data.L_dcup)*1E-5;
            
            
            
            %% LOOP
            % all loop pressure drops and gains
            pd_tot = pd_conc_in_sg + pd_distr_sg + pd_conc_out_sg + pd_distr_riser + pd_conc_riser + pd_distr_dc + pd_conc_in_dc;
            pg_tot = pg_grav_dc;
            
            
            %% FILLING RATIO
            heated_sg_mass = sum(heated_sg_data.discretized_volume*heated_sg_data.density);
            unheated_sg_mass = sum(unheated_sg_data.discretized_volume*unheated_sg_data.density);
            horizontal_riser_mass = sum(horizontal_riser_data.discretized_volume*horizontal_riser_data.density);
            vertical_riser_mass = sum(vertical_riser_data.discretized_volume*vertical_riser_data.density);
            condenser_mass = condenser_data.mass;
            downcomer_mass = 0.5*(rho_1+rho_6)*data.V_dc;
            
            total_mass = heated_sg_mass + unheated_sg_mass + horizontal_riser_mass + vertical_riser_mass + condenser_mass + downcomer_mass;
            M0 = rho_1 * data.V_tot;
            
            FR_computed = total_mass/M0;
            FR_theory = find_FR_theory(p_1, Q);
            
            delta_subcooling_condenser = XSteam('Tsat_p', p_1) - T_6;

            % save data to output struct
            output(i).p_1 = p_1;                                            % p_1
            output(i).Q = Q;                                                % Q
            output(i).FR_theory = input_to_model(3,i);                      % FR_theory
            output(i).converged = 1;                                        % converged
            output(i).FR_computed = FR_computed;                            % FR computed
            output(i).m_dot = m_dot;                                        % m_dot
            output(i).dT_subcooling = delta_subcooling_condenser;           % dT_subcooling
            output(i).x_4 = x_4;                                            % quality at inlet of condenser
            output(i).heated_sg_data = heated_sg_data;                      % steam generator heated data (list of qualities, pressure, enthalpy)
            output(i).unheated_sg_data = unheated_sg_data;                  % steam generator unheated data (list of qualities, pressure, enthalpy)
            output(i).horizontal_riser_data = horizontal_riser_data;        % riser horizontal data (list of qualities, pressure, enthalpy)
            output(i).vertical_riser_data = vertical_riser_data;            % riser vertical data (list of qualities, pressure, enthalpy)
            output(i).condenser_data = condenser_data;                      % condenser data 
            output(i).T_1 = T_1;                                            % T_1
            output(i).h_1 = h_1;                                            % coordinates for 1
            output(i).p_2 = p_2;                                            % coordinates for 2
            output(i).h_2 = h_2;                                            % coordinates for 2
            output(i).p_3 = p_3;                                            % coordinates for 3
            output(i).h_3 = h_3;                                            % coordinates for 3
            output(i).p_4 = p_4;                                            % coordinates for 4
            output(i).h_4 = h_4;                                            % coordinates for 4
            output(i).p_6 = p_6;                                            % coordinates for 6
                
        
        else
            % MOMENTUM DID NOT CONVERGE
            output(i).p_1 = p_1;                                            % p_1
            output(i).Q = Q;                                                % Q
            output(i).FR_theory = input_to_model(3,i);                      % FR_theory
            output(i).converged = 0;                                        % converged
        end

        % aggiorna il file .mat
        save('output_model.mat', 'output');
        clc
end