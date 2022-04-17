function [T_1] = find_temperature(p_1, Q)

display('Looking for T_1');

data = initialize_data();

T_1 = 100;

diff_T = 100;

while diff_T > 0.5
[m_dot, converged, pd_tot_sg, pd_tot_riser, h_4] = momentum_loop(p_1, T_1, Q);

    if ~converged
        T_1
        display('momentum loop did not converge in find_temperature')
        return
    else
       %compute T_5
       p_4 = p_1 + pd_tot_sg + pd_tot_riser;
       T_6 = condensation(p_4, h_4, m_dot);
       diff_T = abs(T_1-T_6);
       display(diff_T, 'Current difference ')
       if diff_T < 10
           T_1 = T_1 + 1
       else
           T_1 = T_1 + 10
       end 
    end
end

T_1 = T_1 - 1

end
