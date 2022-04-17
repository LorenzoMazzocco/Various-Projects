function [L] = compute_L(input_data, core_powers, outer_diameters)
    W = core_powers;
    p = input_data.Psg;
    alpha = input_data.alpha;
    T_in = input_data.Tout;
    T_out = input_data.Tin;
    T_sg = input_data.Tsg;
    delta_T_log = ((T_in - T_sg) - (T_out - T_sg))/log((T_in - T_sg)/(T_out - T_sg));
    D = outer_diameters;
    A_anu = pi*(input_data.RPV/2)^2 - pi*(input_data.B/2)^2;

    L = (W*(p^2))/(alpha*delta_T_log*pi*D*0.8*A_anu);
end