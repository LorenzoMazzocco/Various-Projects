function [FR_theory] = find_FR_theory(p,Q)

T_sat = XSteam('Tsat_p', p) + 273.15;
h_evap = XSteam('hV_p', p) - XSteam('hL_p', p);
p_kpa = p*100;
Q_kw = Q/1000;

LHS = (p_kpa^0.7)*h_evap/(T_sat*(Q_kw^0.86));

%Solve quadratic equation

a = 136;
b = -78;
c = 49 - LHS;

FR_theory = (-b + sqrt(b^2 - 4*a*c))/(2*a);

end