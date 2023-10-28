import equations
gamma_g = 0.7 # Specific gravity of gas
T = 90 # Temperature in C
P1 = 1 # Pressure in mPa
P2 = 5 # Pressure in mPa
P3 = 10 # Pressure in mPa
P4 = 50 # Pressure in mPa
T = equations.C_to_F(T)
P1 = equations.mPa_to_psia(P1)
P2 = equations.mPa_to_psia(P2)
P3 = equations.mPa_to_psia(P3)
P4 = equations.mPa_to_psia(P4)

Ppc = equations.calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc = equations.calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr = (T + 459.67) / Tpc
Ppr1 = P1 / Ppc
Ppr2 = P2 / Ppc
Ppr3 = P3 / Ppc
Ppr4 = P4 / Ppc

mu_g1 = equations.calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr1, Tpr)
mu_g2 = equations.calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr2, Tpr)
mu_g3 = equations.calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr3, Tpr)
mu_g4 = equations.calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr4, Tpr)

print("gas viscosity at 1 mPa is %f cp" %mu_g1)
print("gas viscosity at 5 mPa is %f cp" %mu_g2)
print("gas viscosity at 10 mPa is %f cp" %mu_g3)
print("gas viscosity at 50 mPa is %f cp" %mu_g4)
