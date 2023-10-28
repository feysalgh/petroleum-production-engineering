import equations
gamma_g = 0.7 # Specific gravity of gas
T = 200 # Temperature in F
P1 = 100 # Pressure in psia
P2 = 1000 # Pressure in psia
P3 = 5000 # Pressure in psia
P4 = 10000 # Pressure in psia
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

print("gas viscosity at 100 psia is %f cp" %mu_g1)
print("gas viscosity at 1000 psia is %f cp" %mu_g2)
print("gas viscosity at 5000 psia is %f cp" %mu_g3)
print("gas viscosity at 10000 psia is %f cp" %mu_g4)
