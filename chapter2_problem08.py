import equations
T = 50 # Temperature in C
Rs = 4000 # Standard gas-oil ratio in sm^3/m^3
P = 20 # Pressure in mPa
Pb = 15 # Bubble point pressure in mPa
gamma_g = 0.77 # Specific gravity of gas
gamma_o = 0.8 # Specific gravity of oil

T = equations.C_to_F(T)
P = equations.mPa_to_psia(P)
Pb = equations.mPa_to_psia(Pb)
Rs = Rs / 5.615 # convert sm^3/m^3 to scf/STB

api = equations.calculate_api(gamma_o)
rho_o = equations.calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)
mu_o, mu_ob = equations.calculate_mu_o_Standing_method(api, Rs, P, Pb, T)

print("Above bubble point pressure hence Rs = GOR and since Rs for both is the same rho_o Would be the same as well.")
print("The density of oil is: %f lbm/ft^3" %rho_o)
print("The viscosity of oil at bubble point pressure is: %f cp" %mu_ob)
print("The viscosity of oil at %d psia is: %f cp" %(P, mu_o))