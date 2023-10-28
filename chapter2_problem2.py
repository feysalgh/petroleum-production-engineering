import equations
api = 35 # API gravity
T = 120 # Temperature in F
Rs = 800 # Standard gas-oil ratio in scf/STB
P = 3000 # Pressure in psia
Pb = 2500 # Bubble point pressure in psia

gamma_g = 0.77 # Specific gravity of gas
gamma_o = equations.calculate_gamma_o_using_api(api)
rho_o = equations.calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)
mu_o, mu_ob = equations.calculate_mu_o_Standing_method(api, Rs, P, Pb, T)

print("Above bubble point pressure hence Rs = GOR and since Rs for both is the same rho_o Would be the same as well.")
print("The density of oil is: %f lbm/ft^3" %rho_o)
print("The viscosity of oil at bubble point pressure is: %f cp" %mu_ob)
print("The viscosity of oil at %d psia is: %f cp" %(P, mu_o))