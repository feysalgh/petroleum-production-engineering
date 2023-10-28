import equations
gamma_o = 0.8 # specific gravity of oil
T = 40 # temperature in C
Rs = 0 # Standard gas-oil ratio in scf/STB
gamma_g = 0 # Specific gravity of gas
T = equations.C_to_F(T)
rho_o = equations.calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)
print(rho_o, 'lbm/ft^3')