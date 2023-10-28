import equations
api = 25 # API gravity
T = 100 # Temperature in F
Rs = 0 # Standard gas-oil ratio in scf/STB
gamma_g = 0 # Specific gravity of gas
gamma_o = equations.calculate_gamma_o_using_api(api)
rho_o = equations.calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)

print(rho_o, 'lbm/ft^3')