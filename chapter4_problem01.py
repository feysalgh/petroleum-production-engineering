import equations

q = 1000 # bbl/day
api = 16 # degree API
mu = 5 # cp
alpha = 3 # degree
L = 1000 # ft
d = 2.259 # in
D = d / 12 # ft
delta = 0.001

gamma_o = equations.calculate_gamma_o_using_api(api)
rho = equations.calculate_rho_using_gamma_o(gamma_o)
dz = equations.dz_inclinced_wellbore(L, alpha)
u = equations.u_in_tubing(q, D)
Re = equations.reynolds_number(rho, q, d, mu)
f_F = equations.friction_factor(Re, delta, d)
delta_p = equations.pressure_drop_in_single_phase_flow(rho, dz, 0, u, L, D, f_F)

print("Pressure drop over 1000 ft of tubing is {:.2f} psi".format(delta_p))
