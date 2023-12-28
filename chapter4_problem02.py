import equations

D = 1.66
P_wellhead = 300 # psia
production_rate = 2000 # STB/day
gas_liquid_ratio = 800 # scf/STB
wc = 0.30
api = 40 # api degree
gamma_w = 1.05
gamma_g = 0.7
bw = 1.2 # bbl/STB
wellhead_temperature = 100 # F
L = 8000 # ft
bottomhole_temperature = 170 # F

gor = gas_liquid_ratio / (1 - wc)
wor = wc / (1 - wc)
qo = production_rate * (1 - wc)

gamma_o = equations.calculate_gamma_o_using_api(api)
M = equations.mass_for_1_stb(gamma_o, gamma_w, gamma_g, wor, gor)
Rs_wellhead = equations.calculate_Rs(gamma_g, api, P_wellhead, wellhead_temperature)
Bo_wellhead = equations.calculate_Bo_Standing_method(Rs_wellhead, gamma_o, gamma_g, wellhead_temperature)
Ppc = equations.calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc = equations.calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr = (wellhead_temperature + 460) / Tpc
Ppr = P_wellhead / Ppc
z_wellhead = equations.calculate_zfactor_brill_and_beggs_method(Tpr, Ppr)
V_wellhead = equations.volume_for_1_stb(Bo_wellhead, wor, bw, gor, Rs_wellhead, P_wellhead, wellhead_temperature, z_wellhead)
rho_wellhead = M / V_wellhead
d_rho_v = equations.d_rho_v(M, D, qo)
f2F = equations.f2F(d_rho_v)
k = equations.friction_term(f2F, qo, M, D / 12)
bhp = equations.finding_bhp(P_wellhead, rho_wellhead, k, L, api, bottomhole_temperature, gamma_g, wor, bw, gor, M)

print("The bottomhole using Poettmann-Carpenter method is: {:.2f} psia".format(bhp))
    