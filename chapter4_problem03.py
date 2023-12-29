import equations

L = 8000 # ft
theta = 50 # degree
D =  1.995 / 12 # ft
qg = 500000 # scf/d
gamma_g = 0.75
qo = 2000 # stb/d
gamma_o = 0.85
qw = 500 # stb/d
gamma_w = 1.05
qs = 4 # ft3/d
gamma_s = 2.65
T_head = 100 # F
T_bottom = 170 # F
P_head = 500 # psia

T_avg = (T_head + T_bottom) / 2 + 460 # R
A = 3.14159 * (D * 12) ** 2 / 4
M = equations.mass_for_1_stb(gamma_o, gamma_w, gamma_g, qw / qo, qg / qo)
d_rho_v = equations.d_rho_v(M, D * 12, qo)
fm = equations.f2F(d_rho_v) * 4
p = equations.guo_ghalambor(theta, L, D, qg, gamma_g, qo, gamma_o, qw, gamma_w, qs, gamma_s, T_avg, A, fm, P_head)

print("The pressure at the bottom of the well is %.2f psia" % p)