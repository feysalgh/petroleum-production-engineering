import equations

L = 7000 # ft
theta = 20 # degree
D =  1.995 / 12 # ft
qg = 1000000 # scf/d
gamma_g = 0.7
qo = 1000 # stb/d
gamma_o = 0.85
qw = 300 # stb/d
gamma_w = 1.05
qs = 1 # ft3/d
gamma_s = 2.65
T_head = 100 # F
T_bottom = 224 # F
P_head = 300 # psia

T_avg = (T_head + T_bottom) / 2 + 460 # R
A = 3.14159 * (D * 12) ** 2 / 4
M = equations.mass_for_1_stb(gamma_o, gamma_w, gamma_g, qw / qo, qg / qo)
d_rho_v = equations.d_rho_v(M, D * 12, qo)
fm = equations.f2F(d_rho_v) * 4
p = equations.guo_ghalambor(theta, L, D, qg, gamma_g, qo, gamma_o, qw, gamma_w, qs, gamma_s, T_avg, A, fm, P_head)

print("The pressure at the bottom of the well is %.2f psia" % p)