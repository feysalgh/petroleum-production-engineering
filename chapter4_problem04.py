import equations
import matplotlib.pyplot as plt

L = 6000 # ft
D = 1.995 / 12 # ft
api = 30 # API
mu = 2 # cp
glr = 500 # scf/bbl
gamma_g = 0.65
p_head = 100 # psia
T_head = 80 # F
T_bot = 140 # F
liquid_rate = 1500 # bbl/day
water_cut = 0.2 # fraction
ift = 30 # dynes/cm
gamma_w = 1.05 

qw = liquid_rate * water_cut
gamma_o = equations.calculate_gamma_o_using_api(api)
gamma_l = ((liquid_rate - qw) * gamma_o + qw * gamma_w) / liquid_rate
qg = glr * liquid_rate
A = 3.141592653 * (D) ** 2 / 4
dz = 100
u_sl = liquid_rate * 5.615 / 86400 / A
mu_l = (mu * (liquid_rate - qw) + 0.5 * qw) / liquid_rate
m_t = gamma_l * 62.4 * liquid_rate * 5.615 + 0.0765 * gamma_g * qg


p , depth = equations.Hagedorn_Brown_Correlation(L, p_head, T_head, T_bot, dz, gamma_g, gamma_l, ift, D, m_t, qg, A, mu_l, u_sl)

plt.plot(p, depth)
plt.gca().invert_yaxis()
plt.xlabel('Pressure (psia)')
plt.ylabel('Depth (ft)')
plt.grid(True)
plt.ylim(max(depth), 0)
plt.show()

