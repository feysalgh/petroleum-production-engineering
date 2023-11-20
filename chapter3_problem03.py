import equations
import matplotlib.pyplot as plt

phi = 0.25
k = 100 # md
h = 55 # ft
P = 5000 # psia
pb = 3000 # psia
Bo = 1.2 # rb/stb
mu_o = 1.8 # cp
ct = 0.000013 # psi^-1
A = 640 # acres
re = 2980 # ft
rw = 0.328 # ft
s = 5.5

J_star = equations.IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, Bo, mu_o, rw, re, s)
pwf = list(range(P + 1))
q = []
for i in pwf:
    q.append(equations.vogels_equation_for_q(J_star, P, pb, i))
plt.plot(q, pwf)
plt.xlabel('q (stb/d)')
plt.ylabel('Pwf (psi)')
plt.show()
