import equations
import matplotlib.pyplot as plt

Pp = 2200 # psia
Pf = 1500 # psia
J_star_p = 1.25 # stb/d/psi
mu_o_p = 3.55 # cp
mu_o_f = 3.85 # cp
Bo_p = 1.2 # rb/stb
Bo_f = 1.15 # rb/stb
kr_p = 0.82 # md
kr_f = 0.65 # md

J_star_f = equations.future_IPR_vogels_method(J_star_p, Bo_p, mu_o_p, kr_p, Bo_f, mu_o_f, kr_f)

pwf_p = list(range(Pp + 1))
pwf_f = list(range(Pf + 1))
q_p = []
q_f = []

for i in pwf_p:
    q_p.append(equations.vogels_equation_for_q(J_star_p, Pp, Pp, i))
for i in pwf_f:
    q_f.append(equations.vogels_equation_for_q_future(J_star_f, Pf, i))

plt.plot(q_p, pwf_p, label='present')
plt.plot(q_f, pwf_f, label='future')
plt.xlabel('q (stb/d)')
plt.ylabel('Pwf (psi)')
plt.legend()
plt.show()