import equations
import matplotlib.pyplot as plt

phi = 0.25
k = 10 # md
h = 50 # ft
P = 5000 # psia
pb = 100 # psia
bo = 1.2 # rb/stb
mu_o = 1.5 # cp
ct = 0.0000125 # psi^-1
rw = 0.328 # ft
S = 5
t = 30 * 24

J_star1 = equations.IPR_for_single_liquid_phase_radial_transient(k, h, phi, ct, bo, mu_o, rw, t, S)

pwf1 = list(range(P + 1))
q1 = []
for i in pwf1:
    q1.append(equations.vogels_equation_for_q(J_star1, P, pb, i))

J_star2 = equations.IPR_for_single_liquid_phase_radial_steady_state(k, h, bo, mu_o, rw, P, S)

pwf2 = list(range(P + 1))
q2 = []
for i in pwf2:
    q2.append(equations.vogels_equation_for_q(J_star2, P, pb, i))

J_star3 = equations.IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, bo, mu_o, rw, P, S)

pwf3 = list(range(P + 1))
q3 = []
for i in pwf3:
    q3.append(equations.vogels_equation_for_q(J_star3, P, pb, i))

fig, axs = plt.subplots(1, 3, figsize=(15,5))

axs[0].plot(q1, pwf1)
axs[0].set_xlabel('q (stb/d))')
axs[0].set_ylabel('Pwf (psi)')
axs[0].set_title('transient flow')

axs[1].plot(q2, pwf2)
axs[1].set_xlabel('q (stb/d))')
axs[1].set_ylabel('Pwf (psi)')
axs[1].set_title('steady state flow')

axs[2].plot(q3, pwf3)
axs[2].set_xlabel('q (stb/d))')
axs[2].set_ylabel('Pwf (psi)')
axs[2].set_title('pseudo steady state flow')

plt.show()
