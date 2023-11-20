import equations
import matplotlib.pyplot as plt

P = 3500 # psia
pwf1 = 2500 # psia
q1 = 600 # stb/d
pwf2 = 1500 # psia
q2 = 900 # stb/d

qmax = equations.qmax_for_test_points_two_phase(q1, pwf1, P)
J_star = equations.J_star_using_qmax(qmax, P)

pwf = list(range(P + 1))
q_fetkovich = []
q_vogel = []
for i in pwf:
    q_fetkovich.append(equations.fetkovichs_equation(q1, pwf1, q2, pwf2, P, i))
    q_vogel.append(equations.vogels_equation_for_q(J_star, P, P, i))

plt.plot(q_fetkovich, pwf, label='fetkovich')
plt.plot(q_vogel, pwf, label='vogel')
plt.xlabel('q (stb/d)')
plt.ylabel('Pwf (psi)')
plt.legend()
plt.show()
