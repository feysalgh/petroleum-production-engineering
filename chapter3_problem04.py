import equations
import matplotlib.pyplot as plt

P = 5500 # psia
pb = 3500 # psia
pwf1_A = 4000 # psia
q1_A = 400 # stb/d
pwf1_B = 2000 # psia
q1_B = 1000 # stb/d

J_star_A = equations.IPR_using_test_points(P, pb, pwf1_A, q1_A)
J_star_B = equations.IPR_using_test_points(P, pb, pwf1_B, q1_B)

pwf = list(range(P + 1))
q_A = []
q_B = []
for i in pwf:
    q_A.append(equations.vogels_equation_for_q(J_star_A, P, pb, i))
    q_B.append(equations.vogels_equation_for_q(J_star_B, P, pb, i))

plt.plot(q_A, pwf, label='well A')
plt.plot(q_B, pwf, label='well B')
plt.xlabel('q (stb/d))')
plt.ylabel('Pwf (psi)')
plt.legend()
plt.show()
