import equations
import matplotlib.pyplot as plt

pi = 3000 # psia
J_prime_0 = 4 * 10 ** - 4 # stb/d/psi^2
pe_1 = 2500 # psia
pe_2 = 2000 # psia
pe_3 = 1500 # psia
pe_4 = 1000 # psia

pwf_1 = list(range(pe_1 + 1))
pwf_2 = list(range(pe_2 + 1))
pwf_3 = list(range(pe_3 + 1))
pwf_4 = list(range(pe_4 + 1))

J_prime_1 = equations.J_prime(J_prime_0, pi, pe_1)
J_prime_2 = equations.J_prime(J_prime_0, pi, pe_2)
J_prime_3 = equations.J_prime(J_prime_0, pi, pe_3)
J_prime_4 = equations.J_prime(J_prime_0, pi, pe_4)

q_1 = []
q_2 = []
q_3 = []
q_4 = []

for i in pwf_1:
    q_1.append(equations.fetkovichs_method_for_future_well_performance(pe_1, i, J_prime_1))
for i in pwf_2:
    q_2.append(equations.fetkovichs_method_for_future_well_performance(pe_2, i, J_prime_2))
for i in pwf_3:
    q_3.append(equations.fetkovichs_method_for_future_well_performance(pe_3, i, J_prime_3))
for i in pwf_4:
    q_4.append(equations.fetkovichs_method_for_future_well_performance(pe_4, i, J_prime_4))


plt.plot(q_1, pwf_1, label='pe = 2500')
plt.plot(q_2, pwf_2, label='pe = 2000')
plt.plot(q_3, pwf_3, label='pe = 1500')
plt.plot(q_4, pwf_4, label='pe = 1000')
plt.xlabel('q (stb/d)')
plt.ylabel('Pwf (psi)')
plt.legend()
plt.show()