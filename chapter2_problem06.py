import equations
import matplotlib.pyplot as plt

gamma_g = 0.65 # specific gravity of gas
T = 250 # temperature in F
Ppc = equations.calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Ppr = [P / Ppc for P in range(14, 8000)]

plt.plot(range(14, 8000), Ppr)
plt.xlabel('Pressure (psi)')
plt.ylabel('Ppr')
plt.title('Pseudo-Reduced Pressure vs Pressure')
plt.show()
