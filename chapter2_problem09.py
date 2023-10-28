import equations

y = [0.755, 0.073, 0.011, 0.006, 0.002, 0.003, 0.008, 0.001, 0.001, 0.07, 0.05, 0.02] # Mole fraction of each component
mWi = [16.04, 30.07, 44.10, 58.12, 58.12, 72.15, 72.15, 86.18, 114.23, 28.02, 44.01, 34.08] # Molecular weight of each component

mWa = sum([yi*mwi for yi, mwi in zip(y, mWi)])
gamma_g = equations.calculate_gamma_g_using_apparent_molecular_weight(mWa)
Ppc = equations.calculate_Ppc_ahmed_method(gamma_g, y[9], y[10], y[11])
Tpc = equations.calculate_Tpc_ahmed_method(gamma_g, y[9], y[10], y[11])

print("The apparent molecular weight of the gas is: %f" %mWa)
print("The specific gravity of the gas is: %f" %gamma_g)
print("The pseudo critical pressure of the gas is: %f psia" %Ppc)
print("The pseudo critical temperature of the gas is: %f R" %Tpc)
