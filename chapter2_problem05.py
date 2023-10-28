import equations
gamma_g = 0.65 # specific gravity of gas
T = 150 # temperature in deg F
P1 = 50 # pressure in psia
P2 = 500 # pressure in psia
P3 = 5000 # pressure in psia
Ppc = equations.calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc = equations.calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr = (T + 459.67) / Tpc
Ppr1 = P1 / Ppc
Ppr2 = P2 / Ppc
Ppr3 = P3 / Ppc

z1_1 = equations.calculate_zfactor_brill_and_beggs_method(Tpr, Ppr1)
z1_2 = equations.calculate_zfactor_hall_yarborough_method(Tpr, Ppr1)

z2_1 = equations.calculate_zfactor_brill_and_beggs_method(Tpr, Ppr2)
z2_2 = equations.calculate_zfactor_hall_yarborough_method(Tpr, Ppr2)

z3_1 = equations.calculate_zfactor_brill_and_beggs_method(Tpr, Ppr3)
z3_2 = equations.calculate_zfactor_hall_yarborough_method(Tpr, Ppr3)

print("at P = 50 psia, z for Brill and Beggs method = %f" %z1_1)
print("at P = 50 psia, z for Hall Yarborough method = %f" %z1_2)
print("at P = 500 psia, z for Brill and Beggs method = %f" %z2_1)
print("at P = 500 psia, z for Hall Yarborough method = %f" %z2_2)
print("at P = 5000 psia, z for Brill and Beggs method = %f" %z3_1)
print("at P = 5000 psia, z for Hall Yarborough method = %f" %z3_2)
