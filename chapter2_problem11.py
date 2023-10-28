import equations
gamma_g = 0.65 # specific gravity of gas
T = 80 # temperature in deg C
P1 = 5 # pressure in mPa
P2 = 10 # pressure in mPa
P3 = 50 # pressure in mPa

T = equations.C_to_F(T)
P1 = equations.mPa_to_psia(P1)
P2 = equations.mPa_to_psia(P2)
P3 = equations.mPa_to_psia(P3)

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

print("at P = 5 mPa, z for Brill and Beggs method = %f" %z1_1)
print("at P = 5 mPa, z for Hall Yarborough method = %f" %z1_2)
print("at P = 10 mPa, z for Brill and Beggs method = %f" %z2_1)
print("at P = 10 mPa, z for Hall Yarborough method = %f" %z2_2)
print("at P = 50 mPa, z for Brill and Beggs method = %f" %z3_1)
print("at P = 50 mPa, z for Hall Yarborough method = %f" %z3_2)
