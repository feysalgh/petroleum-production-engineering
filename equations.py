import math, scipy
def calculate_api(gamma_o):
    api = 141.5/gamma_o - 131.5
    return api
def calculate_Rs(gamma_g, api, p, t):
    Rs = gamma_g * (p * 10**(0.0125*(api))/(18 * 10**(0.00091*t))**1.2048)
    return Rs

def calculate_gamma_o_using_api(api):
    gamma_o = (141.5/(api + 131.5))
    return gamma_o

def calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, t):
    rho_o = (62.4 * gamma_o + 0.0136 * Rs * gamma_g)/(0.972 + 0.000147 * (Rs * ((gamma_g/gamma_o) ** 0.5) + 1.25 * t) ** 1.175)
    return rho_o

def calculate_Bo_Standing_method(Rs, gamma_o, gamma_g, t):
    Bo = 0.9759 + 0.00012 * (Rs * ((gamma_g/gamma_o) ** 0.5) + 1.25 * t) ** 1.2
    return Bo

def calculate_mu_o_Standing_method(api, Rs, p, pb, t):
    A = 10 ** (0.43 + (8.33 / (api)))
    mu_od = (0.32 + 1.8 * 10 ** 7 / (api ** 4.53)) * (360/(t + 200)) ** A
    c = 8.62 * 10 ** (-5) * Rs
    d = 1.1 * 10 ** (-3) * Rs
    e = 3.74 * 10 ** (-3) * Rs
    b = 0.68 / 10 ** c + 0.25 / 10 ** d + 0.062 / 10 ** e
    a = Rs * (2.2 * 10 ** (-7) * Rs - 7.4 * 10 ** (-4))
    mu_ob = 10 ** a * mu_od ** b
    mu_o = mu_ob + 0.001 * (p - pb) * (0.024 * mu_ob ** 1.6 + 0.38 * mu_ob ** 0.56)
    return mu_o, mu_ob

def calculate_gamma_g_using_apparent_molecular_weight(apparent_molecular_weight):
    gamma_g = apparent_molecular_weight / 28.97
    return gamma_g

def calculate_apparent_molecular_weight():
    n = input("Enter the number of components: ")
    for i in range(n):
        y = input("Enter the mole fraction of component %d: " %i)
        m = input("Enter the molecular weight of component %d: " %i)
        apparent_molecular_weight += y * m
    return apparent_molecular_weight

def calculate_Ppc_for_composition():
    n = input("Enter the number of components: ")
    for i in range(n):
        y = input("Enter the mole fraction of component %d: " %i)
        Ppc = input("Enter the critical pressure of component %d: " %i)
        Ppc_for_composition += y * Ppc
    return Ppc_for_composition

def calculate_Tpc_for_composition():
    n = input("Enter the number of components: ")
    for i in range(n):
        y = input("Enter the mole fraction of component %d: " %i)
        Tpc = input("Enter the critical temperature of component %d: " %i)
        Tpc_for_composition += y * Tpc
    return Tpc_for_composition

def calculate_Ppc_ahmed_method(gamma_g, yN2, yCO2, yH2S):
    Ppc = 678 - 50 * (gamma_g - 0.5) - 206.7 * yN2 + 440 * yCO2 + 606.7 * yH2S
    return Ppc

def calculate_Tpc_ahmed_method(gamma_g, yN2, yCO2, yH2S):
    Tpc = 326 + 315.7 * (gamma_g - 0.5) - 240 * yN2 - 83.3 * yCO2 + 133.3 * yH2S
    return Tpc

def calculate_mu_carr_et_al_method(gamma_g, yN2, yCO2, yH2S, t, Ppr, Tpr):
    mu1_hc = 8.188 * 10 ** (-3) - 6.15 * 10 ** (-3) * math.log10(gamma_g) + (1.709 * 10 ** (-5) - 2.062 * 10 ** (-6) * gamma_g) * t
    mu1_n2 = (9.59 * 10 **(-3) + 8.48 * 10 ** (-3) * math.log(gamma_g)) * yN2
    mu1_c2 = (6.24 * 10 ** (-3) + 9.08 * 10 ** (-3) * math.log(gamma_g)) * yCO2
    mu1_h2s = (3.73 * 10 ** (-3) + 8.49 * 10 ** (-3) * math.log(gamma_g)) * yH2S
    mu1 = mu1_hc + mu1_n2 + mu1_c2 + mu1_h2s
    a0 = -2.46211820
    a1 = 2.97054714
    a2 = -0.28626405
    a3 = 0.00805420
    a4 = 2.80860949
    a5 = -3.49803305
    a6 = 0.36037302
    a7 = -0.01044324
    a8 = -0.79338568
    a9 = 1.39643306
    a10 = -0.14914493
    a11 = 0.00441016
    a12 = 0.08393872
    a13 = -0.18640885
    a14 = 0.02033679
    a15 = -0.00060958
    mu_r = a0 + a1 * Ppr + a2 * Ppr ** 2 + a3 * Ppr ** 3 + Tpr *(a4 + a5 * Ppr + a6 * Ppr ** 2 + a7 * Ppr ** 3) + Tpr ** 2 * (a8 + a9 * Ppr + a10 * Ppr ** 2 + a11 * Ppr ** 3) + Tpr ** 3 * (a12 + a13 * Ppr + a14 * Ppr ** 2 + a15 * Ppr ** 3)
    mu_g = mu1 * math.exp(mu_r) / Tpr
    return mu_g

def calculate_zfactor_brill_and_beggs_method(Tpr, Ppr):
    F = 0.3106 - 0.49 * Tpr + 0.1824 * Tpr ** 2
    E = 9 * (Tpr -1)
    D = 10 ** F
    C = 0.132 - 0.32 * math.log10(Tpr)
    B = (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr ** 2 + 0.32 * Ppr ** 6 / 10 ** E
    A = 1.39 * (Tpr - 0.92) ** 0.5 - 0.36 * Tpr - 0.10
    z = A + (1 - A) / math.exp(B) + C * Ppr ** D
    return z

def calculate_zfactor_hall_yarborough_method(Tpr, Ppr):
    tr = 1 / Tpr
    A = 0.06125 * tr * math.exp(-1.2 * (1 - tr) ** 2)
    B = tr * (14.76 - 9.76 * tr + 4.58 * tr ** 2)
    C = tr * (90.7 - 242.2 * tr + 42.4 * tr ** 2)
    D = 2.18 + 2.82 * tr
    def f(Y):
        return (Y + Y ** 2 + Y ** 3 - Y ** 4) / (1 - Y) ** 3 - A * Ppr - B * Y ** 2 + C * Y ** D
    sol = scipy.optimize.root_scalar(f, bracket=[0, 0.99999])
    Y = sol.root
    z = A * Ppr / Y
    return z


def calculate_rho_g(gamma_g, p, z, T):
    rho_g = 2.7 * gamma_g * p / (z * T)
    return rho_g

def calculate_Bg(z, T, p):
    Bg = 0.0283 * z * T / p
    return Bg

def C_to_F(C):
    F = (C * 9/5) + 32
    return F

def mPa_to_psia(mPa):
    psia = mPa * 145.038
    return psia

