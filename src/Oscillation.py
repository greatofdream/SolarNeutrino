import numpy as np
from config import oscParas
from config import constant
IOpara = oscParas['nufit5.2IO_woSK']
def V_e(n_es):
    return constant.sqrt2_G_F_N_A * n_es
def MSW_Adiabatic(Es, V_es, delta_m21_2=IOpara['delta_m21_2'], delta_m3l_2=IOpara['delta_m3l_2'], s12_2=IOpara['s12_2'], s13_2=IOpara['s13_2'],):
    # equ2.8 in 10.1016/j.ppnp.2023.104043
    c12_2, c13_2 = 1 - s12_2, 1 - s13_2
    Es, V_es = Es.reshape((-1,1)), V_es.reshape((1,-1))
    V_eEs = V_es * Es
    beta12 = 2 * c13_2 * V_eEs / delta_m21_2
    beta13 = 2 * V_eEs / delta_m3l_2
    s13m_2 = s13_2 * (1 + 2 * beta13)
    c2_12 = 1 - 2*s12_2
    c2_12m = (c2_12 - beta12)/np.sqrt(1 + beta12**2 - 2*c2_12*beta12)
    return c13_2*(1-s13m_2)*(0.5+0.5*c2_12m*c2_12) + s13_2 * s13m_2
