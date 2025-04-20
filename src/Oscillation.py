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

class P_3nu:
    def __init__(self, delta_m21_2=IOpara['delta_m21_2'], delta_m31_2=IOpara['delta_m3l_2'], s12_2=IOpara['s13_2'], s13_2=IOpara['s13_2']):
        # init the ocillation parameters
        self.delta_m21_2, self.delta_m31_2 = delta_m21_2, delta_m31_2
        self.s12_2, self.c12_2 = s12_2, 1 - s12_2
        self.s13_2, self.c13_2 = s13_2, 1 - s13_2
        self.s_2_12, self.s_2_13 = 2 * np.sqrt(self.s12_2 * self.c12_2), 2 * np.sqrt(self.s13_2 * self.c13_2)
        self.s_2_12, self.s_2_13 = self.c12_2 - self.s12_2, self.c13_2 - self.s13_2
        # default E [MeV], E/delta_m21_2 [MeV/eV^2]
        self.Es = np.concatenate([np.arange(0.1, 1, 0.1), np.arange(1, 101), np.arange(500, 10500, 500)])
        self.E_delta_m21_2 = 10**np.arange(3, 11, 0.1)

    def MatterPotential(self, n_e):
        # n_e is in mol with shape[,:,:]
        # Following A_CC, K and Sigma are divided by delta_m21_2
        A_CC = 2 * constant.sqrt2_G_F_N_A * n_e * self.E_delta_m21_2[:, np.newaxis]
        K = A_CC * self.c13_2 - self.c_2_12
        Sigma = self.s_2_12
        Delta_M_2 = np.sqrt(K**2 + Sigma**2)
        Theta_M = np.arccos(K / Delta_M_2) / 2
        return K, Delta_M_2, Theta_M

    def U_M(self, Theta_M):
        # The matrix transfer nu_e, nu_mu base to mass eigenstate
        # Eigenstate vector 2+2 [Re(nu_e),Im(nu_e), Re(nu_mu),Im(nu_mu)]
        # to mass eigenstate U_M.T @ P
        # to flavor eigenstate U_M @ P
        U_M = np.array([
            [np.cos(Theta_M), 0,  np.sin(Theta_M), 0], 
            [0, np.cos(Theta_M), 0,  np.sin(Theta_M)], 
            [-np.sin(Theta_M), 0, np.cos(Theta_M), 0],
            [0, -np.sin(Theta_M), 0, np.cos(Theta_M)],
            ])
        return U_M

    def H_M(self, dx, K, Delta_M_2):
        # H_M @ P
        delta_phase_1 = dx * (K + Delta_M_2) / 4
        delta_phase_2 = dx * (K - Delta_M_2) / 4
        H_M = np.array([
            [np.cos(delta_phase_1), -np.sin(delta_phase_1), 0, 0],
            [np.sin(delta_phase_1), np.cos(delta_phase_1), 0, 0],
            [0, 0, np.cos(delta_phase_2), -np.sin(delta_phase_2)],
            [0, 0, np.sin(delta_phase_2), np.cos(delta_phase_2)],
            ])
        return H_M
 
    def Evolute(self, v_k, U_M, H_M):
        # calculate the effective probability
        # equ2.10 - 2.12 in 10.1007/JHEP09(2013)152 without NSI
        # equ16 in arxiv 2203.11772 without NSI
        # calculate the vector in k+1 step
        v_k1 = U_M @ H_M @ U_M.T @ v_k
        return v_k1
       
    def P_e_3nu(self, v_k):
        # return the probability when current v_k is in vaccum
        return self.c_13_2**2 * (v_k[0]**2 + v_k[1]**2) + self.s_13_2**2

class CrossSection():
    # Equ.A1 in 10.1103/PhysRevD.51.6146
    # currently omit the QED correction in Appendix B in 10.1103/PhysRevD.51.6146
    def __init__(self):
        self.m_e = constant.m_e
        self.s_theta_W_2 = constant.sin_theta_W_2
        self.G_F = constant.G_F
        # due to the value is small, store the power value in self.power
        # -6 comes from G_F GeV->MeV
        # natural unit lead to the (hbar*c)**{-2}, -13 comes from fm->cm
        self.coeff = constant.hbar_c**2
        self.power = -6*2 - 13*2

    def I(self, T):
        x = np.sqrt(1 + 2 * self.m_e / T)
        return (1 / 3 + (3 - x**2) * (x * np.log((x + 1) / (x - 1)) / 2 - 1)) / 6

    def dif_T(self, T_e, E_nu, g_1, g_2):
        z = T_e / E_nu
        T_max = self.T_max(E_nu)
        sigma = 2 * self.G_F**2 * self.m_e / np.pi * (
                g_1**2 + g_2**2 * (1 - z**2) - g_1 * g_2 * self.m_e * z / E_nu
                )
        sigma[T_e>T_max] = 0
        return sigma * self.coeff

    def dif_T_nu_e(self, T_e, E_nu):
        # nu_e
        rho_NC = 1.0126
        kappa = 0.9791 + 0.0097 * self.I(T_e)
        g_1, g_2 = 1 - rho_NC * (0.5 - kappa * self.s_theta_W_2), rho_NC * kappa * self.s_theta_W_2
        return self.dif_T(T_e, E_nu, g_1, g_2)

    def dif_T_nu_mu(self, T_e, E_nu):
        # nu_mu
        rho_NC = 1.0126
        kappa = 0.9970 - 0.00037 * self.I(T)
        g_1, g_2 = rho_NC * (-0.5 + kappa * self.s_theta_W_2), rho_NC * kappa * self.s_theta_W_2
        return self.dif_T(T_e, E_nu, g_1, g_2)

    def T_c_theta(self, T):
        # Equ 5.27 in fudamental neutrino physics
        c_theta_2 = c_theta**2 
        E_nu_2 = E_nu**2
        return np.sqrt(
                T * (self.m_e + E_nu)**2 / E_nu_2 / (T + 2 * self.m_e)
                )

    def c_theta_T(self, E_nu, c_theta):
        # Equ 5.27 in fudamental neutrino physics
        c_theta_2 = c_theta**2 
        E_nu_2 = E_nu**2
        return 2 * self.m_e * E_nu_2 * c_theta_2 / ((self.m_e + E_nu)**2 - E_nu_2 * c_theta_2)

    def T_max(self, E_nu):
        return self.c_theta_T(E_nu, -1)

