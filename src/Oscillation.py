import numpy as np
from config import oscParas
from config import constant
from tqdm import tqdm
import h5py
import concurrent.futures

IOpara = oscParas['nufit5.2IO_woSK']
NOpara = oscParas['nufit5.2NO_woSK']
def V_e(n_es):
    # Attention, there is a 1E-6 factor in the return value
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

def interpolate(ys, v_x):
    # v_x is the interpolating xs
    v_x_f, v_x_i = modf(v_x)
    return ys[v_x_i] * (1 - v_x_f) + ys[v_x_i + 1] * v_x_f

class P_3nu:
    def __init__(self, s12_2=IOpara['s13_2'], s13_2=IOpara['s13_2'], E_Delta_m21_2=1 / IOpara['delta_m21_2']):
        # init the ocillation parameters
        # E_Delta_m21_2 [MeV/eV^2]
        self.s12_2, self.c12_2 = s12_2, 1 - s12_2
        self.s13_2, self.c13_2 = s13_2, 1 - s13_2
        self.s_2_12, self.s_2_13 = 2 * np.sqrt(self.s12_2 * self.c12_2), 2 * np.sqrt(self.s13_2 * self.c13_2)
        self.c_2_12, self.c_2_13 = self.c12_2 - self.s12_2, self.c13_2 - self.s13_2
        # default E [MeV], E/delta_m21_2 [MeV/eV^2], self.E_Delta_m12_2 [eV^-1]
        '''
        self.Es = np.concatenate([np.arange(0.1, 1, 0.1), np.arange(1, 101), np.arange(500, 10500, 500)])
        self.E_delta_m21_2 = 10**np.arange(3, 11, 0.1)
        '''
        self.E_Delta_m12_2 = E_Delta_m21_2 * 1E6
        # The 1E-45 comes from the eV unit transfer; Final unit is eV cm^3
        self.sqrt2_G_F_N_A = np.sqrt(2) * constant.G_F * constant.hbar_c**3 * (constant.N_A * 1E-39)

    def MatterPotential(self, n_e):
        # n_e is in mol with shape[,:,:]
        # Following A_CC, K, Sigma and Delta_M_2 are divide by E
        # A_CC [eV]
        A_CC = 2 * self.sqrt2_G_F_N_A * n_e
        K = self.c_2_12[np.newaxis, np.newaxis, :, np.newaxis] - A_CC[:, np.newaxis, np.newaxis, np.newaxis] * self.c13_2[np.newaxis, np.newaxis, np.newaxis, :] * self.E_Delta_m12_2[np.newaxis, :, np.newaxis, np.newaxis]
        Sigma = self.s_2_12[np.newaxis, np.newaxis, :, np.newaxis]
        Delta_M_2 = np.sqrt(K**2 + Sigma**2)
        # arccos [0, pi]
        Theta_M = np.arccos(K / Delta_M_2) / 2
        # output shape is [n_e, E_delta_m21_2, s_12, s_13]
        # the unit is eV
        return K / self.E_Delta_m12_2[np.newaxis, :, np.newaxis, np.newaxis], Delta_M_2 / self.E_Delta_m12_2[np.newaxis, :, np.newaxis, np.newaxis], Theta_M

    def U_M(self, Theta_M):
        # The matrix transfer nu_e, nu_mu base to mass eigenstate
        # Eigenstate vector 2+2 [Re(nu_e),Im(nu_e), Re(nu_mu),Im(nu_mu)]
        # to mass eigenstate U_M @ P
        # to flavor eigenstate U_M.T @ P
        U_Ms = np.zeros((*(Theta_M.shape), 4, 4))
        U_Ms[..., 0, 0] = U_Ms[..., 1, 1] = U_Ms[..., 2, 2] = U_Ms[..., 3, 3] = np.cos(Theta_M)
        U_Ms[..., 2, 0] = U_Ms[..., 3, 1] = -np.sin(Theta_M)
        U_Ms[..., 0, 2] = U_Ms[..., 1, 3] = -U_Ms[..., 2, 0]
        # output shape is [n_e, E_delta_m21_2, s_12, s_13, 4, 4]
        return U_Ms

    def H_M(self, dx, K, Delta_M_2):
        '''
        H_M @ P
        np.array([
            [np.cos(delta_phase_1), -np.sin(delta_phase_1), 0, 0],
            [np.sin(delta_phase_1), np.cos(delta_phase_1), 0, 0],
            [0, 0, np.cos(delta_phase_2), -np.sin(delta_phase_2)],
            [0, 0, np.sin(delta_phase_2), np.cos(delta_phase_2)],
            ])
        '''
        delta_phase_1 = dx * (K + Delta_M_2) / 4
        delta_phase_2 = dx * (K - Delta_M_2) / 4
        # matmal use the last 2 axis dot
        H_Ms = np.zeros((*(Delta_M_2.shape), 4, 4))
        H_Ms[..., 0, 0] = H_Ms[..., 1, 1] = np.cos(delta_phase_1)
        H_Ms[..., 1, 0] = np.sin(delta_phase_1)
        H_Ms[..., 0, 1] = -H_Ms[..., 1, 0]
        H_Ms[..., 2, 2] = H_Ms[..., 3, 3] = np.cos(delta_phase_2)
        H_Ms[..., 3, 2] = np.sin(delta_phase_2)
        H_Ms[..., 2, 3] = -H_Ms[..., 3, 2]
        # output shape is [n_e, E_delta_m21_2, s_12, s_13, 4, 4]
        return H_Ms
 
    def Evolute(self, v_0, U_Ms, H_Ms):
        # calculate the effective probability
        # equ2.10 - 2.12 in 10.1007/JHEP09(2013)152 without NSI
        # equ16 in arxiv 2203.11772 without NSI
        # calculate the vector in k+1 step
        U_f = np.eye(4)
        for i in range(U_Ms.shape[0]):
            U_M, H_M = U_Ms[i], H_Ms[i]
            U_f = np.matmul(np.matmul(np.matmul(U_M, H_M), np.transpose(U_M, [0, 1, 2, 4, 3])), U_f)
        # output shape is [E_delta_m21_2, s_12, s_13, 4]
        v_k = U_f @ v_0
        v_mass = np.matmul(np.transpose(U_M, [0, 1, 2, 4, 3]), U_f) @ v_0
        return v_k, v_mass

    def P_e_2nu(self, v_k, v_mass):
        # v_k shape is [E_delta_m21_2, s_12, s_13, 4]
        # U_M shape is [E_delta_m21_2, s_12, s_13, 4, 4]
        # Ks, Delta_M_2s, Theta_Ms = self.MatterPotential(np.array([0]))
        # U_M = self.U_M(Theta_Ms)[0]
        # return the probability when current v_k is in vaccum
        # v_mass shape is [E_delta_m21_2, s_12, s_13, 4]
        # v_mass = np.matmul(np.transpose(U_M, [0, 1, 2, 4, 3]), v_k[:, :, :, :, np.newaxis])[..., 0]
        # output shape is [4, E_delta_m21_2, s_12, s_13]
        P2 = v_mass[..., 2]**2 + v_mass[..., 3]**2
        Peff = self.c12_2[np.newaxis, :, np.newaxis] - P2 * (2*self.c12_2 - 1)[np.newaxis, :, np.newaxis]
        return np.stack([v_k[..., 0]**2 + v_k[..., 1]**2, v_k[..., 2]**2 + v_k[..., 3]**2, P2, Peff])

    def P_e_3nu(self, p_2nu):
        # return the probability when current v_k is in vaccum
        return self.c13_2**2 * p_2nu + self.s13_2**2

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
        # 1E5 used for adjusting the final unit as 1E-43, suitable for drawing figure
        self.coeff = constant.hbar_c**2 * 1E5
        self.power = -6*2 - 13*2 - 5

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

    def dif_T_nu_e(self, T_e, E_nu, corr=True):
        # nu_e
        if corr:
            rho_NC = 1.0126
            kappa = 0.9791 + 0.0097 * self.I(T_e)
        else:
            rho_NC = 1
            kappa = 1
        g_1, g_2 = 1 - rho_NC * (0.5 - kappa * self.s_theta_W_2), rho_NC * kappa * self.s_theta_W_2
        return self.dif_T(T_e, E_nu, g_1, g_2)

    def dif_T_nu_mu(self, T_e, E_nu, corr=True):
        # nu_mu
        if corr:
            rho_NC = 1.0126
            kappa = 0.9970 - 0.00037 * self.I(T_e)
        else:
            rho_NC = 1
            kappa = 1
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

    def c_theta_E(self, T, c_theta):
        return self.m_e / (np.sqrt(2 * self.m_e / T + 1) * c_theta - 1)

    def c_theta_min(self, T):
        return 1 / np.sqrt(2 * self.m_e / T_e + 1)

class MSW():
    def __init__(self, radius, rho_e, R):
        self.model = {'Radius': radius, 'rho_e': rho_e}
        # the unit should transfer from cm to nm to eV^-1
        self.R = R * 1E7 / constant.hbar_c

    def setOsc(self, p_3nu):
        # This is necessary to load an oscillation calculator
        self.p_3nu = p_3nu
    
    def evoluteParallel(self, r, c_phi, nu_0):
        p_y = np.sqrt(1 - c_phi**2) * r
        p_x_max = np.sqrt(1 - p_y**2)
        p_xs = np.arange(r * c_phi, p_x_max, self.delta_x)
        p_rs = np.sqrt(p_xs**2 + r**2 * (1 - c_phi**2))
        n_es = np.interp(p_rs, self.model['Radius'], self.model['rho_e'])
        nu_k, nu_mass = self.evolutePath(n_es, nu_0)
        return self.p_3nu.P_e_2nu(nu_k, nu_mass)

    def evoluteCos(self, cos_phi):
        L_half = np.sqrt(1 - cos_phi**2)
        p_l = np.arange(0, 2 * L_half, self.delta_x)
        p_rs = np.sqrt(p_l - L_half)

    def evolutePath(self, n_es, nu_0):
        # calculate the evolution along a specific path with the electron density n_e
        Ks, Delta_M_2s, Theta_Ms = self.p_3nu.MatterPotential(n_es)
        U_Ms = self.p_3nu.U_M(Theta_Ms)
        H_Ms = self.p_3nu.H_M(self.delta_x * self.R, Ks, Delta_M_2s)
        nu_k, nu_mass = self.p_3nu.Evolute(nu_0, U_Ms, H_Ms)
        return nu_k, nu_mass

class Sun(MSW):
    # MSW effect in the sun
    def __init__(self, ssm, flux):
        # the radius of sun is 6.96*1E10 cm
        MSW.__init__(self, ssm['Radius'], ssm['rho_e'], 6.96*1E10)
        self.flux = flux
        # load the model raw data from the output of SSMLoader
        self.ssm, self.flux = ssm, flux
        # The raidus of sun is 1
        # The sampling radius; default is 16 points
        self.coreR_s, self.coreR_e, self.coreR_nbin = 0, 0.32, 16
        self.coreR_bins = np.linspace(self.coreR_s, self.coreR_e, self.coreR_nbin + 1)
        self.coreR_c = (self.coreR_bins[:-1] + self.coreR_bins[1:]) / 2
        # The sampling angle; default is 20 directions
        self.coreCPhi_s, self.coreCPhi_e, self.coreCPhi_nbin = -1, 1, 21
        self.coreCPhi_bins = np.linspace(self.coreCPhi_s, self.coreCPhi_e, self.coreCPhi_nbin)
        # path step, normalized as R=1
        self.delta_x = 2E-4
        # thread pool
        # self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=32)
        # self.executor = concurrent.futures.ProcessPoolExecutor(max_workers=32)

    def merge(self, fs):
        probs_2nu = np.empty((self.coreR_nbin, 4, len(self.p_3nu.E_Delta_m12_2), len(self.p_3nu.s12_2), len(self.p_3nu.s13_2)))
        for i, f in enumerate(fs):
            with h5py.File(f, 'r') as ipt:
                probs_2nu[i] = ipt['sun/prob_2nu'][:]
        return probs_2nu

    def getProbR(self, probs_2nu):
        return np.mean(probs_2nu, axis=0)

    def getProb(self, reaction, probs_2nu_R):
        flux = np.interp(self.coreR_c, self.flux['R'], self.flux[reaction])
        return np.sum(flux[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis] * probs_2nu_R, axis=0) / np.sum(flux)

    def evolute(self, nu_0):
        probs_2nu = np.empty((self.coreR_nbin, 4, len(self.p_3nu.E_Delta_m12_2), len(self.p_3nu.s12_2), len(self.p_3nu.s13_2)))
        for i, r in tqdm(enumerate(self.coreR_c)):
            probs_2nu[i] = self.evoluteR(self, nu_0, r)
        # output shape: 4, len(self.p_3nu.E_Delta_m12_2), len(self.p_3nu.s12_2), len(self.p_3nu.s13_2)
        return probs_2nu

    def evoluteR(self, nu_0, r):
        # probs_2nu = np.empty((self.coreR_nbin, self.coreCPhi_nbin), dtype=[('nu_e', np.float64), ('nu_mu', np.float64), ('nu_1', np.float64), ('nu_2', np.float64)])
        probs_2nu = np.empty((self.coreCPhi_nbin, 4, len(self.p_3nu.E_Delta_m12_2), len(self.p_3nu.s12_2), len(self.p_3nu.s13_2)))
        # sample the start point
        # for j, prob in enumerate(self.executor.map(self.evoluteThread, np.repeat(r, self.coreCPhi_nbin), self.coreCPhi_bins, np.tile(nu_0, (self.coreCPhi_nbin, 1)))):
        for j, c_phi in enumerate(self.coreCPhi_bins):
            # sample the points (p_xs, p_y) in the different path
            probs_2nu[j] = self.evoluteParallel(r, c_phi, nu_0)
        return self.getProbR(probs_2nu)

class Earth(MSW):
    # MSW effect in Earth
    def __init__(self, prem):
        MSW.__init__(self, prem['Radius'], prem['rho_e'])
        # load PREM
        self.prem = prem
        # The raidus of earth is 1
        # The sampling angle; default is 100 directions
        self.c_theta_s, self.c_theta_e, self.c_theta_nbin = 0, 1, 1000
        self.c_theta_bins = np.linspace(self.c_theta_s, self.c_theta_e, self.c_theta_nbin + 1)
        self.c_theta_c = (self.c_theta_bins[:-1] + self.c_theta_bins[1:]) / 2

        self.delta_x = 2E-4

    def evolute(self):
        probs_2nu = np.empty((self.c_theta_nbin, 4, len(self.p_3nu.E_Delta_m12_2), len(self.p_3nu.s12_2), len(self.p_3nu.s13_2)))
        # sample the angle
        for i, c_theta in enumerate(self.c_theta_c):
            n_es = np.interp(p_rs, self.ssm['Radius'], self.ssm['rho_e'])
            nu_k = self.evolutePath(n_es, nu_0)
            probs_2nu[i, j] = self.p_3nu.P_e_2nu(nu_k)
        return probs_2nu

    def getProb(self, reaction, probs_2nu):
        flux = interpolate(self.flux[reaction], self.coreR_c)
        return np.sum(flux * np.mean(probs_2nu, axis=1)) / np.sum(flux)

    def evolutePath(self, n_es, nu_0):
        # calculate the evolution along a specific path with the electron density n_e
        Ks, Delta_M_2s, Theta_Ms = self.p_3nu.MatterPotential(n_es)
        U_Ms = self.p_3nu.U_M(Theta_Ms)
        H_Ms = self.p_3nu.H_M(self.delta_x, Ks, Delta_M_2s)
        nu_k = self.p_3nu.Evolute(nu_0, U_Ms, H_Ms)
        return nu_k

class Livetime():
    # The livetime distribution
    def __init__(self, latitude):
        pass

