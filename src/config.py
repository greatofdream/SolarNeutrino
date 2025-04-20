import matplotlib.colors as mcolors
from scipy.constants import physical_constants
# spectra line style
spectras = ['N13', 'O15', 'F17', 'pp', 'B8', 'hep', 'Be7', 'pep']
labels = ['$^{13}$N', '$^{15}$O', '$^{17}$F', 'pp', '$^{8}$B', 'hep', '$^{7}$Be', 'pep']
errors = [15, 17, 20, 0.6, 12, 30, 6, 1]
xs = [0.1, 0.1, 0.1, 0.3, 1, 0.5, 1, 1]
ys = [3, 0.2, 0.1, 0.2, 1.8, 2, 1.1, 1.2]
rotations = [10, 10, 10, 10, 0, 0, 0, 0]
colors = list(mcolors.TABLEAU_COLORS.keys())
spectraConfig = dict([(s, {'label': label, 'color': c, 'x': x, 'y':y, 'rotation': rot, 'error': err}) for s, label, c, x, y, rot, err in zip(spectras, labels, colors, xs, ys, rotations, errors)])
## append discrete spectra

# neutrino oscillation parameters
## deltam_2 eV^2
oscParas = dict([
    ('nufit5.2NO_woSK', {'s12_2': 0.303, 's13_2': 0.02203, 'delta_m21_2': 7.41*1e-5, 'delta_m3l_2': 2.511*1e-3}),
    ('nufit5.2IO_woSK', {'s12_2': 0.303, 's13_2': 0.02219, 'delta_m21_2': 7.41*1e-5, 'delta_m3l_2': -2.498*1e-3}),

])
class constant():
    # 10^{-6} eV cm^3, 10^{-6} will be eliminated by MeV in E_nu
    # reference: fundamental neutrino physics and astrophysics
    sqrt2_G_F_N_A = 7.63*1e-8
    m_e = physical_constants['electron mass energy equivalent in MeV'][0]
    sin_theta_W_2 = physical_constants['weak mixing angle'][0]
    G_F = physical_constants['Fermi coupling constant'][0]
    hbar_c = physical_constants['reduced Planck constant times c in MeV fm'][0]
