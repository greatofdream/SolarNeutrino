from Oscillation import CrossSection
import argparse
import numpy as np
import h5py

import matplotlib.pyplot as plt
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import ScalarFormatter, LogFormatterSciNotation
from matplotlib import cm
import matplotlib.colors as mcolors

psr = argparse.ArgumentParser()
psr.add_argument('-o', dest='opt', help='output summary file')
args = psr.parse_args()

cs = CrossSection()
E_nu = np.arange(1, 21, 0.1)
T_e = np.arange(0.1, 21, 0.1)
T_max = cs.T_max(E_nu)
d_sigma_nu_e = cs.dif_T_nu_e(T_e[np.newaxis, :], E_nu[:, np.newaxis])

X, Y = np.meshgrid(E_nu, T_e)
with PdfPages(args.opt) as pdf:
    fig, ax = plt.subplots()
    ax.plot(E_nu, T_max, label=r'$T_{e}^{\mathrm{max}}$[MeV]')
    ax.plot(E_nu, E_nu, ls='--')
    ax.set_xlabel(r'$E_{\nu}$[MeV]')
    ax.set_ylabel(r'$T_{e}^{\mathrm{max}}$[MeV]')
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.legend()
    pdf.savefig(fig)
    plt.close()

    fig, ax = plt.subplots()
    surf = ax.pcolormesh(X, Y, d_sigma_nu_e.T, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0)
    ax.set_xlabel(r'$E_{\nu}$[MeV]')
    ax.set_ylabel(r'$T_{e}$[MeV]')
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    # It seems that ScalarFormatter(useMathText=False) can not fix the format of the colorbar
    cbar = fig.colorbar(surf)
    cbar.set_label('$\sigma$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$]', fontsize=15)
    cbar.ax.tick_params(size=1, labelsize=1)
    pdf.savefig(fig)
    plt.close()
