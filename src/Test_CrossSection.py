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
binwidth = 0.1
E_nu = np.arange(1, 21, binwidth)
T_e = np.arange(0.1, 21, binwidth)
T_max = cs.T_max(E_nu)
d_sigma_nu_e_corr = cs.dif_T_nu_e(T_e[np.newaxis, :], E_nu[:, np.newaxis])
d_sigma_nu_mu_corr = cs.dif_T_nu_mu(T_e[np.newaxis, :], E_nu[:, np.newaxis])

d_sigma_nu_e = cs.dif_T_nu_e(T_e[np.newaxis, :], E_nu[:, np.newaxis], corr=False)
d_sigma_nu_mu = cs.dif_T_nu_mu(T_e[np.newaxis, :], E_nu[:, np.newaxis], corr=False)

X, Y = np.meshgrid(E_nu, T_e)

def E_nu_T(fig, ax, X, Y, d_sigma_nu):
    surf = ax.pcolormesh(X, Y, d_sigma_nu.T, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0)
    ax.set_ylabel(r'$T_{e}$[MeV]')
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    # It seems that ScalarFormatter(useMathText=False) can not fix the format of the colorbar
    cbar = fig.colorbar(surf)
    cbar.set_label('$\mathrm{d}\sigma/\mathrm{d}T_e$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]', fontsize=15)
    cbar.ax.tick_params(size=1, labelsize=1)

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
    ax.plot(E_nu, np.sum(d_sigma_nu_e_corr, axis=1) * binwidth, color='b', label=r'$\nu_e$ radiative correction')
    ax.plot(E_nu, np.sum(d_sigma_nu_mu_corr, axis=1) * binwidth, color='g', label=r'$\nu_\mu$ radiative correction')
    ax.plot(E_nu, np.sum(d_sigma_nu_e, axis=1) * binwidth, color='b', ls='--', label=r'$\nu_e$')
    ax.plot(E_nu, np.sum(d_sigma_nu_mu, axis=1) * binwidth, color='g', ls='--', label=r'$\nu_\mu$')
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.set_xlabel(r'$E_{\nu}$[MeV]')
    ax.set_ylabel('$\sigma$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$]')
    ax.legend()
    pdf.savefig(fig)
    plt.close()

    fig, ax = plt.subplots()
    E_nu_T(fig, ax, X, Y, d_sigma_nu_e_corr)
    ax.set_xlabel(r'$E_{\nu_e}$[MeV]')
    pdf.savefig(fig)
    plt.close()

    fig, ax = plt.subplots()
    E_nu_T(fig, ax, X, Y, d_sigma_nu_mu_corr)
    ax.set_xlabel(r'$E_{\nu_\mu}$[MeV]')
    pdf.savefig(fig)
    plt.close()

