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
# check the cross section
binwidth = 0.1
E_nu = np.arange(1, 21, binwidth)
T_e = np.arange(0.1, 21, binwidth)
T_max = cs.T_max(E_nu)
d_sigma_nu_e_corr = cs.dif_T_nu_e(T_e[np.newaxis, :], E_nu[:, np.newaxis])
d_sigma_nu_mu_corr = cs.dif_T_nu_mu(T_e[np.newaxis, :], E_nu[:, np.newaxis])

d_sigma_nu_e = cs.dif_T_nu_e(T_e[np.newaxis, :], E_nu[:, np.newaxis], corr=False)
d_sigma_nu_mu = cs.dif_T_nu_mu(T_e[np.newaxis, :], E_nu[:, np.newaxis], corr=False)

# check the solar angle distribution
# T_e=0.1MeV min(cos_theta_e)=0.299, it is suitable to select 0.25 as the cut edge
theta_e_edges = np.linspace(0.25, 1, 2500)
theta_es = (theta_e_edges[1:] + theta_e_edges[:-1]) / 2
E_nus, mask = cs.c_theta_E(T_e[:, np.newaxis], theta_es[np.newaxis, :])
# E_nus[E_nus<0] = 0
E_nus_mask = np.ma.array(E_nus, mask=mask)
print('negative values number: {}'.format(np.sum(E_nus_mask<=0)))
f_E_e_E_nu = cs.dif_T_nu_e(T_e[:, np.newaxis], E_nus_mask)
f_E_e_theta_e = f_E_e_E_nu * E_nus_mask * (1 + E_nus_mask / cs.m_e)


def E_nu_T(fig, ax, X, Y, d_sigma_nu):
    surf = ax.pcolormesh(X, Y, d_sigma_nu.T, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0, rasterized=True)
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
    X, Y = np.meshgrid(E_nu, T_e)

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

    # check the solar angle
    X, Y = np.meshgrid(T_e, theta_es)

    fig, ax = plt.subplots()
    surf = ax.pcolormesh(X, Y, E_nus_mask.T, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0, rasterized=True)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    cbar = fig.colorbar(surf)
    cbar.set_label(r'$E_\nu$[MeV]', fontsize=15)
    # cbar.ax.tick_params(size=1, labelsize=1)
    ax.set_xlabel(r'$T_{e}$[MeV]')
    ax.set_ylabel(r'$\cos\theta_e$')
    pdf.savefig(fig)
    plt.close()

    fig, ax = plt.subplots()
    surf = ax.pcolormesh(X, Y, f_E_e_theta_e.T, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0, rasterized=True)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    cbar = fig.colorbar(surf)
    cbar.set_label(r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}\frac{\mathrm{d}E_\nu}{\mathrm{d}\cos\theta_e}$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]', fontsize=15)
    # cbar.ax.tick_params(size=1, labelsize=1)
    ax.set_xlabel(r'$T_{e}$[MeV]')
    ax.set_ylabel(r'$\cos\theta_e$')
    pdf.savefig(fig)
    plt.close()
