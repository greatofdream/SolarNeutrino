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
E_nu_s, E_max = 1, 21
T_e_s = 0.1

E_nu = np.arange(E_nu_s, E_max, binwidth)
T_e = np.arange(T_e_s, E_max, binwidth)
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

def get_f_E_e_theta_e(mask):
    E_nus_mask = np.ma.array(E_nus, mask=mask)
    print('negative values number: {}'.format(np.sum(E_nus_mask<=0)))
    f_E_e_E_nu = cs.dif_T_nu_e(T_e[:, np.newaxis], E_nus_mask)
    return E_nus_mask, f_E_e_E_nu * E_nus_mask * (1 + E_nus_mask / cs.m_e)

E_nus_mask, f_E_e_theta_e = get_f_E_e_theta_e(mask)

E_nu_select_1 = [6, 8]
E_nu_select_2 = [16, 18]

E_nus_mask_1, f_E_e_theta_e_1 = get_f_E_e_theta_e(mask | (E_nus <E_nu_select_1[0]) | (E_nus >E_nu_select_1[1]))

E_nus_mask_2, f_E_e_theta_e_2 = get_f_E_e_theta_e(mask | (E_nus <E_nu_select_2[0])| (E_nus >E_nu_select_2[1]))

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
    ax.set_yscale('log')
    pdf.savefig(fig)
    plt.close()

    # select the E_nu_check=10MeV for check
    E_nu_check, scale = 10, 100
    E_nu_i = int((E_nu_check - E_nu_s) / binwidth)
    fig, ax = plt.subplots()
    ax.plot(T_e, scale * d_sigma_nu_e_corr[E_nu_i], color='b', label=r'$\nu_e$ radiative correction')
    ax.plot(T_e, scale * d_sigma_nu_mu_corr[E_nu_i], color='g', label=r'$\nu_\mu$ radiative correction')
    ax.plot(T_e, scale * d_sigma_nu_e[E_nu_i], color='b', ls='--', label=r'$\nu_e$')
    ax.plot(T_e, scale * d_sigma_nu_mu[E_nu_i], color='g', ls='--', label=r'$\nu_\mu$')
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_xlim([0, E_nu_check*1.1])
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$T_{e}$[MeV]')
    ax.set_ylabel(r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}$[$10^{' + '{}'.format(cs.power - 2) + '}\mathrm{cm}^2$]')
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

    def check2d(arr, cbar_label):
        fig, ax = plt.subplots()
        surf = ax.pcolormesh(X, Y, arr, norm=mcolors.LogNorm(), cmap=cm.jet, linewidth=0, rasterized=True)
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        cbar = fig.colorbar(surf)
        cbar.set_label(cbar_label, fontsize=15)
        ax.set_xlabel(r'$T_{e}$[MeV]')
        ax.set_ylabel(r'$\cos\theta_e$')
        pdf.savefig(fig)
        plt.close()

    check2d(E_nus_mask.T, r'$E_\nu$[MeV]')

    check2d(f_E_e_theta_e.T, r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}\frac{\mathrm{d}E_\nu}{\mathrm{d}\cos\theta_e}$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]')

    f_E_e_theta_e[E_nus_mask > E_max] = 0
    check2d(f_E_e_theta_e.T, r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}\frac{\mathrm{d}E_\nu}{\mathrm{d}\cos\theta_e}$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]')

    check2d(f_E_e_theta_e_1.T, r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}\frac{\mathrm{d}E_\nu}{\mathrm{d}\cos\theta_e}$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]')

    check2d(f_E_e_theta_e_2.T, r'$\frac{\mathrm{d}\sigma}{\mathrm{d}T_e}\frac{\mathrm{d}E_\nu}{\mathrm{d}\cos\theta_e}$[$10^{' + '{}'.format(cs.power) + '}\mathrm{cm}^2$/MeV]')
