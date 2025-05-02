'''
Calculate the survival probability of solar neutrino
'''
import numpy as np, pandas as pd, h5py
from SSMLoader import SpectraReaderFactory
import argparse
import matplotlib.pyplot as plt
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from config import spectraConfig
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap
from matplotlib import cm
from config import oscParas

newcolors = cm.jet(np.linspace(0, 1, 8192*4))
cmap = ListedColormap(newcolors)

psr = argparse.ArgumentParser()
psr.add_argument('--sun', dest='sun', help='survival probability hdf5 file for sun')
psr.add_argument('-o', dest='opt', help='output PDF file')
args = psr.parse_args()

# The grid of theta_12, theta_13, E/Delta_m12
s_theta_12_2s = np.arange(0.05, 1, 0.01)
s_theta_13_2s = np.arange(0.015, 0.03, 0.005)
log_E_Delta_m12s = np.arange(3, 11, 0.05)

c_theta_12_2s = 1 - s_theta_12_2s
# load MSW results of sun
with h5py.File(args.sun, 'r') as ipt:
    prob_sun_B8 = ipt['sun/B8'][:]
    prob_sun_hep = ipt['sun/hep'][:]
'''
# load MSW results of earth
with h5py.File(args.earth, 'r') as ipt:
    prob_earth = ipt['earth/prob'][:]

# load MSW results of night
with h5py.File(args.night, 'r') as ipt:
    prob_night_B8 = ipt['night/B8'][:]
    prob_night_hep = ipt['night/hep'][:]
'''
def prob_2d(X, Y, prob, fig, ax, label):
    c = ax.pcolormesh(X, Y, prob, shading='auto', cmap=cmap, rasterized=True)
    cbar = fig.colorbar(c)
    cbar.set_label(label)
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax.set_ylabel(r'$\sin^2{\theta_{12}}$')

with PdfPages(args.opt) as pdf:
    X, Y = np.meshgrid(10**log_E_Delta_m12s, s_theta_12_2s)
    X1, Y = np.meshgrid(10**log_E_Delta_m12s * oscParas['nufit5.2NO_woSK']['delta_m21_2'], s_theta_12_2s)
    for prob, label in zip([prob_sun_B8, prob_sun_hep], [r'$P_{ee}^\mathrm{B8,day}$', r'$P_{ee}^{hep,\mathrm{day}}$']):
    # for prob, label in zip([prob_sun_B8, prob_sun_hep, prob_night_B8, prob_night_hep], prob_earth):
        fig, ax = plt.subplots()
        # select P_e2, using s_theta_13_2 = 0.020
        prob_2d(X, Y, prob[3, :, :, 1].T, fig, ax, label)
        ax.set_xscale('log')
        ax.set_xlabel(r'$E_\nu/\Delta m_{12}^2$')
        pdf.savefig(fig)
        plt.close()


        fig, ax = plt.subplots()
        prob_2d(X1, Y, prob[0, :, :, 1].T, fig, ax, label)
        ax.set_xlim([1, 30])
        ax.set_xlabel(r'$E_\nu$[MeV]')
        pdf.savefig(fig)
        plt.close()


