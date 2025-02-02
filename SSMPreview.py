'''
Preview and summary the SSM results
1. SSM: model.dat, electron density vs radius
2. Flux: flux.dat, flux vs radius
3. TotalFlux: the total flux value
'''
import numpy as np, h5py
from SSMLoader import SSMReaderFactory
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from scipy.constants import physical_constants

psr = argparse.ArgumentParser()
psr.add_argument('-i', dest='ipt', nargs='+', help='SSM model and flux files')
psr.add_argument('-o', dest='opt', help='output summary file')
psr.add_argument('--ssm', dest='ssm', help='SSM format')
psr.add_argument('--plot', default=False, action='store_true')
args = psr.parse_args()
if not args.plot:
    ssmReader = SSMReaderFactory(args.ssm)
    ssmReader.read(args.ipt)
    ssmReader.calcQDensity()
    ssm = ssmReader.SSM
    flux = ssmReader.Flux
    totalFlux = ssmReader.TotalFlux
    with h5py.File(args.opt, 'w') as opt:
        opt.attrs['type'] = args.ssm
        # opt.attrs['Luminosity'] = ssmReader.Luminosity
        # opt.attrs['Radius'] = ssmReader.Radius
        opt.attrs['power'] = ssmReader.FluxPower.to_records()
        opt.create_dataset('SSM', data=ssm.to_records(), compression='gzip')
        opt.create_dataset('Flux', data=flux.to_records(), compression='gzip')
        opt.create_dataset('TotalFlux', data=ssmReader.TotalFlux.to_records(), compression='gzip')
else:
    with h5py.File(args.ipt[0], 'r') as ipt:
       ssm, flux = ipt['SSM'][:], ipt['Flux'][:]
    with PdfPages(args.opt) as pdf:
        # Model Preview
        fig, axs = plt.subplots(5, 1, figsize=(12, 15), sharex=True)
        fig.subplots_adjust(hspace=0)
        axs[0].plot(ssm['Radius'], ssm['Mass'])
        axs[0].set_ylabel('Mass')
        axs[1].plot(ssm['Radius'], ssm['Rho'])
        axs[1].set_ylabel(r'$\rho[\mathrm{g/cm}^3]$')
        axs[2].plot(ssm['Radius'], ssm['Temp'])
        axs[2].set_ylabel('T/K')
        axs[3].plot(ssm['Radius'], ssm['Pres'])
        axs[3].set_ylabel(r'$P[\mathrm{dyn/cm}^2]$')
        axs[4].plot(flux['R'], 10**flux['Log10_e_rho'], label='$e^-$')# electron density distribution
        axs[4].set_ylabel('#[$\mathrm{N}_\mathrm{A}/\mathrm{cm}^3$]')
        axs[-1].set_xlabel('R[$R_{sun}$]')
        axs[-1].xaxis.set_major_locator(MultipleLocator(0.1))
        axs[-1].xaxis.set_minor_locator(MultipleLocator(0.02))
        axs[-1].set_xlim([0, 1])
        pdf.savefig(fig)
        for i in range(1, 4):
            axs[i].set_yscale('log')
        plt.close()

        ## the elements fraction
        fig, ax = plt.subplots(figsize=(12, 6))
        for ele_name in ['H1', 'He4', 'He3', 'C12', 'N14', 'O16']:
            ax.plot(ssm['Radius'], ssm[ele_name], label=ele_name)
        ax.set_xlabel('R[$R_{sun}$]')
        ax.set_ylabel('Mass fraction')
        ax.set_yscale('log')
        ax.set_xlim([0, 1])
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.legend()
        pdf.savefig(fig)

        ## the electron, u, d quark density
        fig, ax = plt.subplots(figsize=(12, 6))
        for ele_name, label_name in zip(['rho_e', 'rho_u', 'rho_d'], ['$e^-$', 'u', 'd']):
            ax.plot(ssm['Radius'], ssm[ele_name], label=label_name)
        ax.set_xlabel('R[$R_{sun}$]')
        ax.set_ylabel('#[$\mathrm{N}_\mathrm{A}/\mathrm{cm}^3$]')
        ax.set_xlim([0, 1])
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.legend()
        pdf.savefig(fig)

        fig, ax = plt.subplots(figsize=(12, 6))
        for ele_name, label_name in zip(['rho_e', 'rho_u', 'rho_d'], ['$e^-$', 'u', 'd']):
            ax.plot(ssm['Radius'], ssm[ele_name] * physical_constants['Avogadro constant'][0], label=label_name)
        ax.set_xlabel('R[$R_{sun}$]')
        ax.set_ylabel('#[$\mathrm{cm}^{-3}$]')
        ax.set_xlim([0, 1])
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.legend()
        pdf.savefig(fig)

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(ssm['Radius'], ssm['rho_u'] / ssm['rho_e'], label='u/$e^-$')
        ax.plot(ssm['Radius'], ssm['rho_d'] / ssm['rho_e'], label='d/$e^-$')
        ax.set_xlabel('R[$R_{sun}$]')
        ax.set_ylabel(r'$\rho_q/\rho_e$')
        ax.set_xlim([0, 1])
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.legend()
        pdf.savefig(fig)

        # Flux Preview
        fig, ax = plt.subplots()
        for nu, label in zip(['pp', 'B8', 'Be7', 'pep', 'hep', 'N13', 'O15', 'F17'], ['pp', '8B', '7Be', 'pep', 'hep', '13N', '15O', '17F']):
            ax.plot(flux['R'], flux[nu], label=label)
        ax.set_ylabel('Fraction')
        ax.set_xlabel('R[$R_{sun}$]')
        ax.legend()
        pdf.savefig(fig)
        plt.close()

