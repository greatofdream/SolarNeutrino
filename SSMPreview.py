'''
Preview and summary the SSM results
'''
import numpy as np, h5py
from SSMLoader import SSMReaderFactory
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages

psr = argparse.ArgumentParser()
psr.add_argument('-i', dest='ipt', nargs='+', help='SSM model and flux files')
psr.add_argument('-o', dest='opt', help='output summary file')
psr.add_argument('--ssm', dest='ssm', help='SSM format')
psr.add_argument('--plot', default=False, action='store_true')
args = psr.parse_args()
if not args.plot:
    ssmReader = SSMReaderFactory(args.ssm)
    ssmReader.read(args.ipt)
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
        fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
        fig.subplots_adjust(hspace=0)
        axs[0].plot(ssm['Radius'], ssm['Mass'])
        axs[0].set_ylabel('Mass')
        axs[1].plot(ssm['Radius'], ssm['Rho'])
        axs[1].set_ylabel(r'$\rho[\mathrm{g/cm}^3]$')
        axs[2].plot(ssm['Radius'], ssm['Temp'])
        axs[2].set_ylabel('T/K')
        axs[3].plot(ssm['Radius'], ssm['Pres'])
        axs[3].set_ylabel(r'$P[\mathrm{dyn/cm}^2]$')
        axs[3].set_xlabel('R[$R_{sun}$]')
        axs[-1].xaxis.set_major_locator(MultipleLocator(0.1))
        axs[-1].xaxis.set_minor_locator(MultipleLocator(0.01))
        axs[-1].set_xlim([0, 1])
        pdf.savefig(fig)
        for i in range(1, 4):
            axs[i].set_yscale('log')
        plt.close()

        # Flux Preview
        fig, ax = plt.subplots()
        for nu, label in zip(['pp', 'B8', 'Be7', 'pep', 'hep', 'N13', 'O15', 'F17'], ['pp', '8B', '7Be', 'pep', 'hep', '13N', '15O', '17F']):
            ax.plot(flux['R'], flux[nu], label=label)
        ax.set_ylabel('Fraction')
        ax.set_xlabel('R[$R_{sum}$]')
        ax.legend()
        pdf.savefig(fig)
        plt.close()
