'''
Preview and summary the SSM results
'''
import numpy as np, h5py
from SSMLoader import SSMReaderFactory
import argparse
import matplotlib.pyplot as plt
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
        opt.attrs['Luminosityy'] = ssmReader.Luminosity
        opt.attrs['Radius'] = ssmReader.Radius
        opt.attrs['power'] = ssmReader.FluxPower
        opt.create_dataset('SSM', data=ssm.to_records(), compression='gzip')
        opt.create_dataset('Flux', data=flux.to_records(), compression='gzip')
        opt.create_dataset('TotalFlux', data=totalFlux.to_records(), compression='gzip')
else:
    with h5py.File(args.ipt, 'r') as ipt:
       ssm, flux = ipt['SSM'][:], ipt['Flux'][:]
    with PdfPages(args.opt + '.pdf') as pdf:
        # Model Preview
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(ssm['R'], ssm['M'])
        axs[0].set_ylabel('Mass')
        axs[1].plot(ssm['R'], ssm['T'])
        axs[1].set_ylabel('T/K')
        axs[2].plot(ssm['R'], ssm['Rho'])
        axs[2].set_ylabel(r'$\rho/(\mathrm{g/cm}^3)$')
        axs[2].set_xlabel('R($R_{sum}$)')
        axs[0].tick_params('x', labelbottom=False)
        axs[1].tick_params('x', labelbottom=False)
        pdf.savefig(fig)
        plt.close()
        # Flux Preview
        fig, ax = plt.subplots()
        for nu in ['pp', 'B8', 'Be7', 'pep', 'hep', 'N13', 'O15', 'F17']:
            ax.plot(flux['R'], flux[nu], label=nu)
        ax.set_ylabel('fraction')
        ax.set_xlabel('R($R_{sum}$)')
        ax.legend()
        pdf.savefig(fig)
        plt.close()
