'''
neutrino evolution in Solar and Vaccum, and Earth
'''
import numpy as np, pandas as pd, h5py
from SSMLoader import SpectraReaderFactory
import argparse
import matplotlib.pyplot as plt
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from Oscillation import V_e, MSW_Adiabatic
from config import spectraConfig
psr = argparse.ArgumentParser()
psr.add_argument('-i', dest='ipt', help='SSM model and flux hdf5 file')
psr.add_argument('-o', dest='opt', help='output summary file')
psr.add_argument('--spectra', nargs='+', help='spectra files')
psr.add_argument('--reactions', nargs='+', help='reaction names')
psr.add_argument('--plot', default=False, action='store_true')
args = psr.parse_args()
if not args.plot:
    with h5py.File(args.ipt, 'r') as ipt:
        fluxData = {
           'TotalFlux': ipt['TotalFlux'][:],
           'power': ipt.attrs['power'],
           'Flux': ipt['Flux'][:]
           }
    files, reactions = args.spectra, args.reactions
    spectraDatas = []
    for reaction, f in zip(reactions, files):
        reactionReader = SpectraReaderFactory(reaction, f)
        spectraDatas.append(reactionReader)
    with h5py.File(args.opt, 'w') as opt:
        for spectraData, reaction in zip(spectraDatas, reactions):
            Pee = MSW_Adiabatic(spectraData.spectra['E'].to_numpy(), V_e(10**fluxData['Flux']['Log10_e_rho']))
            spectra = spectraData.spectra['P'].to_numpy() * np.sum(Pee * fluxData['Flux'][reaction], axis=1) * fluxData['TotalFlux'][reaction] * 10**fluxData['power'][reaction]
            df = pd.DataFrame({ 'E': spectraData.spectra['E'], 'nu': spectra})
            df_origin = pd.DataFrame({ 'E': spectraData.spectra['E'], 'nu': spectraData.spectra['P'].to_numpy() * fluxData['TotalFlux'][reaction] * 10**fluxData['power'][reaction]})
            opt.create_dataset('Earth/'+reaction, data=df.to_records(), compression='gzip')
            opt.create_dataset('Origin/'+reaction, data=df_origin.to_records(), compression='gzip')
else:
    spectraEarthDatas = []
    spectraOriginDatas = []
    with h5py.File(args.ipt, 'r') as ipt:
        for reaction in args.reactions:
            spectraEarthDatas.append(ipt['Earth/'+reaction][:])
            spectraOriginDatas.append(ipt['Origin/'+reaction][:])
    with PdfPages(args.opt) as pdf:
        fig, ax = plt.subplots()
        for spectra, reaction in zip(spectraEarthDatas, args.reactions):
            ax.plot(spectra['E'], spectra['nu'], label=reaction, color=spectraConfig[reaction]['color'])
        ax.set_xlabel(r'$E_{\nu}$[MeV]')
        ax.set_ylabel('Flux[cm$^{-2}$s$^{-1}$MeV$^{-1}$]')
        ax.set_xlim([0.1, 20])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        pdf.savefig(fig)
        plt.close()

        for spectra, spectra_origin, reaction in zip(spectraEarthDatas, spectraOriginDatas, args.reactions):
            fig, ax = plt.subplots()
            ax.plot(spectra['E'], spectra['nu'], label='Earth')
            ax.plot(spectra_origin['E'], spectra_origin['nu'], label='Origin')
            ax.set_xlim([0.1, 20])
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'$E_{\nu}$[MeV]')
            ax.set_ylabel('Flux[cm$^{-2}$s$^{-1}$MeV$^{-1}$]')
            ax.set_title(reaction)
            ax.legend()
            pdf.savefig(fig)
            plt.close()
            if reaction=='B8':
                fig, ax = plt.subplots()
                ax.plot(spectra['E'], spectra['nu']/spectra_origin['nu'])
                ax.set_xscale('log')
                ax.set_xlabel(r'$E_{\nu}$[MeV]')
                ax.set_ylabel('$P_{ee}$')
                ax.set_title(reaction)
                pdf.savefig(fig)
                plt.close()
