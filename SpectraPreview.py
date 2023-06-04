'''
Preview and summary the spectra results
'''
import numpy as np, h5py
from SSMLoader import SpectraReaderFactory
import argparse
import matplotlib.pyplot as plt
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from config import spectraConfig
from Oscillation import V_e, MSW_Adiabatic
psr = argparse.ArgumentParser()
psr.add_argument('-i', dest='ipt', nargs='+', help='spectra files')
psr.add_argument('-o', dest='opt', help='output summary file')
psr.add_argument('--reactions', nargs='+', help='reaction names')
psr.add_argument('--models', nargs='+', help='neutrino flux of SSM models')
args = psr.parse_args()
reactions_CNO = ['N13', 'O15', 'F17']
reactions_pp = ['B8', 'Be7', 'pp', 'pep', 'hep']
files, reactions = args.ipt, args.reactions
figureDatas = []
fluxDatas = []
for model in args.models:
    with h5py.File(model, 'r') as ipt:
       fluxDatas.append({
           'TotalFlux': ipt['TotalFlux'][:],
           'power': ipt.attrs['power'],
           'Flux': ipt['Flux'][:]
           })

if True:
    for reaction, f in zip(reactions, files):
        reactionReader = SpectraReaderFactory(reaction, f)
        figureDatas.append(reactionReader)
    '''
    with h5py.File(args.opt, 'w') as opt:
        opt.attrs['type'] = args.ssm
        opt.attrs['Luminosityy'] = ssmReader.Luminosity
        opt.attrs['Radius'] = ssmReader.Radius
        opt.create_dataset('SSM', data=ssm.to_records(), compression='gzip')
        opt.create_dataset('Flux', data=flux.to_records(), compression='gzip')
    '''
    with PdfPages(args.opt + '.pdf') as pdf:
        # Spectra Preview
        fig, ax = plt.subplots()
        for reactionData, reaction in zip(figureDatas, reactions):
            if reactionData.continous:
                ls = '-' if reaction in reactions_pp else '--'
                ys = reactionData.spectra['P']#/reactionData.binwidth
                ax.plot(reactionData.spectra['E'], ys, label=reaction, color=spectraConfig[reaction]['color'], ls = ls)
                xi = int(np.argmax(reactionData.spectra['P'])*spectraConfig[reaction]['x'])
                ax.text(s=spectraConfig[reaction]['label'], x=reactionData.spectra['E'][xi], y=spectraConfig[reaction]['y'] * ys[xi], color=spectraConfig[reaction]['color'])
            else:
                for i in range(reactionData.spectraNum):
                    ax.axline(reactionData.Es[i])
        ax.set_ylabel('P/MeV')
        ax.set_xlabel(r'$E_{\nu}$/MeV')
        ax.set_xlim([0.03, 20])
        ax.set_yscale('log')
        pdf.savefig(fig)
        ax.set_xscale('log')
        pdf.savefig(fig)
        plt.close()
        # Rate Preview
        for flux in fluxDatas:
            fig, ax = plt.subplots()
            for reactionData, reaction in zip(figureDatas, reactions):
                if reactionData.continous:
                    ls = '-' if reaction in reactions_pp else '--'
                    ys = reactionData.spectra['P'] * (flux['TotalFlux'][reaction]*10**flux['power'])
                    ax.plot(reactionData.spectra['E'], ys, label=reaction, color=spectraConfig[reaction]['color'], ls = ls)
                    xi = int(np.argmax(reactionData.spectra['P'])*spectraConfig[reaction]['x'])
                    ax.text(s=spectraConfig[reaction]['label'], x=reactionData.spectra['E'][xi], y=spectraConfig[reaction]['y'] * ys[xi], color=spectraConfig[reaction]['color'])
                else:
                    for i in range(reactionData.spectraNum):
                        ax.axline(reactionData.Es[i])
            ax.set_ylabel('Flux[cm$^{-2}$s$^{-1}$MeV$^{-1}$]')
            ax.set_xlabel(r'$E_{\nu}$[MeV]')
            ax.set_xlim([0.03, 20])
            ax.set_yscale('log')
            ax.set_xscale('log')
            pdf.savefig(fig)
            plt.close()
            # Upturn Preview
            fig, ax = plt.subplots()
            Es = np.arange(0.1, 20, 0.1)
            Pee = MSW_Adiabatic(Es, V_e(10**flux['Flux']['Log10_e_rho'][0]))
            ax.plot(Es, Pee, label='MSW-LMA')
            ax.set_ylabel('$P_{ee}$')
            ax.set_xlabel(r'$E_{\nu}$[MeV]')
            ax.set_xlim([0.1, 20])
            ax.set_ylim([0.1, 0.8])
            ax.set_xscale('log')
            ax.legend()
            pdf.savefig(fig)
            plt.close()
            fig, ax = plt.subplots()
            V_es = V_e(10**flux['Flux']['Log10_e_rho'])
            X, Y = np.meshgrid(flux['Flux']['R'], Es)
            Pee = MSW_Adiabatic(Es, V_es)
            CS = ax.contour(X, Y, Pee, origin='lower', levels=np.arange(0.2,0.7,0.02), cmap='flag')
            ax.clabel(CS, CS.levels, inline=True, fmt='%.2f', fontsize=12)
            ax.set_xlabel('$R$')
            ax.set_ylabel(r'$E_{\nu}$[MeV]')
            ax.set_ylim([0.1, 20])
            ax.set_xlim([0, 0.5])
            ax.legend()
            pdf.savefig(fig)
            plt.close()

