'''
Preview and summary the spectra results
'''
import numpy as np, h5py
from SSMLoader import SpectraReaderFactory
import argparse
import matplotlib.pyplot as plt
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from config import spectraConfig
from Oscillation import V_e, MSW_Adiabatic
from config import oscParas
NOpara = oscParas['nufit5.2NO_woSK']

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
            ls = '-' if reaction in reactions_pp else '--'
            ys = reactionData.spectra['P']#/reactionData.binwidth
            ax.plot(reactionData.spectra['E'], ys, label=reaction, color=spectraConfig[reaction]['color'], ls = ls)
            xi = int(np.argmax(reactionData.spectra['P'])*spectraConfig[reaction]['x'])
            ax.text(s=spectraConfig[reaction]['label'], x=reactionData.spectra['E'][xi], y=spectraConfig[reaction]['y'] * ys[xi], color=spectraConfig[reaction]['color'])
        ax.set_ylabel('PDF/MeV')
        ax.set_xlabel(r'$E_{\nu}$[MeV]')
        ax.set_xlim([0.03, 20])
        ax.set_ylim(bottom=1E-6)
        ax.set_yscale('log')
        pdf.savefig(fig)
        ax.set_xscale('log')
        pdf.savefig(fig)
        plt.close()

        # Rate Preview
        for flux in fluxDatas:
            fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 6]})
            axs[0].set_frame_on(False)
            axs[0].xaxis.set_visible(False)
            axs[0].yaxis.set_visible(False)
            for reactionData, reaction in zip(figureDatas, reactions):
                ls = '-' if reaction in reactions_pp else '--'
                ys = reactionData.spectra['P'] * (flux['TotalFlux'][reaction] * 10 ** flux['power'][reaction])
                axs[1].plot(reactionData.spectra['E'], ys, label=reaction, color=spectraConfig[reaction]['color'], ls = ls)
                xi = int(np.argmax(reactionData.spectra['P'])*spectraConfig[reaction]['x'])
                # text of the reaction and the uncertainty
                axs[1].text(s=spectraConfig[reaction]['label'] + r'[$\pm$' + '{}%]'.format(spectraConfig[reaction]['error']), x=reactionData.spectra['E'][xi], y=spectraConfig[reaction]['y'] * ys[xi], rotation=spectraConfig[reaction]['rotation'], color=spectraConfig[reaction]['color'])
            axs[0].annotate("Ga", xy=(0.233, 0.1), xycoords='data', xytext=(10, 0.1),
                arrowprops=dict(arrowstyle="<|-, widthA=0.5, widthB=2",
                                connectionstyle="arc3", linewidth=3, color='g'), fontsize=10, color='g')
            axs[0].annotate("Cl", xy=(0.814, 0.3), xycoords='data', xytext=(10, 0.3),
                arrowprops=dict(arrowstyle="<|-, widthA=0.5, widthB=2",
                connectionstyle="arc3", linewidth=3, color='r'), fontsize=10, color='r')
            axs[0].annotate("Cerenkov", xy=(3.5, 0.5), xycoords='data', xytext=(10, 0.5),
                arrowprops=dict(arrowstyle="<|-, widthA=0.5, widthB=2",
                connectionstyle="arc3", linewidth=3, color='b'), fontsize=10, color='b')
            axs[0].annotate("LS", xy=(0.2, 0.7), xycoords='data', xytext=(10, 0.7),
                arrowprops=dict(arrowstyle="<|-, widthA=0.5, widthB=2",
                connectionstyle="arc3", linewidth=3), fontsize=10)
            axs[1].set_ylabel('Flux[cm$^{-2}$s$^{-1}$MeV$^{-1}$]')
            axs[1].set_xlabel(r'$E_{\nu}$[MeV]')
            axs[1].set_xlim([0.03, 20])
            axs[1].set_ylim(bottom=0.1)
            axs[1].set_yscale('log')
            axs[1].set_xscale('log')
            fig.subplots_adjust(hspace=0)
            pdf.savefig(fig)
            plt.close()

            # Upturn Preview
            Es = np.arange(0.1, 20, 0.1)
            V_CC = V_e(10**flux['Flux']['Log10_e_rho'][0])
            Pee = MSW_Adiabatic(Es, V_CC, delta_m21_2=NOpara['delta_m21_2'], delta_m3l_2=NOpara['delta_m3l_2'], s12_2=NOpara['s12_2'], s13_2=NOpara['s13_2'],)
            P_0 = (1 - NOpara['s13_2'])**2 * (1 - 2 * NOpara['s12_2'] * (1 - NOpara['s12_2'])) + NOpara['s13_2']**2
            P_1 = (1 - NOpara['s13_2'])**2 * NOpara['s12_2'] + NOpara['s13_2']**2
            fig, ax = plt.subplots()
            ax.plot(Es, Pee, label='MSW-LMA')
            ax.axhline(P_0, ls='--', c='r', label=r'$c_{13}^4(1-\frac{1}{2}\sin^2{2\theta_{12}})+s_{13}^4$')
            ax.axhline(P_1, ls='--', c='g', label=r'$c_{13}^4c_{12}^2+s_{13}^4$')
            ax.set_ylabel('$P_{ee}$')
            ax.set_xlabel(r'$E_{\nu}$[MeV]')
            ax.set_xlim([0.1, 20])
            ax.set_ylim([0.1, 0.8])
            ax.yaxis.set_minor_locator(MultipleLocator(0.02))
            ax.legend()
            pdf.savefig(fig)
            ax.set_xscale('log')
            pdf.savefig(fig)
            plt.close()

            V_es = V_e(10**flux['Flux']['Log10_e_rho'])
            fig, ax = plt.subplots()
            ax2 = ax.twinx()
            X, Y = np.meshgrid(flux['Flux']['R'], Es)
            Pee = MSW_Adiabatic(Es, V_es)
            CS = ax.contour(X, Y, Pee, origin='lower', levels=np.arange(0.2,0.7,0.02), cmap='flag')
            ax.clabel(CS, CS.levels, inline=True, fmt='%.2f', fontsize=12)
            ax.set_xlabel('$R$')
            ax.set_ylabel(r'$E_{\nu}$[MeV]')
            ax.set_ylim([0.1, 20])
            ax.set_xlim([0, 0.5])
            ax2.plot(flux['Flux']['R'], V_es*1E6, ls='--', color='g')
            ax2.set_ylabel(r'$V_{e}^0[10^{-12}\mathrm{eV}]$', color='g')
            pdf.savefig(fig)
            plt.close()

        # Spectra Preview
        fig, ax = plt.subplots()
        for reactionData, reaction in zip(figureDatas, reactions):
            if reaction in ['B8', 'hep']:
                ys = reactionData.spectra['P']#/reactionData.binwidth
                ax.plot(reactionData.spectra['E'], ys, label=reaction, color=spectraConfig[reaction]['color'], ls = ls)
        ax.set_ylabel('PDF/MeV')
        ax.set_xlabel(r'$E_{\nu}$[MeV]')
        ax.set_xlim([0.03, 20])
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.legend()
        pdf.savefig(fig)
        plt.close()

