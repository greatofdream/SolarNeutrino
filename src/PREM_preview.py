'''
PREM model preview
Modified from https://www.fatiando.org/rockhound/dev/gallery/prem.html#sphx-glr-gallery-prem-py
'''
import rockhound as rh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages

R_core = 0.546 # 10.1103/PhysRevD.56.1792
N_e_to_rho = np.array([0.468, 0.497]) # SK
N_n_to_rho = 1 - N_e_to_rho
N_u_to_rho = 2 * N_e_to_rho + N_n_to_rho
N_d_to_rho = 2 * N_n_to_rho + N_e_to_rho

# Load PREM into a DataFrame
prem = rh.fetch_prem()
R_earth = np.max(prem["radius"])
R_earth_norm = prem["radius"] / R_earth
N_earth = len(prem['radius'])
N_core = np.sum((prem['radius'] / R_earth) <= R_core)
N_mantle = N_earth - N_core

with PdfPages('prem_preview.pdf') as pdf:
    # Plot density and velocities
    fig, axs = plt.subplots(2, 1, figsize=(12, 12), sharex=True)
    fig.subplots_adjust(hspace=0)
    axs[0].plot(R_earth_norm, prem["density"])
    axs[0].axvline(R_core, ls='--', c='k')
    axs[0].set_ylabel(r"$\rho$ [g/cm³]")
    axs[0].grid()
    for velocity in ["Vpv", "Vph", "Vsv", "Vsh"]:
        axs[1].plot(R_earth_norm, prem[velocity], label=velocity)
    axs[1].axvline(R_core, ls='--', c='k')
    axs[1].grid()
    axs[1].legend()
    axs[1].set_ylabel("Velocity [km/s]")
    axs[-1].set_xlabel("Radius [km]")
    axs[-1].xaxis.set_major_locator(MultipleLocator(0.1))
    axs[-1].xaxis.set_minor_locator(MultipleLocator(0.02))
    pdf.savefig(fig)

    # The electron, u, d quark density
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(R_earth_norm, prem['density'] * np.repeat(N_n_to_rho, [N_core, N_mantle]), label='$e^-$')
    ax.plot(R_earth_norm, prem['density'] * np.repeat(N_u_to_rho, [N_core, N_mantle]), label='u')
    ax.plot(R_earth_norm, prem['density'] * np.repeat(N_d_to_rho, [N_core, N_mantle]), label='d')
    ax.legend()
    ax.set_xlabel("Radius [km]")
    ax.set_ylabel('#[$\mathrm{N}_\mathrm{A}/\mathrm{cm}^3$]')
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    pdf.savefig(fig)

