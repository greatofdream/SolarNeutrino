'''
PREM model preview
Modified from https://www.fatiando.org/rockhound/dev/gallery/prem.html#sphx-glr-gallery-prem-py
'''
import rockhound as rh
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.style.use('./journal.mplstyle')
from matplotlib.backends.backend_pdf import PdfPages

# Load PREM into a DataFrame
prem = rh.fetch_prem()

with PdfPages('prem_preview.pdf') as pdf:
    # Plot density and velocities
    fig, axs = plt.subplots(2, 1, figsize=(12, 12), sharex=True)
    fig.subplots_adjust(hspace=0)
    prem.plot("radius", "density", legend=False, ax=axs[0])
    axs[0].set_ylabel("Density [g/cm³]")
    axs[0].grid()
    for velocity in ["Vpv", "Vph", "Vsv", "Vsh"]:
        prem.plot("radius", velocity, legend=False, ax=axs[1], label=velocity)
    axs[1].grid()
    axs[1].legend()
    axs[1].set_ylabel("Velocity [km/s]")
    axs[-1].set_xlabel("Radius [km]")
    pdf.savefig(fig)
