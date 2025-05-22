# Introduction

# Data directory Sturcture
## Data
+ The pattern of model data: `data/[Model Name]/[Abundance Name]/[model|flux].dat`
  + `Model Name`: BP2004, BP2005, B16, ...
  + `Abundance Name`: gs98, ags05, agss09, mb22
    + BP2004 use `gs98`, [data location](https://www.sns.ias.edu/~jnb/SNdata/Export/)
    + BP2005 use `gs98` and `ags05`, [data location](https://www.sns.ias.edu/~jnb/SNdata/Export/)
    + B16 use `gs98` and `agss09`, [data location](https://aliga.ice.csic.es/personal/aldos/Solar_Data.html), [origin data location provided in paper](http://www.ice.cat/personal/aldos/Solar_Data.html) is not accessible.
    + ...: Maybe X25?
  + `model.dat`: the structure (density, temperature, pressure et al.) of the solar; `flux.dat`: flux versus the solar radius.
+ The pattern of spectra data: `SPECTRA/[reaction|element].dat`
  + `reaction|element`: pp, hep, B8, Be7, N13, O15, F17 neutrino energy spectra.
  + The shape of spectra is uncorrelated with the flux. The final spectra is the product of the shape and the flux of the model with different [Abundance]
+ The pattern of crosssection data: `CROSSSECTION/[reaction|element].dat`
+ Preview data and figure
  + `preview.h5`: summary of data
  + `preview.h5.pdf`: preview pictures of data

## Predict
The pattern of the directory is same as `Data` directory.
+ `fluxSolar.h5`: neutrino flux on the solar surface after oscillation in Sun and Vaccum
+ `fluxEarth.h5`: neutrino flux after oscillation in Earth

# Code
The code are executed as `Makefile`. Following code scripts are in the execution order.
+ `SSMPreview.py`: Check the SSM
+ `SpectraPreview.py`: Check the spectra
+ `NuOscApproximate.py`: LMA approximation and simple calculation example
