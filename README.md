# Introduction

# Data directory Sturcture
## Data
+ The pattern of model data: `[Model Name]/[Abundance Name]/[model|flux].dat`
  + `Model Name`: BP2004, BP2005, B16
  + `Abundance Name`: gs98, ags05, ags09, mb22
+ The pattern of spectra data: `SPECTRA/[reaction|element].dat`
  + `reaction|element`: pp, hep, B8, Be7, N13, O15, F17
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
+ `SSMPreview.py`
+ `SpectraPreview.py`
+ `NeutrinoEvolution.py`
