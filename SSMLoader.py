import numpy as np
import pandas as pd
def SSMReaderFactory(SSMType):
    if SSMType == 'BSB':
        reader = BSBReader()
    elif SSMType == 'MB22':
        reader = MB22Reader()
    else:
        return 0
    return reader
def SpectraReaderFactory(reaction, f):
    if reaction == 'pp':
        reader = PPReader(f)
    elif reaction in ['hep', 'N13', 'O15', 'F17']:
        reader = SpectraBareTwoColReader(f)
    elif reaction == 'B8':
        reader = B8Reader(f)
    else:
        return 0
    return reader
class BSBReader():
    def __init__(self):
        self.skipSSMHeader = 23
        self.skipSSMFooter = 3
        self.skipFluxHeader = 27
        self.modelNames = ['M', 'R', 'T', 'Rho', 'P', 'L', 'X', 'Y', 'He3', 'C12', 'N14', 'O16']
        self.fluxNames = ['R', 'T', 'Log10_e_rho', 'M', 'Be7_M', 'pp', 'B8', 'N13', 'O15', 'F17', 'Be7', 'pep', 'hep']
        self.totalFluxNames = ['pp', 'pep', 'hep', 'Be7', 'B8', 'N13', 'O15', 'F17']
    def read(self, files):
        self.readSSM(files[0])
        self.readFlux(files[1])
    def readSSM(self, file):
        self.SSM = pd.read_csv(file, skiprows=self.skipSSMHeader, skipfooter=self.skipSSMFooter, sep='\s+', names=self.modelNames, engine='python')
        self.N = self.SSM.shape[0]
        L_R = np.loadtxt(file, skiprows=self.skipSSMHeader + self.N + 1, delimiter='= ', dtype=[('label', 'U12'), ('value', float)])
        self.Luminosity = L_R[0]['value']
        self.Radius = L_R[1]['value']
    def readFlux(self, file):
        self.Flux = pd.read_csv(file, skiprows=self.skipFluxHeader, sep='\s+', names=self.fluxNames)
        self.TotalFlux = pd.read_csv(file, skiprows=6, sep='\s+', nrows=1, names=self.totalFluxNames)
        self.FluxPower = 10
        # self.FluxUncertainty = 
class MB22Reader():
    def read(self, file):
       pass 
class B8Reader():
    def __init__(self, file, continous=True):
        self.skipFooter = 2
        self.skipHeader = 16
        self.header = ['E', 'P', 'sigmaPlus', 'sigmaMinus']
        self.continous = continous
        self.spectra = pd.read_csv(file, skiprows=self.skipHeader, skipfooter=self.skipFooter, sep='\s+', names=self.header, engine='python')
        self.binwidth = self.spectra.iloc[1]['E']-self.spectra.iloc[0]['E']
class PPReader():
    def __init__(self, file, continous=True):
        self.header = ['E', 'P']
        self.continous = continous
        rawdata = np.loadtxt(file, skiprows=3)
        self.spectra = pd.DataFrame({ self.header[0]: np.ravel(rawdata[:,::2], order='F'), self.header[1]: np.ravel(rawdata[:,1::2], order='F')})
        self.binwidth = self.spectra.iloc[1]['E']-self.spectra.iloc[0]['E']
class SpectraBareTwoColReader():
    def __init__(self, file, continous=True):
        self.header = ['E', 'P']
        self.continous = continous
        self.spectra = pd.read_csv(file, sep='\s+', names=self.header, engine='python')
        self.binwidth = self.spectra.iloc[1]['E']-self.spectra.iloc[0]['E']
class Be7Reader():
    def __init__(self, file, continous=False):
        self.header = ['E', 'P']
        self.N = 92
        self.continous = continous
        self.spectraNum = 2
        self.Es = np.array([np.loadtxt(file, skiprows=1, delimiter=' keV)', max_rows=1), np.loadtxt(file, skiprows=self.N+3+2, delimiter=' keV)', max_rows=1)])
        self.spectra_1 = pd.read_csv(file, sep='\s+', skiprows=3, nrows=self.N, names=self.header, engine='python')
        self.spectra_2 = pd.read_csv(file, sep='\s+', skiprows=3*2+1+self.N, nrows=self.N, names=self.header, engine='python')

