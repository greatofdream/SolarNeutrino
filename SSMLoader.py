import numpy as np
import pandas as pd
from mendeleev import element
import re
def SSMReaderFactory(SSMType):
    model, abundance = SSMType.split('/')
    if model == 'BP2004' or model == 'BP2005':
        reader = BSBReader(model, abundance)
    elif model == 'B16':
        reader = B16Reader(model, abundance)
    elif model == 'MB22':
        reader = MB22Reader(model, abundance)
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
    elif reaction == 'Be7':
        reader = Be7Reader(f)
    elif reaction == 'pep':
        reader = PepReader(f)
    else:
        return 0
    return reader

class Reader():
    def __init__(self, model, abundance):
        self.model = model
        self.abundance = abundance
    def calcQDensity(self):
        # calculate the electron, u, d quark density
        ## load the element information
        Zs = np.empty(len(self.modelNames) - self.ele_offset, dtype=int)# proton
        Ns = np.empty(len(self.modelNames) - self.ele_offset, dtype=int)# neutron
        As = np.empty(len(self.modelNames) - self.ele_offset) # mass
        for i, element_str in enumerate(self.modelNames[self.ele_offset:]):
            ele = Element(element_str)
            Zs[i], Ns[i], As[i] = ele.Z, ele.N, ele.mass

        ## calculation
        fraction = self.SSM.iloc[:, self.ele_offset:].to_numpy()
        self.SSM['rho_e'] = self.SSM['Rho'] * np.sum(fraction * (Zs / As), axis=1)
        self.SSM['rho_u'] = self.SSM['Rho'] * np.sum(fraction * ((2 * Zs + Ns) / As), axis=1)
        self.SSM['rho_d'] = self.SSM['Rho'] * np.sum(fraction * ((Zs + 2 * Ns) / As), axis=1)

class BSBReader(Reader):
    def __init__(self, model, abundance):
        Reader.__init__(self, model, abundance)
        self.skipSSMHeader = 23
        self.skipSSMFooter = 3
        self.skipFluxHeader = 27
        self.modelNames = ['Mass', 'Radius', 'Temp', 'Rho', 'Pres', 'Lumi', 'H1', 'He4', 'He3', 'C12', 'N14', 'O16']# X: H, Y:He4
        self.fluxNames = ['R', 'T', 'Log10_e_rho', 'M', 'Be7_M', 'pp', 'B8', 'N13', 'O15', 'F17', 'Be7', 'pep', 'hep']
        self.totalFluxNames = ['pp', 'pep', 'hep', 'Be7', 'B8', 'N13', 'O15', 'F17']
        self.ele_offset = 6
    def read(self, files):
        self.readSSM(files[0])
        self.readFlux(files[1])
    def readSSM(self, file):
        self.SSM = pd.read_csv(file, skiprows=self.skipSSMHeader, skipfooter=self.skipSSMFooter, sep='\s+', names=self.modelNames, engine='python')
        self.N = self.SSM.shape[0]
        if self.model == "BP2004":
            L_R = np.genfromtxt(file, skip_header=self.skipSSMHeader + self.N + 1, delimiter=' ', dtype=[('label', 'U12'), ('value', float)])
        else:
            L_R = np.genfromtxt(file, skip_header=self.skipSSMHeader + self.N + 1, delimiter='= ', dtype=[('label', 'U12'), ('value', float)])
        self.Luminosity = L_R[0]['value']
        self.Radius = L_R[1]['value']
    def readFlux(self, file):
        self.Flux = pd.read_csv(file, skiprows=self.skipFluxHeader, sep='\s+', names=self.fluxNames)
        self.TotalFlux = pd.read_csv(file, skiprows=6, sep='\s+', nrows=1, names=self.totalFluxNames)
        # power is referred from astro-ph/0412440
        self.FluxPower = pd.DataFrame(columns=self.totalFluxNames, data=np.repeat(10, self.TotalFlux.shape[1])[np.newaxis,:])
        # self.FluxUncertainty = 

class B16Reader(Reader):
    def __init__(self, model, abundance):
        Reader.__init__(self, model, abundance)
        self.skipSSMHeader = 9
        self.skipSSMFooter = 0
        self.skipFluxHeader = 22
        # "Mass     Radius     Temp      Rho       Pres       Lumi      H1       He4      He3       C12       C13       N14        N15      O16       O17       O18        Ne        Na       Mg         Al        Si         P       S         Cl        Ar        K         Ca        Sc         Ti        V         Cr        Mn        Fe        Co        Ni".split()
        self.modelNames = "Mass     Radius     Temp      Rho       Pres       Lumi      H1       He4      He3       C12       C13       N14        N15      O16       O17       O18        Ne        Na       Mg         Al        Si         P       S         Cl        Ar        K         Ca        Sc         Ti        V         Cr        Mn        Fe        Co        Ni".split()
        self.fluxNames = ['R', 'T', 'Log10_e_rho', 'M', 'pp', 'pep', 'hep', 'Be7', 'B8', 'N13', 'O15', 'F17']
        self.totalFluxNames = ['pp', 'pep', 'hep', 'Be7', 'B8', 'N13', 'O15', 'F17']
        self.ele_offset = 6
    def read(self, files):
        self.readSSM(files[0])
        self.readFlux(files[1])
    def readSSM(self, file):
        self.SSM = pd.read_csv(file, skiprows=self.skipSSMHeader, skipfooter=self.skipSSMFooter, sep='\s+', names=self.modelNames, engine='python')
        self.N = self.SSM.shape[0]
    def readFlux(self, file):
        self.Flux = pd.read_csv(file, skiprows=self.skipFluxHeader, sep='\s+', names=self.fluxNames)
        # total flux come from [https://arxiv.org/abs/1611.09867]
        # it is also store in data/B16/flux.dat
        if self.abundance == "gs98":
            self.TotalFlux = pd.DataFrame(columns=self.totalFluxNames, data=np.array([5.98, 1.44, 7.98, 4.93, 5.46, 2.78, 2.05, 5.29])[np.newaxis,:])
        else:
            # agss09
            self.TotalFlux = pd.DataFrame(columns=self.totalFluxNames, data=np.array([6.03, 1.46, 8.25, 4.50, 4.50, 2.04, 1.44, 3.26])[np.newaxis,:])
        self.FluxPower = pd.DataFrame(columns=self.totalFluxNames, data=np.array([10, 8, 3, 9, 6, 8, 8, 6])[np.newaxis,:])
        
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
        self.header = ['dE', 'P']
        self.N = 92
        self.continous = continous
        self.spectraNum = 2
        self.Es = np.array([861.3, 384.3]) #np.array([np.loadtxt(file, skiprows=1, delimiter=' keV)', max_rows=1), np.loadtxt(file, skiprows=self.N+3+2, delimiter=' keV)', max_rows=1)])
        self.spectra_1 = pd.read_csv(file, sep='\s+', skiprows=3, nrows=self.N, names=self.header, engine='python')
        self.spectra_2 = pd.read_csv(file, sep='\s+', skiprows=3*2+1+self.N, nrows=self.N, names=self.header, engine='python')
        self.spectra_1['E'] = (self.Es[0] + self.spectra_1['dE']) * 1E-3
        self.spectra_2['E'] = (self.Es[1] + self.spectra_2['dE']) * 1E-3
        # For the monochromatic spectra, just use the probability
        # self.spectra_1['P'] = self.spectra_1['P'] * 1E3
        # self.spectra_2['P'] = self.spectra_2['P'] * 1E3
        self.branchRatio = np.array([0.897, 0.103]) # 10.1103/PhysRevD.49.3923, ground-state and excited-state
        self.spectra = pd.DataFrame({'E': [self.Es[0]*1E-3, self.Es[0]*1E-3, self.Es[1]*1E-3, self.Es[1]*1E-3], 'P': [self.branchRatio[0], 1E-11, 1E-11, self.branchRatio[1]]})

class PepReader():
    def __init__(self, file, continous=False):
        self.spectra = pd.DataFrame({'E': [1.44, 1.44], 'P': [1, 1E-10]})

class Element():
    def __init__(self, ele_string):
        # split the mass number
        res  = re.match('([^0-9]+)(\d*)', ele_string)
        self.ele_name = res.group(1)
        self.element = element(self.ele_name)
        self.Z = self.element.protons
        if res.group(2)!="":
            self.mass_num = int(res.group(2))
            self.N = self.mass_num - self.Z
            for iso in self.element.isotopes:
                if iso.mass_number == self.mass_num:
                    self.mass = iso.mass
        else:
            self.mass_num = self.element.mass_number
            self.N = self.element.neutrons
            self.mass = self.element.mass
