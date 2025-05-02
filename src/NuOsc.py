'''
Calculate the survival probability of solar neutrino
Due to the electron density is unifrom in different layer in the earth, it is easy to calculate,
Argument: --mode earth
while we split the calculation in different radius in the sun. (The NuOscMerge.py should used as postprocess)
Argument: --mode sun --split
If the oscillation parameters is fixed, it is also quick to calculate the sun directly
Argument: --mode sun
'''
import numpy as np, h5py
from SSMLoader import SpectraReaderFactory
import argparse
from Oscillation import Sun, Earth, P_3nu
from SSMLoader import SSMReaderFactory

psr = argparse.ArgumentParser()
psr.add_argument('-i', dest='ipt', help='SSM model and flux hdf5 file')
psr.add_argument('-o', dest='opt', help='output summary file')
psr.add_argument('--mode', dest='mode', default='sun', help='calculate the evolution in some radius in the sun or in the whole earth')
# The choice of Sun is [0.01, 0.03,..., 0.31]
psr.add_argument('--R', dest='R', type=float, default=-1, help='The radius the calculation in the evolution, negative means use the whole core')
psr.add_argument('--merge', dest='merge', default=False, action='store_true', help='merge the results')
psr.add_argument('--merge_f', nargs='+', dest='merge_f', help='merge the results')

args = psr.parse_args()

# The grid of theta_12, theta_13, E/Delta_m12 [MeV/eV^2]
s_theta_12_2s = np.arange(0.05, 1, 0.01)
s_theta_13_2s = np.arange(0.015, 0.03, 0.005)
log_E_Delta_m12s = np.arange(3, 11, 0.05)
# construct the oscillation calculator
p_3nu = P_3nu(s_theta_12_2s, s_theta_13_2s, 10 ** log_E_Delta_m12s)
# evolution results with the initial neutrino without phase=0
nu_e = np.array([1, 0, 0, 0])
nu_mu = np.array([0, 0, 1, 0])

if args.mode=='sun':
    # load the model of sun
    with h5py.File(args.ipt, 'r') as ipt:
        ssm, flux = ipt['SSM'][:], ipt['Flux'][:]

    sun = Sun(ssm, flux)
    sun.setOsc(p_3nu)
    if args.R<0:
        # evolution from different radius of the core of Sun
        if not args.merge:
            prob_sun_2nu = sun.evolute(nu_e)
        else:
            prob_sun_2nu = sun.merge(args.merge_f)
        prob_sun_b8_2nu = sun.getProb('B8', prob_sun_2nu)
        prob_sun_hep_2nu = sun.getProb('hep', prob_sun_2nu)
        prob_sun_b8 = p_3nu.P_e_3nu(prob_sun_b8_2nu)
        prob_sun_hep = p_3nu.P_e_3nu(prob_sun_hep_2nu)
    else:
        prob_sun_2nu = sun.evoluteR(nu_e, args.R)

    with h5py.File(args.opt, 'w') as opt:
        opt.create_dataset('sun/prob_2nu', data=prob_sun_2nu, compression='gzip')
        if args.R<0:
            opt.create_dataset('sun/B8_2nu', data=prob_sun_b8_2nu, compression='gzip')
            opt.create_dataset('sun/hep_2nu', data=prob_sun_hep_2nu, compression='gzip')
            opt.create_dataset('sun/B8', data=prob_sun_b8, compression='gzip')
            opt.create_dataset('sun/hep', data=prob_sun_hep, compression='gzip')
else:
    # load the model of earth (PREM)
    earth = Earth()
    earth.setOsc(p_3nu)
    # evolution from different of zenith angle of the earth
    prob_earth_2nu = earth.evolute(nu_e)
    prob_earth = p_3nu.P_e_3nu(prob_earth_2nu)
    with h5py.File(args.opt, 'w') as opt:
        opt.create_dataset('earth/prob_2nu', data=prob_earth_2nu, compression='gzip')
        opt.create_dataset('earth/prob', data=prob_earth, compression='gzip')

# the probability in the night is the multiplication of the sun (day) and earth
# prob_night = prob_sun * prob_earth

