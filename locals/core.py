"""
Low-mass Object Characterization by AnaLyzing Slitless Spectroscopy (LOCALS) is a pure-Python software package which ingests JWST pipeline reduced NIRISS WFSS exposures and outputs a detailed catalog of each detected point source. For each point source in the exposure the software will:
- extract a 1D spectrum from the WFSS trace,
- identify the coordinates of the point source from the undispersed image,
- search Simbad and Vizier for supplemental data (photometry, astrometry, spectral type, etc.) or flag as new object candidate,
- construct an SED from the NIR spectrum and available photometry,
- perform MCMC model fit of the SED to estimate Teff, log(g), and metallicity,
- if distance is known or can be estimated from spectral type, calculate fundamental parameters (Lbol, Teff, mass)
- add point source to the output catalog of all collected data and derived fundamental parameters.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from SEDkit import seds
from astrodbkit import astrodb, astrocat
