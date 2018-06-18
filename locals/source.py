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
from astrodbkit import astrodb
from SEDkit import sed
from astropy.io import fits
import astropy.table as at
from astroquery.vizier import Vizier
import glob
import astropy.coordinates as coord
import pkg_resources
import astropy.units as q
import h5py

PARAMS = {'II/328/allwise': {'col':'band', 'magnitude':['W1mag','W2mag','W3mag','W4mag'], 'errs':['e_W1mag','e_W2mag','e_W3mag','e_W4mag']},
          'II/246/out': {'col':'band', 'magnitude':['Jmag','Hmag','Kmag'], 'errs':['e_Jmag','e_Hmag','e_Kmag']}
      }
      
VIZ_SVO = {'Jmag':'2MASS.J', 'Hmag':'2MASS.H', 'Kmag':'2MASS.Ks', 'W1mag':'WISE.W1', 'W2mag':'WISE.W2', 'W3mag':'WISE.W3', 'W4mag':'WISE.W4'}

class Source(sed.SED):
    """
    A class to store parsed grizli output
    """
    def __init__(self, ra, dec, filepath=None, name='My Target', **kwargs):
        """
        Initialize the Source object with a file path
        """
        # Inherit from ArraySpectrum
        super().__init__()
        
        self.name = name
        self.filepath = filepath or '/Users/jfilippazzo/Desktop/dataset.hdf5'
        self.search_radius = 15*q.arcsec
        self.sky_coords = ra, dec
        
        # Store metadata from source_list
        for k,v in kwargs.items():
            setattr(self, k, v)
        
        
    @property
    def info(self):
        """
        Print all the source info
        """
        for attr in dir(self):
            if not attr.startswith('_') and attr not in ['info','find_data']:
                print('{0: <10}: {1}'.format(attr, getattr(self, attr)))


    def save(filepath=''):
        """
        Make an HDF5 file with all the data needed to MakeSED
        """
        # Override stored filepath
        self.filepath = filepath or self.filepath
        
        # Create the file if necessary
        if not os.path.isfile(self.filepath):
            f = h5py.File(self.filepath, "w")

