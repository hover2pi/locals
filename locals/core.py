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
from SEDkit import sed
from astrodbkit import astrodb, astrocat
from astropy.io import fits
import astropy.table as at
from astroquery.vizier import Vizier
import glob
import astropy.coordinates as coord

def test():
    """
    A function to test SourceCatalog creation
    """
    
    # Initialize empty source catclog
    cat = SourceCatalog('/grp/jwst/ssb/test_build7.1/examples_for_dms/level2/nrc_grism/')
    
    # Ingest the source data
    x1d_file = glob.glob(os.path.join(cat.dirpath,'*_x1d.fits'))[0]
    cat.ingest_sources(x1d_file)
    
    return cat

class Source(object):
    """
    A class to store parsed grizli output
    """
    def __init__(self):
        """
        Initialize the Source object with a file path
        """
        self.id = ''

class SourceCatalog(object):
    """
    A class to ingest a JWST pipeline output to produce a source catalog
    """
    
    def __init__(self, dirpath):
        """
        Initialize the SourceCatalog object
        
        Parameters
        ----------
        dirpath: str
            The path to the JWST pipeline output
        """
        # The path to the pipeline output directory
        self.dirpath = dirpath
        
        # Get the source catalog (_cat.ecsv)
        self.cat_file = glob.glob(os.path.join(self.dirpath,'*.ecsv'))[0]
        self.source_list = at.Table.read(self.cat_file, format='ascii.ecsv')
        
        self.sources = []
        self.x1d_files = []
        
    def ingest_sources(self, x1d_file):
        """
        Ingest the FITS files in the given directory and parse
        """
        # Get the calibrated 2D spectroscopic data (_cal.fits)
        # DO WE NEED THIS?
        # self.cal_file = glob.glob(os.path.join(self.dirpath,'*_cal.fits'))[0]
        
        # Get the 1D extracted spectra (_x1d.fits)
        self.x1d_files.append(x1d_file)
        
        x1d_hdu = fits.open(x1d_file)
        for n in range(len(x1d_hdu)):

            # Initialize source object
            source = Source()

            # Add the FITS data to the source object
            if x1d_hdu[n].name=='EXTRACT1D':

                # Parse the data
                for col in x1d_hdu[n].data.dtype.names:
                    setattr(source, col, x1d_hdu[n].data[col])

                # Parse the header
                for card in x1d_hdu[n].header.cards:
                    setattr(source, card[0], card[1])

                self.sources.append(source)

        x1d_hdu.close()
            
    def xmatch_sources(self, radius=0.00001):
        """
        Run astroquery to get ancillary data for each source and add to it's Source instance
        """
        pass
        # crds = 
        # data = Vizier.query_region()
