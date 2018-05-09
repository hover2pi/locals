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
from .source import Source
from astropy.io import fits
import astropy.table as at
from astroquery.vizier import Vizier
import glob
import astropy.coordinates as coord
import pkg_resources
import astropy.units as q
import h5py


class SourceCatalog(object):
    """
    A class to ingest a JWST pipeline output to produce a source catalog
    """
    
    def __init__(self, dirpath, dummy=True):
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
        self.source_ids = []
        self.x1d_file = glob.glob(os.path.join(self.dirpath,'*_x1d.fits'))[0]
        
        # Make a Source object for each row in the source_list
        for n,row in enumerate(self.source_list):
            ra = row['icrs_centroid'].ra
            dec = row['icrs_centroid'].dec
            source = Source(ra=ra, dec=dec, **{k:row[k] for k in row.colnames})
            
            # Look for photometry
            source.find_photometry()

            # Look for distance
            source.find_parallax()
            
            # Add observed JWST photometry to the source
            # TODO
            
            # Add observed WFSS spectrum to the source
            wave_units = q.um
            flux_units = q.erg/q.s/q.cm**2/q.AA
            source.add_spectrum_file(self.x1d_file, wave_units, flux_units, ext=('EXTRACT1D', n+1))
            
            # Add the source to the catalog
            self.source_ids.append(int(source.id))
            self.sources.append(source)


def make_dummy_data(x1d_file):
    """Put real spectra into the x1d file for testing"""
    import copy
    
    hdu = fits.open(x1d_file)
    
    # LHS 2924
    dummy_file = '/Users/jfilippazzo/Dropbox/BDNYC_spectra/SpeX/IRTF Library (Prism+LXD)/M9V_LHS2924.fits'
    dummy_data = fits.getdata(dummy_file)
    dummy_rec = np.rec.array(list(zip(*dummy_data)), formats='float32,float32,float32', names='WAVELENGTH,FLUX,ERROR')
    hdu[2].data = dummy_rec
    
    # APMP 0559 (subdwarf!)
    dummy_file = '/Users/jfilippazzo/Dropbox/BDNYC_spectra/SpeX/Prism/! Subdwarf spectra that need source_ids (not in db yet!)/spex-prism_APMPM0559-2903_20051231_BUR06A.txt'
    dummy_data = np.genfromtxt(dummy_file, unpack=True)
    dummy_rec = np.rec.array(list(zip(*dummy_data)), formats='float32,float32,float32', names='WAVELENGTH,FLUX,ERROR')
    hdu[1].data = dummy_rec
    
    # Write the new file
    hdu.writeto(x1d_file.replace('_x1d','_dummy_x1d'), overwrite=True)
    
    hdu.close()
    
