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
        self.x1d_files = glob.glob(os.path.join(self.dirpath,'*_x1d.fits'))
        
        # Make a Source object for each row in the source_list
        for row in self.source_list:
            ra = row['icrs_centroid'].ra
            dec = row['icrs_centroid'].dec
            source = Source(ra=ra, dec=dec, **{k:row[k] for k in row.colnames})
            
            # Look for photometry
            source.find_photometry()

            # Look for distance
            source.find_parallax()
            
            self.source_ids.append(int(source.id))
            self.sources.append(source)
            
        # Ping Julia and Kevin to see if there are aperture phot results with WFSS output to do color cuts here
        
        # Generate JWST colors for BDs 
            
        # Open up x1d files and add spectra to the sources
        # for x1d_file in self.x1d_files:
        #
        #     x1d_hdu = fits.open(x1d_file)
        #     for n in range(len(x1d_hdu)):
        #
        #         # Add the FITS data to the source object
        #         if x1d_hdu[n].name=='EXTRACT1D':
        #
        #             # Get the source_id for this spectrum
        #             source_id = int(x1d_hdu[n].header['SOURCEID'])
        #             source_idx = self.source_ids.index(source_id)
        #             source = self.sources[source_idx]
        #
        #             # Put in dummy data
        #             if dummy:
        #                 fake_data = np.genfromtxt(pkg_resources.resource_filename('locals', 'data/STSci_Vega.txt'), unpack=True)
        #                 w = fake_data[0]*q.um
        #                 f = fake_data[1]*q.erg/q.s/q.cm**2/q.AA
        #
        #             source.add_spectrum(w, f)
        #
        #             self.sources[source_idx] = source
        #
        #     x1d_hdu.close()
        
    # def ingest_sources(self, x1d_file, dummy=True):
    #     """
    #     Ingest the FITS files in the given directory and parse
    #
    #     Parameters
    #     ----------
    #     x1d_file: str
    #         The path to the Level 2 x1d file
    #     dummy: bool
    #         Add dummy flux from file
    #     """
    #     # Get the calibrated 2D spectroscopic data (_cal.fits)
    #     # DO WE NEED THIS?
    #     # self.cal_file = glob.glob(os.path.join(self.dirpath,'*_cal.fits'))[0]
    #
    #     # Get the 1D extracted spectra (_x1d.fits)
    #     self.x1d_files.append(x1d_file)
    #
    #     x1d_hdu = fits.open(x1d_file)
    #     for n in range(len(x1d_hdu)):
    #
    #         # Initialize source object
    #         source = Source()
    #
    #         # Add the data from the source_list
    #         # for col in self.source_list.colnames:
    #         #     print(self.source_list.loc['id'])
    #             # setattr(self, col, self.source_list[col][n])
    #
    #         # Add the FITS data to the source object
    #         if x1d_hdu[n].name=='EXTRACT1D':
    #
    #             # Parse the data
    #             for col in x1d_hdu[n].data.dtype.names:
    #                 setattr(source, col, x1d_hdu[n].data[col])
    #
    #             # Parse the header
    #             for card in x1d_hdu[n].header.cards:
    #                 setattr(source, card[0], card[1])
    #
    #             # Put in dummy data
    #             if dummy:
    #                 fake_data = np.genfromtxt(pkg_resources.resource_filename('locals', 'data/STSci_Vega.txt'), unpack=True)
    #                 source.FLUX = np.interp(source.WAVELENGTH, fake_data[0]/10000., fake_data[1])
    #
    #             self.sources.append(source)
    #
    #     x1d_hdu.close()
    #
    # def xmatch_sources(self, radius=10*q.arcsec, catalogs=['II/246/out','II/328/allwise','V/147/sdss12','II/243/denis']):
    #     """
    #     Run astroquery to get ancillary data for each source and add to it's Source instance
    #     """
    #     for source in self.sources:
    #
    #         # Get the ICRS coordinates from the table and perform Vizier query
    #         result = Vizier.query_region(source['icrs_centroid'], radius=radius, catalog=catalogs)
    #
    #         self.sources[n].raw_data = result
