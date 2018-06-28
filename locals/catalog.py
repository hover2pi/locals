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
import glob
import pkg_resources
import h5py
import astropy.units as q
import astropy.table as at
import astropy.coordinates as coord
import astropy.io.ascii as ii
from astropy.io import fits
from astroquery.vizier import Vizier
from SEDkit import sed, catalog
from . import colors, make_data, source
    

class SourceCatalog(catalog.SEDCatalog):
    """
    A class to ingest a JWST pipeline output to produce a source catalog
    """
    
    def __init__(self, dirpath, color_cut=None, verbose=False):
        """
        Initialize the SourceCatalog object
        
        Parameters
        ----------
        dirpath: str
            The path to the JWST pipeline output
        color_cut: str
            The name of the color cuts to use
        """
        # Inherit from SEDkit.catalog.SEDCatalog
        super().__init__()
        
        # The path to the pipeline output directory
        self.dirpath = dirpath
        self.verbose = verbose
        
        # Get the source catalog (_cat.ecsv)
        self.cat_file = glob.glob(os.path.join(self.dirpath,'*.ecsv'))[0]
        self.source_list = at.Table.read(self.cat_file, format='ascii.ecsv')
        self.x1d_files = glob.glob(os.path.join(self.dirpath,'*_x1d.fits'))
        self.phot_files = glob.glob(os.path.join(self.dirpath,'*_phot.csv'))
        
        # Put all photometry into one table
        self.photometry = at.vstack([ii.read(f) for f in self.phot_files])
        
        # Make a Source object for each row in the source_list
        for n,row in enumerate(self.source_list):
            ra = row['icrs_centroid'].ra
            dec = row['icrs_centroid'].dec
            name = 'Source {}'.format(row['id'])
            src = source.Source(ra=ra, dec=dec, name=name,
                                verbose=self.verbose,
                                **{k:row[k] for k in row.colnames})
            
            # Add the JWST photometry for this source
            for phot in self.photometry:
                if phot['id']==row['id']:
                    src.add_photometry(phot['band'], phot['magnitude'],
                                       phot['magnitude_unc'])
            
            # # Look for photometry (Need real coordinates for this)
            # src.find_SDSS()
            # src.find_2MASS()
            # src.find_WISE()
            # src.find_PanSTARRS()
            
            # Look for distance
            # src.find_Gaia()
            
            # Check to see if the source makes the color cut
            src.photometry.add_index('band')
            keep = colors.in_color_range(src.photometry, color_cut)
            if keep:
                
                # Add observed WFSS spectra to the source
                for x1d in self.x1d_files:
                    funit = q.erg/q.s/q.cm**2/q.AA
                    src.add_spectrum_file(x1d, q.um, funit, ext=n+1)
                    
                # Fit a blackbody
                src.fit_blackbody()
                    
                # Add the source to the catalog
                self.add_SED(src)
                
        print("{}/{} sources added to catalog{}"\
              .format(len(self.results), len(self.source_list),
              " after applying '{}' color cuts".format(color_cut)
              if color_cut is not None else ''))
            
            
    # @property
    # def catalog(self):
    #     """
    #     Generate a table of the results for every source
    #     """
    #     # Make a table for each result
    #     tables = []
    #     for src in self.sources:
    #
    #         # Get the values
    #         res = src.results
    #         names = ['{} [{}]'.format(i,j) for i,j in zip(list(res['param']),list(res['units']))]
    #         values = at.Column(['{} +/- {}'.format(i,j) if isinstance(i,(float,int)) else i for i,j in zip(list(res['value']),list(res['unc']))])
    #
    #         # Fix some values
    #         names[0] = 'name'
    #
    #         # Make the table
    #         r_table = at.Table(values, names=names)
    #
    #         tables.append(r_table)
    #
    #     # Make master table
    #     final = at.vstack(tables)
    #
    #     return final
    
