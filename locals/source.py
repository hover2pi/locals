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
        self.ra = ra
        self.dec = dec
        self.icrs_centroid = None
        
        # Store metadata from source_list
        for k,v in kwargs.items():
            setattr(self, k, v)
        
        # Set ra and dec
        if isinstance(ra, q.quantity.Quantity) and isinstance(dec, q.quantity.Quantity):
            self.icrs_centroid = coord.SkyCoord(ra=ra, dec=dec, frame='icrs')
        
        
    @property
    def info(self):
        """
        Print all the source info
        """
        for attr in dir(self):
            if not attr.startswith('_') and attr not in ['info','find_data']:
                print('{0: <10}: {1}'.format(attr, getattr(self, attr)))
                
            
    def find_WISE(self, search_radius=None, catalogs=['II/328/allwise']):
        """
        Search for source parallax
        """
        # Make sure there are coordinates
        if isinstance(self.icrs_centroid, coord.sky_coordinate.SkyCoord):
            
            viz_cat = Vizier.query_region(self.icrs_centroid, radius=search_radius or self.search_radius, catalog=catalogs)
            
            if len(viz_cat)>0:
                wise = viz_cat[0][0]
            
                for band in ['W1mag','W2mag','W3mag','W4mag']:
                    try:
                        mag, unc = list(wise[[band,'e_'+band]])
                        self.add_photometry(VIZ_SVO.get(band), float(mag), float(unc))
                    except IOError:
                        pass
            
        else:
            print("Can't find WISE photometry without coordinates!")
            
    def find_2MASS(self, search_radius=None, catalogs=['II/246/out']):
        """
        Search for source parallax
        """
        # Make sure there are coordinates
        if isinstance(self.icrs_centroid, coord.sky_coordinate.SkyCoord):
            
            viz_cat = Vizier.query_region(self.icrs_centroid, radius=search_radius or self.search_radius, catalog=catalogs)
            
            if len(viz_cat)>0:
                tmass = viz_cat[0][0]
            
                for band in ['Jmag','Hmag','Kmag']:
                    try:
                        mag, unc = list(tmass[[band,'e_'+band]])
                        self.add_photometry(VIZ_SVO.get(band), float(mag), float(unc))
                    except IOError:
                        pass
            
        else:
            print("Can't find 2MASS photometry without coordinates!")
            
            
    def find_photometry(self, search_radius=None, catalogs=['II/246/out','II/328/allwise']):#,'V/147/sdss12','II/243/denis']):
        """
        Search for source photometry
        """
        # Make sure there are coordinates
        if isinstance(self.icrs_centroid, coord.sky_coordinate.SkyCoord):
            
            # Get the ICRS coordinates from the table and perform Vizier query
            self.raw_photometry = Vizier.query_region(self.icrs_centroid, radius=search_radius or self.search_radius, catalog=catalogs)
            keys = self.raw_photometry.keys()
        
            # Process the photometry
            new_tables = []
            for name,table in zip(keys,self.raw_photometry):
            
                new_tables.append(self._process_table(name, table))
            
            # Put the data together
            if new_tables:
                photometry = at.vstack(new_tables)
            
                # Get the SVO names so the Bandpass object can be created
                photometry['band'] = [VIZ_SVO.get(b) for b in photometry['band']]
            
                # Add them to the photometry table
                for row in photometry:
                    data = list(row)
                    self.add_photometry(data[0], float(data[1]), float(data[2]))
            
        else:
            print("Can't find photometry without coordinates!")
            
            
    def find_parallax(self, search_radius=None, catalogs=['I/345/gaia2']):
        """
        Search for source parallax
        """
        # Make sure there are coordinates
        if isinstance(self.icrs_centroid, coord.sky_coordinate.SkyCoord):
            parallaxes = Vizier.query_region(self.icrs_centroid, radius=search_radius or self.search_radius, catalog=catalogs)
            
            if parallaxes:
                # Get the ICRS coordinates from the table and perform Vizier query
                self.raw_parallaxes = parallaxes
            
                parallax = list(self.raw_parallaxes[0][0][['Plx','e_Plx']])
            
                self.parallax = parallax[0]*q.mas, parallax[1]*q.mas
            
                # Get Gband while we're here
                try:
                    mag, unc = list(self.raw_parallaxes[0][0][['Gmag','e_Gmag']])
                    self.add_photometry('Gaia.G', mag, unc)
                except:
                    pass
            
        else:
            print("Can't find parallaxes without coordinates!")
        
        
    @staticmethod
    def _process_table(name, table):
        """Split a row of multiple measurements into a table
        
        Parameters
        ----------
        name: str
            The name of the table
        table: astropy.table.Table
            The table of data
        
        Returns
        -------
        astropy.table.Table
            The new table
        """
        # Get the table params
        params = PARAMS.get(name)
        
        if params:
        
            # Make empty table
            new_table = at.Table()
        
            # Get the values to seed
            for k,v in params.items():
                if k not in ['col','errs']:
                
                    # Get the measurement column
                    new_col = k
                
                    # Add the data
                    new_col_keys = params[new_col]
                    new_col_vals = list(table[new_col_keys][0])
                    new_col_errs = list(table[params['errs']][0])
        
            # Add a new column with the keys to split
            new_table[params['col']] = new_col_keys
        
            # Add a new column with the values and errors
            new_table[new_col] = new_col_vals
            new_table[new_col+'_unc'] = new_col_errs
        
            return new_table
            
        else:
            print('Please add {} to PARAMS.'.format(name))


    def save(filepath=''):
        """
        Make an HDF5 file with all the data needed to MakeSED
        """
        # Override stored filepath
        self.filepath = filepath or self.filepath
        
        # Create the file if necessary
        if not os.path.isfile(self.filepath):
            f = h5py.File(self.filepath, "w")

