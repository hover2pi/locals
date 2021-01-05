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
import glob

import astropy.units as q
import astropy.table as at
from astropy.io import fits
from sedkit import SED, Catalog, BTSettl

from . import colors


class SourceCatalog(Catalog):
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
        self.cat_file = glob.glob(os.path.join(self.dirpath, '*.ecsv'))[0]
        self.source_list = at.Table.read(self.cat_file, format='ascii.ecsv')
        self.x1d_files = glob.glob(os.path.join(self.dirpath, '*_x1d.fits'))

        # Load BT Settl grid
        bt = BTSettl(resolution=1000, trim=(0.3 * q.um, 5 * q.um))

        # Make a Source object for each row in the source_list
        for n, row in enumerate(self.source_list):

            # Only process stars
            if row['is_star']:

                # Metadata
                ra = row['sky_centroid'].ra
                dec = row['sky_centroid'].dec
                name = 'Source {}'.format(row['id'])
                src = SED(ra=ra, dec=dec, name=name, verbose=self.verbose, **{k: row[k] for k in row.colnames})
                print("SED created for {}".format(name))

                # JWST photometry
                mag = row['aper_total_vegamag']
                err = row['aper_total_vegamag_err']
                src.add_photometry('NIRCam.F356W', mag, err)

                # Check for all-sky catalog photometry
                src.find_SDSS()
                src.find_2MASS()
                src.find_WISE()

                # Look for distance
                src.find_Gaia()

                # Check to see if the source makes the color cut
                src.photometry.add_index('band')
                keep = colors.in_color_range(src.photometry, color_cut)
                if keep:

                    # Add observed WFSS spectra to the source
                    for x1d in self.x1d_files:
                        header = fits.getheader(x1d, ext=n+1)
                        funit = q.erg / q.s / q.cm**2 / q.AA

                        # Get extracted spectrum from x1d file
                        data = fits.getdata(x1d)
                        src.add_spectrum([data['WAVELENGTH'] * q.um, data['FLUX'] * funit, data['ERROR'] * funit])

                        # # Save the params for verification
                        # src.Teff_model = header['TEFF']
                        # src.logg_model = header['LOGG']
                        # src.FeH_model = header['FEH']

                    # # Fit a blackbody
                    # # src.fit_blackbody()
                    #
                    # # Fit spectral type
                    # src.fit_spectral_type()

                    # Fit model grid
                    src.fit_modelgrid(bt)

                    # Add the source to the catalog
                    self.add_SED(src)

        print("{}/{} sources added to catalog {}".format(len(self.results), len(self.source_list), " after applying '{}' color cuts".format(color_cut) if color_cut is not None else ''))
