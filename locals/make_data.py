"""A module to generate fake spectra and photometry for WFSS"""
from copy import copy
import glob
import os
from pkg_resources import resource_filename
from random import randrange

from astroquery.simbad import Simbad
import astropy.table as at
import astropy.units as q
import astropy.io.ascii as ii
from astropy.io import fits
from astropy.coordinates import SkyCoord
import numpy as np
from sedkit import Spectrum, SED, Catalog
from sedkit.relations import DwarfSequence
from svo_filters import Filter


# List the NIRISS bands available for WFSS mode and their wavelength (min,max) in microns
NIRISS_bands = ['NIRISS.F090W', 'NIRISS.F115W', 'NIRISS.F140M', 'NIRISS.F150W', 'NIRISS.F158M', 'NIRISS.F200W']

# Customize Simbad query
Simbad.ROW_LIMIT = -1
Simbad.TIMEOUT = 300
Simbad.add_votable_fields('otype', 'sptype', 'id', 'ubv', 'fluxdata(V)', 'fluxdata(J)')


def field_simulation(target, instrument='NIRISS', radius=5*q.arcmin):
    """
    Simulate a JWST pipeline catalog of WFSS data for all the stars in the vicinity
    of the provided target

    Parmaters
    ---------
    target: str
        The name of the central target
    instrument: str
        The instrument to simulate, ['NIRISS', 'NIRCam']
    radius: astropy.units.quantity.Quantity
        The radius of the catalog
    """
    # Query sources
    results = Simbad.query_region(target, radius=radius)
    print('{} sources found within {} of {}'.format(len(results), radius, target))

    # Make a Catalog
    cat = Catalog()

    # Make an SED for each and add it to the catalog
    for star in results:

        # Omit planets
        if 'planet' not in star['OTYPE'].decode("utf-8").lower():

            # Make the SED from the MAIN_ID and find photometry in 2MASS
            x = SED(name=star['MAIN_ID'], method_list=['find_2MASS'])
            phot = x.photometry

            # Make sure there is 2MASS photometry
            if phot is not None:

                # If no spectral type, infer from 2MASS
                if x.spectral_type is None:

                    # 2MASS photometry
                    jmag, junc = list(p[p['band'] == '2MASS.J'][0][['app_magnitude', 'app_magnitude_unc']])
                    hmag, hunc = list(p[p['band'] == '2MASS.H'][0][['app_magnitude', 'app_magnitude_unc']])
                    kmag, kunc = list(p[p['band'] == '2MASS.Ks'][0][['app_magnitude', 'app_magnitude_unc']])

                    # Colors
                    j_h, j_h_unc = jmag - hmag, np.sqrt(junc**2 + hunc**2)
                    j_k, j_k_unc = jmag - kmag, np.sqrt(junc**2 + kunc**2)
                    h_k, h_k_unc = hmag - kmag, np.sqrt(hunc**2 + kunc**2)

                    # Infer spectral type from color
                    ds = DwarfSequence('')

                # Require spectral type and 2MASS photometry
            

                # Make a spectrum for the measured (or randomized) spectral type
                # spt = x.spectral_type[0]

                # Resample to the GR150 wavelengths
            

                # Add it to the SED so it gets normalized to the photometry
            

                # Calculate the synthetic photometry for the NIRISS bands
            

                # Add the SED to the catalog
                cat.add_SED(x)

    return results, cat

def generate_data(path=resource_filename('locals', 'data/fake/'), mag_range=(11.13,18)):
    """Generate a fake JWST pipeline catalog replete with photometry and spectra
    using a range of model atmospheres

    Parameters
    ----------
    path: str
        The path to the target directory
    """
    # Get some random spectra
    try:
        files = glob.glob('/user/jfilippazzo/Models/ACES/default/*.fits')[::50]
    except:
        files = glob.glob('/Users/jfilippazzo/Documents/Modules/_DEPRECATED/limb_dark_jeff/limb/specint/*.fits')[::20]

    # Make a fake source catalog (with only essential columns for now)
    catpath = os.path.join(path,'fake_source_catalog.ecsv')
    ids = list(range(len(files)))
    coords = SkyCoord([89.7455]*len(ids), [-29.05744]*len(ids), unit='deg', frame='icrs')
    cat = at.QTable([ids,coords], names=('id','icrs_centroid'))
    cat.write(catpath)

    # Open the x1d file
    header = fits.getheader(resource_filename('locals', 'data/template_x1d.fits'))

    # Make Spectrum objects from models at R=150
    wavelength = np.arange(0.05,2.6,0.0001)[::66]*q.um

    # Normalize the spectra to a random F200W magnitude
    spectra = []
    f200w = Filter('NIRISS.F200W')
    f200w.wave_units = q.um
    for file in files:

        # Create Spectrum
        flux = fits.getdata(file)[-1][::66]*q.erg/q.s/q.cm**2/q.AA
        unc = flux/50.
        spec = Spectrum(wavelength, flux, unc)

        # Normalize to F200W
        mag = np.random.uniform(*mag_range)
        norm_spec = spec.renormalize(mag, f200w)
        spectra.append(norm_spec)

    # Make a separate x1d file and photometry file for each bandpass
    # containing data for each source
    for band in NIRISS_bands:

        try:

            # Get the Filter object
            bp = Filter(band)
            bp.wave_units = q.um

            # Make x1d file for spectra
            x1d_file = os.path.join(path,'{}_x1d.fits'.format(band))
            x1d_hdu = fits.HDUList(fits.PrimaryHDU(header=header))

            # Make csv file for photometry
            phot_file = os.path.join(path,'{}_phot.csv'.format(band))
            phot_data = at.Table(names=('id','band','magnitude','magnitude_unc'), dtype=(int,'S20',float,float))

            # Iterate over spectra
            for id,(f,spec) in enumerate(zip(files,spectra)):

                # Trim spectrum to bandpass for x1d file
                spec = Spectrum(*spec.spectrum, trim=[(0*q.um,bp.WavelengthMin*1E-4*q.um),(bp.WavelengthMax*1E-4*q.um,10*q.um)])

                # Calculate magnitude and add to photometry table
                mag, mag_unc = spec.synthetic_magnitude(bp, force=True)
                phot_data.add_row([id, band, mag, mag_unc])

                # Add source spectrum params for verification
                params = f.split('/')[-1].split('-')
                header['TEFF'] = int(params[0].replace('lte',''))
                header['LOGG'] = float(params[1][:4])
                header['FEH'] = float(params[-6][:-8].split('+')[-1])
                header['FILEPATH'] = f
                header['PUPIL'] = band

                # Put spectrum in x1d fits file
                data = fits.BinTableHDU(data=np.rec.array(list(zip(*spec.data)),
                                             formats='float32,float32,float32',
                                             names='WAVELENGTH,FLUX,ERROR'),
                                        header=header)
                data.name = 'EXTRACT1D'

                x1d_hdu.append(data)

            # Write the photometry file
            phot_data.write(phot_file, format='ascii.csv')
            del phot_data

            # Write the x1d file
            x1d_hdu.writeto(x1d_file, overwrite=True)
            del x1d_hdu

        except IOError:
            pass
            