"""A module to generate fake spectra and photometry for WFSS"""
from copy import copy
import glob
import os
from pkg_resources import resource_filename

import astropy.table as at
import astropy.units as q
import astropy.io.ascii as ii
from astropy.io import fits
from astropy.coordinates import SkyCoord
import numpy as np
from sedkit.spectrum import Spectrum
from svo_filters import Filter


# List the NIRISS bands available for WFSS mode and their wavelength (min,max) in microns
NIRISS_bands = ['NIRISS.F090W', 'NIRISS.F115W', 'NIRISS.F140M', 'NIRISS.F150W', 'NIRISS.F158M', 'NIRISS.F200W']


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
            