"""A module to generate fake spectra and photometry for WFSS"""
import glob
import os
from pkg_resources import resource_filename
import random
import time

from astroquery.simbad import Simbad
import astropy.units as q
from astropy.io import fits
import numpy as np
from sedkit import Spectrum, SED, Catalog, BTSettl, ModelSpectrum
from sedkit.relations import DwarfSequence
from sedkit.modelgrid import ModelGrid
import sedkit.utilities as u
from svo_filters import Filter
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column

# Customize Simbad query
Simbad.ROW_LIMIT = -1
Simbad.TIMEOUT = 300
Simbad.add_votable_fields('otype', 'sptype', 'id', 'ubv', 'fluxdata(V)', 'fluxdata(J)')


def generate_data(catalog_obj, outdir, photometry=['NIRCam.F115W', 'NIRCam.F150W', 'NIRCam.F277W', 'NIRCam.F356W', 'NIRCam.F444W']):
    """
    Generate a fake JWST pipeline catalog replete with photometry and spectra
    using a range of model atmospheres

    Parameters
    ----------
    catalog_obj: sedkit.catalog.Catalog
        The catalog to convert into JWST mock data
    outdir: str
        The path for the new files
    """
    # Get the SED Catalog results
    data = catalog_obj.results

    # Generate new source catalog for each bandpass
    for band in photometry:

        # Populate the table
        ids = Column(np.arange(len(data)))
        coords = SkyCoord(ra=data['ra'] * q.deg, dec=data['dec'] * q.deg, frame='icrs')
        is_star = np.ones_like(ids)
        aper_total_vegamag = data[band].astype(np.float64)
        aper_total_vegamag_err = data['{}_unc'.format(band)].astype(np.float64)

        # Write the file
        source_catalog = QTable([ids, coords, is_star, aper_total_vegamag, aper_total_vegamag_err], names=['id', 'sky_centroid', 'is_star', 'aper_total_vegamag', 'aper_total_vegamag_err'])
        catpath = os.path.join(outdir, 'jw01345005001_01101_00012_{}_cat.ecsv'.format(band.replace('.', '_')))
        source_catalog.write(catpath, format='ascii.ecsv', overwrite=True)

    # Open the x1d file
    header = fits.getheader(resource_filename('locals', 'data/template_x1d.fits'))

    # Make x1d file for spectra
    x1d_file = os.path.join(outdir, '{}_x1d.fits'.format(band))
    x1d_hdu = fits.HDUList(fits.PrimaryHDU(header=header))

    # Put spectra in x1d fits file
    for id, source in enumerate(data):

        # Get the spectrum from the SED
        spec = source['SED'].spectra['spectrum'][0]

        # Add the data to the BinTable
        data = fits.BinTableHDU(data=np.rec.array(list(zip(*spec.data)), formats='float32,float32,float32', names='WAVELENGTH,FLUX,ERROR'), header=header)
        data.name = 'EXTRACT1D'

        # Add the BinTable to the HDUList
        x1d_hdu.append(data)

    # Write the x1d file
    x1d_hdu.writeto(x1d_file, overwrite=True)


def field_simulation(ra='14 20 26.6772', dec='+52 57 19.02', bandpass='NIRCam.F356W', photometry=['NIRCam.F115W', 'NIRCam.F150W', 'NIRCam.F200W', 'NIRCam.F277W', 'NIRCam.F356W', 'NIRCam.F444W'], radius=None, outdir=resource_filename('locals', 'data/fake/'), model_dir='/Users/jfilippazzo/Documents/Data/Models/bt-settl_400_10k/'):
    """
    Simulate a JWST pipeline catalog of WFSS data for all the stars in the vicinity
    of the provided target

    Parmaters
    ---------
    target: str
        The name of the central target
    bandpass: str
        The bandpass to simulate
    radius: astropy.units.quantity.Quantity
        The radius of the catalog
    """
    # Get instrument
    instrument = bandpass.split('.')[0]

    # Set search radius
    if radius is None:
        if instrument == 'NIRCam':
            radius = 4.4 * q.arcmin
        if instrument == 'NIRISS':
            radius = 2.2 * q.arcmin

    # Query sources
    coord = SkyCoord(ra=ra, dec=dec, unit=(q.degree, q.degree), frame='icrs')
    results = Gaia.query_object_async(coordinate=coord, width=radius, height=radius)
    results = results[[isinstance(row['teff_val'], np.float32) for row in results]]
    print('{} sources found within {} of ra={}, dec={}'.format(len(results), radius, ra, dec))

    # Add two brown dwarfs to the results
    results.insert_row(random.randint(0, len(results)), {'designation': 'BD1', 'teff_val': 800, 'ra': 14.35, 'dec': 52.960})
    results.insert_row(random.randint(0, len(results)), {'designation': 'BD2', 'teff_val': 1600, 'ra': 14.34, 'dec': 52.962})

    # Make a Catalog
    cat = Catalog()

    # Get the bandpass used with the disperser and the Gaia G-band
    disp_bp = Filter(bandpass)
    gband = Filter('Gaia.G')

    # Get SpT-Teff relation
    spt_teff_rel = DwarfSequence()

    # Get the modelgrid
    mg = ModelGrid('BT-Settl', ['alpha', 'logg', 'teff', 'meta'], q.AA, q.erg/q.s/q.cm**2/q.AA)
    mg.load(model_dir)

    # Make an SED for each and add it to the catalog
    for star in results:

        # Make the SED
        new_sed = SED(name=star['designation'])
        new_sed.ra = star['ra'] * q.deg
        new_sed.dec = star['dec'] * q.deg
        new_sed.evo_model = 'parsec12_solar'

        # Add the Gaia G band magnitude
        new_sed.find_Gaia(include=['photometry'])

        if star['parallax'] > 0:
            new_sed.parallax = star['parallax'] * q.mas, star['parallax_error'] * q.mas
        new_sed._calibrate_photometry()

        # Get Teff and SpT
        teff = star['teff_val']
        spt_teff_rel.derive(xparam='Teff', yparam='spt', order=1, xrange=(teff - 200, teff + 200))
        new_sed.spectral_type = spt_teff_rel.evaluate(teff)[0].round(1)

        # Make a spectrum for the given effective temperature
        new_sed.Teff_Gaia = teff
        mod_spec = mg.get_spectrum(teff=teff, closest=True)
        wav, flx = mod_spec.spectrum

        # Add uncertainty
        snr = 10
        unc = flx / snr
        flx = np.random.normal(loc=flx, scale=unc) * flx.unit
        raw_spec = Spectrum(wav, flx, unc)

        # Resample the raw spectrum to the bandpass wavelengths (with wings)
        bp_wave_ext = disp_bp.wave[0].to(raw_spec.wave_units)
        bp_wave_ext[0] *= 0.999
        bp_wave_ext[-1] *= 1.0001
        new_spec = raw_spec.resamp(bp_wave_ext, name='NIRCam WFSS (F356W)')

        # Add the resampled spectrum to the SED
        new_sed.add_spectrum(new_spec)

        # Normalize raw model spectrum to photometry and
        # get synthetic mags as the direct image aperture photometry
        norm = u.mag2flux(gband, new_sed.photometry[0]['app_magnitude'])[0] / raw_spec.synthetic_flux(gband, force=True)[0]
        norm_spec = Spectrum(wav, flx * norm, unc * norm)
        for band in photometry:
            bp = Filter(band)
            mag, mag_unc = norm_spec.synthetic_magnitude(bp)
            new_sed.add_photometry(bp, mag, mag_unc)

        # Do the model fit
        new_sed.fit_modelgrid(mg)

        # Drop the Gaia band
        new_sed.drop_photometry(0)

        # Add the SED to the catalog
        cat.add_SED(new_sed)
        new_sed.make_sed()

        # Pause so astroquery doesn't abort (no more than 6 calls per second)
        time.sleep(0.1)

    if outdir is None:
        return results, cat

    else:
        generate_data(cat, outdir, photometry=photometry)


# def field_simulation_Simbad(target, bandpass, radius=None):
#     """
#     Simulate a JWST pipeline catalog of WFSS data for all the stars in the vicinity
#     of the provided target
#
#     Parmaters
#     ---------
#     target: str
#         The name of the central target
#     bandpass: str
#         The bandpass to simulate
#     radius: astropy.units.quantity.Quantity
#         The radius of the catalog
#     """
#     # Get instrument
#     instrument = bandpass.split('.')[0]
#
#     # Set search radius
#     if radius is None:
#         if instrument == 'NIRCam':
#             radius = 4.4 * q.arcmin
#         if instrument == 'NIRISS':
#             radius = 2.2 * q.arcmin
#
#     # Query sources
#     results = Simbad.query_region(target, radius=radius)
#     results = results[[i['OTYPE'].decode("utf-8") == 'Star' for i in results]]
#     print('{} sources found within {} of {}'.format(len(results), radius, target))
#
#     # Make a Catalog
#     cat = Catalog()
#
#     # Get the bandpass
#     bp = Filter(bandpass)
#
#     # Get SpT-Teff relation
#     spt_teff_rel = DwarfSequence()
#     spt_teff_rel.derive(xparam='spt', yparam='Teff', order=9)
#
#     # Get the modelgrid
#     mg = ModelGrid('BT-Settl', ['alpha', 'logg', 'teff', 'meta'], q.AA, q.erg/q.s/q.cm**2/q.AA)
#     mg.load('/Users/jfilippazzo/Desktop/bt-settl/')
#
#     # Make an SED for each and add it to the catalog
#     for star in results:
#
#         # Make the SED from the MAIN_ID and find photometry in 2MASS
#         new_sed = SED(name=star['MAIN_ID'], method_list=['find_Gaia', 'find_2MASS'])
#         new_sed._calibrate_photometry()
#         phot = new_sed.photometry
#
#         # Determine spectral type
#         if new_sed.spectral_type is None:
#
#             # If no photometry, pick one at random
#             if phot is None:
#
#                 # Set random spectral type (A0 - M9)
#                 new_sed.spectral_type = random.randint(20, 69)
#
#             # Otherwise, infer from color-color relations
#             else:
#
#                 # # 2MASS photometry
#                 # jmag, junc = list(phot[phot['band'] == '2MASS.J'][0][['app_magnitude', 'app_magnitude_unc']])
#                 # hmag, hunc = list(phot[phot['band'] == '2MASS.H'][0][['app_magnitude', 'app_magnitude_unc']])
#                 # kmag, kunc = list(phot[phot['band'] == '2MASS.Ks'][0][['app_magnitude', 'app_magnitude_unc']])
#                 #
#                 # # Colors
#                 # j_h, j_h_unc = jmag - hmag, np.sqrt(junc**2 + hunc**2)
#                 # j_k, j_k_unc = jmag - kmag, np.sqrt(junc**2 + kunc**2)
#                 # h_k, h_k_unc = hmag - kmag, np.sqrt(hunc**2 + kunc**2)
#                 #
#                 # # Infer effective temperature from color
#                 # ds = DwarfSequence()
#                 # ds.derive(xparam='J-H', yparam='spt', order=3)
#                 # new_sed.spectral_type = ds.evaluate(j_h)
#                 new_sed.spectral_type = random.randint(20, 69)
#
#         # Make a spectrum for the measured (or randomized) spectral type
#         teff = spt_teff_rel.evaluate(new_sed.spectral_type[0])[0]
#         mod_spec = ModelSpectrum(teff=teff)
#         wav, flx = mod_spec.spectrum
#
#         # Add uncertainty
#         snr = 10
#         unc = flx / snr
#         flx = np.random.normal(loc=flx, scale=unc) * flx.unit
#         raw_spec = Spectrum(wav, flx, unc)
#
#         # Normalize raw model spectrum to photometry and
#         # get synthetic mag as the direct image aperture photometry
#         norm_spec = raw_spec.norm_to_mags(phot)
#         mag, mag_unc = norm_spec.synthetic_magnitude(bp, force=True)
#         new_sed.add_photometry(bp, mag, mag_unc)
#
#         # Resample the raw spectrum to the bandpass wavelengths (with wings)
#         bp_wave = bp.wave[0].to(raw_spec.wave_units)
#         bp_wave[0] *= 0.999
#         bp_wave[-1] *= 1.0001
#         new_spec = raw_spec.resamp(bp_wave)
#
#         # Add the resampled spectrum to the SED
#         new_sed.add_spectrum(new_spec)
#
#         # Fit model to SED
#         new_sed.fit_modelgrid(mg)
#
#         # Add the SED to the catalog
#         cat.add_SED(new_sed)
#
#         # Pause so astroquery doesn't abort (no more than 6 calls per second)
#         time.sleep(0.1)
#
#     return results, cat


# def generate_data(path=resource_filename('locals', 'data/fake/'), mag_range=(11.13,18)):
#     """Generate a fake JWST pipeline catalog replete with photometry and spectra
#     using a range of model atmospheres
#
#     Parameters
#     ----------
#     path: str
#         The path to the target directory
#     """
#     # Get some random spectra
#     try:
#         files = glob.glob('/user/jfilippazzo/Models/ACES/default/*.fits')[::50]
#     except:
#         files = glob.glob('/Users/jfilippazzo/Documents/Modules/_DEPRECATED/limb_dark_jeff/limb/specint/*.fits')[::20]
#
#     # Make a fake source catalog (with only essential columns for now)
#     catpath = os.path.join(path,'fake_source_catalog.ecsv')
#     ids = list(range(len(files)))
#     coords = SkyCoord([89.7455]*len(ids), [-29.05744]*len(ids), unit='deg', frame='icrs')
#     cat = at.QTable([ids,coords], names=('id','icrs_centroid'))
#     cat.write(catpath)
#
#     # Open the x1d file
#     header = fits.getheader(resource_filename('locals', 'data/template_x1d.fits'))
#
#     # Make Spectrum objects from models at R=150
#     wavelength = np.arange(0.05,2.6,0.0001)[::66]*q.um
#
#     # Normalize the spectra to a random F200W magnitude
#     spectra = []
#     f200w = Filter('NIRISS.F200W')
#     f200w.wave_units = q.um
#     for file in files:
#
#         # Create Spectrum
#         flux = fits.getdata(file)[-1][::66]*q.erg/q.s/q.cm**2/q.AA
#         unc = flux/50.
#         spec = Spectrum(wavelength, flux, unc)
#
#         # Normalize to F200W
#         mag = np.random.uniform(*mag_range)
#         norm_spec = spec.renormalize(mag, f200w)
#         spectra.append(norm_spec)
#
#     # Make a separate x1d file and photometry file for each bandpass
#     # containing data for each source
#     for band in NIRISS_bands:
#
#         try:
#
#             # Get the Filter object
#             bp = Filter(band)
#             bp.wave_units = q.um
#
#             # Make x1d file for spectra
#             x1d_file = os.path.join(path,'{}_x1d.fits'.format(band))
#             x1d_hdu = fits.HDUList(fits.PrimaryHDU(header=header))
#
#             # Make csv file for photometry
#             phot_file = os.path.join(path,'{}_phot.csv'.format(band))
#             phot_data = at.Table(names=('id','band','magnitude','magnitude_unc'), dtype=(int,'S20',float,float))
#
#             # Iterate over spectra
#             for id,(f,spec) in enumerate(zip(files,spectra)):
#
#                 # Trim spectrum to bandpass for x1d file
#                 spec = Spectrum(*spec.spectrum, trim=[(0*q.um,bp.WavelengthMin*1E-4*q.um),(bp.WavelengthMax*1E-4*q.um,10*q.um)])
#
#                 # Calculate magnitude and add to photometry table
#                 mag, mag_unc = spec.synthetic_magnitude(bp, force=True)
#                 phot_data.add_row([id, band, mag, mag_unc])
#
#                 # Add source spectrum params for verification
#                 params = f.split('/')[-1].split('-')
#                 header['TEFF'] = int(params[0].replace('lte',''))
#                 header['LOGG'] = float(params[1][:4])
#                 header['FEH'] = float(params[-6][:-8].split('+')[-1])
#                 header['FILEPATH'] = f
#                 header['PUPIL'] = band
#
#                 # Put spectrum in x1d fits file
#                 data = fits.BinTableHDU(data=np.rec.array(list(zip(*spec.data)),
#                                              formats='float32,float32,float32',
#                                              names='WAVELENGTH,FLUX,ERROR'),
#                                         header=header)
#                 data.name = 'EXTRACT1D'
#
#                 x1d_hdu.append(data)
#
#             # Write the photometry file
#             phot_data.write(phot_file, format='ascii.csv')
#             del phot_data
#
#             # Write the x1d file
#             x1d_hdu.writeto(x1d_file, overwrite=True)
#             del x1d_hdu
#
#         except IOError:
#             pass
