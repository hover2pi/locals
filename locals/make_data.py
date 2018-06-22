"""A module to generate fake spectra and photometry for WFSS"""
import os
import numpy as np
import astropy.table as at
import astropy.units as q
import astropy.io.ascii as ii
import glob
from copy import copy
from astropy.io import fits
from astropy.coordinates import SkyCoord
from pkg_resources import resource_filename
from SEDkit.spectrum import Spectrum
from SEDkit.synphot import Bandpass

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
    # files = glob.glob('/user/jfilippazzo/Models/ACES/default/*.fits')[::50]
    files = glob.glob('/Users/jfilippazzo/Documents/Modules/_DEPRECATED/limb_dark_jeff/limb/specint/*.fits')[::20]
     
    # Make a fake source catalog (with only essential columns for now)
    catpath = os.path.join(path,'fake_source_catalog.ecsv')
    ids = list(range(len(files)))
    coords = [SkyCoord(ra=89.7455*q.deg, dec=-29.05744*q.deg, frame='icrs') for i in ids]
    cat = at.QTable([ids,coords], names=('id','icrs_centroid'))
    cat.write(catpath)
    
    # Open the x1d file
    hdu = fits.open(resource_filename('locals', 'data/template_x1d.fits'))
    
    # Make Spectrum objects from models at R=150
    wavelength = np.arange(0.05,2.6,0.0001)[::66]*q.um
    
    # Normalize the spectra to a random F200W magnitude
    spectra = []
    f200w = Bandpass('NIRISS.F200W')
    f200w.wave_units = q.um
    
    # Iterate over sources
    for file in files[:1]:
        
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
        
            # Get the Bandpass object
            bp = Bandpass(band)
            bp.wave_units = q.um
        
            # Make x1d file for spectra
            x1d_file = os.path.join(path,'{}_x1d.fits'.format(band))
            x1d_hdu = fits.HDUList()
            x1d_hdu.append(copy(hdu[0]))
        
            # Make csv file for photometry
            phot_file = os.path.join(path,'{}_phot.csv'.format(band))
            phot_data = at.Table(names=('id','band','magnitude','magnitude_unc'), dtype=(int,'S20',float,float))
        
            # Iterate over spectra
            for id,spec in enumerate(spectra):
            
                # Trim spectrum to bandpass for x1d file
                spec = Spectrum(*spec.spectrum, trim=[(0*q.um,bp.WavelengthMin*1E-4*q.um),(bp.WavelengthMax*1E-4*q.um,10*q.um)])
        
                # Calculate magnitude and add to photometry table
                mag, mag_unc = spec.synthetic_magnitude(bp, force=True)
                phot_data.add_row([id, band, mag, mag_unc])

                # Put spectrum in x1d fits file
                data = fits.BinTableHDU(data=np.rec.array(list(zip(*spec.data)), formats='float32,float32,float32', names='WAVELENGTH,FLUX,ERROR'))
                data.name = 'EXTRACT1D'
                
                x1d_hdu.append(data)
            
            # Write the photometry file
            phot_data.write(phot_file, format='ascii.csv')
            del phot_data
            
            # Write the x1d file
            x1d_hdu.append(hdu[-1])
            x1d_hdu.writeto(x1d_file, overwrite=True)
            del x1d_hdu
            
        except IOError:
            pass
            