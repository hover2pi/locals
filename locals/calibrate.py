#! /usr/bin/env python

"""
Calibrate WFSS _uncal.fits files through to _x1d.fits files for analysis with locals
"""
__author__ = "Joe Filippazzo"
__date__ = "11-18-2020"
__maintainer__ = "Joe Filippazzo"
__email__ = "jfilippazzo@stsci.edu"

import glob
import os
from configobj import ConfigObj

from astropy.io import fits
from sedkit.utilities import isnumber
from jwst.pipeline import collect_pipeline_cfgs, Detector1Pipeline, Image2Pipeline, Image3Pipeline, Spec2Pipeline

TERM = os.get_terminal_size()
SYM = "~"


def printnice(level):
    """
    Print the current pipeline level and make it nice

    Parameters
    ----------
    level: str
        The thing you want to print nicely.
    """
    print("\n", SYM*TERM.columns, "\n", SYM, level.center(TERM.columns-2), SYM, "\n", SYM*TERM.columns, "\n", sep="")


def process_image(direct, configdir='.', outdir=None, params={}):
    """
    Calibrate a direct JWST image

    Parameters
    ----------
    direct: str
        Path to direct image
    configdir: str
        Path to directory with config files
    outdir: str
        Path to directory where output files should be written
    params: dict
        The key/value pairs of pipeline/step to override in the calibration
    """
    # Make sure cfg file are there
    collect_pipeline_cfgs.collect_pipeline_cfgs(configdir)

    if outdir is None:
        outdir = os.path.dirname(direct)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("Created directory {}".format(outdir))

    # Ensure that the source_catalog step has save_results set to True
    sc = ConfigObj(os.path.join(configdir, "source_catalog.cfg"))
    if "save_results" not in sc:
        print("Modifying {} to add save_results option".format(os.path.join(configdir, "source_catalog.cfg")))
        sc["save_results"] = True
        sc.write()

    # Run the direct image through DETECTOR1
    printnice("DETECTOR1")
    image1 = Detector1Pipeline(save_results=True, config_file=os.path.join(configdir, "calwebb_detector1.cfg"), output_dir=outdir)
    image1 = mod_pipe_steps(image1, params.get('DETECTOR1', []))
    image1.run(direct)
    image1_file = direct.replace('_uncal.fits', '_rate.fits')

    # Run the direct image through IMAGE2
    printnice("IMAGE2")
    image2 = Image2Pipeline(save_results=True, config_file=os.path.join(configdir, "calwebb_image2.cfg"), output_dir=outdir)
    image2 = mod_pipe_steps(image2, params.get('IMAGE2', []))
    image2.run(image1_file)
    image2_file = direct.replace('_uncal.fits', '_cal.fits')

    # Run the direct image through IMAGE3
    printnice("IMAGE3")
    image3 = Image3Pipeline(save_results=True, config_file=os.path.join(configdir, "calwebb_image3.cfg"), output_dir=outdir)
    image3 = mod_pipe_steps(image3, params.get('IMAGE3', []))
    image3.run(image2_file)


def process_wfss(direct, dispersed, configdir='.', outdir=None, params={}):
    """
    Calibrate a combination of direct and dispersed WFSS images

    Parameters
    ----------
    direct: str
        Path to direct image
    dispersed: str
        Path to dispersed image
    configdir: str
        Path to directory with config files
    outdir: str
        Path to directory where output files should be written
    params: dict
        The key/value pairs of pipeline/step to override in the calibration
    """
    # Make sure cfg file are there
    collect_pipeline_cfgs.collect_pipeline_cfgs(configdir)

    if outdir is None:
        outdir = os.path.dirname(direct)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("Created directory {}".format(outdir))

    # Ensure that the source_catalog step has save_results set to True
    sc = ConfigObj(os.path.join(configdir, "source_catalog.cfg"))
    if "save_results" not in sc:
        print("Modifying {} to add save_results option".format(os.path.join(configdir, "source_catalog.cfg")))
        sc["save_results"] = True
        sc.write()

    # Process the direct image
    process_image(direct, configdir=configdir, outdir=outdir, params=params)

    # Run the dispersed image through DETECTOR1
    printnice("DETECTOR1")
    spec1 = Detector1Pipeline(save_results=True, config_file=os.path.join(configdir, "calwebb_detector1.cfg"), output_dir=outdir)
    spec1 = mod_pipe_steps(spec1, params.get('DETECTOR1', []))
    spec1.run(dispersed)
    spec1_file = dispersed.replace('_uncal.fits', '_rate.fits')

    # Add the source catalog file to the spectral rate file
    scatfile = glob.glob(os.path.join(outdir, "*.ecsv"))[0]
    with fits.open(spec1_file, mode="update") as hdulist:
        hdr0 = hdulist[0].header
        hdr0.set("SCATFILE", scatfile)

    # Run the dispersed image through SPEC2
    printnice("SPEC2")
    spec2 = Spec2Pipeline(save_results=True, config_file=os.path.join(configdir, "calwebb_spec2.cfg"), output_dir=outdir)
    spec2 = mod_pipe_steps(spec2, params.get('SPEC2', []))
    spec2.run(spec1_file)


def mod_pipe_steps(pipe, params):
    """
    Update pipeline object attributes

    Parameters
    ----------
    pipe:
        The pipeline object
    params: str
        The list of (step.param = value) of the parameter to modify
    """
    for param in params:

        try:

            # Split the parameter
            step, val = param.split('=')
            step, par = step.split('.')
            step, par, val = step.strip(), par.strip(), val.strip()

            # Check if it's a number
            if isnumber(val):
                if val.isnumeric():
                    val = int(val)
                else:
                    val = float(val)

            # Update the attribute
            setattr(getattr(pipe, step), par, val)

        except IOError:

            print('{} not applied'.format(param))

    return pipe
