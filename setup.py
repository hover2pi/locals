#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup, find_packages
    setup
except ImportError:
    from distutils.core import setup
    setup

from codecs import open
from os import path

setup(
    name='locals',
    version='0.2.0',
    description='Low-mass Object Characterization by AnaLyzing Slitless Spectroscopy',
    url='https://github.com/hover2pi/locals',
    author='Joe Filippazzo',
    author_email='jfilippazzo@stsci.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords='astrophysics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy', 'astropy', 'bokeh', 'sedkit', 'svo_filters', 'h5py', 'astroquery'],
    include_package_data=True,

)