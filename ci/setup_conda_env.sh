#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n lcl python=$PYTHON_VERSION || exit 1
source activate lcl

echo "Installing packages..."
conda install flake8 beautifulsoup4 lxml numpy astropy h5py
pip install astroquery sedkit svo_filters pytest pytest-cov coveralls