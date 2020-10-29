========
`locals`
========

.. image:: https://img.shields.io/github/v/release/hover2pi/locals?label=locals
        :alt: GitHub release
        :target: https://github.com/hover2pi/locals/releases

.. image:: https://img.shields.io/travis/hover2pi/locals.svg
        :target: https://travis-ci.org/hover2pi/locals.svg?branch=master

.. image:: https://readthedocs.org/projects/locals/badge/?version=latest
        :target: https://locals.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/hover2pi/locals/badge.svg
        :target: https://coveralls.io/github/hover2pi/locals


`locals` is a pure Python package that ingests JWST Wide-Field Slitless Spectroscopy data and returns a source catalog of all the low-mass stars in the field along with their calculated fundamental and secondary parameters.


Installation
############

.. code-block:: bash

   pip install locals

Usage
#####

For each spectrum in our WFSS data products, a `Source` object is created. `locals` then searches Vizier to find ancillary photometry, parallaxes, and spectral types. Then it assembles an SED from the available data and calculates fundamental and atmospheric parameters.

.. code-block:: python

   from locals import SourceCatalog
   from pkg_resources import resource_filename

   # Path to the JWST pipeline output directory
   cat_path = resource_filename('locals', 'data/wfss/')

   # Initialize a source catalog from the directory of pipeline products
   cat = SourceCatalog(cat_path)


Let's see the source list of all objects in this catalog.

.. code-block:: python

   cat.source_list

All the sources are stored as `Source` instances within the `SourceCatalog.sources` attribute.

Let's inspect an individual source and plot its SED.

.. code-block:: python

   # Get the source instance
   targ = cat.sources[0]

   # Calculate the results
   print(targ.results)

.. code-block:: html

    <i>Table length=15</i>
    <table id="table4767599808" class="table-striped table-bordered table-condensed">
    <thead><tr><th>param</th><th>value</th><th>unc</th><th>units</th></tr></thead>
    <thead><tr><th>object</th><th>object</th><th>object</th><th>object</th></tr></thead>
    <tr><td>name</td><td>LHS 2924</td><td>--</td><td>--</td></tr>
    <tr><td>age</td><td>6.0</td><td>4.0</td><td>Gyr</td></tr>
    <tr><td>distance</td><td>10.99</td><td>0.02</td><td>pc</td></tr>
    <tr><td>parallax</td><td>90.9962</td><td>0.12710000574588776</td><td>mas</td></tr>
    <tr><td>radius</td><td>0.8786209573091851</td><td>0.06782337214316517</td><td>solRad</td></tr>
    <tr><td>spectral_type</td><td>--</td><td>--</td><td>--</td></tr>
    <tr><td>membership</td><td>--</td><td>--</td><td>--</td></tr>
    <tr><td>fbol</td><td>6.95e-11</td><td>7.88e-13</td><td>erg / (cm2 s)</td></tr>
    <tr><td>mbol</td><td>13.913</td><td>0.012</td><td>--</td></tr>
    <tr><td>Lbol</td><td>1e+30</td><td>1.2e+28</td><td>erg / s</td></tr>
    <tr><td>Lbol_sun</td><td>-3.58</td><td>0.005</td><td>--</td></tr>
    <tr><td>Mbol</td><td>13.708</td><td>0.012</td><td>--</td></tr>
    <tr><td>logg</td><td>4.5</td><td>0.07</td><td>--</td></tr>
    <tr><td>mass</td><td>0.8896913720017506</td><td>0.0</td><td>solMass</td></tr>
    <tr><td>Teff</td><td>783</td><td>30</td><td>K</td></tr>
    </table>

Hooray!

Let's see a table of all the results.

.. code-block:: python

   # Get the source instance
   targ = cat.sources[0]

.. code-block:: html

    <i>Table length=2</i>
    <table id="table4629905304" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>age [Gyr]</th><th>distance [pc]</th><th>parallax [mas]</th><th>radius [solRad]</th><th>spectral_type [--]</th><th>membership [--]</th><th>fbol [erg / (cm2 s)]</th><th>mbol [--]</th><th>Lbol [erg / s]</th><th>Lbol_sun [--]</th><th>Mbol [--]</th><th>logg [--]</th><th>mass [solMass]</th><th>Teff [K]</th></tr></thead>
    <thead><tr><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th></tr></thead>
    <tr><td>Source 10</td><td>6.0 +/- 4.0</td><td>79.19 +/- 0.68</td><td>12.6276 +/- 0.10779999941587448</td><td>nan +/- nan</td><td>--</td><td>--</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>-9.22e+18 +/- -9.22e+18</td></tr>
    <tr><td>Source 13</td><td>6.0 +/- 4.0</td><td>10.99 +/- 0.02</td><td>90.9962 +/- 0.12710000574588776</td><td>0.8786209573091851 +/- 0.06782337214316517</td><td>--</td><td>--</td><td>6.1e-11 +/- 3.41e-13</td><td>14.055 +/- 0.006</td><td>8.82e+29 +/- 5.89e+27</td><td>-3.64 +/- 0.003</td><td>13.85 +/- 0.006</td><td>4.5 +/- 0.07</td><td>0.8896913720017506 +/- 0.0</td><td>758 +/- 29</td></tr>
    </table>

Fantastic!
