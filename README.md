
# How to use `locals`
`locals` is a pure Python package that ingests JWST Wide-Field Slitless Spectroscopy data and returns a source catalog of all the low-mass stars in the field along with their calculated fundamental and secondary parameters.

### Requirements
- pip install SEDkit
- pip install bokeh


```python
# Imports
from locals import source, catalog
from SEDkit import sed, spectrum, synphot
import astropy.units as q
import astropy.table as at
import numpy as np
from pkg_resources import resource_filename
from bokeh.io import output_notebook, show
output_notebook()
```


## Making an SED for a single source (no WFSS spectra)

For each spectrum in our WFSS data products, a `Source` object is created. `locals` then searches Vizier to find ancillary photometry, parallaxes, and spectral types. Then it assembles an SED from the available data and calculates fundamental and atmospheric parameters.

Let's pretend one of the brown dwarfs in our field is [LHS2924](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=LHS2924&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id). Here is what `locals` is doing for each source under the hood.


```python
# Coordinates of LHS 2924 (Usually taken from the WFSS source_list)
ra = 217.180137*q.deg
dec = 33.177539*q.deg

# Make the source object
src = source.Source(ra=ra, dec=dec, name='LHS 2924')

# Find some photometry
src.find_photometry()

# Find a distance from Gaia DR2
src.find_parallax()

# Take a look
print('\n',src.photometry,'\n')
```

    Setting age to (<Quantity 6.0 Gyr>, <Quantity 4.0 Gyr>)
    Make this handle asymmetric uncertainties!
    
       band          eff         app_magnitude ...      abs_flux_unc      bandpass
                     um                       ... erg / (Angstrom cm2 s)         
    -------- ------------------ ------------- ... ---------------------- --------
      Gaia.G 0.6597013110557894       16.6299 ...                    0.0   Gaia.G
     2MASS.J 1.2393093155660664 11.9899997711 ... 2.6782632688510406e-17  2MASS.J
     2MASS.H 1.6494947091246095 11.2250003815 ... 2.3139101841398836e-17  2MASS.H
    2MASS.Ks  2.163860596453316 10.7440004349 ... 1.3684444634813076e-17 2MASS.Ks
     WISE.W1 3.3897048577485647 10.4280004501 ...  3.543992004006416e-18  WISE.W1
     WISE.W2    4.6406356615375 10.1619997025 ... 1.3561312900080802e-18  WISE.W2
     WISE.W3 12.567593607521216 9.64599990845 ...  7.026240830733788e-20  WISE.W3
     WISE.W4 22.314225592308713 9.43700027466 ...   7.62677792182739e-20  WISE.W4 
    


Now we can calculate the fundamental parameters and plot the results.


```python
src.results
```

    Setting radius to (<Quantity 0.8786209573091851 solRad>, <Quantity 0.06782337214316517 solRad>)





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


## Make SEDs for all sources in the WFSS field
Now that we see the output of a single source, we can run this on the entire `SourceCatalog` generated from the WFSS output. The difference here is that a WFSS spectrum will be included in each `Source`.


```python
# Path to the JWST pipeline output directory
cat_path = resource_filename('locals', 'data/wfss/')

# Initialize empty source catalog
cat = catalog.SourceCatalog(cat_path)
```

    Setting age to (<Quantity 6.0 Gyr>, <Quantity 4.0 Gyr>)
    Make this handle asymmetric uncertainties!
    Setting age to (<Quantity 6.0 Gyr>, <Quantity 4.0 Gyr>)
    Make this handle asymmetric uncertainties!
    Warning, 1 of 546 bins contained negative fluxes; they have been set to zero.


Let's see the source list of all objects in this catalog.


```python
cat.source_list
```




<i>Table length=2</i>
<table id="table4704793376" class="table-striped table-bordered table-condensed">
<thead><tr><th>id</th><th>xcentroid</th><th>ycentroid</th><th>ra_icrs_centroid</th><th>dec_icrs_centroid</th><th>xmin</th><th>xmax</th><th>ymin</th><th>ymax</th><th>abmag</th><th>abmag_error</th><th>sky_bbox_ll</th><th>sky_bbox_ur</th><th>sky_bbox_lr</th><th>sky_bbox_ul</th><th>icrs_centroid</th></tr></thead>
<thead><tr><th></th><th>pix</th><th>pix</th><th></th><th></th><th>pix</th><th>pix</th><th>pix</th><th>pix</th><th></th><th></th><th>deg,deg</th><th>deg,deg</th><th>deg,deg</th><th>deg,deg</th><th>deg,deg</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>object</th><th>object</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>object</th><th>object</th><th>object</th><th>object</th><th>object</th></tr></thead>
<tr><td>10</td><td>1023.89938168</td><td>1023.19367348</td><td>-3.4283850489866588e-06</td><td>-2.694720849093827e-06</td><td>940.0</td><td>1107.0</td><td>943.0</td><td>1106.0</td><td>9</td><td>0.0816191841269</td><td>0.000729834375136,-0.000694686445181</td><td>359.999270901,0.000711962615737</td><td>0.000729834375136,-0.000694686445181</td><td>359.999270901,0.000711962615737</td><td>89.7455,-29.05744</td></tr>
<tr><td>13</td><td>1023.89938168</td><td>1023.19367348</td><td>-3.4283850489866588e-06</td><td>-2.694720849093827e-06</td><td>940.0</td><td>1107.0</td><td>943.0</td><td>1106.0</td><td>9</td><td>0.0816191841269</td><td>0.000729834375136,-0.000694686445181</td><td>359.999270901,0.000711962615737</td><td>0.000729834375136,-0.000694686445181</td><td>359.999270901,0.000711962615737</td><td>217.180137,33.177539</td></tr>
</table>



All the sources are stored as `Source` instances within the `SourceCatalog.sources` attribute. 

For example, LHS 2924 is part of this catalog so let's make sure we get the same results as when it is run individually.


```python
# Get the source instance
lhs2924 = cat.sources[1]

# Calculate the results
lhs2924.results

```
    Setting radius to (<Quantity 0.8786209573091851 solRad>, <Quantity 0.06782337214316517 solRad>)
```

Hooray!

Let's take a look at the other one, 2MASS 0558-2903, an extreme subdwarf M7.


```python
# Get the source instance
subdwarf = cat.sources[0]

# Calculate the results
print(subdwarf.results)
```

    Setting radius to (<Quantity nan solRad>, <Quantity nan solRad>)
        param       value           unc             units    
    ------------- --------- ------------------- -------------
             name Source 10                  --            --
              age       6.0                 4.0           Gyr
         distance     79.19                0.68            pc
         parallax   12.6276 0.10779999941587448           mas
           radius       nan                 nan        solRad
    spectral_type        --                  --            --
       membership        --                  --            --
             fbol       nan                 nan erg / (cm2 s)
             mbol       nan                 nan            --
             Lbol       nan                 nan       erg / s
         Lbol_sun       nan                 nan            --
             Mbol       nan                 nan            --
             logg       nan                 nan            --
             mass       nan                 nan       solMass
             Teff -9.22e+18           -9.22e+18             K


Fantastic!

Let's see a table of the results.


```python
cat.results
```




<i>Table length=2</i>
<table id="table4629905304" class="table-striped table-bordered table-condensed">
<thead><tr><th>name</th><th>age [Gyr]</th><th>distance [pc]</th><th>parallax [mas]</th><th>radius [solRad]</th><th>spectral_type [--]</th><th>membership [--]</th><th>fbol [erg / (cm2 s)]</th><th>mbol [--]</th><th>Lbol [erg / s]</th><th>Lbol_sun [--]</th><th>Mbol [--]</th><th>logg [--]</th><th>mass [solMass]</th><th>Teff [K]</th></tr></thead>
<thead><tr><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th><th>str42</th></tr></thead>
<tr><td>Source 10</td><td>6.0 +/- 4.0</td><td>79.19 +/- 0.68</td><td>12.6276 +/- 0.10779999941587448</td><td>nan +/- nan</td><td>--</td><td>--</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>nan +/- nan</td><td>-9.22e+18 +/- -9.22e+18</td></tr>
<tr><td>Source 13</td><td>6.0 +/- 4.0</td><td>10.99 +/- 0.02</td><td>90.9962 +/- 0.12710000574588776</td><td>0.8786209573091851 +/- 0.06782337214316517</td><td>--</td><td>--</td><td>6.1e-11 +/- 3.41e-13</td><td>14.055 +/- 0.006</td><td>8.82e+29 +/- 5.89e+27</td><td>-3.64 +/- 0.003</td><td>13.85 +/- 0.006</td><td>4.5 +/- 0.07</td><td>0.8896913720017506 +/- 0.0</td><td>758 +/- 29</td></tr>
</table>




