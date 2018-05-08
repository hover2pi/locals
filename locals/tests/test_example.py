def test_primes():
    from ..example_mod import primes
    assert primes(10) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]


def testSourceCatalog():
    """
    Test to see if SourceCatalog can be assembled
    """
    # Cat path
    cat_path = '/Users/jfilippazzo/Desktop/LOCALS_data/sources/'
    
    # Initialize empty source catalog
    cat = SourceCatalog(cat_path)
    
    return cat
    
def testSource():
    """
    Test to see that Vega source can be assembled
    """
    ra = 279.23473479*q.deg
    dec = 38.78368896*q.deg
    
    # Make source
    source = src.Source(ra=ra, dec=dec, name='Vega')
    
    # Find photometry
    source.find_photometry()
    
    return source
    
def testSED():
    
    # Run test catalog and grab Vega
    cat = testSourceCatalog()
    source = cat.sources[1]
    
    # Pass data to MakeSED
    db = astrodb.Database('/Users/jfilippazzo/Documents/Modules/BDNYCdevdb/bdnycdev.db')
    s = sed.MakeSED(source.id, db, from_dict=source)
    s.plot()