COLOR_CUTS = {'brown_dwarfs': {'NIRISS.F150W-NIRISS.F200W':(0.3,0.5), 'NIRISS.F158M-NIRISS.F150W':(-0.2,-0.1)}}

def in_color_range(photometry, color_cut):
    """Test to see if the source should be excluded given the color cut
    
    Parameters
    ----------
    photometry: astropy.tables.Table
        The table of photometry to check
    color_cut: str, dict
        The name of the built-in color cut to apply or
        a dictionary of the 
    
    Returns
    -------
    bool
        The 
    """
    # Convert name into dict if a valid string is provided
    if isinstance(color_cut, str):
        name = color_cut
        color_cut = COLOR_CUTS.get(color_cut)
        if not isinstance(color_cut, dict):
            print("Could not find color cut named", name)

    # Return True if null color cut or bad name
    if color_cut is None:
        return True

    # Apply color cut if it is a dictionary
    elif isinstance(color_cut, dict):

        # Apply all criteria
        keep = True
        for bands, bounds in color_cut.items():
            bounds = sorted(bounds)
            b1, b2 = bands.split('-')
            if b1 in photometry['band'] and b2 in photometry['band']:

                # Calculate color
                mag1 = photometry[photometry['band']==b1][0]['app_magnitude']
                mag2 = photometry[photometry['band']==b2][0]['app_magnitude']

                # Apply criterion
                if not bounds[0] <= mag1-mag2 <= bounds[1]:
                    keep = False

        return keep

    else:
        raise ValueError("I do not understand the color_cut input {}.\
                          Please enter a string or dictionary."\
                          .format(color_cut))