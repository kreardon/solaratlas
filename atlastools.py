###############
# ATLAS TOOLS #
###############

''' Intended for the use alongside Atlas object instantiating & use.'''

# 
import numpy as np
from astropy.io import fits
from specutils import Spectrum1D
from astropy import units as u
from astropy.coordinates import earth
import re

class observatory:
  def __init__(observatoryobj, obsname, obscoord, instrument):
    observatoryobj.name       = obsname
    #observatoryobj.location   = earth.EarthLocation.of_site(obsname)
    observatoryobj.location   = obscoord
    observatoryobj.instrument = instrument

# Issue #3 - very rudimentary atlas object - needs to be fleshed out 
class atlas:
  def __init__(atlasobject, specobj_sun, specobj_telluric, specobj_all, target, source, observatory):
    atlasobject.target      = target
    atlasobject.source      = source
    atlasobject.observatory = observatory
    # atlas data here
    atlasobject.sun         = specobj_sun           # when populated, should be a spec1D object
    atlasobject.atm         = specobj_telluric      # when populated, should be a spec1D object
    atlasobject.components  = specobj_all           # when populated, should be a dictionary of spec1D object
        
#
def make_dictionary(filename, extension):
    '''
    Generates a dictionary based on file specified 
    extension keywords and values.

    PARAMETERS:
        filename: path of file to generate
                  dictionary from (must be a string)

        extension: specific extension of file to 
                   generate dictionary from (integer)

    RETURNS:
        dictionary: the dictionary version of file
                    extension's keywords and values.
    '''
    
    f = fits.open(filename)
    
    # creating Primary Header dictionary
    dictionary = {}
    
    for kw in f[0].header:
            dictionary[kw] = f[0].header[kw]
            
    # adding in all extension-specified dictionary values
    for kw in f[extension].header:
        dictionary[kw] = f[extension].header[kw]

    units        = search_key('TUNIT', dictionary)
    type_columns = search_key('TTYPE', dictionary)
    print(type_columns)
     
    return dictionary

#
def search_key(kw, dictionary):
    '''
    Searches through a dictionary (generated
    for file via function 'make_dictionary')
    for specified keyword and returns a list
    with all values matching that keyword.

    PARAMETERS:
        kw: keyword the user wants to search the
            dictionary for (must be a string)

        dictionary: a python dictionary containing
                    all keywords and values for 
                    specified atlas file (generated
                    via function 'make_dictionary')
    RETURNS:
        res: list containing the value results of the search
    '''

    search_key = kw
    res = [val for key, val in dictionary.items() if search_key in key]

    return res

#
def split_column_ttype(ttype_value):

    ttype_info     = {}
    matchidx       = re.search(r'[A-Za-z]',ttype_value[::-1])
    coltype_base   = (ttype_value[:-matchidx.start(0)]).strip()
    coltype_idx    = (ttype_value[-matchidx.start(0):]).strip()
    ttype_info['Column Type']        = coltype_base
    ttype_info['Column Type Number'] = int(coltype_idx)
    print(coltype_base,"  --  ",coltype_idx)

    return ttype_info

#
def column_information(dictionary):
    
    columns = search_key('TFORM', dictionary)
    
    column_info = []
    
    for i in range(len(columns)):
        comptype = search_key('TTYPE', dictionary)
        column_info.append(comptype)
        
        units = search_key('TUNIT', dictionary)
        column_info.append(units)

        axislabels = search_key('TDESC', dictionary)
        column_info.append(axislabels)

        target = search_key('TOBJC', dictionary)
        column_info.append(target)

        deriv_methd = search_key('TMTHD', dictionary)
        column_info.append(deriv_methd)

        tel = search_key('TWATM', dictionary)
        column_info.append(tel)
    
    for a, b, c, d, e, f in zip(column_info[0],column_info[1],column_info[2],column_info[3],column_info[4], column_info[5]):
        print("-----COLUMN-----\n\nComponent type: {}\nColumn data units: {}\nAxis Labels: {}\nTarget: {}\nDerivation Method: {}\nIncludes Telluric Absorption: {}\n".format(a,b,c,d,e,f))

#
def find_column_index(dictionary, column_key):
    # extract the FITS keyword column index for the column being plotted
    # (as defined by the column_key input)
    # i.e. the keywords have names like TTITLn - we need to find the value 
    # of "n" corresponding to the column name provided by column_key

    kw = list(dictionary.keys())[list(dictionary.values()).index(column_key)]
    #integer_search = re.compile(r'\d+(?:\.\d+)?')
    integer_search = re.compile(r'[A-Z]+(\d+)')
    col_index   = integer_search.findall(kw)
    
    return col_index
#
def telluric_info(dictionary):
    '''
    Return series of print statements notifying whether 
    telluric absoprtion is included in stored data.

    PARAMETERS:
        dictionary: a python dictionary containing
                    all keywords and values for 
                    specified atlas file (generated
                    via function 'make_dictionary')

    RETURNS:
        *NO VALUE RETURNED* - prints information
    '''
    types = []
    present = []
    
    for value in search_key('TTYPE', dictionary):
        types.append(value)
    for value in search_key('TWATM', dictionary):
        present.append(value)
        
    for i in present:
        str(i)
    for t, p in zip(types, present):
        print("Data type: {} - Includes telluric absorption?: {}\n".format(t, p))

#
def filecontent_map(filename):
    '''
    Quick file content map of general file information
    & extensions.
    
    PARAMETERS:
        filename: path of file to provide information on
                  (must be a string)
    RETURNS:
        **NO RETURN VALUE, just printed information**
    '''
    
    file = fits.open(filename, memmap=True)
    
    
    print("General file information:\n")
    
    print("Object: {}".format(file[0].header['OBJECT']))
    print("Atlas source: {}".format(file[0].header['ATL_SOUR']))
    print("Atlas acquisition site: {}".format(file[0].header['ATL_OBS']))
    print("Atlas wavelength coverage: {:.2f} {} – {:.2f} {}\n".format(file[0].header['WAVEMIN'], file[0].header['CUNIT1'], file[0].header['WAVEMAX'], file[0].header['CUNIT1']))

# skips the primary HDU because not a binary table of data and info
    for i in range(1,len(file)):
        print("EXTENSION {}:\n".format(i))
        
        atlasdict = make_dictionary(filename,i)
        
        columns = search_key('COLNUM', atlasdict)
        count = 1
            
        minimums = search_key('TDMIN', atlasdict)
        maximums = search_key('TDMAX', atlasdict)

        waverange = [x for y in zip(minimums, maximums) for x in y]
        res = list(zip(waverange, waverange[1:] + waverange[:1]))
        finalranges = res[::2]

        units = search_key('TUNIT', atlasdict)
        column_type = search_key('TTYPE', atlasdict)
        #matchidx = re.search(r'[A-Za-z]',column_type[::-1])
        #coltype_base = (column_type[:-matchidx.start(0)]).strip()
        #coltype_idx  = (column_type[-matchidx.start(0):]).strip()

        newcount = 0
        for waveinterval in finalranges:
            column_type_in = column_type[newcount]
            matchidx       = re.search(r'[A-Za-z]',column_type_in[::-1])
            coltype_base   = (column_type_in[:-matchidx.start(0)]).strip()
            coltype_idx    = (column_type_in[-matchidx.start(0):]).strip()

            print("Column {} - {} : value range = {:.2f} – {:.2f} {}".format(newcount, coltype_base, waveinterval[0], waveinterval[1], units[newcount]))
            newcount+=1
        print("\n\n")
    count += 1
            
#
def store_data(filename, extension, startwave=1*u.nm, endwave=1*u.nm):
    '''
    Function that stores all data within specified extension
    and returns two dictionaries: one of the data from the file
    with units and one of that same data that has been pushed
    through Spectrum1D
    
    PARAMETERS:
        filename: full name of the file (must be a string)
        extension: specific extension to pull and store
                   data from
    
    RETURNS:
        file_data: data with proper astropy units straight 
                   from the file associated with the data 
                   type- essentially the metadata straight 
                   from file(dictionary with {data type 
                   label : value})
        spec_data: metadata that has been fed through Spectrum1D
                   to generate usable atlas data (dictionary
                   with data type label : Spec1D object)
    '''
    
    import numpy as np
    from astropy.io import fits
    import atlastools
    from astropy import units as u
    from specutils import Spectrum1D
    
    # opening the .fits file for access
    f = fits.open(filename)
    
    # generating an extension-specific dictionary
    atlasdict = atlastools.make_dictionary(filename, extension)
    
    # searching the extension dictionary for keywords tagged TTYPE which are the data
    # and storing them in a list
    column_types = atlastools.search_key('TTYPE', atlasdict)
    print("checking ",column_types)

    # now we need to find out which column of the table contains the wavelength
    # scale to apply to the data
    # so we will first search for all columns where the TTYPE keyword contains the
    # 'Wavelength Scale" identifier
    #
    # since there could technically be multiple columns with wavelength data,
    # we will search for the column with the lowest unique identifier index.
    # normally this should be "1", but you never know what people will do...
    wavelength_col_idx = 9999
    wavescale_identifier = 'Wavelength Scale'
    # loop through all the TTYPEn values and find those with a wavelength scale
    for ttype_value in column_types:
        column_type = split_column_ttype(ttype_value)
        if column_type['Column Type'] == wavescale_identifier:
            # now check to see if we have found the column with the lowest 
            # wavelength scale identifier
            if column_type['Column Type Number'] < wavelength_col_idx:
                wavelength_col_idx = column_type['Column Type Number']
                wavelength_primary_id = ttype_value
    wavelength_scale_ttype = wavelength_primary_id

    # this has the wavelength units hardwired, but it really should be determined from
    # the units given in TUNIT for this column
    wavelength_unit  = u.nm
    wavelength_scale = f[1].data['Wavelength Scale   1'] * wavelength_unit

    atlas_min = np.min(wavelength_scale)
    atlas_max = np.max(wavelength_scale)

    if ((startwave < atlas_min) or 
        (startwave > atlas_max)) : 
         startwave = atlas_min
    if ((endwave <= startwave)  or 
        (endwave >= atlas_max))  : 
         endwave = atlas_max

    mask = ((wavelength_scale > startwave) & (wavelength_scale < endwave))
    print(np.sum(mask))
    
    # reading in the actual data per TTYPE keyword, and storing that data in a new
    # dictionary associated with the type of data it is
    stored = {}
    for value in column_types:
        full_column   = f[1].data[value]
        selected_range = full_column[mask]
        stored[value] = selected_range
        
    # searching the extension dictionary for keywords tagged TUNIT which are the units
    # associated with the data and storing them in a list
    units_to_store = atlastools.search_key('TUNIT', atlasdict)
    
    # associating the data type with the units it required and storing in a new dictionary
    unit_tags = {}
    for key, unit in zip(column_types, units_to_store):
        unit_tags[key] = unit
        
    # searching through the unit_tags dictionary and changing the value to the proper
    # astropy unit quantity (required for reading into Spectrum1D)
    for key, value in unit_tags.items():
        if value == 'Relint':
            unit_tags[key] = u.dimensionless_unscaled
        if value == 'Nanometers':
            unit_tags[key] = u.nm
        if unit_tags[key] == 'cm^-1':
            unit_tags[key] = u.k     # kayser: CGS unit of wavenumber
        if unit_tags[key] == 'W m^(-2) sr^(-1) Angstrom^(-1)':
            unit_tags[key] = u.W / u.m / u.m / u.sr / u.AA     # kayser: CGS unit of wavenumber
    
    # generating a new dictionary that has the data type keyword now associated with the
    # actual data with the astropy units attached (this is the metadata dictionary returned)
    file_data = {}
    for key in column_types:
        data = stored[key]
        unit = unit_tags[key]
        combined = data * unit
        file_data[key] = combined
    
    # Finally, pushing the file metadata through Spectrum1D to create a dictionary of Spec1D objects
    # associated with the type of data they are (this is the Spec1D dictionary returned)
    # ( I intend to come back and fix this to use a regex string pattern of 'Wavelength Scale' rather
    # than nested if statements)
    spec_data = {}
    for key in file_data:
        if re.search(r'Wavelength Scale*', key) is None:
#        if key != 'Wavelength Scale   1':
#            if key != 'Wavelength Scale   2':
                col_index = find_column_index(atlasdict, key)
                finalized = Spectrum1D(spectral_axis=file_data['Wavelength Scale   1'], flux=file_data[key])
                finalized.meta.update(unit            = atlastools.search_key('TUNIT' + col_index[0], atlasdict))
                finalized.meta.update(has_telluric    = atlastools.search_key('TWATM' + col_index[0], atlasdict))
                finalized.meta.update(has_solar       = atlastools.search_key('TWSUN' + col_index[0], atlasdict))
                finalized.meta.update(observed_object = atlastools.search_key('TOBJC' + col_index[0], atlasdict))
                finalized.meta.update(observed_target = atlastools.search_key('TTRGT' + col_index[0], atlasdict))
                finalized.meta.update(col_title       = atlastools.search_key('TTITL' + col_index[0], atlasdict))
                finalized.meta.update(col_label       = atlastools.search_key('TLABL' + col_index[0], atlasdict))
                finalized.meta.update(col_description = atlastools.search_key('TDESC' + col_index[0], atlasdict))
                spec_data[key] = finalized
                
    return file_data, spec_data

#
def make_atlas(filename, extension, startwave=1*u.nm, endwave=1*u.nm, loaddata=0):
    '''
    Generates an atlas object based on file and extension specified 
    using extension keywords and values.

    PARAMETERS:
        filename: path of file to generate
                  dictionary from (must be a string)

        extension: specific extension of file to 
                   generate dictionary from (integer)

    RETURNS:
        dictionary: the dictionary version of file
                    extension's keywords and values.
    '''
    from astropy.coordinates import EarthLocation
    import atlastools

    atlasdict = make_dictionary(filename, extension)

    obslocation = EarthLocation(lat=atlastools.search_key('ATL_LAT', atlasdict), 
                                lon=atlastools.search_key('ATL_LONG', atlasdict), 
                                height=atlastools.search_key('ATL_ALT', atlasdict))
    
    atlas_obs = observatory(atlastools.search_key('ATL_OBS', atlasdict), obslocation, 'FTS')
    #atlas_obs = observatory('Kitt Peak', EarthLocation.of_site('Kitt Peak'), 'FTS')
    atlas = atlastools.atlas(1, 1, 1, atlastools.search_key('OBJECT', atlasdict), atlastools.search_key('ATL_SOUR', atlasdict), atlas_obs)
    
    # Issue #1 - hardcoded extraction of column information into spec1D objects
    if loaddata is not 0:
    
        (file_data, spec_data) = store_data(filename, extension, startwave=startwave, endwave=endwave)
        search_key = 'Local Intensity   1'
        spect_sun  = spec_data.get(search_key)
        search_key = 'Telluric Spectrum   1'
        spect_atm  = spec_data.get(search_key)

        atlas.sun = spect_sun
        atlas.atm = spect_atm
        atlas.components = spec_data
    
        # Issue #2 - need more column specific information associated with Spec1D object
        # is putting a range of individual items in the meta dictionary the best way to do this?
        atlas.sun.meta.update(title=atlastools.search_key('TTITL2', atlasdict))
        atlas.atm.meta.update(title=atlastools.search_key('TTITL4', atlasdict))
        
    return atlas

#
def atlas_spectrum_plot(column_key, atlas_obj, startwave=1*u.nm, endwave=1*u.nm, plot_unit='AA'):
    '''
    Function that returns a single plot of the 
    provided data given  associated dictionary.
    
    PARAMETERS:
    
        column_key:    the keyword value used to index the Spectrum1D
                    object within the dictionary              (TYPE: string)

        atlas_obj:  atlas object containing spec1D object dictionary with
                    atlas components to be plotted (generated from the atlastools
                    'make_atlas' function, indexed as
                    atlas_obj.components['column_key'], where the column_key
                    is a string key value)                    (TYPE: class object)

        dictionary : dictionary associated with extension data being plotted
                     (generated from 'make_dictionary')       (TYPE: dictionary)

        startwave:  starting wavelength of plot range         (TYPE: integer * astropy.units)

        endwave:    ending wavelength of plot range           (TYPE: integer * astropy.units)
    
        plot_unit:  wavelength units in which to make the plot
                    typical values are ['AA', 'nm', 'micron'] (TYPE: string)

    RETURNS:
    
        **NO VALUE(S) RETURNED:** generates plot
    '''
    
    import matplotlib.pyplot as plt
    from astropy import units as u
    import re
    import time
    
    plt.figure(figsize=(11,8))
    plt.rc('font', family='serif',size=14)

    # convert the spectral axis and the spectral plot ranges
    # to the requested plotting units
    component_to_plot = atlas_obj.components[column_key]
    spectral_axis     = component_to_plot.spectral_axis
    spectral_data     = component_to_plot.flux
    
    plot_spectral_axis = spectral_axis.to(plot_unit)
    startwave = startwave.to(plot_unit)
    endwave   = endwave.to(plot_unit)

    plt.plot(plot_spectral_axis, spectral_data)
    
    # extract the FITS keyword column index for the column being plotted
    # (as by the column_key input)
    # i.e. the keywords have names like TTITLn - we need to find the value 
    # of "n" corresponding to the column name provided by column_key
    #kw = list(dictionary.keys())[list(dictionary.values()).index(column_key)]
    #integer_search = re.compile(r'\d+(?:\.\d+)?')
    #col_index   = integer_search.findall(kw)
    
    # define the keyword names that contain information for the plotted spectrum
    # old way of doing this, pulling keyword values from the header dictionary
    #col_index = find_column_index(dictionary, column_key)
    #title   = dictionary['TTITL' + col_index[0]]
    #ylabel  = dictionary['TUNIT' + col_index[0]]
    #descrip = dictionary['TDESC' + col_index[0]]
    # new way using values pulled from spec1D object metadata 
    title    = (component_to_plot.meta['col_title'])[0]
    ylabel   = (component_to_plot.meta['unit'])[0]
    descrip  = (component_to_plot.meta['col_description'])[0]
    
    # explicitly searching for the minimum and maximum in the 
    # full, large atlas array can be very slow,
    # but since the wavelength scale _should_ be arrayed in 
    # montonically increasing order, we can just select
    # the first and last element instead.
    full_extrema_search = 0
    if full_extrema_search == 1:
        plot_spectral_axis_min = min(plot_spectral_axis)
        plot_spectral_axis_max = max(plot_spectral_axis)
    else:
        plot_spectral_axis_min = plot_spectral_axis[0]
        plot_spectral_axis_max = plot_spectral_axis[-1]

    # check to make sure the starting and ending wavelengths for
    # the plotting actually fall within the range of the atlas
    # to be plotted
    if ((startwave < plot_spectral_axis_min) or 
         (startwave > plot_spectral_axis_max)) : 
             startwave = plot_spectral_axis_min
    if ((endwave <= startwave)                 or 
        (endwave >= plot_spectral_axis_max))    : 
             endwave = plot_spectral_axis_max
    plt.xlim(startwave.value, endwave.value)
    # this should be dynamically set and/or be an optional user input
    plt.ylim(0.01,1.02)
    
    plt.title(title)
    #plt.xlabel(startwave.unit)
    xlabel_str = "Wavelength [ " + "{0.unit:s}".format(startwave) + " ]"
    plt.xlabel(xlabel_str)
    plt.ylabel(ylabel)
    plt.figtext(0.35,-0.05,descrip);
