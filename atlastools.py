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
from operator import itemgetter 
import re

class atlas_file_extension_info:
    def __init__(atlas_file_exten_obj, wavemin, wavemax, wavenum, numcols, coltypes, has_solar, has_telluric):
        atlas_file_obj.wavemin        = wavemin   # 400 * u.nm
        atlas_file_obj.wavemax        = wavemax   # 700 * u.nm
        atlas_file_obj.wavenum        = wavenum   # column length (e.g. 324561) 
        atlas_file_obj.numcols        = numcols   # number of columns in extension
        atlas_file_obj.coltypes       = coltypes  # list ['Wavelength Scale 1', 'Local Intensity 1', 'Local Intensity 2',...]
        atlas_file_obj.has_solar      = has_solar #    list [-1, 1, 1, 0]
        atlas_file_obj.has_telluric   = has_telluric # list [-1, 0, 1, 1]

class atlas_file_content:
    def __init__(atlas_file_obj, filename, filepath, source, observatory, version, obsobject, extension_info):
        atlas_file_obj.filename       = filename
        atlas_file_obj.filepath       = filepath
        atlas_file_obj.source         = source
        atlas_file_obj.observatory    = observatory
        atlas_file_obj.version        = version
        atlas_file_obj.obsobject      = obsobject
        atlas_file_obj.extension_info = extension_info
        
class observatory:
  def __init__(observatoryobj, obsname, obscoord, instrument):
    observatoryobj.name       = obsname
    #observatoryobj.location   = earth.EarthLocation.of_site(obsname)
    observatoryobj.location   = obscoord
    observatoryobj.instrument = instrument

class target:
  def __init__(target, object_name, target_name, solar_mu, solar_rad_distance, integration_area):
    targetobj.name            = object_name
    targetobj.feature_type    = target_name
    targetobj.mu              = solar_mu
    targetobj.raddist         = solar_rad_distance
    targetobj.area            = integration_area
    
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
     
    return dictionary

#
def search_key(keyword_string, dictionary):
    '''
    Searches through a dictionary (generated
    for file via function 'make_dictionary')
    for specified keyword and returns a list
    with all values matching that keyword.

    PARAMETERS:
        keyword_string: keyword the user wants to search the
            dictionary for (must be a string)

        dictionary: a python dictionary containing
                    all keywords and values for 
                    specified atlas file (generated
                    via function 'make_dictionary')
    RETURNS:
        res: (n,2) list containing the matched keyword and 
             the value associated with that keyword
    '''

    # initialize output lists
    keyword_list   = []
    keyword_values = []

    count = 0
    for keyword, value in dictionary.items():
        if keyword_string in keyword:
            keyword_list.append(keyword)
            keyword_values.append(value)
            count += 1
    
    # if no matching keywords were found, return default values
    if count == 0:
        keyword_list   = ['No Keyword Match']
        keyword_values = [np.nan]

    return (keyword_values, keyword_list)
  
#
def filesearch():
    '''
    Directory and file searching function that prints lists 
    of atlas-specific files in user's current directory + 
    subdirectories, directory where atlastools.py resides + 
    subdirectories, and the atlasfiles/ directory + 
    subdirectories if it exists.
    
    PARAMETERS:
    ----------
    
        **NO PARAMETER(S) REQUIRED**
        
    RETURNS:
    -------
    
        **NO RETURN VALUE, just printed information**
    '''

    import os
    from fnmatch import fnmatch

    def path_leaf(path):
        import ntpath
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)

    root = os.getcwd()

    # finding where atlastools.py is
    module = 'atlastools.py'
    for path, subdirs, files in os.walk(root):
        for name in files:
            if fnmatch(name, module):
                atlastools_location = os.path.join(path, module)

    # keeping this in case the module needs to be manually loaded in- 
    # necessary for using 'atlastools.__file__'
    #from importlib.machinery import SourceFileLoader
    #atlastools = SourceFileLoader("atlastools", atlastools_location).load_module()

    # generating path source for atlastools.py location
    module_source = os.path.dirname(atlastools_location)

    # likely types of atlas file extensions to look for
    pattern = "*.fits"
    pattern2 = "*.gz"

    # locating atlas files in the current directory
    print("Current directory: {}\n".format(root))
    print("Available atlas files:\n")

    for path, subdirs, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                print(os.path.join(path_leaf(path), name))
            if fnmatch(name, pattern2):
                print(os.path.join(path_leaf(path), name))
    print("______________________________")

    # locating atlas files in the atlastools.py directory        
    print("\natlastools.py directory: {}\n".format(module_source))
    print("Available atlas files:\n")

    for path, subdirs, files in os.walk(module_source):
        for name in files:
            if fnmatch(name, pattern):
                print(os.path.join(path_leaf(path), name))
            if fnmatch(name, pattern2):
                print(os.path.join(path_leaf(path), name))
    print("______________________________")

    # locating atlas files in the atlasfiles/ directory
    atlasfiles_source = module_source + '/atlasfiles/'

    if os.path.isdir(atlasfiles_source) == True:
        print("\natlasfiles/ directory: {}\n".format(module_source))
        print("Available atlas files:\n")
        for path, subdirs, files in os.walk(atlasfiles_source):
            for name in files:
                if fnmatch(name, pattern):
                    print(os.path.join(path_leaf(path), name))
                if fnmatch(name, pattern2):
                    print(os.path.join(path_leaf(path), name))
        print("______________________________")
        
#
def split_column_ttype(ttype_value):
    '''
    Splits a value provided in the TTYPEn FITS header
    into its base data type (e.g. 'Wavelength Scale') 
    and the unique identifier (an incrementing integer) 
    appended on the type string.
    Since the data type may contain spaces and the number 
    of spaces between the datatype and identifier is not
    fixed, we need to split the string relying on the 
    assumption that the data type contains only letters 
    (upper or lower case) and spaces.

    PARAMETERS:
        ttype_value: input string containing composite
                     column identifier string

    RETURNS:
        ttype_info: dictionary contain data type string
                    and indentifier value (as an integer)
    '''

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
    
    columns = (search_key('TFORM', dictionary))[0]
    
    column_info = []
    
    for i in range(len(columns)):
        coltypes     = (search_key('TTYPE', dictionary))[0]
        column_info.append(coltypes[0])
        
        units        = (search_key('TUNIT', dictionary))[0]
        column_info.append(units)

        axislabels   = (search_key('TDESC', dictionary))[0]
        column_info.append(axislabels)

        targets      = (search_key('TOBJC', dictionary))[0]
        column_info.append(targets)

        deriv_methds = (search_key('TMTHD', dictionary))[0]
        column_info.append(deriv_methds)

        has_telluric = (search_key('TWATM', dictionary))[0]
        column_info.append(has_telluric)
    
    for a, b, c, d, e, f in zip(column_info[0],column_info[1],column_info[2],column_info[3],column_info[4], column_info[5]):
        print("-----COLUMN-----\n\nComponent type: {}\nColumn data units: {}\nAxis Labels: {}\nTarget: {}\nDerivation Method: {}\nIncludes Telluric Absorption: {}\n".format(a,b,c,d,e,f))

#
def find_column_index(dictionary, column_key):
    # extract the FITS keyword column index for the column being plotted
    # (as defined by the column_key input)
    # i.e. the keywords have names like TTITLn - we need to find the value 
    # of "n" corresponding to the column name provided by column_key

    #if any(word in keywords for word in text)
    
    keyword_match = []
    for keyword, value in dictionary.items():
        if value == column_key:
            keyword_match.append(keyword)
    print(type(keyword_match))
    
    #kw = list(dictionary.keys())[list(dictionary.values()).index(column_key)]
    #integer_search = re.compile(r'\d+(?:\.\d+)?')
    column_indices = []
    for match in keyword_match:
        print(match)
        integer_search = re.compile(r'[A-Z]+(\d+)')
        column_indices.append(integer_search.findall(match))
    
    return (column_indices, keyword_match)

def tunit_str_to_unit(fits_unit_text):
    # convert the (allowed) string values for the TUNITn keyword
    # into a proper astropy unit
    # this is required for reading the data into a Spectrum1D
    
    if fits_unit_text.lower() == 'relint' or fits_unit_text.lower() == 'flux' :
        unit_type = u.dimensionless_unscaled
    elif fits_unit_text.lower() == 'nanometers' or fits_unit_text.lower() == 'nm':
        unit_type = u.nm
    elif fits_unit_text.lower() == 'angstrom' or fits_unit_text.upper() == 'A':
        unit_type = u.nm
    elif fits_unit_text.lower() == 'cm^-1':
        unit_type = u.k     # kayser: CGS unit of wavenumber
    elif fits_unit_text.lower() == 'W m^(-2) sr^(-1) Angstrom^(-1)':
        unit_type = u.W / u.m / u.m / u.sr / u.AA     # kayser: CGS unit of wavenumber
    else:
        unit_type = u.dimensionless_unscaled
    
    return unit_type
               
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
    
    for value in (search_key('TTYPE', dictionary))[0]:
        types.append(value)
    for value in (search_key('TWATM', dictionary))[0]:
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
    print("Atlas wavelength coverage: {:.2f} {} – {:.2f} {}\n".format(file[0].header['WAVEMIN'], 
                                                                      file[0].header['CUNIT1'], file[0].header['WAVEMAX'],
                                                                      file[0].header['CUNIT1']))

# skips the primary HDU because not a binary table of data and info
    for i in range(1,len(file)):
        print("EXTENSION {}:\n".format(i))
        
        atlasdict = make_dictionary(filename,i)
        
        count = 1
            
        minimums = (search_key('TDMIN', atlasdict))[0]
        maximums = (search_key('TDMAX', atlasdict))[0]

        waverange = [x for y in zip(minimums, maximums) for x in y]
        res = list(zip(waverange, waverange[1:] + waverange[:1]))
        finalranges = res[::2]

        units       = (search_key('TUNIT', atlasdict))[0]
        column_type = (search_key('TTYPE', atlasdict))[0]
        #matchidx = re.search(r'[A-Za-z]',column_type[::-1])
        #coltype_base = (column_type[:-matchidx.start(0)]).strip()
        #coltype_idx  = (column_type[-matchidx.start(0):]).strip()

        newcount = 0
        for waveinterval in finalranges:
            column_type_in = column_type[newcount]
            matchidx       = re.search(r'[A-Za-z]',column_type_in[::-1])
            coltype_base   = (column_type_in[:-matchidx.start(0)]).strip()
            coltype_idx    = (column_type_in[-matchidx.start(0):]).strip()

            print("Column {} - {} : value range = {:.2f} – {:.2f} {}".format(newcount, 
                                    coltype_base, waveinterval[0], waveinterval[1], units[newcount]))
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
    column_types,column_types_kwrds = (atlastools.search_key('TTYPE', atlasdict))
    print("checking ",column_types)
    print("keywords ",column_types_kwrds)

    # now we need to find out which column of the table contains the wavelength
    # scale to apply to the data
    # so we will first search for all columns where the TTYPE keyword contains the
    # 'Wavelength Scale" identifier
    #
    # since there could technically be multiple columns with wavelength data 
    # (e.g. wavelength and wavenumber), we will search for the column with 
    # the lowest unique identifier index.
    # normally this should be "1", but you never know what people will do...
    wavelength_col_idx       = 9999
    wavescale_identifier     = 'Wavelength Scale'
    wavelength_primary_id    = ''
    wavelength_primary_ttype = ''
    # loop through all the TTYPEn values and find those with a wavelength scale
    for ttype_value in column_types:
        print(ttype_value)
        # this splits the TTYPE value into its (non-unique) type and the incrementing index value
        column_type = split_column_ttype(ttype_value)
        if column_type['Column Type'] == wavescale_identifier:
            # now check to see if we have found the column with a lower
            # wavelength scale identifier
            # if so, repopulate values of wavelength scale identifiers
            if column_type['Column Type Number'] < wavelength_col_idx:
                wavelength_col_idx       = column_type['Column Type Number']
                wavelength_primary_id    = ttype_value
                wavelength_primary_ttype = column_types_kwrds[column_types.index(wavelength_primary_id)]
                print("Primary Wavelength Scale",wavelength_primary_ttype,wavelength_primary_id)

    wavelength_primary_colindex = (find_column_index(atlasdict, wavelength_primary_id))[0]
    unit             = (atlastools.search_key('TUNIT' + wavelength_primary_colindex[0].pop(), atlasdict))[0]
    print(type(unit[0]),unit[0])
    wavelength_unit  = tunit_str_to_unit(unit[0])
    wavelength_scale = f[1].data[wavelength_primary_id] * wavelength_unit

    atlas_min = np.min(wavelength_scale)
    atlas_max = np.max(wavelength_scale)

    if ((startwave < atlas_min) or 
        (startwave > atlas_max)) : 
         startwave = atlas_min
    if ((endwave  <= startwave)  or 
        (endwave  >= atlas_max))  : 
         endwave   = atlas_max

    mask = ((wavelength_scale > startwave) & (wavelength_scale < endwave))
    print(np.sum(mask))
    
    # reading in the actual data per TTYPE keyword, and storing that data in a new
    # dictionary associated with the type of data it is
    stored = {}
    for value in column_types:
        full_column    = f[1].data[value]
        selected_range = full_column[mask]
        stored[value]  = selected_range
        
    # searching the extension dictionary for keywords tagged TUNIT which are the units
    # associated with the data and storing them in a list
    units_to_store = (atlastools.search_key('TUNIT', atlasdict))[0]
    
    # associating the data type with the units it required and storing in a new dictionary
    unit_tags = {}
    for key, unit in zip(column_types, units_to_store):
        unit_tags[key] = unit
        
    # searching through the unit_tags dictionary and changing the value to the proper
    # astropy unit quantity (required for reading into Spectrum1D)
    for key, value in unit_tags.items():
        unit_tags[key] = tunit_str_to_unit(value)

    # generating a new dictionary that has the data type keyword now associated with the
    # actual data with the astropy units attached (this is the metadata dictionary returned)
    file_data = {}
    for key in column_types:
        key_value            = key
        data                 = stored[key_value]
        unit                 = unit_tags[key_value]
        combined             = data * unit
        file_data[key_value] = combined
    
    # Finally, pushing the file metadata through Spectrum1D to create a dictionary of Spec1D objects
    # associated with the type of data they are (this is the Spec1D dictionary returned)
    # ( I intend to come back and fix this to use a regex string pattern of 'Wavelength Scale' rather
    # than nested if statements)
    spec_data = {}
    for key in file_data:
        if re.search(r'Wavelength Scale*', key) is None:
                col_index = (find_column_index(atlasdict, key))[0].pop()
                finalized = Spectrum1D(spectral_axis=file_data['Wavelength Scale   1'], flux=file_data[key])
                
                finalized.meta.update(unit            = (atlastools.search_key('TUNIT' + col_index[0], atlasdict))[0])
                finalized.meta.update(has_telluric    = (atlastools.search_key('TWATM' + col_index[0], atlasdict))[0])
#                finalized.meta.update(has_solar       = (atlastools.search_key('TWSUN' + col_index[0], atlasdict))[0])
                finalized.meta.update(observed_object = (atlastools.search_key('TOBJC' + col_index[0], atlasdict))[0])
                finalized.meta.update(observed_target = (atlastools.search_key('TTRGT' + col_index[0], atlasdict))[0])
                finalized.meta.update(col_title       = (atlastools.search_key('TTITL' + col_index[0], atlasdict))[0])
                finalized.meta.update(col_label       = (atlastools.search_key('TLABL' + col_index[0], atlasdict))[0])
                finalized.meta.update(col_description = (atlastools.search_key('TDESC' + col_index[0], atlasdict))[0])
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

    obslocation = EarthLocation(lat   = (atlastools.search_key('ATL_LAT', atlasdict))[0], 
                                lon   = (atlastools.search_key('ATL_LONG', atlasdict))[0], 
                                height= (atlastools.search_key('ATL_ALT', atlasdict))[0])
    
    atlas_obs = observatory((atlastools.search_key('ATL_OBS', atlasdict))[0], obslocation, 'FTS')
    #atlas_obs = observatory('Kitt Peak', EarthLocation.of_site('Kitt Peak'), 'FTS')
    atlas = atlastools.atlas(0, 0, 1, 
                             (atlastools.search_key('OBJECT', atlasdict))[0], 
                             (atlastools.search_key('ATL_SOUR', atlasdict))[0], atlas_obs)
    
    # Issue #1 - hardcoded extraction of column information into spec1D objects
    if loaddata is not 0:
    
        (file_data, spec_data) = store_data(filename, extension, startwave=startwave, endwave=endwave)
        search_key = 'Local Intensity   1'
        spect_sun  = spec_data.get(search_key)
        if spect_sun is not None:
            print("populating solar atlas object")
            atlas.sun = spect_sun
            spect_sun_index = (find_column_index(atlasdict, search_key))[0].pop()
            atlas.sun.meta.update(title=(atlastools.search_key('TTITL' + spect_sun_index[0], atlasdict))[0])

        search_key = 'Telluric Spectrum   1'
        spect_atm  = spec_data.get(search_key)
        if spect_atm is not None:
            print("populating telluric atlas object")
            spect_atm_index = (find_column_index(atlasdict, search_key))[0].pop()
            atlas.atm = spect_atm
            atlas.atm.meta.update(title=(atlastools.search_key('TTITL' + spect_atm_index[0], atlasdict))[0])

        atlas.components = spec_data
    
        # Issue #2 - need more column specific information associated with Spec1D object
        # is putting a range of individual items in the meta dictionary the best way to do this?
        
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
