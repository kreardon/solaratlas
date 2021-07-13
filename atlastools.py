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
from pathlib import Path

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

class atlas_file_extension_info:
    def __init__(atlas_file_exten_obj, wavemin, wavemax, col_length, numcols, coltypes, coltypes_idx, has_solar, has_telluric, colvalue_min, colvalue_max, colvalue_units, col_targets, col_radialdists, col_titles):
        atlas_file_exten_obj.wavemin        = wavemin         # unit aware minimum wavelength in extension (e.g. 400 * u.nm)
        atlas_file_exten_obj.wavemax        = wavemax         # unit aware maximum wavelength in extension (e.g. 700 * u.nm)
        atlas_file_exten_obj.col_length     = col_length      # length of extension columns (e.g. 324561)
        atlas_file_exten_obj.numcols        = numcols         # number of columns in extension
        atlas_file_exten_obj.coltypes       = coltypes        # list of column types from TTYPEn (e.g. ['Wavelength Scale', 'Local Intensity', etc. 2'])
        atlas_file_exten_obj.coltypes_idx   = coltypes_idx    # list of column type indices from TTYPEn (e.g. ['1', '1', '2'])
        atlas_file_exten_obj.has_solar      = has_solar       # list of values of TWSOLn indicating presence of solar component in column (e.g. [-1, 1, 1, 0])
        atlas_file_exten_obj.has_telluric   = has_telluric    # list of values of TWATMn indicating presence of solar component in column (e.g. [-1, 0, 1, 1])
        atlas_file_exten_obj.colmins        = colvalue_min    # list of minimum value in each column (e.g. [300, -0.05, 0.0])
        atlas_file_exten_obj.colmaxs        = colvalue_max    # list of maximum value in each column (e.g. [700, 1.05, 1.0])
        atlas_file_exten_obj.colunits       = colvalue_units  # list of unit type in each column from TUNITn (e.g. ['nanometers', 'relative intensity', 'flux'])
        atlas_file_exten_obj.targets        = col_targets     # list of column targets from TTRGTn (e.g. ['Quiet Sun', 'Quiet Sun', 'Quiet Sun'])
        atlas_file_exten_obj.radists        = col_radialdists # list of radial distance of target region from TRADn (e.g. [0.0, 0.0, 0.989000])
        atlas_file_exten_obj.titles         = col_titles      # list of column titles from TDESCn or TTITLn (e.g. [ 'Stenflo/SS3 Disk Center Atlas', ..., 'Stenflo/SS3 Limb Atlas' ])

class atlas_file_content:
    def __init__(atlas_file_obj, filename, filepath, source, observatory, version, obsobject, extension_info):
        atlas_file_obj.filename       = filename              # name of source file being described
        atlas_file_obj.filepath       = filepath              # full or relative path to source file being described
        atlas_file_obj.source         = source                # description of atlas source from ATL_SOUR
        atlas_file_obj.observatory    = observatory           # name of observatory where atlas was acquired from ATL_OBS
        atlas_file_obj.version        = version               # version number of atlas file
        atlas_file_obj.obsobject      = obsobject             # observed object for atlas
        atlas_file_obj.extension_info = extension_info        # list of atlas_file_extension_info objects contain details on individual extensions in file

    def file_info(atlas_filemap_obj):
        '''
        Lists file top-line information, details, number of extensions
        NOTE: this method assumes a filemap object had been created

        PARAMETERS:
        ----------
            [None]

        RETURNS:
        -------
            [None, printed information]

        EXAMPLE CALL:
        ------------
            [filemap_object_variable].file_info()
        '''
        filename = atlas_filemap_obj.filename
        file = fits.open(filename)
        print("General file information:\n")

        print("Object: {}".format(file[0].header['OBJECT']))
        print("Atlas source: {}".format(file[0].header['ATL_SOUR']))
        print("Atlas acquisition site: {}".format(file[0].header['ATL_OBS']))
        print("Atlas wavelength coverage: {:.2f} {} – {:.2f} {}".format(file[0].header['WAVEMIN'],
                                                                          file[0].header['CUNIT1'], file[0].header['WAVEMAX'],
                                                                          file[0].header['CUNIT1']))
        print("Number of Extensions: {}\n".format(len(file)))

    def list_extensions(atlas_filemap_obj):
        '''
        Loops through atlas file extensions and prints column information e.g.,
        waveranges & components.
        NOTE: this method assumes a filemap object had been created

        PARAMETERS:
        ----------
            [None]

        RETURNS:
        -------
            [None, printed information]

        EXAMPLE CALL:
        ------------
            [filemap_object_variable].list_extensions()
        '''
        filename = atlas_filemap_obj.filename
        file = fits.open(filename)

        ext_count = 1 # creating an extension count to increment
        # starting the extension loop (skips the primary HDU because not a binary table of data and info)
        for i in range(1,len(file)):

            atlasdict = make_dictionary(filename,i)

            fullminimums = search_key('TDMIN', atlasdict)
            minimums = fullminimums[0] # indexing the first list that the function 'search_key' returns ('search_key' returns a list of values and a list of keywords)
            fullmaximums = search_key('TDMAX', atlasdict)
            maximums = fullmaximums[0] # indexing the first list that the function 'search_key' returns
            waverange = [x for y in zip(minimums, maximums) for x in y]
            res = list(zip(waverange, waverange[1:] + waverange[:1]))
            finalranges = res[::2]

            fullunits       = search_key('TUNIT', atlasdict)
            units = fullunits[0] # indexing the first list that the function 'search_key' returns
            column_type = search_key('TTYPE', atlasdict)

            print("EXTENSION {}:\n".format(ext_count)) # printing extension number per column info

            # looping through each column within the bintable for waveranges:
            column_count = 0 # creating a column count to increment
            print("\tWaveranges per column:\n") # labeling the following printouts for clarity
            for i in finalranges:
                print("\tColumn {} wavelength range: {} {}".format(column_count, i, units[column_count]))
                column_count+=1
            print("\n")

            # looping through each column within the bintable for component types:
            ttype_count = 0 # making a ttype_count to increment
            ttype_search = search_key('TTYPE', atlasdict)
            split = []
            print("\tSpecific column components:\n") # labeling the printout of 'split_column_ttype' function for clarity
            for string in ttype_search[0]: # indexing the first list that the function 'search_key' returns and looping through it
                print("\tColumn {}:".format(ttype_count)) # associating column number with component
                value = split_column_ttype(string, printout=True)
                split.append(value)
                ttype_count+=1 # incrementing the ttype count

            ext_count+=1 # incrementing extension count
            print("\n")


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
    Generates a dictionary based on file specified extension keywords and
    values.

    PARAMETERS:
    ----------
        filename: path of file to generate dictionary from (must be a string)

        extension: specific extension of file to generate dictionary from
                   (integer)

    RETURNS:
    -------
        dictionary: the dictionary version of file extension's keywords and
                    values.
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
    Searches through a dictionary (generated for file via function
    'make_dictionary') for specified keyword and returns a list with all values
    matching that keyword.

    PARAMETERS:
    ----------
        keyword_string: keyword the user wants to search the dictionary for
                        (must be a string)

        dictionary: a python dictionary containing all keywords and values for
                    specified atlas file (generated via function
                    'make_dictionary')
    RETURNS:
    -------
        res: (n,2) list containing the matched keyword and the value associated
             with that keyword
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
def filesearch(additional_directory=False):
    '''
    Directory and file searching function that prints lists of atlas-specific
    files in user's current directory + subdirectories, directory where
    atlastools.py resides + subdirectories, and the atlasfiles/ directory +
    subdirectories if it exists.

    PARAMETERS:
    ----------
        **NO PARAMETER(S) REQUIRED**

        multi_directory: optional parameter to have the function search through
                         user-specified directories for atlas files
                         (TYPE: string)

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

    # code for searching additional directories provided by user outside of auto searches
    if additional_directory == True:
        directories = []
        print("To exit additional directory input, type 'exit'")
        while True:
            prompt = input("Additional directory to search: ")
            input_directory = str(prompt)
            directories.append(input_directory)

            if prompt == 'exit':
                break
        for directory in directories:
            if directory == 'exit':
                directories.remove('exit')
        print("\nDirectories to search:\n\t{}".format(directories))

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

    # going through the final additional_directory search:

    # defining functions needed
    def directory_walk(directory):
        while True:
            if os.path.isdir(directory) == True:
                os.chdir(directory)
                break
            else:
                os.chdir('..')

    def find_atlases(directory):
        import glob
        # marking starting directory to go back to
        start = os.getcwd()
        # walking up to specified directory
        current = directory_walk(directory)
        pattern = "*.fits"
        pattern2 = "*.gz"

        atlases = []
        fits = glob.glob(pattern)
        gz = glob.glob(pattern2)
        # only storing the file name if glob search is successful
        if fits:
            atlases.append(fits)
        if gz:
            atlases.append(gz)

        # moving back to the starting place for user
        os.chdir(start)
        print("\nADDITIONAL SEARCH---")
        print("\nAtlases found in '{}/':\n".format(directory))
        for i in atlases:
            for j in i:
                print(j)

    # the additional directory search itself:
    if additional_directory == True:
        for i in directories:
            find_atlases(i)

#
def split_column_ttype(ttype_value, printout=False):
    '''
    Splits a value provided in the TTYPEn FITS header into its base data type
    (e.g. 'Wavelength Scale') and the unique identifier (an incrementing
    integer) appended on the type string. Since the data type may contain spaces
    and the number of spaces between the datatype and identifier is not fixed,
    we need to split the string relying on the assumption that the data type
    contains only letters (upper or lower case) and spaces.

    PARAMETERS:
        ttype_value: input string containing composite column identifier string

        printout:    optional keyword argument; if True, will print out column
                     components

    RETURNS:
        ttype_info: dictionary contain data type string and indentifier value
                    (as an integer)
    '''

    ttype_info     = {}
    matchidx       = re.search(r'[A-Za-z]',ttype_value[::-1])
    coltype_base   = (ttype_value[:-matchidx.start(0)]).strip()
    coltype_idx    = (ttype_value[-matchidx.start(0):]).strip()
    ttype_info['Column Type']        = coltype_base
    ttype_info['Column Type Number'] = int(coltype_idx)
    if printout == True:
        print("\t",coltype_base,"  --  ",coltype_idx)

    return ttype_info

#
def column_information(dictionary):
    '''
    Basic printout of various column information within the fits files

    PARAMETERS:
    ----------
        dictionary: dictionary (created via 'make_dictionary' function) of the
                    file being used

    RETURNS:
    -------
        [*NO VALUE RETURNED* - prints information]
    '''

    columns = (search_key('TFORM', dictionary))[0]

    column_info = []

    for i in range(len(columns)):
        colnum      = (search_key('COLNUM', dictionary))[0]
        column_info.append(colnum)

        coltypes     = (search_key('TTYPE', dictionary))[0]
        column_info.append(coltypes)

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

    for a, b, c, d, e, f, g in zip(column_info[0],column_info[1],column_info[2],column_info[3],column_info[4], column_info[5], column_info[6]):
        print("{}\n\nComponent type: {}\nColumn data units: {}\nAxis Labels: {}\nTarget: {}\nDerivation Method: {}\nIncludes Telluric Absorption: {}\n".format(a,b,c,d,e,f,g))

#
def find_column_index(dictionary, column_key, keyword_base=False):
    # extract the FITS keyword column index for the column being plotted
    # (as defined by the column_key input)
    # i.e. the keywords have names like TTITLn - we need to find the value
    # of "n" corresponding to the column name provided by column_key

    #if any(word in keywords for word in text)
    if keyword_base == False:
        keyword_base = '[0-9]$'

    keyword_match = []
    for keyword, value in dictionary.items():
        if value == column_key and re.search(keyword_base,keyword):
            keyword_match.append(keyword)
    #print("type - ",type(keyword_match))

    #kw = list(dictionary.keys())[list(dictionary.values()).index(column_key)]
    #integer_search = re.compile(r'\d+(?:\.\d+)?')
    column_indices = []
    for match in keyword_match:
    #    print("match - ", match)
        integer_search = re.compile(r'[A-Z]+(\d+)')
        column_indices.append(integer_search.findall(match))

    return (column_indices, keyword_match)

def tunit_str_to_unit(fits_unit_text):
    # convert the (allowed) string values for the TUNITn keyword
    # into a proper astropy unit
    # this is required for reading the data into a Spectrum1D

    if fits_unit_text.lower() == 'relint' or fits_unit_text.lower() == 'flux' :
        unit_type = u.dimensionless_unscaled
    elif fits_unit_text.lower() == 'nanometers' or fits_unit_text.lower() == 'nm' or fits_unit_text.lower() == 'Nanometers':
        unit_type = u.nm
    elif fits_unit_text.lower() == 'angstrom' or fits_unit_text.upper() == 'A' or fits_unit_text.lower() == 'Angstrom':
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
    Return series of print statements notifying whether telluric absoprtion is
    included in stored data.

    PARAMETERS:
        dictionary: a python dictionary containing all keywords and values for
                    specified atlas file (generated via function
                    'make_dictionary')

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
    Content map of the atlas file, including general information and extension
    information, returned as an atlas object.

    PARAMETERS:
    ----------
        filename: path of file to provide information on (must be a string)

    RETURNS:
    -------
        atlas_filemap: atlas object containing both atlas_file_content &
                       atlas_extention_info (objectname.extension_info[n] for n
                       extensions)
    '''

    file = fits.open(filename, memmap=True)

    atlas_extension_all = []         # creating the end-all list that will have all extension objects in it
    ext_count = 1                    # creating an extension count to increment
    for i in range(1,len(file)):     # skips the primary HDU because not a binary table of data and info

        atlasdict = make_dictionary(filename,i)

        fullminimums = search_key('TDMIN', atlasdict)
        minimums = fullminimums[0]
        fullmaximums = search_key('TDMAX', atlasdict)
        maximums = fullmaximums[0]
        waverange = [x for y in zip(minimums, maximums) for x in y]
        res = list(zip(waverange, waverange[1:] + waverange[:1]))
        finalranges = res[::2]

        fullunits       = search_key('TUNIT', atlasdict)
        units = fullunits[0]
        column_type = search_key('TTYPE', atlasdict)


        # finding file-specific header unit to tag onto integers-- will add in a lot of other if/elif statements
        # like what's found in the tunit_str_to_unit function (Unless there's a good way to strip single quotes
        # off of a string that I'm unaware of?)
        header_unit = file[0].header['CUNIT1']

        tunit_str_to_unit
        if header_unit == 'nm' or 'nanometers':
            unit = u.nm

        #wave_col_id  = file[0].header['ATLWVCOL']
        wave_col_id  = (search_key('ATLWVCOL',atlasdict))[0]
        waveref_col  = (find_column_index(atlasdict, wave_col_id[0]))[0].pop()
        wavemin      = search_key('TDMIN' + waveref_col[0], atlasdict)
        wavemax      = search_key('TDMAX' + waveref_col[0], atlasdict)
        waveunit_str = (search_key('TCUNI' + waveref_col[0], atlasdict))[0]
        waveunit     = tunit_str_to_unit(waveunit_str[0])
        wavemin      = wavemin[0] * waveunit
        wavemax      = wavemax[0] * waveunit
        
        #wavemin      = file[0].header['WAVEMIN'] * unit
        #wavemax      = file[0].header['WAVEMAX'] * unit
        col_length   = (search_key('NAXIS2',atlasdict))[0]

        # using the atlastools.split_column_ttype function to separate list of TTYPE strings from TTYPE integer:
        ttype_search = search_key('TTYPE', atlasdict)
        numcols      = len(ttype_search[1])
        split = []
        types = []
        indicies = []

        for string in ttype_search[0]:
            value = split_column_ttype(string)
            split.append(value)
        for dictionary in split:
            types.append(dictionary['Column Type'])
        for dictionary in split:
            indicies.append(dictionary['Column Type Number'])

        
        coltypes     = types
        coltypes_idx = indicies
        has_solar    = (search_key('TWSOL', atlasdict))[0]
        has_telluric = (search_key('TWATM', atlasdict))[0]
        colmins      = (search_key('TDMIN', atlasdict))[0]
        colmaxs      = (search_key('TDMAX', atlasdict))[0]
        colunits     = (search_key('TUNIT', atlasdict))[0]
        targets      = (search_key('TTRGT', atlasdict))[0]
        radists      = (search_key('TRAD', atlasdict))[0]
        titles       = (search_key('TTITL', atlasdict))[0]


        atlas_extension = atlas_file_extension_info(wavemin, wavemax, col_length, numcols, coltypes, coltypes_idx, has_solar, has_telluric, colmins, colmaxs, colunits, targets, radists, titles)

        # adding the finalized extension object to the large list of objects:
        atlas_extension_all.append(atlas_extension)
        ext_count+=1 # incrementing extension count


    inputfilename_fullpath   = (Path(filename).resolve())
    inputfilename_directory  = inputfilename_fullpath.parent
    source                   = file[0].header['ATL_SOUR']
    observatory              = file[0].header['ATL_OBS']
    version                  = '1.1'
    obsobject                = file[0].header['OBJECT']
    extension_info           = atlas_extension_all

    atlas_filemap = atlas_file_content(filename, inputfilename_directory, source, observatory, version, obsobject, atlas_extension_all)

    return atlas_filemap

#
def store_data(filename, extension, startwave=1*u.nm, endwave=1*u.nm):
    '''
    Function that stores all data within specified extension and returns two
    dictionaries: one of the data from the file with units and one of that same
    data that has been pushed through Spectrum1D.

    PARAMETERS:
    ----------
        filename:   full name of the file (must be a string)
        extension:  specific extension to pull and store
                    data from (integer)
        startwave:  Wavelength start for data to be generated at (integer must
                    have astropy unit quantity)
        endwave:    Wavelength end for data to be generated to (integer must
                    have astropy unit quantity)

    RETURNS:
    -------
        file_data: data with proper astropy units straight from the file
                   associated with the data type- essentially the metadata
                   straight from file(dictionary with {data type label : value})
        spec_data: metadata that has been fed through Spectrum1D to generate
                   usable atlas data (dictionary with data type label : Spec1D
                   object)
    '''

    import numpy as np
    from astropy.io import fits
    import atlastools
    from astropy import units as u
    from specutils import Spectrum1D

    # opening the .fits file for access
    f = fits.open(filename)

    # generating an extension-specific dictionary
    atlasdict = make_dictionary(filename, extension) #%#% removed atlastools.make_dictionary

    # searching the extension dictionary for keywords tagged TTYPE which are the data
    # and storing them in a list
    column_types,column_types_kwrds = (search_key('TTYPE', atlasdict)) #%#% removed atlastools.search_key
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

    wavelength_primary_colindex = (find_column_index(atlasdict, wavelength_primary_id, keyword_base='TTYPE'))[0]
    print(wavelength_primary_id,wavelength_primary_colindex)
    unit             = (search_key('TUNIT' + wavelength_primary_colindex[0].pop(), atlasdict))[0] #%#% removed atlastools.search_key
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
    units_to_store = (search_key('TUNIT', atlasdict))[0] #%#% removed atlastools.search_key

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

                finalized.meta.update(col_type        = (search_key('TTYPE' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(unit            = (search_key('TUNIT' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(has_telluric    = (search_key('TWATM' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(has_solar       = (search_key('TWSOL' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(observed_object = (search_key('TOBJC' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(observed_target = (search_key('TTRGT' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(col_title       = (search_key('TTITL' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(col_label       = (search_key('TLABL' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(col_description = (search_key('TDESC' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                finalized.meta.update(col_method      = (search_key('TMTHD' + col_index[0], atlasdict))[0]) #%#% removed atlastools.search_key
                spec_data[key] = finalized

    return file_data, spec_data

#
def make_atlas(filename, extension, startwave=1*u.nm, endwave=1*u.nm, loaddata=0):
    '''
    Generates an atlas object based on file and extension specified using
    extension keywords and values.

    PARAMETERS:
        filename: path of file to generate dictionary from (must be a string)

        extension: specific extension of file to generate dictionary from
                   (integer)

    RETURNS:
        dictionary: the dictionary version of file extension's keywords and
                    values.
    '''
    from astropy.coordinates import EarthLocation
    import atlastools

    atlasdict = make_dictionary(filename, extension) #%#% removed atlastools.make_dictionary

    obslocation = EarthLocation(lat   = (atlastools.search_key('ATL_LAT', atlasdict))[0],
                                lon   = (atlastools.search_key('ATL_LONG', atlasdict))[0],
                                height= (atlastools.search_key('ATL_ALT', atlasdict))[0])

    atlas_obs = observatory((atlastools.search_key('ATL_OBS', atlasdict))[0], obslocation, 'FTS')
    #atlas_obs = observatory('Kitt Peak', EarthLocation.of_site('Kitt Peak'), 'FTS')
    atlas = atlastools.atlas(0, 0, 1, (atlastools.search_key('OBJECT', atlasdict))[0],
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

#        spec1d_all = []
#        for key in spec_data:
#            spec1d_all.append(spec_data[key])
        atlas.components = spec_data

        # Issue #2 - need more column specific information associated with Spec1D object
        # is putting a range of individual items in the meta dictionary the best way to do this?

    return atlas

#
def atlas_spectrum_plot(column_key, atlas_obj, startwave=1*u.nm, endwave=1*u.nm, plot_unit='AA'):
    '''
    Function that returns a single plot of the provided data given associated
    dictionary.

    PARAMETERS:
    ----------
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
    -------
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
