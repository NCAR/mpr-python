
import sys, os
import warnings

from datetime import datetime
import numpy as np
import xarray as xr

FILL_VALUE = -9999.0


class NoneError(Exception):
    pass


def progress(r):
    try:
        from tqdm import tqdm
        return tqdm(r)
    except ImportError:
        return r


def load_geophysical_attributes(io_cfg,  **kwargs):
    print('\n Loading geophysical attributes...', flush=True)

    root=io_cfg['INPUT']['DIRE']
    geotype=io_cfg['GEO_ATTR_TYPE']

    if isinstance(geotype, dict):
        geolist = geotype.keys()
    elif isinstance(geotype, list):
        geolist = geotype
    else:
        raise NoneError("geotype must be provided in dictionary key or list")

    attr_data = xr.Dataset()
    for categ in geolist:
        attr_data = attr_data.merge(xr.open_dataset(os.path.join(root, f'{categ}_data.nc'), **kwargs))

    return attr_data


def load_mapping_data(io_cfg, var_list=None, **kwargs):
    print('\n Loading mapping weight data...', flush=True)

    root=io_cfg['INPUT']['DIRE']
    file=io_cfg['INPUT']['MAP_FILE']

    with xr.open_dataset(os.path.join(root, file), **kwargs) as ds:

        drop_var = []
        if var_list is not None:

            for var in var_list:
                if not var in ds.variables:
                    warnings.warn('variable: "%s" not exist' % var)

            drop_var = [var for var in ds.variables if not var in var_list]

        mat_data = process_mapping_data(ds.drop_vars(drop_var), io_cfg, var_list=var_list)

    return mat_data


def process_mapping_data(mapping_data, io_cfg, var_list=None):
    print("\n Pre-process mapping data arrays")

    # input:  xarray dataset,  mapping_data
    #         list,            mapping variable list to be processed
    # return: dictionary,      {variable name: numpy array}

    # data dimension variable reformat
    #   mat_weight      [nHRU x maxOverlaps]
    #   mat_i_index     [nHRU x maxOverlaps]
    #   mat_j_index     [nHRU x maxOverlaps]
    #   mat_intersector [nHRU x maxOverlaps]

    # enforce to include ['polyid', 'overlaps', 'weight']
    if not all(v in var_list for v in ['polyid', 'overlaps', 'weight']):
        raise Exception('Sorry, you need to include "polyid" AND "overlaps" AND "weight" in var_list')

    # remove variables in var_list NOT exist in MAPPING_VARS_META
    var_list = [var for var in var_list if var in io_cfg['MAPPING_VARS_META'].keys()]

    var_dic = {}
    for var in var_list:
        var_dic[var] = mapping_data[var].values
        if var in ['i_index', 'j_index']:
             # j_index is 1-based and starts at south ... make 0-based
             # i_index is 1-based and starts at West ... make 0-based
            var_dic[var] = var_dic[var] - 1

    # check number of target HRUs, maximum number of overlapping source polygons
    nHRU = len(var_dic['polyid'])
    maxOverlaps = var_dic['overlaps'].max()

    # Check sum of overlapping polygon number is equal to data dimension
    flag = True
    if var_dic['overlaps'].sum() == len(var_dic['weight']):
        flag = False

    # convert 1D data dimension array to a matrix (nHRU x maxOverlaps)
    mat_dic = {}
    for var in var_list:
        if io_cfg['MAPPING_VARS_META'][var]['dim'] == 'data':
            mat_dic[var] = np.zeros((nHRU, maxOverlaps), dtype=io_cfg['MAPPING_VARS_META'][var]['type'])

    ix2=0;
    for p in range(0, nHRU):
        if var_dic['overlaps'][p]>0:
            ix1 = ix2
            ix2 = ix1+var_dic['overlaps'][p]
        elif var_dic['overlaps'][p]==0:
            if flag: # NOTE: grid2poly.py put overlaps=0 output skip data in data dimension, but poly2poly put overlaps=1 and missing values in data dimension
                ix1 = ix2
                ix2 = ix2+1
            matWgts[p, 0] = 1.0
            continue

        for var in var_list:
            if io_cfg['MAPPING_VARS_META'][var]['dim'] == 'data':
                if var == 'weight': # will normalize in case sum of weight is not 1
                    mat_dic[var][p, 0:var_dic['overlaps'][p]] = var_dic[var][ix1:ix2]/var_dic[var][ix1:ix2].sum()
                else:
                    mat_dic[var][p, 0:var_dic['overlaps'][p]] = var_dic[var][ix1:ix2]
            else: # add polyid dimension
                mat_dic[var] = var_dic[var]

    return mat_dic


def write_param(param_data, param_meta, hru_data, hru_meta, io_cfg, par_list=None, output=None):
    print("\n Write parameters", flush=True)

    if output is None:
        output = os.path.join(io_cfg['OUTPUT']['DIRE'], io_cfg['OUTPUT']['PARAM_FILE_NAME'])

    ds = xr.Dataset()
    ds[hru_meta['name']]=((hru_meta['dim']), hru_data)

    if par_list is None:
        par_list = []
        for par, meta in param_meta.items():
            if meta['compute']:
                par_list.append(par)
        #par_list = [*param_data]

    else:  # Make sure that parameter in [par_list] exist in data, remove if not exist
        for par in par_list:
            if not par in param_data.keys():
                warnings.warn(f'parameter: "{par}" not exist')
                par_list.remove(par)

    for name in par_list:
        ds[name] = ((param_meta[name]['dim']), param_data[name])

        ds[name].attrs['long_name'] = param_meta[name]['long_name']
        ds[name].attrs['units'] = param_meta[name]['units']
        ds[name].encoding['_FillValue'] = FILL_VALUE

    history = '{}: {}\n'.format(datetime.now().strftime('%c'),' '.join(sys.argv))
    ds.attrs={'Conventions':'xxxx', 'title':'Hydrologic model parameter', 'history':history}

    ds.to_netcdf(output)
