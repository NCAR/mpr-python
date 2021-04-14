
import os
import warnings

import xarray as xr

# TODO:

GEO_ATTR_ROOT_DIR = '/glade/p/ral/hap/mizukami/pnw-extrems/geospatial_data/geophysical/data_mpr'

#  var_type: {var_name: var_name_in_file}
#  file name is f{var_type}_data.nc
GEO_ATTR_TYPE = {'climate':    {'precipitation':'prec', 'tempeature':'tavg', 'wind':'wind', 'humidity':'rh', 'aridity':'ai'}, 
                 'soil':       {'sand':'sand_pct', 'clay':'clay_pct', 'silt':'silt_pct', 'bulk_density':'bulk_density', 'organic_carbon':'soc'},
                 'vegetation': {'IGPG':'vegclass', 'canopy_height':'ch', 'lai':'lai'}, 
                 'topography': {'ele':'ele_mean', 'roughness':'ele_std_dev', 'slope':'slope_mean'}
                 }

mapping_file = 'spatialweights_grid_600m_to_HUC12.nc'


class NoneError(Exception):
    pass


def progress(r):
    try:
        from tqdm import tqdm
        return tqdm(r)
    except ImportError:
        return r


def load_geophysical_attributes(root=GEO_ATTR_ROOT_DIR, geotype=GEO_ATTR_TYPE,  **kwargs):
    print('loading geophysical attributes...', flush=True)

    if isinstance(geotype, dict):
        geolist = GEO_ATTR_TYPE.keys()
    elif isinstance(geotype, list):
        geolist = GEO_ATTR_TYPE
    else:
        raise NoneError("geotype must be provided in dictionary key or list")

    attr_data = xr.Dataset()
    for categ in geolist:
        attr_data = attr_data.merge(xr.open_dataset(os.path.join(root, f'{categ}_data.nc'), **kwargs))

    return attr_data


def load_mapping_data(root=GEO_ATTR_ROOT_DIR, file=mapping_file, var_list=None, **kwargs):
    print('loading mapping weight data...', flush=True)
    
    with xr.open_dataset(os.path.join(root, file), **kwargs) as ds:
       
        drop_var = []
        if var_list is not None:
            
            for var in var_list:
                if not var in ds.variables:
                    warnings.warn('variable: "%s" not exist' % var)
                    
            drop_var = [var for var in ds.variables if not var in var_list]
            
        ds = ds.drop_vars(drop_var)
    
    return ds