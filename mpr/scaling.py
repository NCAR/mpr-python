
import os
import warnings
import numpy as np

# TODO:

FILL_VALUE = -9999.0

def horizontal_weighted_mean(mapping_data, origArrays, pvalue, default=FILL_VALUE):
    """ Brief: Compute areal weighted generalized mean value of for all target HRUs, given pvalue

        Details:
        input:  mapping_data, reformated mapping data,        dictionary {mapping variable name: numpy array}
                ogirArray,    fine resolution parameter grid, numpy array (2D or 3D see below)
                pvalue,       parameter in generalized mean operator
        return: wgtedVal,     remapped parameter array,

        ogirArray: 2D [lat, lon] -      -> wgtedVal: 1D [hru]
        ogirArray: 3D [Month, lat, lon] -> wgtedVal: 2D [Month, hru]
                      [lyr, lat, lon]   -> wgtedVal: 2D [lyr, hru]

    """
    # numpy broadcasting rule
    # https://numpy.org/doc/stable/user/basics.broadcasting.html

    array_shape = origArrays.shape
    nDims       = len(array_shape)

    # TODO
    # move these out of function
    nOutHRUs    = len(mapping_data['overlaps'])
    maxOverlaps = mapping_data['overlaps'].max()

    if nDims == 2:
        wgtedVals   = np.zeros((nOutHRUs), dtype='float32')
        matDataVals = np.zeros((nOutHRUs, maxOverlaps), dtype='float32')
    elif nDims == 3:
        wgtedVals   = np.zeros((array_shape[0], nOutHRUs), dtype='float32')
        matDataVals = np.zeros((array_shape[0], nOutHRUs, maxOverlaps), dtype='float32')
    else:
        pass # add error check - array with 4 or larger dimension is not supported.

    # reformat var data into regular matrix matching weights format (nOutPolygons, maxOverlaps)
    #   used advanced indexing to extract matching input grid indices
    for iHru in range(0, nOutHRUs):
        if mapping_data['overlaps'][iHru]>0:
            if nDims == 2:
                matDataVals[iHru, 0:mapping_data['overlaps'][iHru]] = \
                    origArrays[mapping_data['j_index'][iHru, 0:mapping_data['overlaps'][iHru]], mapping_data['i_index'][iHru, 0:mapping_data['overlaps'][iHru]] ]
            elif nDims == 3:
                matDataVals[:, iHru, 0:mapping_data['overlaps'][iHru]] = \
                    origArrays[:, mapping_data['j_index'][iHru, 0:mapping_data['overlaps'][iHru]], mapping_data['i_index'][iHru, 0:mapping_data['overlaps'][iHru]] ]
        else:
            if nDims == 2:
                matDataVals[iHru, 0] = default
            elif nDims == 3:
                matDataVals[:, iHru, 0] = default

    # adjust weight matrix again (put zero where parameter is nan and re-scale weight)
    if nDims == 2:
        idx = np.where(np.isnan(matDataVals[:,:]))
    elif nDims == 3:
        idx = np.where(np.isnan(matDataVals[0,:,:]))
    mapping_data['weight'][idx] = 0
    sum_weight = np.tile(1.0/np.sum(mapping_data['weight'], axis=1), (mapping_data['weight'].shape[1],1)).T
    weight = sum_weight * mapping_data['weight']

    # Compute mean at all the target hrus
    if abs(pvalue) < 0.00001: # geometric mean
        wgtedVals = exp(np.nansum(log(matDataVals)* weight, axis=nDims-1))
    else:
        wgtedVals = np.nansum((matDataVals**pvalue) * weight, axis=nDims-1) **(1.0/pvalue)   # produces vector of weighted values

    return wgtedVals


def vertical_weighted_mean(mapping_data, origArrays, pvalue, default=FILL_VALUE):
    """Compute areal weighted generalized mean value of for vertical layer
       """
    # numpy broadcasting rule
    # https://numpy.org/doc/stable/user/basics.broadcasting.html

    # ogirinal grid: 3D [lat, lon, soil_lyr] -> target grid: 3D [lat, lon, model_lyr]

    orgArrays_reshaped = np.moveaxis(origArrays, 0, -1)
    array_shape        = orgArrays_reshaped.shape  # original array [soil_lyr, lat, lon] -> [lat, lon, soil_lyr]
    nDims              = len(array_shape)

    # TODO
    # move these out of function
    nMlyr       = len(mapping_data['overlaps'])
    maxOverlaps = mapping_data['overlaps'].max()

    if nDims == 3:
        wgtedVals   = np.zeros((array_shape[0], array_shape[1], nMlyr), dtype='float32')
        matDataVals = np.full((array_shape[0], array_shape[1], nMlyr, maxOverlaps), np.nan, dtype='float32')
    else:
        pass # add error check - array with other dimension is not supported.

    # reformat var data into regular matrix matching weights format (nOutPolygons, maxOverlaps)
    #   used advanced indexing to extract matching input grid indices
    for ixm in range(0, nMlyr):
        if mapping_data['overlaps'][ixm]>0:
            if nDims == 3:
                matDataVals[:, :, ixm, 0:mapping_data['overlaps'][ixm]] = \
                    orgArrays_reshaped[:, :, mapping_data['soil_index'][ixm,0:mapping_data['overlaps'][ixm]]]
        else:
            if nDims == 3:
                matDataVals[:, :, ixm, 0] = default

    weight = np.broadcast_to(mapping_data['weight'], (array_shape[0], array_shape[1], *mapping_data['weight'].shape ))
    if abs(pvalue) < 0.00001: # geometric mean
        wgtedVals = exp(_nansum(log(matDataVals)* weight, axis=nDims))
    else:
        wgtedVals = _nansum(matDataVals**pvalue * weight, axis=nDims) **(1.0/pvalue)   # produces vector of weighted values

    return np.moveaxis(wgtedVals, -1, 0)


def _nansum(a, **kwargs):
    mx = np.isnan(a).all(**kwargs)
    res = np.nansum(a, **kwargs)
    res[mx] = np.nan
    return res


def get_index_array(a_array, b_array):
    '''
    Get index array where each index points to locataion in a_array. The order of index array corresponds to b_array

      e.g.,
      a_array = [2, 4, 1, 8, 3, 10, 5, 9, 7, 6]
      b_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
      result  = [2, 0, 4, 1, 6, 9, 8, 3, 7, 5]

    https://stackoverflow.com/questions/8251541/numpy-for-every-element-in-one-array-find-the-index-in-another-array
    '''
    index = np.argsort(a_array)
    sorted_a_array = a_array[index]
    sorted_index = np.searchsorted(sorted_a_array, b_array)

    yindex = np.take(index, sorted_index, mode="clip")
    mask = a_array[yindex] != b_array

    result = np.ma.array(yindex, mask=mask)

    return result

