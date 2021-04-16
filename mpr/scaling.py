
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

    # Compute mean at all the target hrus
    if abs(pvalue) < 0.00001: # geometric mean
        wgtedVals = exp(np.nansum(log(matDataVals)* mapping_data['weight'], axis=nDims-1))
    else:
        wgtedVals = np.nansum(matDataVals**pvalue * mapping_data['weight'], axis=nDims-1) **(1.0/pvalue)   # produces vector of weighted values

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
        matDataVals = np.zeros((array_shape[0], array_shape[1], nMlyr, maxOverlaps), dtype='float32')
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
        wgtedVals = exp(np.nansum(log(matDataVals)* weight, axis=nDims))
    else:
        wgtedVals = np.nansum(matDataVals**pvalue * weight, axis=nDims) **(1.0/pvalue)   # produces vector of weighted values

    return np.moveaxis(wgtedVals, -1, 0)
