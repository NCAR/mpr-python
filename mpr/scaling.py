
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
