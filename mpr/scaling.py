
import os
import warnings
import numpy as np

import xarray as xr

# TODO:

FILL_VALUE = -9999.0

def horizontal_weighted_mean(matWgts, matIndex_i, matIndex_j, targetPolyIDs, overlaps, origArrays, pvalue, default=FILL_VALUE):
    """Compute areal weighted generalized mean value of for all output polygons
       """
    
    # numpy broadcasting rule
    # https://numpy.org/doc/stable/user/basics.broadcasting.html
    
    # ogirinal grid: 2D [lat, lon] -       > target grid: 1D [hru]
    # ogirinal grid: 3D [Month, lat, lon] -> target grid: 2D [Month, hru]
    #                   [lyr, lat, lon] -> target grid: 2D [lyr, hru]
    
    array_shape = origArrays.shape
    nDims       = len(array_shape)
    
    # TODO
    # move these out of function
    nOutHRUs    = len(overlaps)
    maxOverlaps = overlaps.max()
    
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
        if overlaps[iHru]>0:
            if nDims == 2:
                matDataVals[iHru, 0:overlaps[iHru]] = \
                    origArrays[matIndex_j[iHru, 0:overlaps[iHru]], matIndex_i[iHru, 0:overlaps[iHru]] ]
            elif nDims == 3:
                matDataVals[:, iHru, 0:overlaps[iHru]] = \
                    origArrays[:, matIndex_j[iHru, 0:overlaps[iHru]], matIndex_i[iHru, 0:overlaps[iHru]] ]
        else:
            if nDims == 2:
                matDataVals[iHru, 0] = default
            elif nDims == 3:
                matDataVals[:, iHru, 0] = default
 
    if abs(pvalue) < 0.00001: # geometric mean
        wgtedVals = exp(np.nansum(log(matDataVals)* matWgts, axis=nDims-1))
    else:
        wgtedVals = np.nansum(matDataVals**pvalue * matWgts, axis=nDims-1) **(1.0/pvalue)   # produces vector of weighted values
        
    return wgtedVals
