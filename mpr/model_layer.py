
import numpy as np

def comp_layer_weight(soil_thickness, model_thickness):

    #-- Compute for depths to 1)top and 2) bottom of soil and model layer
    
    nMlyr = len(model_thickness)
    
    soil_bot = np.cumsum(soil_thickness)
    model_bot = np.cumsum(model_thickness)

    model_top = shift(model_bot, 1, fill_value=0)
    soil_top  = shift(soil_bot, 1, fill_value=0)

    #-- Compute index of soil layer where top of model layer is located
    idxTop = np.zeros(len(model_thickness), dtype='int32')
    idxBot = np.zeros(len(model_thickness), dtype='int32')

    for ixm in range(nMlyr):
        for ixs in range(len(soil_thickness)):
            if model_top[ixm] < soil_bot[ixs]:
                idxTop[ixm] = ixs
                break

    #-- Compute index of soil layer where the bottom of model layer is located
    for ixm in range(nMlyr):
        for ixs in range(len(soil_thickness)):
            if model_bot[ixm] < soil_bot[ixs]:
                idxBot[ixm] = ixs
                break

    maxSlyr = (idxBot-idxTop+1).max()
    mapping_data = {'overlaps':   np.zeros(nMlyr, dtype='int32'),
                    'soil_index': np.zeros((nMlyr, maxSlyr), dtype='int32'),
                    'weight' :    np.zeros((nMlyr, maxSlyr), dtype='float32'),
                   }

    for ixm in range(len(model_thickness)):

        if idxTop[ixm] == idxBot[ixm]: # if model layer is completely within soil layer
            mapping_data['overlaps'][ixm]      = 1
            mapping_data['soil_index'][ixm, 0] = ixs
            mapping_data['weight'][ixm, 0]     = 1.0
            continue

        # loop frm the upper most soil layer to the lowest soil layer that intersect current model layer
        counter = 0
        for ixs in range(idxTop[ixm], idxBot[ixm]+1):
            # if model layer contains multiple soil layers
            if ixs == idxTop[ixm]:           # for the upper most soil layer that intersect model layer
                mapping_data['weight'][ixm, counter] = (soil_bot[ixs] - model_top[ixm])/model_thickness[ixm]
            elif ixs == idxBot[ixm]:         # for the lowest soil layer that intersect model layer
                mapping_data['weight'][ixm, counter] = (model_bot[ixm] - soil_top[ixs])/model_thickness[ixm]
            else:                            # for soil layers that completely in model layer
                mapping_data['weight'][ixm, counter] = soil_thickness[ixs]/model_thickness[ixm]
            mapping_data['soil_index'][ixm, counter] = ixs
            counter += 1
        mapping_data['overlaps'][ixm]      = counter
        
    return mapping_data


def shift(arr, num, fill_value=np.nan):
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result
