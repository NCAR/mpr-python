# Transfer function definitions
# Add or experiment transfer fucntion as needed
#
# ---- convention
# def <transfer_function_name>(attr_data, param_data, param_cfg, coef):
#
#   <transfer_function_name> should match with ['param']['<param_name>']['tf'] in param_meta.yml
#
#   -- arguments
#   attr_data:  dictionary, geophysical data (native resolution)
#   param_data: dictionary, parameter data (native resolution)
#   param_cfg:  dictionary, parameter config (see param_meta.yml)
#   coef:       list,       transfer function coefficient (see tf.yml)

import numpy as np
import mpr.constant as const


class transfer_function():

    def k_macropore_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def k_soil_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = 10.0**(-6.0 + 0.0126* attr_data['sand_pct'] - 0.0064* attr_data['clay_pct'])*const.INCH2M/const*HR2SEC
        return coef[0]*a

    def theta_sat_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = 0.788 + 0.001*attr_data['clay_pct']- 0.263*attr_data['bulk_density']*const.KGCM2GCM #+ coef[1] * attr_data['precip']
        return coef[0]*a

    def qSurfScale_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def aquiferBaseflowExp_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def aquiferBaseflowRate_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def Fcapil_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def summerLAI_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = attr_data['lai'].sel(month=[6,7,8]).mean(dim='month')/10
        return coef[0]*a

    def heightCanopyTop_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        return coef[0]*attr_data['ch']

    def heightCanopyBottom_tf1(attr_data=None, param_data=None, coef=None, param_cfg=None):
        return coef[0]*param_data['heightCanopyTop']

    def frozenPrecipMultip_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = logistic_fun(attr_data['wind'].sel(month=[12,1,2]).mean(dim='month'), A=1.0, L=param_cfg['max'], x0=coef[0], k=coef[1])
        return a


def logistic_fun(x, L=1, k=1, x0=0, A=0):
    return A + (L-A) / (1 + np.exp(-k*(x-x0)))
