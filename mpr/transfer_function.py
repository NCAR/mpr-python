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

    def norm_prec_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # normalized precipitation
        a = attr_data['prec'].mean(dim='month')
        b = a/a.mean()
        return b

    def retention_slope_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # cosby et al., 1984
        a = 3.1 + 0.157*attr_data['clay_pct']- 0.003*attr_data['sand_pct']
        return coef[0]*a

    def matric_potential_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # cosby et al., 1984. unit is kPa
        a = -1* (10.0**(1.54- 0.0095* attr_data['sand_pct']+ 0.0063* attr_data['silt_pct']))*const.cmH2O2kPa # cosby equation give cm-H2O
        return coef[0]*a

    def theta_sat_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # Zacharias & Wessolek 200
        a = coef[0]*(0.788 + 0.001*attr_data['clay_pct']- 0.263*attr_data['bulk_density']*const.KGCM2GCCM) #+ coef[1] * attr_data['precip']
        b = logistic_fun(a, A=param_cfg['min'], L=param_cfg['max'], x0=(param_cfg['max']-param_cfg['min'])/2, k=12.5)
        return b

    def theta_sat_tf2(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # cosby et al., 1984
        a = 50.5 - 0.142*attr_data['sand_pct']- 0.037*attr_data['clay_pct']
        return coef[0]*a/100.0

    def fieldCapacity_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        psi_fc=-10.0 #kPa
        a = param_data['theta_sat']*(psi_fc/param_data['matric_potential'])**(-1.0/param_data['retention_slope'])
        return coef[0]*a

    def critSoilWilting_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        psi_wp = -1500.0 #kPa
        a = param_data['theta_sat']*(psi_wp/param_data['matric_potential'])**(-1.0/param_data['retention_slope'])
        return coef[0]*a

    def critSoilTranspire_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        return coef[0]*param_data['theta_sat'] + (1.0 - coef[0])* param_data['critSoilWilting']

    def theta_res_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        return coef[0]*param_data['critSoilWilting']

    def k_soil_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # cosby et al., 1984 given in inch/hr
        a = coef[0]*10.0**(-0.6 + 0.0126* attr_data['sand_pct'] - 0.0064* attr_data['clay_pct'])*const.INCH2M/const.HR2SEC
        b = logistic_fun(param_data['norm_prec'], A=0.5, L=1.5, x0=1.0, k=coef[1])
        return a*b

    def k_macropore_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = logistic_fun(param_data['norm_prec'], A=0.0005, L=0.09, x0=1.0, k=coef[0])
        return a

    def qSurfScale_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = logistic_fun(attr_data['slope_mean'], A=param_cfg['max'], L=param_cfg['min'], x0=coef[0], k=coef[1])
        return a

    def aquiferBaseflowExp_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = logistic_fun(attr_data['slope_mean'], A=1.0, L=param_cfg['max'], x0=coef[0], k=coef[1])
        return a

    def aquiferBaseflowRate_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        return  (10**coef[0])*param_data['k_soil']

    def Fcapil_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = param_data['qSurfScale'].copy()
        a = a*0.0 + 1.0
        return coef[0]*a

    def summerLAI_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = attr_data['lai'].sel(month=[6,7,8]).mean(dim='month')/10
        return coef[0]*a

    def heightCanopyTop_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = attr_data['ch'].where(attr_data['ch']>0, 0.1)
        return coef[0]*a

    def heightCanopyBottom_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        return coef[0]*param_data['heightCanopyTop']

    def frozenPrecipMultip_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        a = logistic_fun(attr_data['wind'].sel(month=[12,1,2]).mean(dim='month'), A=1.0, L=param_cfg['max'], x0=coef[0], k=coef[1])
        return a

    def vGn_alpha_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # Vereecken et al., 1989: Estimating the soil moisture retention characteristic from texture, bulk density, and carbon content. Soil Sci. 148
        a = -1.0* np.exp(-2.486 + 0.025* attr_data['sand_pct'] - 0.351* attr_data['soc']*0.1 -2.617* attr_data['bulk_density']*const.KGCM2GCCM -
                         0.023* attr_data['clay_pct'])/const.CM2M
        return coef[0]*a

    def vGn_n_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # Vereecken et al., 1989: Estimating the soil moisture retention characteristic from texture, bulk density, and carbon content. Soil Sci. 148
        a = coef[0]* np.exp(0.053 + 0.009* attr_data['sand_pct'] - 0.013* attr_data['clay_pct'] + 0.00015* attr_data['sand_pct']**2)
        a = a.where(a>=1.0, 1.01)
        return a

    def wettingFrontSuction_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        # Rawls et al., Green‐ampt Infiltration Parameters from Soils Data 1982.   or Clapp and Horgberger 1978
        # parameter is pore-size distribution index 0.1 - 3.0
        a = (2.0*param_data['retention_slope'] + 3.0)/(param_data['retention_slope'] + 3.0)*param_data['matric_potential'] * const.kPa2mH2O * (-1.0)
        return coef[0]*a

    def zScale_TOPMODEL_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def routingGammaScale_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass

    def routingGammaShape_tf1(attr_data=None, param_data=None, param_cfg=None, coef=None):
        pass


def logistic_fun(x, L=1, k=1, x0=0, A=0):
    return A + (L-A) / (1 + np.exp(-k*(x-x0)))

# additional van Genuchten parameters
# Zacharias, S. and Wessolek, G. (2007), Excluding Organic Matter Content from Pedotransfer Predictors of Soil Water Retention. Soil Sci. Soc. Am. J., 71: 43-50. https://doi.org/10.2136/sssaj2006.0098
