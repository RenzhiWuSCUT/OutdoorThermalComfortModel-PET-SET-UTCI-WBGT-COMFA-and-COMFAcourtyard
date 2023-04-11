import math
from typing import Union

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module Name: OutdoorThermalIndex
Module Description: 
    This module contains code for common outdoor thermal comfort indices (OCTs), including:
    - PET (Physiological Equivalent Temperature),
    - SET* (Standard Equivalent Temperature),
    - UTCI (Universal Thermal Climate Index),
    - WBGToutdoor (Outdoor Wet Bulb Globe Temperature),
    - COMFA,
    - COMFAcourtayrd,
    These indices are widely used to evaluate the thermal environment and predict
    the thermal sensation of a person in outdoor spaces. The code implements various
    formulas and algorithms to calculate the OCT values based on the inputs of weather
    conditions, clothing, and human physiological parameters.
Author: Renzhi WU
Email: renzhiwuscut@gmail.com
Date: Feb 7, 2023
"""


# Reference:
# [1].	Wu, R., et al., The COMFA model for assessing courtyard thermal comfort in hot and humid regions: A comparative study with existing models. Building and Environment, 2023. 234: p. 110150.
class OCT:

    # Reference:
    # [1] K. Blazejczyk, J. Baranowski, A. Blazejczyk, Heat stress and occupational health and safety - Spatial and temporal differentiation, Misc. Geogr. 18 (2014) 61–67. https://doi.org/10.2478/mgrsd-2014-0011.
    # [2] Chriswmackey, ladybug-legacy, (2018). https://github.com/ladybug-tools/ladybug-legacy.
    @staticmethod
    def cal_WBGToutdoor(
            tAir_C: float = 18.081,
            tMrt_C: float = 28.94,
            vAir: float = 2.13,
            rh: float = 43.912) -> float:

        """
        Parameters 参数:
        tAir_C : float : Air temperature [℃]. 空气温度[℃]。
        tMrt_C : float : Mean radiant temperature [℃]. 平均辐射温度[℃]。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。

        Returns 返回值:
        float : Wet Bulb Globe Temperature for outdoor conditions [℃]. 户外湿球温度[℃]。
        """

        coefL = [-2836.5744, -6028.076559, 19.54263612, -0.02737830188, 0.000016261698, (7.0229056e-10),
                 (-1.8680009e-13)]
        tAir_K = tAir_C + 273.15
        qAir_hPa = 2.7150305 * math.log(tAir_K)
        for iCoef, coef in enumerate(coefL):
            qAir_hPa = qAir_hPa + (coef * (tAir_K ** (iCoef - 2)))
        qAir_hPa = math.exp(qAir_hPa) * rh * 0.01 / 100
        tWet_C = 1.885 + 0.3704 * tAir_C + 0.4492 * qAir_hPa
        tGlobe_C = 2.098 - 2.561 * vAir + 0.5957 * tAir_C + 0.4017 * tMrt_C
        WBGTout = 0.7 * tWet_C + 0.2 * tGlobe_C + 0.1 * tAir_C
        return WBGTout

    # Reference:
    # [1] D. Fiala, Dynamic Simulation of Human Heat Transfer and Thermal Comfort, Sustain. Dev. 45 (1998) 1. https://www.dora.dmu.ac.uk/xmlui/handle/2086/4129.
    # [2] D. Fiala, K.J. Lomas, M. Stohrer, A computer model of human thermoregulation for a wide range of environmental conditions: The passive system, J. Appl. Physiol. 87 (1999) 1957–1972. https://doi.org/10.1152/jappl.1999.87.5.1957.
    # [3] D. Fiala, G. Havenith, P. Bröde, B. Kampmann, G. Jendritzky, UTCI-Fiala multi-node model of human heat transfer and temperature regulation, Int. J. Biometeorol. 56 (2012) 429–441. https://doi.org/10.1007/s00484-011-0424-7.
    # [4] ISB Commission 6, C.A. 730, UTCI Universal Thermal Climate Index Documents, (n.d.). http://www.utci.org/utci_doku.php.
    # [5] B. P, F. D, B. K, H. I, J. G, K. B, T. B, H. G, Deriving the operational procedure for the Universal Thermal Climate Index (UTCI)., Int. J. Biometeorol. 56 (2012) 481–494. https://pubmed.ncbi.nlm.nih.gov/21626294/.
    # [6] Chriswmackey, ladybug-legacy, (2018). https://github.com/ladybug-tools/ladybug-legacy.
    @staticmethod
    def cal_UTCI(tAir_C: float = 18.081,
                 tMrt_C: float = 28.94,
                 vAir: float = 2.13,
                 rh: float = 43.912) -> float:
        """
        Parameters 参数:
        tAir_C : float : Air temperature [℃]. 空气温度[℃]。
        tMrt_C : float : Mean radiant temperature [℃]. 平均辐射温度[℃]。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。

        Returns 返回值:
        float : Universal Thermal Climate Index [℃]. 通用热气候指数UTCI [℃]。
        """

        def cal_svp(tAir_C):
            coefL = [-2836.5744, -6028.076559, 19.54263612, -0.02737830188, 0.000016261698, (7.0229056 * (10 ** (-10))),
                     (-1.8680009 * (10 ** (-13)))]
            tAir_K = tAir_C + 273.15
            qAirSaturated_hPa = 2.7150305 * math.log(tAir_K)
            for iCoef, coef in enumerate(coefL):
                qAirSaturated_hPa = qAirSaturated_hPa + (coef * (tAir_K ** (iCoef - 2)))
            qAirSaturated_hPa = math.exp(qAirSaturated_hPa) * 0.01
            return qAirSaturated_hPa

        vAir = min(max(vAir, 0.5), 17)
        qAir_hPa = cal_svp(tAir_C) * (rh / 100.0)
        d_mrt = tMrt_C - tAir_C
        qAir_kPa = qAir_hPa / 10.0

        UTCI_approx = tAir_C + \
                      (0.607562052) + \
                      (-0.0227712343) * tAir_C + \
                      (8.06470249 * (10 ** (-4))) * tAir_C * tAir_C + \
                      (-1.54271372 * (10 ** (-4))) * tAir_C * tAir_C * tAir_C + \
                      (-3.24651735 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * tAir_C + \
                      (7.32602852 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C + \
                      (1.35959073 * (10 ** (-9))) * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C + \
                      (-2.25836520) * vAir + \
                      (0.0880326035) * tAir_C * vAir + \
                      (0.00216844454) * tAir_C * tAir_C * vAir + \
                      (-1.53347087 * (10 ** (-5))) * tAir_C * tAir_C * tAir_C * vAir + \
                      (-5.72983704 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * tAir_C * vAir + \
                      (-2.55090145 * (10 ** (-9))) * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C * vAir + \
                      (-0.751269505) * vAir * vAir + \
                      (-0.00408350271) * tAir_C * vAir * vAir + \
                      (-5.21670675 * (10 ** (-5))) * tAir_C * tAir_C * vAir * vAir + \
                      (1.94544667 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * vAir * vAir + \
                      (1.14099531 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * tAir_C * vAir * vAir + \
                      (0.158137256) * vAir * vAir * vAir + \
                      (-6.57263143 * (10 ** (-5))) * tAir_C * vAir * vAir * vAir + \
                      (2.22697524 * (10 ** (-7))) * tAir_C * tAir_C * vAir * vAir * vAir + \
                      (-4.16117031 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * vAir * vAir * vAir + \
                      (-0.0127762753) * vAir * vAir * vAir * vAir + \
                      (9.66891875 * (10 ** (-6))) * tAir_C * vAir * vAir * vAir * vAir + \
                      (2.52785852 * (10 ** (-9))) * tAir_C * tAir_C * vAir * vAir * vAir * vAir + \
                      (4.56306672 * (10 ** (-4))) * vAir * vAir * vAir * vAir * vAir + \
                      (-1.74202546 * (10 ** (-7))) * tAir_C * vAir * vAir * vAir * vAir * vAir + \
                      (-5.91491269 * (10 ** (-6))) * vAir * vAir * vAir * vAir * vAir * vAir + \
                      (0.398374029) * d_mrt + \
                      (1.83945314 * (10 ** (-4))) * tAir_C * d_mrt + \
                      (-1.73754510 * (10 ** (-4))) * tAir_C * tAir_C * d_mrt + \
                      (-7.60781159 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * d_mrt + \
                      (3.77830287 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * tAir_C * d_mrt + \
                      (5.43079673 * (10 ** (-10))) * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C * d_mrt + \
                      (-0.0200518269) * vAir * d_mrt + \
                      (8.92859837 * (10 ** (-4))) * tAir_C * vAir * d_mrt + \
                      (3.45433048 * (10 ** (-6))) * tAir_C * tAir_C * vAir * d_mrt + \
                      (-3.77925774 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * vAir * d_mrt + \
                      (-1.69699377 * (10 ** (-9))) * tAir_C * tAir_C * tAir_C * tAir_C * vAir * d_mrt + \
                      (1.69992415 * (10 ** (-4))) * vAir * vAir * d_mrt + \
                      (-4.99204314 * (10 ** (-5))) * tAir_C * vAir * vAir * d_mrt + \
                      (2.47417178 * (10 ** (-7))) * tAir_C * tAir_C * vAir * vAir * d_mrt + \
                      (1.07596466 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * vAir * vAir * d_mrt + \
                      (8.49242932 * (10 ** (-5))) * vAir * vAir * vAir * d_mrt + \
                      (1.35191328 * (10 ** (-6))) * tAir_C * vAir * vAir * vAir * d_mrt + \
                      (-6.21531254 * (10 ** (-9))) * tAir_C * tAir_C * vAir * vAir * vAir * d_mrt + \
                      (-4.99410301 * (10 ** (-6))) * vAir * vAir * vAir * vAir * d_mrt + \
                      (-1.89489258 * (10 ** (-8))) * tAir_C * vAir * vAir * vAir * vAir * d_mrt + \
                      (8.15300114 * (10 ** (-8))) * vAir * vAir * vAir * vAir * vAir * d_mrt + \
                      (7.55043090 * (10 ** (-4))) * d_mrt * d_mrt + \
                      (-5.65095215 * (10 ** (-5))) * tAir_C * d_mrt * d_mrt + \
                      (-4.52166564 * (10 ** (-7))) * tAir_C * tAir_C * d_mrt * d_mrt + \
                      (2.46688878 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * d_mrt * d_mrt + \
                      (2.42674348 * (10 ** (-10))) * tAir_C * tAir_C * tAir_C * tAir_C * d_mrt * d_mrt + \
                      (1.54547250 * (10 ** (-4))) * vAir * d_mrt * d_mrt + \
                      (5.24110970 * (10 ** (-6))) * tAir_C * vAir * d_mrt * d_mrt + \
                      (-8.75874982 * (10 ** (-8))) * tAir_C * tAir_C * vAir * d_mrt * d_mrt + \
                      (-1.50743064 * (10 ** (-9))) * tAir_C * tAir_C * tAir_C * vAir * d_mrt * d_mrt + \
                      (-1.56236307 * (10 ** (-5))) * vAir * vAir * d_mrt * d_mrt + \
                      (-1.33895614 * (10 ** (-7))) * tAir_C * vAir * vAir * d_mrt * d_mrt + \
                      (2.49709824 * (10 ** (-9))) * tAir_C * tAir_C * vAir * vAir * d_mrt * d_mrt + \
                      (6.51711721 * (10 ** (-7))) * vAir * vAir * vAir * d_mrt * d_mrt + \
                      (1.94960053 * (10 ** (-9))) * tAir_C * vAir * vAir * vAir * d_mrt * d_mrt + \
                      (-1.00361113 * (10 ** (-8))) * vAir * vAir * vAir * vAir * d_mrt * d_mrt + \
                      (-1.21206673 * (10 ** (-5))) * d_mrt * d_mrt * d_mrt + \
                      (-2.18203660 * (10 ** (-7))) * tAir_C * d_mrt * d_mrt * d_mrt + \
                      (7.51269482 * (10 ** (-9))) * tAir_C * tAir_C * d_mrt * d_mrt * d_mrt + \
                      (9.79063848 * (10 ** (-11))) * tAir_C * tAir_C * tAir_C * d_mrt * d_mrt * d_mrt + \
                      (1.25006734 * (10 ** (-6))) * vAir * d_mrt * d_mrt * d_mrt + \
                      (-1.81584736 * (10 ** (-9))) * tAir_C * vAir * d_mrt * d_mrt * d_mrt + \
                      (-3.52197671 * (10 ** (-10))) * tAir_C * tAir_C * vAir * d_mrt * d_mrt * d_mrt + \
                      (-3.36514630 * (10 ** (-8))) * vAir * vAir * d_mrt * d_mrt * d_mrt + \
                      (1.35908359 * (10 ** (-10))) * tAir_C * vAir * vAir * d_mrt * d_mrt * d_mrt + \
                      (4.17032620 * (10 ** (-10))) * vAir * vAir * vAir * d_mrt * d_mrt * d_mrt + \
                      (-1.30369025 * (10 ** (-9))) * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (4.13908461 * (10 ** (-10))) * tAir_C * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (9.22652254 * (10 ** (-12))) * tAir_C * tAir_C * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (-5.08220384 * (10 ** (-9))) * vAir * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (-2.24730961 * (10 ** (-11))) * tAir_C * vAir * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (1.17139133 * (10 ** (-10))) * vAir * vAir * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (6.62154879 * (10 ** (-10))) * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (4.03863260 * (10 ** (-13))) * tAir_C * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (1.95087203 * (10 ** (-12))) * vAir * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (-4.73602469 * (10 ** (-12))) * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt + \
                      (5.12733497) * qAir_kPa + \
                      (-0.312788561) * tAir_C * qAir_kPa + \
                      (-0.0196701861) * tAir_C * tAir_C * qAir_kPa + \
                      (9.99690870 * (10 ** (-4))) * tAir_C * tAir_C * tAir_C * qAir_kPa + \
                      (9.51738512 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * tAir_C * qAir_kPa + \
                      (-4.66426341 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * tAir_C * tAir_C * qAir_kPa + \
                      (0.548050612) * vAir * qAir_kPa + \
                      (-0.00330552823) * tAir_C * vAir * qAir_kPa + \
                      (-0.00164119440) * tAir_C * tAir_C * vAir * qAir_kPa + \
                      (-5.16670694 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * vAir * qAir_kPa + \
                      (9.52692432 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * tAir_C * vAir * qAir_kPa + \
                      (-0.0429223622) * vAir * vAir * qAir_kPa + \
                      (0.00500845667) * tAir_C * vAir * vAir * qAir_kPa + \
                      (1.00601257 * (10 ** (-6))) * tAir_C * tAir_C * vAir * vAir * qAir_kPa + \
                      (-1.81748644 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * vAir * vAir * qAir_kPa + \
                      (-1.25813502 * (10 ** (-3))) * vAir * vAir * vAir * qAir_kPa + \
                      (-1.79330391 * (10 ** (-4))) * tAir_C * vAir * vAir * vAir * qAir_kPa + \
                      (2.34994441 * (10 ** (-6))) * tAir_C * tAir_C * vAir * vAir * vAir * qAir_kPa + \
                      (1.29735808 * (10 ** (-4))) * vAir * vAir * vAir * vAir * qAir_kPa + \
                      (1.29064870 * (10 ** (-6))) * tAir_C * vAir * vAir * vAir * vAir * qAir_kPa + \
                      (-2.28558686 * (10 ** (-6))) * vAir * vAir * vAir * vAir * vAir * qAir_kPa + \
                      (-0.0369476348) * d_mrt * qAir_kPa + \
                      (0.00162325322) * tAir_C * d_mrt * qAir_kPa + \
                      (-3.14279680 * (10 ** (-5))) * tAir_C * tAir_C * d_mrt * qAir_kPa + \
                      (2.59835559 * (10 ** (-6))) * tAir_C * tAir_C * tAir_C * d_mrt * qAir_kPa + \
                      (-4.77136523 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * tAir_C * d_mrt * qAir_kPa + \
                      (8.64203390 * (10 ** (-3))) * vAir * d_mrt * qAir_kPa + \
                      (-6.87405181 * (10 ** (-4))) * tAir_C * vAir * d_mrt * qAir_kPa + \
                      (-9.13863872 * (10 ** (-6))) * tAir_C * tAir_C * vAir * d_mrt * qAir_kPa + \
                      (5.15916806 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * vAir * d_mrt * qAir_kPa + \
                      (-3.59217476 * (10 ** (-5))) * vAir * vAir * d_mrt * qAir_kPa + \
                      (3.28696511 * (10 ** (-5))) * tAir_C * vAir * vAir * d_mrt * qAir_kPa + \
                      (-7.10542454 * (10 ** (-7))) * tAir_C * tAir_C * vAir * vAir * d_mrt * qAir_kPa + \
                      (-1.24382300 * (10 ** (-5))) * vAir * vAir * vAir * d_mrt * qAir_kPa + \
                      (-7.38584400 * (10 ** (-9))) * tAir_C * vAir * vAir * vAir * d_mrt * qAir_kPa + \
                      (2.20609296 * (10 ** (-7))) * vAir * vAir * vAir * vAir * d_mrt * qAir_kPa + \
                      (-7.32469180 * (10 ** (-4))) * d_mrt * d_mrt * qAir_kPa + \
                      (-1.87381964 * (10 ** (-5))) * tAir_C * d_mrt * d_mrt * qAir_kPa + \
                      (4.80925239 * (10 ** (-6))) * tAir_C * tAir_C * d_mrt * d_mrt * qAir_kPa + \
                      (-8.75492040 * (10 ** (-8))) * tAir_C * tAir_C * tAir_C * d_mrt * d_mrt * qAir_kPa + \
                      (2.77862930 * (10 ** (-5))) * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (-5.06004592 * (10 ** (-6))) * tAir_C * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (1.14325367 * (10 ** (-7))) * tAir_C * tAir_C * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (2.53016723 * (10 ** (-6))) * vAir * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (-1.72857035 * (10 ** (-8))) * tAir_C * vAir * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (-3.95079398 * (10 ** (-8))) * vAir * vAir * vAir * d_mrt * d_mrt * qAir_kPa + \
                      (-3.59413173 * (10 ** (-7))) * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (7.04388046 * (10 ** (-7))) * tAir_C * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (-1.89309167 * (10 ** (-8))) * tAir_C * tAir_C * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (-4.79768731 * (10 ** (-7))) * vAir * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (7.96079978 * (10 ** (-9))) * tAir_C * vAir * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (1.62897058 * (10 ** (-9))) * vAir * vAir * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (3.94367674 * (10 ** (-8))) * d_mrt * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (-1.18566247 * (10 ** (-9))) * tAir_C * d_mrt * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (3.34678041 * (10 ** (-10))) * vAir * d_mrt * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (-1.15606447 * (10 ** (-10))) * d_mrt * d_mrt * d_mrt * d_mrt * d_mrt * qAir_kPa + \
                      (-2.80626406) * qAir_kPa * qAir_kPa + \
                      (0.548712484) * tAir_C * qAir_kPa * qAir_kPa + \
                      (-0.00399428410) * tAir_C * tAir_C * qAir_kPa * qAir_kPa + \
                      (-9.54009191 * (10 ** (-4))) * tAir_C * tAir_C * tAir_C * qAir_kPa * qAir_kPa + \
                      (1.93090978 * (10 ** (-5))) * tAir_C * tAir_C * tAir_C * tAir_C * qAir_kPa * qAir_kPa + \
                      (-0.308806365) * vAir * qAir_kPa * qAir_kPa + \
                      (0.0116952364) * tAir_C * vAir * qAir_kPa * qAir_kPa + \
                      (4.95271903 * (10 ** (-4))) * tAir_C * tAir_C * vAir * qAir_kPa * qAir_kPa + \
                      (-1.90710882 * (10 ** (-5))) * tAir_C * tAir_C * tAir_C * vAir * qAir_kPa * qAir_kPa + \
                      (0.00210787756) * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (-6.98445738 * (10 ** (-4))) * tAir_C * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (2.30109073 * (10 ** (-5))) * tAir_C * tAir_C * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (4.17856590 * (10 ** (-4))) * vAir * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (-1.27043871 * (10 ** (-5))) * tAir_C * vAir * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (-3.04620472 * (10 ** (-6))) * vAir * vAir * vAir * vAir * qAir_kPa * qAir_kPa + \
                      (0.0514507424) * d_mrt * qAir_kPa * qAir_kPa + \
                      (-0.00432510997) * tAir_C * d_mrt * qAir_kPa * qAir_kPa + \
                      (8.99281156 * (10 ** (-5))) * tAir_C * tAir_C * d_mrt * qAir_kPa * qAir_kPa + \
                      (-7.14663943 * (10 ** (-7))) * tAir_C * tAir_C * tAir_C * d_mrt * qAir_kPa * qAir_kPa + \
                      (-2.66016305 * (10 ** (-4))) * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (2.63789586 * (10 ** (-4))) * tAir_C * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (-7.01199003 * (10 ** (-6))) * tAir_C * tAir_C * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (-1.06823306 * (10 ** (-4))) * vAir * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (3.61341136 * (10 ** (-6))) * tAir_C * vAir * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (2.29748967 * (10 ** (-7))) * vAir * vAir * vAir * d_mrt * qAir_kPa * qAir_kPa + \
                      (3.04788893 * (10 ** (-4))) * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (-6.42070836 * (10 ** (-5))) * tAir_C * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (1.16257971 * (10 ** (-6))) * tAir_C * tAir_C * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (7.68023384 * (10 ** (-6))) * vAir * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (-5.47446896 * (10 ** (-7))) * tAir_C * vAir * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (-3.59937910 * (10 ** (-8))) * vAir * vAir * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (-4.36497725 * (10 ** (-6))) * d_mrt * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (1.68737969 * (10 ** (-7))) * tAir_C * d_mrt * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (2.67489271 * (10 ** (-8))) * vAir * d_mrt * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (3.23926897 * (10 ** (-9))) * d_mrt * d_mrt * d_mrt * d_mrt * qAir_kPa * qAir_kPa + \
                      (-0.0353874123) * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.221201190) * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.0155126038) * tAir_C * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-2.63917279 * (10 ** (-4))) * tAir_C * tAir_C * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.0453433455) * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.00432943862) * tAir_C * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (1.45389826 * (10 ** (-4))) * tAir_C * tAir_C * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (2.17508610 * (10 ** (-4))) * vAir * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-6.66724702 * (10 ** (-5))) * tAir_C * vAir * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (3.33217140 * (10 ** (-5))) * vAir * vAir * vAir * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.00226921615) * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (3.80261982 * (10 ** (-4))) * tAir_C * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-5.45314314 * (10 ** (-9))) * tAir_C * tAir_C * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-7.96355448 * (10 ** (-4))) * vAir * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (2.53458034 * (10 ** (-5))) * tAir_C * vAir * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-6.31223658 * (10 ** (-6))) * vAir * vAir * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (3.02122035 * (10 ** (-4))) * d_mrt * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-4.77403547 * (10 ** (-6))) * tAir_C * d_mrt * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (1.73825715 * (10 ** (-6))) * vAir * d_mrt * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-4.09087898 * (10 ** (-7))) * d_mrt * d_mrt * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.614155345) * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.0616755931) * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.00133374846) * tAir_C * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.00355375387) * vAir * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-5.13027851 * (10 ** (-4))) * tAir_C * vAir * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (1.02449757 * (10 ** (-4))) * vAir * vAir * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.00148526421) * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-4.11469183 * (10 ** (-5))) * tAir_C * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-6.80434415 * (10 ** (-6))) * vAir * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-9.77675906 * (10 ** (-6))) * d_mrt * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.0882773108) * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (-0.00301859306) * tAir_C * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.00104452989) * vAir * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (2.47090539 * (10 ** (-4))) * d_mrt * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa + \
                      (0.00148348065) * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa * qAir_kPa

        return UTCI_approx

    # Reference:
    # [1] J.A.. Stolwijk, A Mathematical Model of Physiological Temperature Regulation in Man, NASA Contract. Reports. 303 (1971) 14–30.
    # [2] GAGGE AP, STOLWIJK JAJ, NISHI Y, Effective temperature scale based on a simple model of human physiological regulatory response, ASHRAE Trans. 77 (1971) 247–263.
    # [3] A.P. Gagge, A.P. Fobelets, L.G. Berglund, Standard Predictive Index of Human Response To the Thermal Environment., ASHRAE Trans. 92 (1986) 709–731.
    # [4] C. Fountain, Marc;Huizenga, A Thermal Sensation Model For Use By The Engineering Profession, (1995). https://escholarship.org/uc/item/89d5c8k7.
    # [5] Building simulation resources: SET_star, (2002). http://news-sv.aij.or.jp/kankyo/s12/Resource/ap/SET_star/SET_star.htm.
    @staticmethod
    def cal_SETstar(tAir_C: float = 18.081,
                    tMrt_C: float = 28.94,
                    vAir: float = 2.13,
                    rh: float = 43.912,
                    height_m: float = 1.596,
                    weight_kg: float = 55.3,
                    Q_movement: float = 58.15,
                    rCloClo: float = 0.7) -> float:
        """
        Parameters 参数:
        tAir_C : float : Air temperature [℃]. 空气温度[℃]。
        tMrt_C : float : Mean radiant temperature [℃]. 平均辐射温度[℃]。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。
        height_m : float : Height of the person [m]. 身高[m]。
        weight_kg : float : Weight of the person [kg]. 体重[kg]。
        Q_movement : float : Metabolic heat production [W/m2]. 新陈代谢产热[W/m2]。
        rCloClo : float : Clothing insulation [clo]. 服装热阻[clo]。

        Returns 返回值:
        set_star: float : Standard Effective Temperature (star) [℃]. 标准等效温度SET* [℃]。
        et_star: float : Standard Effective Temperature (star) [℃]. 标准等效温度SET* [℃]。
        """

        def cal_svp(tAir_C):
            return math.exp((186686. - 40301830. / (tAir_C + 235.)) / 10000.)

        pBaro_kPa = 101.325
        pAir_MMHG = 0

        areaSrf = 71.84e-4 * ((height_m * 100) ** 0.725) * (weight_kg ** 0.425)
        pBaro_MMHG = pBaro_kPa / 101.325 * 760
        pAirSaturated_MMHG = cal_svp(tAir_C)
        if pAir_MMHG == 0: pAir_MMHG = rh * pAirSaturated_MMHG / 100.

        # Purpose【常数设定】
        SBC = 5.67 * 10 ** (-8)
        fCloAreaCoef = 0.25
        Q_work = 0
        coefSweat = 170
        coefVasodilation = 200
        coefVasoconstriction = 0.5

        # Purpose【初始值设定】
        tSkin_C_setting = 33.7
        tCore_C_setting = 36.8
        tBody_C_setting = 36.49
        flowBloodSkin_setting = 6.3
        tSkin_C = tSkin_C_setting
        tCore_C = tCore_C_setting
        flowBloodSkin = flowBloodSkin_setting
        Q_shiver = 0
        alpha_RatioMass_skin2core = 0.1

        # Purpose：Unit conversion
        Q_latent = 0.1 * Q_movement
        pBaro_atm = pBaro_MMHG / 760
        rClo = 0.155 * rCloClo
        fCloArea = 1 + 0.15 * rCloClo
        coefLewis = 2.2 / pBaro_atm
        Q_meta = Q_movement
        if vAir == 0: vAir = 0.0000

        # Purpose：wet ratio upper limit
        if rCloClo <= 0.:
            ratioWetMax = 0.38 * (vAir ** (-0.29))
            coefCloWetInherent = 1
        else:
            ratioWetMax = 0.59 * (vAir ** (-0.08))
            coefCloWetInherent = 0.45

        # Purpose：convective heat transfer coefficient
        hConv = 3 * (pBaro_atm ** (0.53))
        if ((Q_movement / 58.2) < 0.85):
            hConvActivity = 0
        else:
            hConvActivity = 5.66 * ((((Q_movement / 58.2) - 0.85) * pBaro_atm) ** 0.39)
        hConvVelocity = 8.600001 * ((vAir * pBaro_atm) ** 0.53)
        if (hConv <= hConvActivity):    hConv = hConvActivity
        if (hConv < hConvVelocity):    hConv = hConvVelocity

        # Purpose：initial clothing surface temperature
        hRadLinear = 4.7
        hTotal = hRadLinear + hConv
        rAir = 1 / (fCloArea * hTotal)
        tOperation_C = (hRadLinear * tMrt_C + hConv * tAir_C) / hTotal
        tClo_C = tOperation_C + (tSkin_C - tOperation_C) / (hTotal * (rAir + rClo))

        # 目的：衣服表面温度和线性辐射换热系数联立求解
        # Purpose：solving clothing surface temperature and linear radiation heat transfer coefficient
        convergeFlag = False
        while (convergeFlag == False):
            tClo_C_thisStep = tClo_C
            hRadLinear = 4 * SBC * (((tClo_C + tMrt_C) / 2 + 273.15) ** 3) * 0.72
            hTotal = hRadLinear + hConv
            rAir = 1 / (fCloArea * hTotal)
            tOperation_C = (hRadLinear * tMrt_C + hConv * tAir_C) / hTotal
            tClo_C = (rAir * tSkin_C + rClo * tOperation_C) / (rAir + rClo)
            if (abs(tClo_C - tClo_C_thisStep) > 0.01):
                convergeFlag = False
            else:
                convergeFlag = True

        # 目的：计算步长为一秒钟
        # Purpose: Calculate the step size as one second
        secondsOfOneMin60 = 60
        for thisSecond in range(secondsOfOneMin60):

            # 目的：重新计算衣服表面温度
            # Purpose: Recalculate the clothing surface temperature
            while (abs(tClo_C - tClo_C_thisStep) > 0.01):
                tClo_C_thisStep = tClo_C
                hRadLinear = 4 * SBC * (((tClo_C + tMrt_C) / 2 + 273.15) ** 3) * 0.72
                hTotal = hRadLinear + hConv
                rAir = 1 / (fCloArea * hTotal)
                tOperation_C = (hRadLinear * tMrt_C + hConv * tAir_C) / hTotal
                tClo_C = (rAir * tSkin_C + rClo * tOperation_C) / (rAir + rClo)

            Q_sensible = (tSkin_C - tOperation_C) / (rAir + rClo)
            Q_core2skin = (tCore_C - tSkin_C) * (5.28 + 1.163 * flowBloodSkin)
            Q_res_sensible = 0.0023 * Q_meta * (44 - pAir_MMHG)
            Q_res_latent = 0.0014 * Q_meta * (34 - tAir_C)

            Q_balance_core = Q_meta - Q_core2skin - Q_res_sensible - Q_res_latent - Q_work
            Q_balance_skin = Q_core2skin - Q_sensible - Q_latent

            cp_skin = 0.97 * alpha_RatioMass_skin2core * weight_kg
            cp_core = 0.97 * (1 - alpha_RatioMass_skin2core) * weight_kg

            tSkin_delta = ((Q_balance_skin * areaSrf) / cp_skin) / 60
            tCore_delta = ((Q_balance_core * areaSrf) / cp_core) / 60

            oneSecond1 = 1
            tSkin_C = tSkin_C + tSkin_delta * oneSecond1
            tCore_C = tCore_C + tCore_delta * oneSecond1
            t_body = alpha_RatioMass_skin2core * tSkin_C + (1 - alpha_RatioMass_skin2core) * tCore_C

            # >>>> 目的：冷热信号的计算
            # >>>> Purpose：Calculation of Cold and Heat Signals
            skinSignal = (tSkin_C - tSkin_C_setting)
            warm_skin = skinSignal * (skinSignal > 0)
            cold_skin = -skinSignal * ((-skinSignal) > 0)
            coreSignal = (tCore_C - tCore_C_setting)
            warm_core = coreSignal * (coreSignal > 0)
            cold_core = -coreSignal * ((-coreSignal) > 0)
            bodySignal = (t_body - tBody_C_setting)
            warm_body = bodySignal * (bodySignal > 0)
            cold_body = -bodySignal * ((-bodySignal) > 0)

            # >>>> 目的：皮肤血流量
            # >>>> Purpose: Skin Blood Flow
            flowBloodSkin = (flowBloodSkin_setting + coefVasodilation * warm_core) / (
                    1 + coefVasoconstriction * cold_skin)
            if (flowBloodSkin > 90): flowBloodSkin = 90
            if (flowBloodSkin < 0.5): flowBloodSkin = 0.5

            flowSweatRegulated = coefSweat * warm_body * math.exp(warm_skin / 10.7)
            if (flowSweatRegulated > 500): flowSweatRegulated = 500
            Q_sweat = 0.68 * flowSweatRegulated

            rAirEvap = 1 / (coefLewis * fCloArea * hConv)
            rCloEvap = rClo / (coefLewis * coefCloWetInherent)

            Q_diffusion_max = (cal_svp(tSkin_C) - pAir_MMHG) / (rAirEvap + rCloEvap)
            rateWet_sweat = Q_sweat / Q_diffusion_max
            rateWet = 0.06 + 0.94 * rateWet_sweat
            Q_diffusion = rateWet * Q_diffusion_max - Q_sweat

            Q_latent = Q_sweat + Q_diffusion

            if (rateWet > ratioWetMax):
                rateWet = ratioWetMax
                rateWet_sweat = (ratioWetMax - 0.06) / 0.94
                Q_sweat = rateWet_sweat * Q_diffusion_max
                Q_diffusion = 0.06 * (1 - rateWet_sweat) * Q_diffusion_max
                Q_latent = Q_sweat + Q_diffusion

            if (Q_diffusion_max < 0):
                # Q_diffusion = 0
                # Q_sweat = 0
                rateWet = ratioWetMax
                # rateWet_sweat = ratioWetMax
                Q_latent = Q_diffusion_max

            Q_shiver = 19.4 * cold_skin * cold_core
            Q_meta = Q_movement + Q_shiver

            alpha_RatioMass_skin2core = 0.0417737 + 0.7451833 / (flowBloodSkin + 0.585417)
            if (thisSecond < secondsOfOneMin60):
                tClo_C = (rAir * tSkin_C + rClo * tOperation_C) / (rAir + rClo)

        Q_output_skin = Q_sensible + Q_latent
        hSensible = 1 / (rAir + rClo)
        hLatent = 1 / (rAirEvap + rCloEvap)
        ratioWetting = rateWet
        qSkinSaturated = cal_svp(tSkin_C)

        hRadLinear_SSET = hRadLinear
        if ((Q_movement / 58.2) < 0.85):
            hConv_SSET = 3
        else:
            hConv_SSET = 5.66 * (((Q_movement / 58.2) - 0.85) ** 0.39)
            if (hConv_SSET < 3): hConv_SSET = 3
        hTotal_SSET = hConv_SSET + hRadLinear_SSET

        rClo_clo_SSET = 1.3264 / ((Q_movement - Q_work) / 58.15 + 0.7383) - 0.0953
        rClo_MKnW_SSET = 0.155 * rClo_clo_SSET

        fCloArea_SSET = 1 + fCloAreaCoef * rClo_clo_SSET
        hClo_SSET = 1 / (1 + 0.155 * fCloArea_SSET * hTotal_SSET * rClo_clo_SSET)

        coefIM = 0.45
        coefCloWetInherent_SSET = coefIM * hConv_SSET / hTotal_SSET * (1 - hClo_SSET) / (
                hConv_SSET / hTotal_SSET - hClo_SSET * coefIM)

        rAir_SSET = 1 / (fCloArea_SSET * hTotal_SSET)
        rCloSensible_SSET = 1 / (coefLewis * fCloArea_SSET * hConv_SSET)
        rCloLatent_SSET = rClo_MKnW_SSET / (coefLewis * coefCloWetInherent_SSET)

        hSensible_SSET = 1 / (rAir_SSET + rClo_MKnW_SSET)
        hLatent_SSET = 1 / (rCloSensible_SSET + rCloLatent_SSET)

        delta = 0.0001

        # Purpose【ET*】
        xold = tSkin_C - Q_output_skin / hSensible
        err1 = Q_output_skin - hSensible * (tSkin_C - xold) - ratioWetting * hLatent * (
                qSkinSaturated - 0.5 * cal_svp(xold))
        yold = xold + delta
        err2 = Q_output_skin - hSensible * (tSkin_C - yold) - ratioWetting * hLatent * (
                qSkinSaturated - 0.5 * cal_svp(yold))
        t_evi = xold - delta * err1 / (err2 - err1)
        while (abs(t_evi - xold) > 0.01):
            xold = t_evi
            err1 = Q_output_skin - hSensible * (tSkin_C - xold) - ratioWetting * hLatent * (
                    qSkinSaturated - 0.5 * cal_svp(xold))
            yold = xold + delta
            err2 = Q_output_skin - hSensible * (tSkin_C - yold) - ratioWetting * hLatent * (
                    qSkinSaturated - 0.5 * cal_svp(yold))
            t_evi = xold - delta * err1 / (err2 - err1)
        et_star = t_evi

        # Purpose【SET*】
        xold = tSkin_C - Q_output_skin / hSensible_SSET
        err1 = Q_output_skin - hSensible_SSET * (tSkin_C - xold) - ratioWetting * hLatent_SSET * (
                qSkinSaturated - 0.5 * cal_svp(xold))
        yold = xold + delta
        err2 = Q_output_skin - hSensible_SSET * (tSkin_C - yold) - ratioWetting * hLatent_SSET * (
                qSkinSaturated - 0.5 * cal_svp(yold))
        t_evi = xold - delta * err1 / (err2 - err1)

        while (abs(t_evi - xold) > 0.01):
            xold = t_evi
            err1 = Q_output_skin - hSensible_SSET * (tSkin_C - xold) - ratioWetting * hLatent_SSET * (
                    qSkinSaturated - 0.5 * cal_svp(xold))
            yold = xold + delta
            err2 = Q_output_skin - hSensible_SSET * (tSkin_C - yold) - ratioWetting * hLatent_SSET * (
                    qSkinSaturated - 0.5 * cal_svp(yold))
            t_evi = xold - delta * err1 / (err2 - err1)
        set_star = t_evi
        return set_star, et_star

    # Reference:
    # [1] P. Höppe, Die Energiebilanz des Menschen, Wiss Mitt Meteorol Inst Univ München. 49 (1984) 173.
    # [2] Hoppe PR, Heat balance modelling, Experientia. 49 (1993) 741–746.
    # [3] VDI, Umweltmeteorologie - Methoden zur human-biometeorologischen Bewertung von Klima und Lufthygiene für die Stadt- und Regionalplanung - Teil I: Klima - VDI 3787, (2008) 1–32.
    # [4] Chinese University of Hong Kong department of Architecture, Urban Climatic Map and Standards for Wind Environment - Feasibility Study Technical Input Report No . 1 : Methodologies and Findings of User ’ s Wind Comfort Level Survey Nov 2008, (2008) 151. http://www.pland.gov.hk/pland_en/p_study/prog_s/ucmapweb/ucmap_project/content/reports/Comfort_Level_Survey.pdf.
    # [5] E. Walther, Q. Goestchel, The P.E.T. comfort index: Questioning the model, Build. Environ. 137 (2018) 1–10. https://doi.org/10.1016/j.buildenv.2018.03.054.   @staticmethod
    # [6] Chriswmackey, ladybug-legacy, (2018). https://github.com/ladybug-tools/ladybug-legacy.
    @staticmethod
    def cal_PET(tAir_C: float = 21.0000,
                tMrt_C: float = 22.0000,
                vAir: float = 1.0000,
                rh: float = 50.0000,
                age: int = 35, sex: Union[int, str] = 1,
                height_m: float = 1.7500,
                weight_kg: float = 75.0000,
                bodyPosture: str = 'standing',
                Q_meta: float = 79.152,
                rCloClo: float = 0.9) -> float:

        """
        Parameters 参数:
        tAir_C : float : Air temperature [℃]. 空气温度[℃]。
        tMrt_C : float : Mean radiant temperature [℃]. 平均辐射温度[℃]。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。
        age : int : Age [years]. 年龄[岁]。
        sex : Union[int, str] : Sex of the person. 1 for male, 2 for female, 3 for average sex. 性别，1或'male'为男性，2或 'female'为女性，3或'average sex'为不考虑性别。
        height_m : float : Height of the person [m]. 身高[m]。
        weight_kg : float : Weight of the person [kg]. 体重[kg]。
        bodyPosture : str : Body posture, can be 'standing', 'sitting', 'crouching'. 身体姿势，取值范围为'standing', 'sitting',或 'crouching'。
        Q_meta : float : Metabolic heat production [W/m2]. 新陈代谢产热[W/m2]。
        rCloClo : float : Clothing insulation [clo]. 服装热阻[clo]。

        Returns 返回值:
        float : Physiological Equivalent Temperature [℃]. 生理等效温度PET [℃]。
        """

        qAir_ = rh / 100.0 * 6.105 * math.exp(17.27 * tAir_C / (237.7 + tAir_C))
        age = age  # in years
        if sex == "male":
            sex = 1
        elif sex == "female":
            sex = 2
        elif sex == "average sex":
            sex = 3
        sex = 1 if sex < 1.5 else 2

        rCloClo = 0.02 if rCloClo - 0.02 < 0.01 else rCloClo
        eta = 0.0
        fCloAreaExpansion = 1.15  # Höppe
        # fCloAreaExpansion = 1 + (0.31 * rCloClo)  # Spasic

        if bodyPosture == "sitting":
            fEffRad = 0.696
        elif bodyPosture == "standing":
            fEffRad = 0.725
        elif bodyPosture == "crouching":
            fEffRad = 0.67

        pBaroHPA_setting, pBaroHPA = 1013.25, 1013.25
        densitiyBlood, cpBlood = 1.06, 3640.0
        Q_food, emissivitySkin, emissivityClo = 0.0, 0.99, 0.95
        coefEvap = 2.42 * 10 ** 6
        sigma = 5.67 * 10 ** -8
        areaBody = 0.203 * math.pow(weight_kg, 0.425) * math.pow(height_m, 0.725)
        # inner body energy
        Q_metaBase_female = 3.19 * math.pow(weight_kg, 0.75) * (1.0 + 0.004 * (30.0 - age) + 0.018 * (
                height_m * 100.0 / math.pow(weight_kg, 1.0 / 3.0) - 42.1))
        Q_metaBase_male = 3.45 * math.pow(weight_kg, 0.75) * (1.0 + 0.004 * (30.0 - age) + 0.01 * (
                height_m * 100.0 / math.pow(weight_kg, 1.0 / 3.0) - 43.4))
        Q_meta_male = Q_meta + Q_metaBase_male
        Q_meta_female = Q_meta + Q_metaBase_female
        Q_meta_people = 0.0
        if sex == 1:
            Q_meta_people = Q_meta_male
        elif sex == 2:
            Q_meta_people = Q_meta_female
        elif sex == 3:
            Q_meta_people = (Q_meta_male + Q_meta_female) / 2
        Q_meta_eff = Q_meta_people * (1.0 - eta)
        # sensible respiratory energy
        cAir = 1010.0
        tAirExpired = 0.47 * tAir_C + 21.0
        flowAir2Lung = 1.44 * math.pow(10.0, -6.0) * Q_meta_people
        Q_res_sensible = cAir * (tAir_C - tAirExpired) * flowAir2Lung
        # deferred respiration energy
        qAirRes = 6.11 * math.pow(10.0, 7.45 * tAirExpired / (235.0 + tAirExpired))
        Q_res_latent = 0.623 * coefEvap / pBaroHPA * (qAir_ - qAirRes) * flowAir2Lung
        Q_res = Q_res_sensible + Q_res_latent

        f_clo_area_ = (173.51 * rCloClo - 2.36 - 100.76 * rCloClo * rCloClo + 19.28 * math.pow(
            rCloClo, 3.0)) / 100.0
        rateWetSkin = 0.0
        # c = [None for i in range(11)]
        tcore = [None for i in range(7)]

        hConv = 2.67 + 6.5 * math.pow(vAir, 0.67)
        hConv = hConv * math.pow(pBaroHPA / pBaroHPA_setting, 0.55)

        if f_clo_area_ > 1.0:
            f_clo_area_ = 1.0
        rClo = rCloClo / 6.45 / f_clo_area_
        if rCloClo >= 2.0:
            ratioHeightCylClo = 1.0
        if rCloClo > 0.6 and rCloClo < 2.0:
            ratioHeightCylClo = (height_m - 0.2) / height_m
        if rCloClo <= 0.6 and rCloClo > 0.3:
            ratioHeightCylClo = 0.5
        if rCloClo <= 0.3 and rCloClo > 0.0:
            ratioHeightCylClo = 0.1
        radius_cyl_outer = areaBody * (fCloAreaExpansion - 1.0 + f_clo_area_) / (
                6.28 * height_m * ratioHeightCylClo)
        radius_cyl_inner = f_clo_area_ * areaBody / (6.28 * height_m * ratioHeightCylClo)
        radius_difference = radius_cyl_outer - radius_cyl_inner
        # skin temperatures
        for idxTcoreRoot in range(1, 7):
            tsk = 34.0
            iBit = 0
            tcl = (tAir_C + tMrt_C + tsk) / 3.0
            Q_balance2 = 0.0
            while True:
                for i100 in range(1, 100):
                    areaClo = areaBody * f_clo_area_ + areaBody * (fCloAreaExpansion - 1.0)
                    Q_L = emissivityClo * sigma * (
                            math.pow(tcl + 273.2, 4.0) - math.pow(tMrt_C + 273.2, 4.0)) * fEffRad
                    hCloCond = 6.28 * height_m * ratioHeightCylClo * radius_difference / (
                            rClo * math.log(radius_cyl_outer / radius_cyl_inner) * areaClo)
                    tsk = 1.0 / hCloCond * (hConv * (tcl - tAir_C) + Q_L) + tcl
                    # radiation balance
                    AreaEffRad = areaBody * fEffRad
                    Q_rad_bare = AreaEffRad * (
                            1.0 - f_clo_area_) * emissivitySkin * sigma * (
                                         math.pow(tMrt_C + 273.2, 4.0) - math.pow(tsk + 273.2, 4.0))
                    Q_rad_clo = fEffRad * areaClo * emissivityClo * sigma * (
                            math.pow(tMrt_C + 273.2, 4.0) - math.pow(tcl + 273.2, 4.0))
                    Q_rad = Q_rad_bare + Q_rad_clo
                    # convection
                    Q_conv_bare = hConv * (tAir_C - tsk) * areaBody * (
                            1.0 - f_clo_area_)
                    Q_conv_clo = hConv * (tAir_C - tcl) * areaClo
                    Q_conv = Q_conv_bare + Q_conv_clo
                    # core temperature
                    # 注1：式1的系数为bnCoef_A，bnCoef_B，bnCoef_C；是舒张与收缩公式
                    # Note1: Coefficients for Equation 1 (relaxation and contraction formula) are bnCoef_A, bnCoef_B, bnCoef_C
                    # 注2：式2的系数为bnCoef_A，bnCoef_B2，bnCoef_C2；是舒张公式
                    # Note2: Coefficients for Equation 2 (relaxation formula) are bnCoef_A, bnCoef_B2, bnCoef_C2
                    Q_core = Q_meta_eff + Q_res
                    blood_density_cp_product = areaBody * densitiyBlood * cpBlood
                    intermForEq1_1 = 18.0 - 0.5 * tsk
                    intermForEq1_2 = 5.28 * areaBody * intermForEq1_1
                    bnCoef_A1 = 13.0 / 625.0 * blood_density_cp_product
                    intermForEq2 = 0.76075 * blood_density_cp_product
                    bnCoef_B1 = intermForEq1_2 - intermForEq2 - tsk * bnCoef_A1
                    bnCoef_C1 = -Q_core * intermForEq1_1 - tsk * intermForEq1_2 + tsk * intermForEq2
                    bnCoef_Delta = bnCoef_B1 * bnCoef_B1 - 4.0 * bnCoef_A1 * bnCoef_C1
                    bnCoef_B2 = 5.28 * areaBody - intermForEq2 - bnCoef_A1 * tsk
                    bnCoef_C2 = intermForEq2 * tsk - Q_core - 5.28 * areaBody * tsk
                    bnCoef_Delta2 = bnCoef_B2 * bnCoef_B2 - 4.0 * bnCoef_A1 * bnCoef_C2
                    if tsk == 36.0:  tsk = 36.01

                    # 注：默认血流速率的解
                    # Note：Solution of the setting blood flow rate
                    tcore[6] = Q_core / (5.28 * areaBody + blood_density_cp_product * 6.3 / 3600.0) + tsk
                    # 注：自皮肤的血流量和冷信号的值
                    # Note：Values of blood flow from skin and cold signals
                    tcore[2] = Q_core / (
                            5.28 * areaBody + blood_density_cp_product * 6.3 / 3600.0 / (
                            1.0 + 0.5 * (34.0 - tsk))) + tsk
                    if bnCoef_Delta2 >= 0.0:
                        tcore[0] = (-bnCoef_B2 + math.pow(bnCoef_Delta2, 0.5)) / (2.0 * bnCoef_A1)
                        tcore[5] = (-bnCoef_B2 - math.pow(bnCoef_Delta2, 0.5)) / (2.0 * bnCoef_A1)
                    if bnCoef_Delta >= 0.0:
                        tcore[1] = (-bnCoef_B1 + math.pow(abs(bnCoef_Delta), 0.5)) / (2.0 * bnCoef_A1)
                        tcore[4] = (-bnCoef_B1 - math.pow(abs(bnCoef_Delta), 0.5)) / (2.0 * bnCoef_A1)
                    # 注：最大血流速率解
                    # Note：Solution of the maxing blood flow rate
                    tcore[3] = Q_core / (5.28 * areaBody + blood_density_cp_product * 1.0 / 40.0) + tsk
                    # transpiration
                    tbody = 0.1 * tsk + 0.9 * tcore[idxTcoreRoot - 1]
                    flow_sweat_male = 304.94 * (tbody - 36.6) * areaBody / 3600000.0
                    qSkinSaturated = 6.11 * math.pow(10.0, 7.45 * tsk / (235.0 + tsk))
                    if tbody <= 36.6:
                        flow_sweat_male = 0.0
                    flow_sweat_female = 0.7 * flow_sweat_male
                    if sex == 1:
                        flow_sweat_aveage = flow_sweat_male
                    if sex == 2:
                        flow_sweat_aveage = flow_sweat_female
                    if sex == 3:
                        flow_sweat_aveage = (flow_sweat_male + flow_sweat_female) / 2
                    Q_sweatLatent = -flow_sweat_aveage * coefEvap
                    rVel = 0.633 * hConv / (pBaroHPA * cAir)
                    fCloTransMass = 1.0 / (1.0 + 0.92 * hConv * rClo)
                    Q_sweatSensible = rVel * (
                            qAir_ - qSkinSaturated) * areaBody * coefEvap * fCloTransMass
                    rateWetSkin = Q_sweatLatent / Q_sweatSensible
                    if rateWetSkin > 1.0:
                        rateWetSkin = 1.0
                    Q_sweat_difference = Q_sweatLatent - Q_sweatSensible
                    if Q_sweat_difference <= 0.0:
                        Q_sweat = Q_sweatSensible
                    if Q_sweat_difference > 0.0:
                        Q_sweat = Q_sweatLatent
                    if Q_sweat > 0.0:
                        Q_sweat = 0.0
                    # diffusion
                    rSkinEvap_ = 0.79 * math.pow(10.0, 7.0)
                    rCloEvap_ = 0.0
                    Q_diffusion = coefEvap / (
                            rSkinEvap_ + rCloEvap_) * areaBody * (
                                          1.0 - rateWetSkin) * (
                                          qAir_ - qSkinSaturated)
                    # max vb
                    vb1 = 34.0 - tsk
                    vb2 = tcore[idxTcoreRoot - 1] - 36.6
                    if vb2 < 0.0:
                        vb2 = 0.0
                    if vb1 < 0.0:
                        vb1 = 0.0
                    vb = (6.3 + 75.0 * vb2) / (1.0 + 0.5 * vb1)
                    # energy balance
                    Q_balance = Q_meta_eff + Q_diffusion + Q_res + Q_sweat + Q_conv + Q_rad + Q_food
                    # clothing temperature
                    if iBit == 0:
                        xx = 1.0
                    if iBit == 1:
                        xx = 0.1
                    if iBit == 2:
                        xx = 0.01
                    if iBit == 3:
                        xx = 0.001
                    if Q_balance > 0.0:
                        tcl = tcl + xx
                    if Q_balance < 0.0:
                        tcl = tcl - xx
                    # Purpose【Q_balance和Q_balance2符号不相等时，说明收敛了】
                    if (Q_balance > 0.0 or Q_balance2 <= 0.0) and (
                            Q_balance < 0.0 or Q_balance2 >= 0.0):
                        Q_balance2 = Q_balance
                        i100 += 1
                    else:
                        break
                if iBit == 0.0 or iBit == 1.0 or iBit == 2.0:
                    iBit = iBit + 1
                    Q_balance2 = 0.0
                else:
                    break

            # 疑惑：求了这个根之后，判断根是否位置不对，似乎并不需要迭代20次，迭代1次就行
            # Question: Determine the root after calculation and check if the root position is not correct. It seems that only one iteration is needed, not 20 iterations.
            # 注：根的位置是被排列好的，每个根大概什么位置，什么数值是能被预知的
            # Note：The position of the root is arranged, and the approximate position and value of each root can be predicted.
            for k in range(20):
                # 注：如果是式2解，则bnCoef_Delta2 >= 0.0
                # Note: If it's the solution of equation 2, then bnCoef_Delta2 >= 0.0.
                if iBit == 3.0 and (idxTcoreRoot != 2 and idxTcoreRoot != 5):
                    if idxTcoreRoot != 6 and idxTcoreRoot != 1:
                        if idxTcoreRoot != 3:
                            if idxTcoreRoot != 7:
                                if idxTcoreRoot == 4:
                                    convergeFlag = True
                                    break
                            else:
                                if tcore[idxTcoreRoot - 1] >= 36.6 or tsk <= 34.0:
                                    convergeFlag = False
                                    break
                                convergeFlag = True
                                break
                        else:
                            if tcore[idxTcoreRoot - 1] >= 36.6 or tsk > 34.0:
                                convergeFlag = False
                                break
                            convergeFlag = True
                            break
                    else:
                        if bnCoef_Delta2 < 0.0 or (tcore[idxTcoreRoot - 1] < 36.6 or tsk <= 33.85):
                            convergeFlag = False
                            break
                        convergeFlag = True
                        break
                if bnCoef_Delta < 0.0 or (tcore[idxTcoreRoot - 1] < 36.6 or tsk > 34.05):
                    convergeFlag = False
                    break

            if convergeFlag == False:
                continue
            else:
                if (idxTcoreRoot == 4 or vb < 91.0) and (idxTcoreRoot != 4 or vb >= 89.0):
                    if vb > 90.0:
                        vb = 90.0
                    # water loss
                    ws = flow_sweat_aveage * 3600.0 * 1000.0
                    if ws > 2000.0:
                        ws = 2000.0
                    wd = Q_diffusion / coefEvap * 3600.0 * (-1000.0)
                    wr = Q_res_latent / coefEvap * 3600.0 * (-1000.0)
                    wsum = ws + wr + wd
                    break
                    # return tcore[idxTcoreRoot - 1], Q_rad, Q_conv, Q_diffusion
            # water loss
            ws = flow_sweat_aveage * 3600.0 * 1000.0
            wd = Q_diffusion / coefEvap * 3600.0 * (-1000.0)
            wr = Q_res_latent / coefEvap * 3600.0 * (-1000.0)
            wsum = ws + wr + wd
            if idxTcoreRoot - 3 < 0:
                index = 3
            else:
                index = idxTcoreRoot - 3
            # tcore[index]
            break

        valPET = tAir_C
        Q_balance_PET2 = 0.0
        iBit = 0
        while iBit != 4:
            hConv_PET = 2.67 + 6.5 * math.pow(0.1, 0.67)
            hConv_PET = hConv_PET * math.pow(pBaroHPA / pBaroHPA_setting, 0.55)
            # radiation saldo
            AreaEffRad = areaBody * fEffRad
            Q_rad_bare_PET = AreaEffRad * (1.0 - f_clo_area_) * emissivitySkin * sigma * (
                    math.pow(valPET + 273.2, 4.0) - math.pow(tsk + 273.2, 4.0))
            Q_rad_clo_PET = fEffRad * areaClo * emissivityClo * sigma * (
                    math.pow(valPET + 273.2, 4.0) - math.pow(tcl + 273.2, 4.0))
            Q_rad_PET = Q_rad_bare_PET + Q_rad_clo_PET
            # convection
            Q_conv_bare_PET = hConv_PET * (valPET - tsk) * areaBody * (1.0 - f_clo_area_)
            Q_conv_clo_PET = hConv_PET * (valPET - tcl) * areaClo
            Q_conv_PET = Q_conv_bare_PET + Q_conv_clo_PET
            # diffusion
            Q_diffusion_PET = coefEvap / (rSkinEvap_ + rCloEvap_) * areaBody * (
                    1.0 - rateWetSkin) * (12.0 - qSkinSaturated)
            # breathing
            tAirExpired_PET = 0.47 * valPET + 21.0
            Q_res_sensible_PET = cAir * (valPET - tAirExpired_PET) * flowAir2Lung
            qAirRes_PET = 6.11 * math.pow(10.0, 7.45 * tAirExpired_PET / (235.0 + tAirExpired_PET))
            Q_res_latent_PET = 0.623 * coefEvap / pBaroHPA * (
                    12.0 - qAirRes_PET) * flowAir2Lung
            Q_res_PET = Q_res_sensible_PET + Q_res_latent_PET
            # energy balance
            Q_balance_PET = Q_meta_eff + Q_diffusion_PET + Q_res_PET + Q_sweat + Q_conv_PET + Q_rad_PET
            if iBit == 0:
                xx = 1.0
            if iBit == 1:
                xx = 0.1
            if iBit == 2:
                xx = 0.01
            if iBit == 3:
                xx = 0.001
            if Q_balance_PET > 0.0:
                valPET = valPET - xx
            if Q_balance_PET < 0.0:
                valPET = valPET + xx
            if (Q_balance_PET > 0.0 or Q_balance_PET2 <= 0.0) and (
                    Q_balance_PET < 0.0 or Q_balance_PET2 >= 0.0):
                Q_balance_PET2 = Q_balance_PET
            else:
                iBit = iBit + 1

        return valPET

    # Reference:
    # [1] W. Cheng, R.D. Brown, An energy budget model for estimating the thermal comfort of children, Int. J. Biometeorol. 64 (2020) 1355–1366. https://doi.org/10.1007/s00484-020-01916-x.
    # [2] N.A. Kenny, J.S. Warland, R.D. Brown, T.G. Gillespie, Part A: Assessing the performance of the comfa outdoor thermal comfort model on subjects performing physical activity, Int. J. Biometeorol. 53 (2009) 415–428. https://doi.org/10.1007/s00484-009-0226-3.
    # [3] N.A. Kenny, J.S. Warland, R.D. Brown, T.G. Gillespie, Part B: Revisions to the COMFA outdoor thermal comfort model for application to subjects performing physical activity, Int. J. Biometeorol. 53 (2009) 429–441. https://doi.org/10.1007/s00484-009-0227-2.
    # [4] J.K. Vanos, J.S. Warland, T.J. Gillespie, N.A. Kenny, Improved predictive ability of climate-human-behaviour interactions with modifications to the COMFA outdoor energy budget model, Int. J. Biometeorol. 56 (2012) 1065–1074. https://doi.org/10.1007/s00484-012-0522-1.
    # [5] R.D. Brown, T.J. Gillespie, Estimating outdoor thermal comfort using a cylindrical radiation thermometer and an energy budget model, Int. J. Biometeorol. 30 (1986) 43–52. https://doi.org/10.1007/BF02192058.
    @staticmethod
    def cal_COMFA(tAirC: float = 28.0,
                  ksolar: float = 100.0,
                  vAir: float = 5.0,
                  rh: float = 60.0,
                  height_m: float = 1.8,
                  weight_kg: float = 78,
                  rcl_clo: float = 0.264,
                  Ma: float = 116.0,
                  sex: Union[int, str] = 1,
                  age: int = 38,
                  alt: float = 40.0,
                  beam: float = 1.0,
                  SVF: float = 1.0) -> float:
        """
        Parameters 参数:
        tAirC : float : Air temperature [℃]. 空气温度[℃]。
        ksolar : float : Solar radiation (or Global horizontal irradiance, GHI) [W/m2]. 太阳辐射[W/m2] (或全球水平辐照度GHI)。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。
        height_m : float : Height of the person [m]. 身高[m]。
        weight_kg : float : Weight of the person [kg]. 体重[kg]。
        rcl_clo : float : Clothing insulation [clo]. 服装热阻[clo]。
        Ma : float : Metabolic rate [W/m2]. 新陈代谢率[W/m2]。
        sex : Union[int, str] : Sex of the person. 1 for male, 2 for female. 性别，1或'male'为男性，2或 'female'为女性。
        age : int : Age [years]. 年龄[岁]。
        alt : float : Solar altitude angle [deg]. 太阳高度角[°]。
        beam : float : Solar beam fraction [−]. 直射辐射透过率[−]。
        SVF : float : Sky View Factor [−]. 天空角系数[−]。

        Returns 返回值:
        Q_Balance: float: Energy balance budget (COMFAcourtyard) [W/m2]. COMFA[W/m2]。
        Q_Meta: float: Metabolic heat production [W/m2]. 代谢产热量[W/m2]。
        Q_Rad: float: Radiation absorption heat [W/m2]. 辐射吸收热量[W/m2]。
        Q_Conv: float: Convective heat dissipation [W/m2]. 对流散热[W/m2]。
        Q_Evap: float: Evaporative heat dissipation [W/m2]. 蒸发散热[W/m2]。
        Q_L: float: Longwave emission heat [W/m2]. 长波发射热量 [W/m2]。
        level: int: Thermal comfort level [-]. 热舒适级别 [-]。
        """

        def clo2sm(clo):
            rclKMn2Wn1 = 0.155 * clo
            rSnCM = rclKMn2Wn1 / 0.082
            rSnM = rSnCM * 100.0
            return rSnM

        def calRad(SVF, alt, beam, kSolar, tAirK):
            ## Purpose【albedo】
            albedoCo = 37 / 100;
            difFrac = 10.0
            albedoSky = 0.0
            albedoGnd = 25.0
            DHI = (difFrac * kSolar / 100)
            difFracPercent = difFrac / 100
            RadianConvert = 0.017453293
            altRad = alt * RadianConvert
            K_Flt = ((1 - difFracPercent) * kSolar) / math.tan(altRad)
            K_Cyl = K_Flt / 3.141592654
            albedoGndPercent = albedoGnd / 100
            pi = 3.14
            L_Ref = (((kSolar * beam * albedoGndPercent) / 2) + ((kSolar * beam * albedoSky / 100) / 2 / pi))
            K_Abs = (((K_Cyl * beam) + (SVF * DHI) + DHI * (1 - SVF) * albedoSky) + (L_Ref)) * (1 - albedoCo)
            L_Sky = (1.2 * (5.67E-8 * (tAirK ** 4))) - 171
            L_Wall = 0.98 * (5.67E-8 * (tAirK ** 4)) * (1 - SVF)
            L_Gnd = 0.98 * (5.67E-8 * (tAirK ** 4))
            L_Abs = (((L_Wall + (SVF * L_Sky)) * 0.5) + (L_Gnd * 0.5)) * 0.98
            R_Abs = (K_Abs + L_Abs) * 0.8
            return R_Abs

        vActivity = 0.0
        sex = 1 if sex < 1.5 else 2
        if sex == 1:
            rmr = 10 * weight_kg + 6.25 * height_m * 100 - 5 * age + 5
        else:
            rmr = 9.99 * weight_kg + 6.25 * height_m * 100 - 4.92 * age - 161
        Watt = rmr / 24 * 1.163
        facActivity = Ma / 80.0
        Ma = Watt * facActivity

        # Purpose: temperature and radiance
        tAirK = tAirC + 273.15
        kSolar = ksolar
        # Purpose: clothing
        rCoStatic = clo2sm(rcl_clo);
        rCoVaporStatic = rCoStatic * 2.74524064171123;
        vRelative = 1.92
        rCoVapor = rCoVaporStatic * (-0.8 * (1 - math.exp(-vRelative / 1.095)) + 1);
        rClo = rCoStatic * (-0.37 * (1 - math.exp(-vActivity / 0.72)) + 1)  # 50

        # Rabs
        R_Abs = calRad(SVF, alt, beam, kSolar, tAirK)

        # Metabolism
        rhPercent = rh / 100
        evap = ((0.6108 * (math.exp((17.269 * tAirC) / (tAirC + 237.3)))) * rhPercent)
        f_Latent = (0.15 - (0.0173 * evap) - (0.0014 * tAirC))
        Q_Meta = (1 - f_Latent) * Ma

        # Purpose【Tc】
        tCoreC = 36.5 + ((0.0043) * Q_Meta)

        # Purpose【Tsk】
        ReNumber = 11333 * vAir
        if ReNumber >= 40000:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.0266 * (ReNumber ** 0.805)))
        elif ReNumber >= 4000 and ReNumber < 40000:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.683 * (ReNumber ** 0.466)))
        else:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.193 * (ReNumber ** 0.618)))
        eSweat = 0.42 * (Q_Meta - 58)
        cAir = 1212.
        rTissue = cAir / (0.13 * eSweat + 15)
        tSkinC = (((tCoreC - tAirC) / (rTissue + rClo + rAerodynamic)) * (rAerodynamic + rClo)) + tAirC

        # CONV
        Q_Conv = cAir * ((tSkinC - tAirC) / (rClo + rAerodynamic))
        # Tsf
        tSrfC = ((tSkinC - tAirC) / (rClo + rAerodynamic)) * rAerodynamic + tAirC
        # L
        Q_L = 0.75 * ((0.95 * 5.67E-8) * ((tSrfC + 273) ** 4))
        # qs
        qSkin = 0.6108 * (math.exp((17.269 * tSkinC) / (tSkinC + 237.3)))

        # Evap
        pLv = (1.16 * 2442)
        qAir = evap
        rAerodynamicVapor = (0.92 * rAerodynamic)
        rTissueVapor = 7700
        eDiffusion = pLv * (qSkin - qAir) / (rCoVapor + rAerodynamicVapor + rTissueVapor)
        eMax = pLv * ((qSkin - qAir) / (rCoVapor + rAerodynamicVapor))
        eSweat = 0.42 * (Q_Meta - 58)
        Q_Evap = eSweat + eDiffusion
        if eMax <= Q_Evap:  Q_Evap = eMax

        Q_Rad = R_Abs
        Q_Balance = Q_Meta + Q_Rad - Q_Conv - Q_Evap - Q_L

        def get_q_balance_level(q_balance: float) -> str:
            if q_balance < -150:
                return -2  # "Would prefer to be much warmer"
            elif -150 <= q_balance < -50:
                return -1  # "Would prefer to be warmer"
            elif -50 <= q_balance <= 50:
                return 0  # "No change"
            elif 50 < q_balance <= 150:
                return 1  # "Would prefer to be cooler"
            else:
                return 2  # "Would prefer to be much cooler"

        level = get_q_balance_level(Q_Balance)

        return Q_Balance, Q_Meta, Q_Rad, -Q_Conv, - Q_Evap, - Q_L, level

    # Reference:
    # [1].	Wu, R., et al., The COMFA model for assessing courtyard thermal comfort in hot and humid regions: A comparative study with existing models. Building and Environment, 2023. 234: p. 110150.
    # [2] W. Cheng, R.D. Brown, An energy budget model for estimating the thermal comfort of children, Int. J. Biometeorol. 64 (2020) 1355–1366. https://doi.org/10.1007/s00484-020-01916-x.
    # [3] N.A. Kenny, J.S. Warland, R.D. Brown, T.G. Gillespie, Part A: Assessing the performance of the comfa outdoor thermal comfort model on subjects performing physical activity, Int. J. Biometeorol. 53 (2009) 415–428. https://doi.org/10.1007/s00484-009-0226-3.
    # [4] N.A. Kenny, J.S. Warland, R.D. Brown, T.G. Gillespie, Part B: Revisions to the COMFA outdoor thermal comfort model for application to subjects performing physical activity, Int. J. Biometeorol. 53 (2009) 429–441. https://doi.org/10.1007/s00484-009-0227-2.
    # [5] J.K. Vanos, J.S. Warland, T.J. Gillespie, N.A. Kenny, Improved predictive ability of climate-human-behaviour interactions with modifications to the COMFA outdoor energy budget model, Int. J. Biometeorol. 56 (2012) 1065–1074. https://doi.org/10.1007/s00484-012-0522-1.
    # [6] R.D. Brown, T.J. Gillespie, Estimating outdoor thermal comfort using a cylindrical radiation thermometer and an energy budget model, Int. J. Biometeorol. 30 (1986) 43–52. https://doi.org/10.1007/BF02192058.
    @staticmethod
    def cal_COMFAcourtyard(tAir_C: float = 28.0,
                           vAir: float = 5.0,
                           rh: float = 60.0,
                           height_m: float = 1.8,
                           weight_kg: float = 78,
                           rcl_clo: float = 0.264,
                           Q_movement: float = 116.0,
                           sex: Union[int, str] = 'male',
                           age: int = 38,
                           tMrt_C: float = 30.0,
                           bodyPosture: str = 'standing') -> float:

        """
        Parameters 参数:
        tAir_C : float : Air temperature [℃]. 空气温度[℃]。
        vAir : float : Air velocity [m/s]. 风速[m/s]。
        rh : float : Relative humidity [%]. 相对湿度[%]。
        height_m : float : Height of the person [m]. 身高[m]。
        weight_kg : float : Weight of the person [kg]. 体重[kg]。
        rcl_clo : float : Clothing insulation [clo]. 服装热阻[clo]。
        Q_movement : float : Metabolic rate [W/m2]. 新陈代谢率[W/m2]。
        sex : Union[int, str] : Sex of the person. 'male' for male, 'female' for female. 性别，'male'为男性，'female'为女性。
        age : int : Age [years]. 年龄[岁]。
        tMrt_C : float : Mean radiant temperature [℃]. 平均辐射温度[℃]。
        bodyPosture : str : Body posture, can be 'standing', 'sitting', 'crouching'. 身体姿势，取值范围为'standing', 'sitting',或 'crouching'。

        Returns 返回值:
        Q_Balance: float: Energy balance budget (COMFAcourtyard) [W/m2]. COMFAcourtyard[W/m2]。
        Q_Meta: float: Metabolic heat production [W/m2]. 代谢产热量[W/m2]。
        Q_Rad: float: Radiant heat transfer [W/m2]. 辐射换热[W/m2]。
        Q_Conv: float: Convective heat transfer [W/m2]. 对流换热[W/m2]。
        Q_Evap: float: Evaporative heat transfer [W/m2]. 蒸发换热[W/m2]。
        Q_Cond: float: Conductive heat transfer [W/m2]. 导热换热[W/m2]。
        level: int: thermal comfort level [-]. 热舒适级别 [-]。
        """

        def clo2sm(clo):
            rclKMn2Wn1 = 0.155 * clo
            rSnCM = rclKMn2Wn1 / 0.082
            rSnM = rSnCM * 100.0
            return rSnM

        vActivity = 0.0

        if type(bodyPosture) == int or type(bodyPosture) == float:
            if bodyPosture < 1.5:
                bodyPosture = "sitting"
            elif bodyPosture < 2.5:
                bodyPosture = "standing"
            else:
                bodyPosture = "crouching"

        if type(sex) == int or type(sex) == float:
            sex = 1 if sex < 1.5 else 2
        if sex == 1 or sex == 'male':
            rateMetaRestOneDay_kCal = 10 * weight_kg + 6.25 * height_m * 100 - 5 * age + 5
        else:
            rateMetaRestOneDay_kCal = 9.99 * weight_kg + 6.25 * height_m * 100 - 4.92 * age - 161
        rateMetaRest_WMn2 = rateMetaRestOneDay_kCal / 24 * 1.163
        facActivity = Q_movement / 80.0
        Q_meta_total = rateMetaRest_WMn2 * facActivity

        rCoStatic = clo2sm(rcl_clo);
        rCoVaporStatic = rCoStatic * 2.74524064171123;
        vRelative = 1.92
        # vRelative = vAir
        rCoVapor = rCoVaporStatic * (-0.8 * (1 - math.exp(-vRelative / 1.095)) + 1);

        # Purpose【rCo】
        rClo = rCoStatic * (-0.37 * (1 - math.exp(-vActivity / 0.72)) + 1)

        # Purpose【M】
        rhPercent = rh / 100
        q_air_ = ((0.6108 * (math.exp((17.269 * tAir_C) / (tAir_C + 237.3)))) * rhPercent)
        friction_latent = (0.15 - (0.0173 * q_air_) - (0.0014 * tAir_C))
        Q_Meta = (1 - friction_latent) * Q_meta_total

        # Purpose【Tc】
        tCoreC = 36.5 + ((0.0043) * Q_Meta)

        # Purpose【Tsk】
        ReNumber = 11333 * vAir
        if ReNumber >= 40000:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.0266 * (ReNumber ** 0.805)))
        elif ReNumber >= 4000 and ReNumber < 40000:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.683 * (ReNumber ** 0.466)))
        else:
            rAerodynamic = (0.17 / ((0.71 ** 0.33) * 0.000022 * 0.193 * (ReNumber ** 0.618)))
        Q_Sweat = 0.42 * (Q_Meta - 58)
        rTissue = 1212 / (0.13 * Q_Sweat + 15)
        tSkinC = (((tCoreC - tAir_C) / (rTissue + rClo + rAerodynamic)) * (rAerodynamic + rClo)) + tAir_C

        # Purpose【Rabs】
        areaBody = 0.203 * (weight_kg ** 0.425) * (height_m ** 0.725)  # area of the body
        fCloArea = (173.51 * rcl_clo - 2.36 - 100.76 * rcl_clo * rcl_clo + 19.28 * (rcl_clo ** 3)) / 100  # 服装面积系数
        fCloAreaExpansion = 1.15  # 服装交换面积系数
        factor_area_effective_radiation = 0.696  # 辐射有效面积[m2]
        if bodyPosture == "sitting" or bodyPosture == 1:
            factor_area_effective_radiation = 0.696
        elif bodyPosture == "standing" or bodyPosture == 2:
            factor_area_effective_radiation = 0.725
        elif bodyPosture == "crouching" or bodyPosture == 3:
            factor_area_effective_radiation = 0.67
        areaClo = areaBody * fCloArea + areaBody * (fCloAreaExpansion - 1.)  # 服装面积[m2]
        emissivity_lw_sk = 0.877  # 皮肤长波发射率[-]
        emissivity_lw_cl = 0.591  # 衣服长波发射率[-]
        sigma = 5.67 * (10 ** (-8))  # Stefan-Boltzmann常数[W/(m2*K^(-4))]
        area_effect_radiation = areaBody * factor_area_effective_radiation  # 辐射有效面积

        # >>>> tClo
        tCloC = ((tSkinC - tAir_C) / (rClo + rAerodynamic)) * rAerodynamic + tAir_C
        # >>>> radius_difference
        if fCloArea > 1.0:
            fCloArea = 1.0
        rClo = rcl_clo / 6.45 / fCloArea
        ratioHeightCylClo = 1.0
        if rcl_clo >= 2.0:
            ratioHeightCylClo = 1.0
        if rcl_clo > 0.6 and rcl_clo < 2.0:
            ratioHeightCylClo = (height_m - 0.2) / height_m
        if rcl_clo <= 0.6 and rcl_clo > 0.3:
            ratioHeightCylClo = 0.5
        if rcl_clo <= 0.3 and rcl_clo > 0.0:
            ratioHeightCylClo = 0.1
        radius_cyl_outer = areaBody * (fCloAreaExpansion - 1.0 + fCloArea) / (
                6.28 * height_m * ratioHeightCylClo)
        radius_cyl_inner = fCloArea * areaBody / (6.28 * height_m * ratioHeightCylClo)
        radius_difference = radius_cyl_outer - radius_cyl_inner

        # >>>> hCloCond
        hCloCond = 6.28 * height_m * ratioHeightCylClo * radius_difference / (
                rClo * math.log(radius_cyl_outer / radius_cyl_inner) * areaClo)

        # >>>> SET*:hconv
        pBaroHPA = 1013.25
        pBaroAtm = pBaroHPA / 1013.25
        hConv = 3 * (pBaroAtm ** (0.53))
        if ((Q_movement / 58.2) < 0.85):
            hConvActivity = 0
        else:
            hConvActivity = 5.66 * ((((Q_movement / 58.2) - 0.85) * pBaroAtm) ** 0.39)
        hConvVelocity = 8.600001 * ((vAir * pBaroAtm) ** 0.53)
        if (hConv <= hConvActivity):    hConv = hConvActivity
        if (hConv < hConvVelocity):    hConv = hConvVelocity

        # >>>> tCloDelta
        hConvSkin = 1212 / (rAerodynamic)
        # >>>> >>>> Q_conv_clo
        QQ_conv_clo = hConv * (tAir_C - tCloC) * factor_area_effective_radiation * areaClo
        # >>>> >>>> Q_rad_clo
        QQ_rad_clo = factor_area_effective_radiation * areaClo * emissivity_lw_cl * sigma * (
                ((tMrt_C + 273.2) ** 4) - ((tCloC + 273.2) ** 4))
        hRadLinear = 4 * sigma * (((tCloC + tMrt_C) / 2 + 273.15) ** 3) * 0.72
        tCloDelta = (QQ_rad_clo / hRadLinear + QQ_conv_clo / hConv) / areaBody
        # >>>> Q_Cond
        Q_Cond = hCloCond * ((tCloDelta + tCloC) - tSkinC)
        # >>>> Q_Conv
        QQ_conv_bare = - hConv * (tAir_C - tSkinC) * areaBody * (1. - fCloArea)
        Q_Conv = QQ_conv_bare / areaBody
        # >>>> Q_Rad
        QQ_rad_bare = area_effect_radiation * (1. - fCloArea) * emissivity_lw_sk * sigma * (
                ((tMrt_C + 273.2) ** 4) - ((tSkinC + 273.2) ** 4))
        Q_Rad = QQ_rad_bare / areaBody

        # qs
        qSkin = 0.6108 * (math.exp((17.269 * tSkinC) / (tSkinC + 237.3)))
        # Y
        pLv = (1.16 * 2442)
        qAir = q_air_
        rAerodynamicVapor = (0.92 * rAerodynamic)
        rTissueVapor = 7700
        Q_Diffusion = pLv * (qSkin - qAir) / (rCoVapor + rAerodynamicVapor + rTissueVapor)
        eMax = pLv * ((qSkin - qAir) / (rCoVapor + rAerodynamicVapor))
        Q_Sweat = 0.42 * (Q_Meta - 58)
        Q_Evap = Q_Sweat + Q_Diffusion
        if eMax <= Q_Evap:  Q_Evap = eMax

        Q_Balance = Q_Meta + Q_Rad - Q_Conv - Q_Evap + Q_Cond

        def get_q_balance_level(q_balance: float) -> int:
            if q_balance < -118.85181:
                return -3  # Very Cold
            elif -118.85181 <= q_balance < -62.29072398:
                return -2  # Cold
            elif -62.29072398 <= q_balance < -5.729638009:
                return -1  # Cool (Comfort Zone)
            elif -5.729638009 <= q_balance < 50.83144796:
                return 0  # Neutral (Comfort Zone)
            elif 50.83144796 <= q_balance < 107.3925339:
                return 1  # Warm (Comfort Zone)
            elif 107.3925339 <= q_balance < 163.9536199:
                return 2  # Hot
            else:
                return 3  # Very Hot

        level = get_q_balance_level(Q_Balance)

        return (Q_Balance, Q_Meta, Q_Rad, - Q_Conv, - Q_Evap, Q_Cond, level)


class mainT:

    @staticmethod
    def _main1():
        val_PET = OCT.cal_PET()
        val_SETstar = OCT.cal_SETstar()[0]
        val_UTCI = OCT.cal_UTCI()
        val_WBGT = OCT.cal_WBGToutdoor()
        val_COMFA = OCT.cal_COMFA()[0]
        val_COMFAcourtyard = OCT.cal_COMFAcourtyard()[0]
        print("PET: ", val_PET)
        print("SET*: ", val_SETstar)
        print("UTCI: ", val_UTCI)
        print("WBGT: ", val_WBGT)
        print("COMFA: ", val_COMFA)
        print("COMFAcourtyard: ", val_COMFAcourtyard)


if __name__ == '__main__':
    mainT._main1()
