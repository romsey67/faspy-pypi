import numba
import numpy as np



@numba.njit('float64(float64,float64,float64,float64)')
def bond_price(tenor, yld, cpn, fre):

    total_cpn = np.int32(tenor * fre)
    cpn_values = 0.00
    for period in range(1,total_cpn + 1, 1):
        cpn_value = (cpn / fre) / (1 + (yld * 0.01) / fre) ** (period)
        cpn_values += cpn_value


    pv_fv = 100 / (1 + (yld  * 0.01) / fre ) ** (tenor * fre)
    price = cpn_values + pv_fv
    #if yld == 5 :
    #    print(tenor, yld, cpn, fre, price)
    return price


@numba.njit('float64(float64,float64,float64,float64)')
def bond_duration(tenor, yld, cpn, fre):
    inc = 0.0001
    price0 = bond_price(tenor, yld, cpn, fre)
    price1 = bond_price(tenor, yld + inc, cpn, fre)
    duration = (price1 - price0)/inc
    return duration


@numba.njit('float64(float64,float64,float64,float64)')
def bond_convexity(tenor, yld, cpn, fre):
    inc = 0.0001
    duration0 = bond_duration(tenor, yld, cpn, fre)
    duration1 = bond_duration(tenor, yld + inc, cpn, fre)
    convexity = (duration1 - duration0)/inc
    return convexity
