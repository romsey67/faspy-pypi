import numpy as np
from faspy.interestrate.rmp_stats import rmpcdf, rmpppf
import sympy as sy
from numba import njit

# %%
def value_at_risk(assets, correlations):
    np_assets = np.array(assets)
    np_corr = np.array(correlations)
    value = np_assets.dot(np_corr)
    return float(np.sqrt(value.dot(np.transpose(np_assets))))


def logreturns(prices, period=1):
    """
    Calculate the log returns of prices
        Parameters:
            prices: numpy.array. An array of prices. Variables are column based
            period: int. number of observation period to calulate the returns
        Returns:
            numpy.array. An array of log returns

    """
    rows, columns = prices.shape
    # slicing arrays
    price_current = prices[: rows - period]
    price_last = prices[period: ]
    returns = np.divide(price_current, price_last)
    log_returns = np.log(returns)
    return log_returns


def volatilities(data):
    return np.std(data, axis=0)


def vol_from_rates(rate, period=1):
    rows, columns = rate.shape
    # slicing arrays
    price_current = rate[: rows - period]
    price_last = rate[period:]
    result = price_current - price_last
    return np.std(result, axis=0)


def returns_from_rates(data, period=1):
    rows, columns = data.shape
    price_current = data[: rows - period]
    price_last = data[period:]
    result = price_current - price_last
    return result


def returns_from_prices(data, period=1):
    rows, columns = data.shape
    price_current = data[: rows - period]
    price_last = data[period:]
    returns = np.divide(price_current, price_last)
    log_returns = np.log(returns)
    return log_returns


def correlation(data):
    return np.corrcoef(data, rowvar=False)


def get_cfweight(point1, point2, point):
    weight = sy.Symbol('weight')

    expr = weight * point1 + (1 - weight) * point2
    solved_rates = sy.solveset(expr - point, weight, domain=sy.S.Reals)
    return solved_rates


def get_cfalloc(vol1, vol2, vol, cor1_2):

    alpha = sy.Symbol('alpha')
    # a = sy.Symbol('a')
    # b = sy.Symbol('b')
    # c = sy.Symbol('c')
    a = vol1**2 + vol2**2 - 2*cor1_2*vol1*vol2
    b = 2 * cor1_2 * vol1 * vol2 - 2 * vol2 ** 2
    c = vol2**2 - vol**2

    expr = a*alpha**2 + b*alpha + c
    solved_rates = sy.solveset(expr, alpha, domain=sy.S.Reals)
    return solved_rates


# mapping of cash flow to vertices is only applicable to forward products
# hence the requirement for discount factors
def map_cf_to_var_vertices(cashflows, disfac, var_vertices, corr, vola):
    # cashflows is an array of dictionary [{"times": float, "pv": float},]
    cf_vertices = list(map(lambda x: 0, var_vertices))
    mymax = max(var_vertices)
    mymin = min(var_vertices)
    for cf in cashflows:
        if cf['times'] < mymin:
            index = 0
            vol = vola[index]
            cf_vertices[index] = cf_vertices[index] + cf['pv']
        
        elif cf['times'] > mymax:
            index = len(var_vertices) -1
            vol = vola[index]
            cf_vertices[index] = cf_vertices[index] + cf['pv']
        else:
            booltime = list(map(lambda x: x > cf["times"], var_vertices))
            index = booltime.index(True)
            
            point1, point2 = disfac[index - 1], disfac[index]
            weight = get_cfweight(point1, point2, cf['df'])
            weight = list(weight)
            weight = float(weight[0])

            vol1, vol2 = vola[index-1], vola[index]
            vol = weight * vol1+ (1-weight)* vol2
            cor1_2 = corr[index-1][index]

            cf_alloc = get_cfalloc(vol1,vol2,vol,cor1_2)
            cf_alloc = list(cf_alloc)

            cf_alloc = list(map(lambda x: float(x),cf_alloc))
            cf_alloc = list(filter(lambda x: x>=0 and x <=1,cf_alloc))
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #The following need to be check
            if len(cf_alloc)==0:
                cf_alloc = 1
            else:
                cf_alloc = cf_alloc[0]
            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            cf_vertices[index-1] += cf_alloc * cf['pv']
            cf_vertices[index] += (1-cf_alloc) * cf['pv']
    return cf_vertices


def var_asset_weightbyprice(assets, confidence, vol, fxrate=1):
    np_assets = np.array(assets) * np.array(vol)
    ppf = rmpppf(confidence)
    return np_assets * ppf * fxrate


def var_asset_weightbyrate(assets, confidence, vol, bpvs, fxrate=1):
    # assume rates used for vol calculation is expressed in terms of percentage
    np_assets = np.array(assets) * np.array(vol) * 10000 * np.array(bpvs)
    ppf = rmpppf(confidence)
    return np_assets * ppf * fxrate


def pvbp01_continuousrate(data):
    # data is an array of array [[time, rate]]
    # time is a fraction of a year
    # rate is in terms of percentage
    data_copy = list(data)
    for datum in data_copy:
        pvbp01 = (np.exp(-(datum[0] * datum[1] / 100 + 0.0001)) -
                  np.exp(-(datum[0] * datum[1] / 100)))
        datum.append(pvbp01)

    return data_copy
