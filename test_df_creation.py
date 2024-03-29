#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:36:21 2020

@author: RMS671214
"""

from faspy2.interestrate import discount_curve as dcurve
from numba import types, njit, jit
from numba.typed import Dict
import numpy as np
from faspy2.interestrate.rmp_dates import tenor_to_maturity as ttm, \
day_count_factor as day_cf
import time
from collections import deque

# %%
rates = []
holidays = []
value_date = np.datetime64('2021-01-04')
st_set = Dict.empty(key_type=types.unicode_type, value_type=types.unicode_type,)
st_set['business_day'] = 'Modified Following'
st_set['rate_basis'] = 'Simple'
st_set['day_count'] = 'Actual/365'
st_tenor = np.array(['1W', '1M', '3M', '6M', '12M'])
st_rates = np.array([2.35, 2.45, 2.55, 2.65, 2.75])





#lst_tenor = list(st_tenor)
#njit
def test(valuedate, tenors, settings):
        size = tenors.shape[0]
        dates = np.empty(shape=size, dtype='datetime64[s]')
        dcfs = np.empty(shape=size, dtype='float64')
        dfs = np.empty(shape=size, dtype='float64')
        for i in range(size):
                dates[i] = ttm(valuedate, tenors[i], convention=settings['day_count'],
                                business_day=settings['business_day'],
                                holidays=holidays) 
                dcfs[i] = day_cf(settings['day_count'], valuedate, dates[i])
        dfs = 1/(1 + (np.array(dcfs) * st_rates /100))
        return dfs
    

start0 = time.perf_counter()
result = test(value_date,st_tenor, st_set)
end0 = time.perf_counter()
print(f"Time required for list comprehension: {end0 - start0:0.8f} seconds")
print (result)


lt_set = Dict.empty(key_type=types.unicode_type, value_type=types.unicode_type,)
lt_set['business_day'] = 'No Adjustment'
lt_set['frequency'] = 'Semi_Annual'
lt_set['day_count'] = 'Actual/Actual'
lt_set['date_gen'] = 'Forward from issue date'


rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Simple", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}
rates.append(rate)
rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Simple", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 3.30, '1W': 3.35, '1M': 3.45, '3M': 3.55,
                  '6M': 3.65, '12M': 3.75, '1Y': 3.70, '2Y': 3.80, '3Y': 3.90,
                  '5Y': 4.00, '10Y': 4.10, '30Y': 4.25}}
rates.append(rate)
#print(rates)
start0 = time.perf_counter()
dfs = dcurve.discount_factor_gen(rate, return_type="times")

end0 = time.perf_counter()
print(f"Time required with old method: {end0 - start0:0.8f} seconds")
#print(dfs)

# %%

flat = dcurve.flat_curve("2020-10-30", "2020-12-30", 5.30, rate_basis="Money Market",
                         day_count="Actual/365", bus_day="No Adjustment",
                         tenors=None, return_type="days")
