#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from interestrate import fixincome as fi
import numpy as np
from interestrate import rmp_dates as rd
from interestrate import discount_curve as dcurve
from interestrate import rmp_curves as rcurve
import pandas as pd
# %%

myfl = {}
myfl['issue_date'] = np.datetime64('2018-10-22')
myfl['value_date'] = np.datetime64('2021-10-22')
myfl['maturity'] = np.datetime64('2028-10-22')
myfl['day_count'] = 'Actual/365 Fixed'
myfl['frequency'] = 'Semi-Annual'
myfl['business_day'] = 'No Adjustment'
myfl['date_generation'] = rd.date_gen_method[0]
myfl['face_value'] = 1000000
myfl['coupon'] = 10
myfl["principal_exchange"] = False


# %%

rate = {"value_date": str(myfl['value_date']),
        "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}


dfs = dcurve.discount_factor_gen(rate, return_type="times")
x_axis = [x["times"] for x in dfs]
y_axis = [x["df"] for x in dfs]
df_func = rcurve.interpolation(x_axis, y_axis, float(1/366), is_function=True)
pd_ori = pd.DataFrame(dfs)

# %%
results = fi.fixleg(myfl, val_curve=rate)
#print("FIX LEG STRUCTURE")
#print("=================")
#print(results["structure"])

print("FIX LEG RIKS")
print("============")
print(results["risks"])

pd_struc = pd.DataFrame(results["structure"])


# %%
