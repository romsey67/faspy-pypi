#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from interestrate.fixincome import floatbond_value, \
    date_structures, floatbond
import numpy as np
from interestrate import rmp_dates as rd
from interestrate import discount_curve as dcurve
from interestrate import rmp_curves as rcurve
import pandas as pd

# %%
# create the float structure

mybond = {}
mybond['issue_date'] = np.datetime64('2018-10-22')
mybond['value_date'] = np.datetime64('2018-10-22')
mybond['maturity'] = np.datetime64('2020-10-22')
mybond['day_count'] = 'Actual/365 Fixed'
mybond['frequency'] = 'Semi-Annual'
mybond['business_day'] = 'No Adjustment'
mybond['date_generation'] = rd.date_gen_method[1]
mybond['face_value'] = 100
mybond["current_coupon"] = 3.65
mybond['spread'] = 1.5
mybond["margin"] = 1
mybond["fixing_basis"] = "Same Day"

# %%

rate = {"value_date": str(mybond['value_date']),
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
# Test Date Structure
coupon_dates = date_structures(mybond)

# %%

zcurve = rcurve.discount_factor_from_zspread(mybond['value_date'],
                                             coupon_dates,
                                             mybond['day_count'],
                                             mybond['frequency'],
                                             dfs, mybond['spread'])
pd_z = pd.DataFrame(zcurve)



# %%
fbond = floatbond(mybond, rate)
pd0 = pd.DataFrame(fbond["structure"])
risks = fbond["risks"]
print(fbond["risks"])

#%%

# Using Valuation Curve
val_curve = dict(rate)
val_curve["rates"] = rcurve.shift_curve(val_curve["rates"], bp=mybond["spread"])
fbond2 = floatbond(mybond, rate, val_curve=val_curve)
risks = fbond2["risks"]
print(fbond2["risks"])




