#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:36:21 2020

@author: RMS671214
"""

from faspy.interestrate import discount_curve as dcurve
from faspy.interestrate import rmp_curves as rcurve

# %%

rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}


dfs = dcurve.discount_factor_gen(rate, return_type="times")


# %%

x_axis = [x["times"] for x in dfs]
y_axis = [x["df"] for x in dfs]

start = 0.5
end = 0.75
# use this to calculate a fwd discount factor
fwd1 = rcurve.calc_fwd_df(start, end, x_axis, y_axis)

ifunc = rcurve.interpolation(x_axis, y_axis, start, is_function=True)
print(type(ifunc))
# use this to calculate a long list of forward discount factor
fwd2 = rcurve.calc_fwd_df(start, end, ifunc=ifunc)

# e.g of using list
fwd_period = [[0.1, 0.2], [0.2, 0.3], [0.4, 0.7], [0.7, 0.9]]
fwd_dfs = [rcurve.calc_fwd_df(x[0], x[1], ifunc=ifunc) for x in fwd_period]
