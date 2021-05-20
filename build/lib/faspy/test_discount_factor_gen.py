#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:36:21 2020

@author: RMS671214
"""

from interestrate import discount_curve as dcurve

# %%
rates = []
rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
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
                  '6M': None, '12M': 3.75, '1Y': None, '2Y': 3.80, '3Y': 3.90,
                  '5Y': 4.00, '10Y': 4.10, '30Y': 4.25}}
rates.append(rate)
print(rates)

dfs = dcurve.discount_factor_gen(rate, return_type="time")
print(dfs)

# %%

flat = dcurve.flat_curve("2020-10-30", "2020-12-30", 5.30, rate_basis="Money Market",
                         day_count="Actual/365", bus_day="No Adjustment",
                         tenors=None, return_type="days")