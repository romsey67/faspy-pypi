#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:05:24 2020

@author: RMS671214
"""

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