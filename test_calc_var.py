#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 07:32:58 2020

@author: RMS671214
"""
from faspy.interestrate import discount_curve as dcurve
from numpy import datetime64 as dt64
vdate = dt64("2020-10-30", "D")

# %%
# Historical Discount Factors
rates = []
rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}
rates.append(rate)

rate = {"value_date": "2020-10-29", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.31, '1W': 2.34, '1M': 2.44, '3M': 2.56,
                  '6M': 2.67, '1Y': 2.75, '2Y': 2.85, '3Y': 2.95,
                  '5Y': 3.04, '10Y': 3.11, '30Y': 3.20}}
rates.append(rate)

rate = {"value_date": "2020-10-28", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.40, '1W': 2.44, '1M': 2.54, '3M': 2.66,
                  '6M': 2.77, '1Y': 2.85, '2Y': 2.95, '3Y': 3.05,
                  '5Y': 3.14, '10Y': 3.21, '30Y': 3.30}}
rates.append(rate)

rate = {"value_date": "2020-10-27", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.38, '1W': 2.44, '1M': 2.55, '3M': 2.76,
                  '6M': 2.87, '1Y': 3.05, '2Y': 3.15, '3Y': 3.25,
                  '5Y': 3.34, '10Y': 3.41, '30Y': 3.50}}
rates.append(rate)

rate = {"value_date": "2020-10-27", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.35, '1W': 2.40, '1M': 2.50, '3M': 2.65,
                  '6M': 2.75, '1Y': 2.85, '2Y': 2.95, '3Y': 3.05,
                  '5Y': 3.15, '10Y': 3.25, '30Y': 3.35}}
rates.append(rate)

dfs = dcurve.discount_factor_gen(rates, return_type="time")
# print(dfs)

# %%
# Define VaR Vertices using tenors
# convert the tenors to date then to no of days
from faspy.interestrate.rmp_dates import tenor_to_maturity as ttm, \
    day_count_factor as day_cf
var_vertices = ["6M", "1Y", "3Y", "5Y", "10Y", "20Y"]
var_dates = list(map(lambda tenor: ttm(vdate, tenor,
                                       business_day="Modified Following"),
                     var_vertices))
var_days = list(map(lambda date: (date - vdate).astype("int"), var_dates))
var_time = list(map(lambda date: day_cf("Actual/365", vdate, date),
                    var_dates))

# %%
# get historical discount factors
hdf = []
for df in dfs:
    x_axis = [x["times"] for x in df]
    y_axis = [y["df"] for y in df]
    idf = dcurve.interpolation(x_axis, y_axis, var_time)
    hdf.append(idf)
#print(hdf)

# %%
# 1.Calculate the log normal return
# 2. Calculate volatility
# 3. Calculate Correlation
from faspy.risk import var_utils as vu
from numpy import array as arr

returns = vu.returns_from_prices(arr(hdf))
vols = vu.volatilities(returns)
corr = vu.correlation(returns)
# print(returns, vols, corr)

# %%
# Bond Exposure
from faspy.interestrate.cashflows import fixbond_structures, fixbond_value
import numpy as np
from faspy.interestrate import rmp_dates as rd


mybond = {}
# mybond['issue_date'] = np.datetime64('2018-10-22')
mybond['value_date'] = np.datetime64('2020-10-22')
mybond['maturity'] = np.datetime64('2028-10-22')
mybond['day_count'] = 'Actual/365 Fixed'
mybond['frequency'] = 'Semi-Annual'
mybond['business_day'] = 'No Adjustment'
mybond['date_generation'] = rd.date_gen_method[1]
mybond['face_value'] = 1000000
mybond['coupon'] = 10
# mybond['ytm'] = 12.00
mybond['type'] = 'Fixed Rate Bond'

structures = list(fixbond_structures(mybond))
#print(structures)

# %%
# Reformat bond cash flows array

latest_df = dfs[0]
x_axis = [x["times"] for x in latest_df]
y_axis = [y["df"] for y in latest_df]
idf = dcurve.interpolation(x_axis, y_axis, var_days, is_function=True)
cf = []
for structure in structures:
    datum = {}
    if structure["end_date"] > vdate:
        datum["times"] = day_cf("Actual/365", vdate, structure["end_date"])
        datum["df"] = idf(datum["times"]) 
        datum["pv"] = datum["df"] * structure["cash_flow"]
        cf.append(datum)

# %%
# mapped bond cash flow to var vertices
vdf = hdf[0]
mcf = vu.map_cf_to_var_vertices(cf, vdf, var_time, corr, vols)
confidence = 0.99
# weigh the asset
weighted = vu.var_asset_weightbyprice(mcf, confidence, vols)
# calculate value at risk
va_risk = vu.value_at_risk(weighted, corr)

