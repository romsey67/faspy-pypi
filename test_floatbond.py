#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from faspy.interestrate.cashflows import floatbond_value, \
    coupon_structures, floatbond_structures, floatbond_value2
import numpy as np
from faspy.interestrate import rmp_dates as rd
from faspy.interestrate import discount_curve as dcurve
import pandas as pd

# %%

mybond = {}
mybond['issue_date'] = np.datetime64('2018-10-22')
mybond['value_date'] = np.datetime64('2021-10-22')
mybond['maturity'] = np.datetime64('2028-10-22')
mybond['day_count'] = 'Actual/365 Fixed'
mybond['frequency'] = 'Semi-Annual'
mybond['business_day'] = 'No Adjustment'
mybond['date_generation'] = rd.date_gen_method[1]
mybond['face_value'] = 1000000
mybond["current_coupon"] = 3.40
mybond['spread'] = 0.75
mybond["margin"] = 0.75
mybond["fixing_basis"] = "Same Day"

structures = list(floatbond_structures(mybond))
pd1 = pd.DataFrame(structures)


# %%
rate = {"value_date": str(mybond["value_date"]), "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/365",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}

dfs = dcurve.discount_factor_gen(rate, return_type="time")

# %%
new_structures = floatbond_value(mybond["value_date"], list(structures),
                                 mybond["spread"], mybond["day_count"], dfs=dfs)
pd2 = pd.DataFrame(new_structures)


# %%
value = 0
accrued = new_structures[0]["accrued"]
for structure in new_structures:
    value += structure["pv"]


# %%
spread = 0.75
m_rate = dict(rate)
rates = m_rate["rates"]
for key in rates:
    rates[key] = rates[key] + spread

print(m_rate)

m_dfs = dcurve.discount_factor_gen(m_rate, return_type="time")

# %%
new_structures2 = floatbond_value2(mybond["value_date"], list(structures),
                                   mybond["day_count"],
                                   ref_df=dfs, market_df=m_dfs)
pd3 = pd.DataFrame(new_structures2)
value2 = 0
accrued = new_structures2[0]["accrued"]
for structure in new_structures2:
    value2 += structure["pv"]