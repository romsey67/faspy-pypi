#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from faspy.interestrate.cashflows import loans_structures, fixbond_value, \
    coupon_structures
import numpy as np
from faspy.interestrate import rmp_dates as rd
# %%

myloan = {}
myloan['issue_date'] = np.datetime64('2018-10-22')
myloan['value_date'] = np.datetime64('2021-10-22')
myloan['maturity'] = np.datetime64('2028-10-22')
myloan['day_count'] = 'Actual/365 Fixed'
myloan['frequency'] = 'Semi-Annual'
myloan['business_day'] = 'No Adjustment'
myloan['date_generation'] = rd.date_gen_method[1]
myloan['face_value'] = 1000000
myloan['coupon'] = 10
myloan['ytm'] = 12.00
myloan['type'] = 'Fixed Rate Bond'
myloan["rate_type"] = "fixed"

structures = list(loans_structures(myloan))
print(structures)

# %%
fixbond_value(mybond["value_date"], structures, mybond["ytm"],
              mybond["day_count"], mybond["frequency"])
print(structures)
value = 0
accrued = 0
for structure in structures:
    value += structure["pv"]
    accrued += structure["accrued"]
    
print("value: ", value)
print("accrued: ", accrued)

# %%

coupon_dates = list(coupon_structures(mybond))
print(coupon_dates)

#%%
a = ("a", "b")
b = list(a)