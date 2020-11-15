#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from faspy.interestrate.cashflows import fixbond_structures, fixbond_value, \
    coupon_structures
import numpy as np
from faspy.interestrate import rmp_dates as rd
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
mybond['coupon'] = 10
mybond['ytm'] = 12.00
mybond['type'] = 'Fixed Rate Bond'

structures = list(fixbond_structures(mybond))
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