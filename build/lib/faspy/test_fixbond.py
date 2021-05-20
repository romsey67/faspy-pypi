#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from interestrate.fixincome import fixbond_structures, fixbond_value, \
    date_structures, fixbond
import numpy as np
from numpy import datetime64 as dt64
from interestrate import rmp_dates as rd
# %%

mybond = {}
mybond['issue_date'] = np.datetime64('2018-10-22')
mybond['value_date'] = np.datetime64('2021-10-22')
mybond['maturity'] = np.datetime64('2128-10-22')
mybond['day_count'] = 'Actual/365 Fixed'
mybond['frequency'] = 'Semi-Annual'
mybond['business_day'] = 'No Adjustment'
mybond['date_generation'] = rd.date_gen_method[1]
mybond['face_value'] = 1000000
mybond['coupon'] = 10.00
mybond['ytm'] = 10.00
mybond['type'] = 'Fixed Rate Bond'

structures = list(fixbond_structures(mybond))
print("STRUCTURES")
print("===========")
print(structures)

try:
    import pandas as pd
    pd1 = pd.DataFrame(structures)
except:
    pass


# %%

testdata = {"value_date":dt64("2020-05-01"), "maturity": dt64("2025-10-01"),
            "day_count": "Actual/Actual", "frequency":"Semi-Annual",
            "business_day": "No Adjustment",
            "date_generation": "Backward from maturity date",
            "face_value": 10000000, "coupon": 2.00, "ytm": 2.00}
val = fixbond(testdata)
print(val["risks"])


try:
    pd3 = pd.DataFrame(val["structure"])

except:
    import pandas as pd
    pd3 = pd.DataFrame(val["structure"])
    
    
