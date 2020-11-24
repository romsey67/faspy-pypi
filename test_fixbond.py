#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from faspy.interestrate.fixincome import fixbond_structures, fixbond_value, \
    date_structures
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
mybond['ytm'] = None
mybond['type'] = 'Fixed Rate Bond'

structures = list(fixbond_structures(mybond))
print("STRUCTURES")
print("=============")
print(structures)

try:
    import pandas as pd
    pd1 = pd.DataFrame(structures)
except:
    pass

# %%
print("TESTING fixbond_value")
print("=======================")
print(structures[-1]["end_date"])
structures = fixbond_value(mybond["value_date"], structures, 12,
              mybond["day_count"], mybond["frequency"])
print(structures)
try:
    import pandas as pd
    pd_data = pd.DataFrame(structures)
except:
    pass
value = 0
accrued = 0
for structure in structures:
    value += structure["pv"]
    accrued += structure["accrued"]

print("value: ", value)
print("accrued: ", accrued)

try:
    import pandas as pd
    pd_newdata = pd.DataFrame(structures)
except:
    pass


# %%

coupon_dates = list(date_structures(mybond))
print(coupon_dates)
