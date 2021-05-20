#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from interestrate.fixincome import date_structures, calc_customfix_structures, \
    value_customfix_structures, create_structures_from_dates
import numpy as np
from interestrate import rmp_dates as rd
from interestrate import discount_curve as dcurve
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

dates = list(date_structures(mybond))
print("DATES")
print("=====")
print(dates)

for date in dates:
    date["start_date"] = str(date["start_date"])
    date["end_date"] = str(date["end_date"])
    
print("DATES AS STRING")
print("===============")
print(dates)

try:
    import pandas as pd
    pd1 = pd.DataFrame(dates)
except:
    pass


# %%
# Create Sample structure

lenstruc = len(dates)
face_values = [10_000_000.00  + i * 100_000.00 for i in range(lenstruc)]
print("FACE VALUES")
print("===========")
print(face_values)
coupons = [3.5  + i * 0.20 for i in range(lenstruc)]
print("COUPONS")
print("=======")
print(coupons)

fv_flows = [i * 100_000.00 for i in range(lenstruc)]
print("FV FLOWS")
print("========")
print(fv_flows)
fv_flows[-1] = -face_values[-1]

new_structures = create_structures_from_dates(dates, coupons, face_values, fv_flows)
print("CUSTOM STRUCTURE")
print("================")
print(new_structures)
try:
    import pandas as pd
    pd2 = pd.DataFrame(new_structures)
except:
    pass

# %% value structure using YTM

astructures = value_customfix_structures(mybond["value_date"], new_structures,
                                         mybond["day_count"],
                                         mybond["frequency"],
                                         5.00)

print("VALUATION OF CUSTOM STRUCTURE")
print("=============================")
print(astructures)
try:
    import pandas as pd
    pd3 = pd.DataFrame(astructures)
except:
    pass

# %% value struture using disount factors

rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}

dfs = dcurve.discount_factor_gen(rate, return_type="time")

astructures1 = value_customfix_structures(mybond["value_date"], new_structures,
                                         mybond["day_count"],
                                         mybond["frequency"],
                                         dfs)

print("VALUATION OF CUSTOM STRUCTURE")
print("=============================")
print(astructures)
try:
    import pandas as pd
    pd4 = pd.DataFrame(astructures1)
except:
    pass