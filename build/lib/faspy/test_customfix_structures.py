#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:02:58 2020

@author: RMS671214
"""

from faspy.interestrate.fixincome import date_structures, calc_customfix_structures, \
    value_customfix_structures
import numpy as np
from faspy.interestrate import rmp_dates as rd
from faspy.interestrate import discount_curve as dcurve
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

structures = list(date_structures(mybond))
print("DATES")
print("=====")
print(structures)

try:
    import pandas as pd
    pd1 = pd.DataFrame(structures)
except:
    pass


# %%
# Create Sample structure

lenstruc = len(structures)
for i in range(lenstruc):
    structure = structures[i]
    if i == 0:
        structure["face_value"] = 10_000_000
        structure["fv_flow"] = 1_000_000
        structure["coupon"] = 5.00
    else:
        prevstruc = structures[i-1]
        structure["face_value"] = prevstruc["face_value"] + prevstruc["fv_flow"]
        structure["fv_flow"] = 1_000_000
        structure["coupon"] = prevstruc["coupon"] + 0.25

structures[-1]["fv_flow"] = structures[-1]["face_value"]
new_structures = calc_customfix_structures(structures,
                                           mybond["day_count"],
                                           mybond["frequency"],
                                           mybond["business_day"])
print("CUSTOM STRUCTURE")
print("================")
print(new_structures)
try:
    import pandas as pd
    pd2 = pd.DataFrame(new_structures)
except:
    pass

# %%
# valuing the structures using flat curve

flatcurve = dcurve.flat_curve(mybond['value_date'], mybond['value_date'] + 50,
                              3.5)
print(flatcurve)

newstruct = value_customfix_structures(mybond["value_date"], new_structures,
                                       mybond["day_count"],
                                       mybond["frequency"], flatcurve)
print("VALUING CUSTOM STRUCTURE")
print("========================")
print(newstruct)
try:
    import pandas as pd
    pd3 = pd.DataFrame(newstruct)
except:
    pass


# %%

coupon_dates = list(date_structures(mybond))
print(coupon_dates)
