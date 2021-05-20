#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:33:15 2020

@author: RMS671214
"""


convs = ('Actual/365', 'Actual/365 Fixed', 'Actual/365L', 'Actual/364',
         'Actual/360', 'Eurobond basis (ISDA 2006)', '30A/360', '30/360 US',
         '30/360 ICMA', 'Eurobond basis (ISDA 2000)', 'Actual/Actual',
         'Actual/Actual ISDA', 'Actual/Actual AFB')  # 'Actual/Actual ICMA'
day_counts = convs
all_cons = []
for conv in convs:
    tup_conv = (conv, conv)
    all_cons.append(tup_conv)


business_days = ("No Adjustment", "Following", "Modified Following",
                 "Preceeding", "Modified Preceeding", "End Of Month")

date_gen_method = ("Forward from issue date", "Backward from maturity date")
frequencies = {"Annual": 12, "Semi-Annual": 6, "Quarterly": 3, "Monthly": 1}
start_basis = {'Spot': 2, 'Tom': 1, 'Same Day': 0}
