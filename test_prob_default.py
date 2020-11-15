#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:34:05 2020

@author: RMS671214
"""
from faspy.credit.credit import calc_prob_default, prob_default_interpolation as pdi
import pandas as pd

spreads = []
spreads.append({"tenor": 1, "spread": 50, "rec_rate": 0.4})
spreads.append({"tenor": 2, "spread": 70, "rec_rate": 0.4})
spreads.append({"tenor": 3, "spread": 100, "rec_rate": 0.4})
spreads.append({"tenor": 5, "spread": 150, "rec_rate": 0.4})
spreads.append({"tenor": 10, "spread": 200, "rec_rate": 0.4})
spreads.append({"tenor": 20, "spread": 350, "rec_rate": 0.4})

defaults = calc_prob_default(spreads)
pd_def = pd.DataFrame(defaults)

tenors = list(range(20))
interp = pdi(defaults, tenors)
pd_interp = pd.DataFrame(interp)