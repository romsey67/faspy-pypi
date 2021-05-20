#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:34:05 2020

@author: RMS671214
"""
from faspy.basel.credit.credit import calc_prob_default, \
    prob_default_interpolation as pdi, expected_loss, \
    mean_rr_with_betadist
        
import pandas as pd

# %%

spreads = []
spreads.append({"tenor": 1, "spread": 50, "rec_rate": 0.4})
spreads.append({"tenor": 2, "spread": 70, "rec_rate": 0.4})
spreads.append({"tenor": 3, "spread": 100, "rec_rate": 0.4})
spreads.append({"tenor": 5, "spread": 150, "rec_rate": 0.4})
spreads.append({"tenor": 10, "spread": 200, "rec_rate": 0.4})
spreads.append({"tenor": 20, "spread": 350, "rec_rate": 0.4})

print(spreads)
defaults = calc_prob_default(spreads)
pd_def = pd.DataFrame(defaults)

tenors = list(range(20))
print(tenors)
interp = pdi(defaults, tenors)
pd_interp = pd.DataFrame(interp)

# %%
# Test expected Loss

prob_default = 0.5
exposure_at_default = 1_000_000
rec_rate= 0.40
loss_given_default = None

eloss = expected_loss(prob_default, exposure_at_default, recovery_rate=rec_rate,
                      loss_given_default=loss_given_default)
print(eloss)

# %%
rr = mean_rr_with_betadist(0.5, 1)