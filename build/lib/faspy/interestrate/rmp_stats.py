#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 14:11:12 2020

@author: RMS671214
"""
import math
import numba
from time import time

@numba.njit('float64(float64)')
def rmpcdf(x):
    return 0.5 * (1 + math.erf(x/(2**0.5)))

@numba.njit('float64(float64)')
def rmpppf(x):
    cd = 0.75
    value0 = rmpcdf(cd)
    value1 = rmpcdf(cd + 0.000000001)
    der = (value1 - value0) / 0.000000001
    counter = 0
    while counter <500:
        cd = cd - (value0 - x)/der
        value0 = rmpcdf(cd)
        value1 = rmpcdf(cd + 0.000000001)
        if abs(value0-x) < 0.00000000000000000000001:
            counter = 500

        der = (value1 - value0) / 0.000000001
        counter += 1

    return cd
