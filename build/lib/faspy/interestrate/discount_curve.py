#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:49:32 2020

@author: RMS671214
"""
from faspy.interestrate.rmp_curves import *
from faspy.interestrate.fas_ircls import STRates, LTRates, Rates


# %%

def discount_factor_gen(rates):
    # rates must contain is an array of dictionaries
    # eg [{value_date:, st_busday:, st_ratebasis:, st_daycount:, lt_busday:,
    # lt_frequency:, lt_daycount:, rates:, {}}]
    # rates key in the dictionary is a dictionary of rates
    # eg  {'1W': None, '2W': None, '3W': None, '1M': None, '2M': None,
    # '3M': None, '4M': None, '5M': None, '6M': None, '9M': None,
    # '12M': None}'1Y': None, '2Y': None, '3Y': None, '4Y': None, '5Y': None,
    # '6Y': None, '7Y': None, '10Y': None, '15Y': None, '20Y': None,
    # '30Y': None}

    if isinstance(rates, dict):
        return _discount_factor_generate(rates)

    elif isinstance(rates, list):
        results = []
        for rate in rates:
            if isinstance(rate, dict):
                results.append(_discount_factor_generate(rate))

    return results



def _discount_factor_generate(rate):
    strates = STRates()
    strates.rate_basis = rate.get("st_ratebasis")
    strates.business_day = rate.get("st_busday")
    strates.day_count = rate.get("st_daycount")
    strates.rates = rate["rates"]

    ltrates = LTRates()
    ltrates.frequency = rate.get("lt_frequency")
    ltrates.business_day = rate.get("lt_busday")
    ltrates.day_count = rate.get("lt_daycount")
    ltrates.rates = rate["rates"]

    myrates = Rates()
    myrates.strates = strates
    myrates.ltrates = ltrates
    myrates.date_gen_method = 'Forward from issue date'
    myrates.value_date = rate.get("value_date")
    myrates.calcdf()
    return myrates.df.data
