#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:49:32 2020

@author: RMS671214
"""
from faspy.interestrate.rmp_curves import *
from faspy.interestrate.fas_ircls import STRates, LTRates, Rates


# %%

def discount_factor_gen(rates, return_type="time"):
    """
    Generate the discount factors.

            Parameters:

                rates: a dictionary with the following keys - value_date,
                st_busday, st_ratebasis, st_daycount, lt_busday,
                lt_frequency, lt_daycount and rates. 'rates' key in
                the dictionary is a dictionary of interest rate in
                percentage having the following keys - 1W,  2W, 3W, 1M,
                2M, 3M, 4M, 5M, 6M, 9M, 12M, 1Y, 2Y, 3Y, 4Y, 5Y, 6Y, 7Y,
                10Y, 15Y, 20Y, 30Y.

            Returns:

                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """
    if isinstance(rates, dict):
        dfs = _discount_factor_generate(rates)
        if return_type == "days":
            df = [{"days": x["days"], "df": x["df"]} for x in dfs]
            df.insert(0, {"days": 0.00, "df": 1.00})
        else:
            df = [{"times": x["time"], "df": x["df"]} for x in dfs]
            df.insert(0, {"times": 0.00, "df": 1.00})
        results = df

    elif isinstance(rates, list):
        results = []
        for rate in rates:
            if isinstance(rate, dict):
                dfs = _discount_factor_generate(rate)
                if return_type == "days":
                    df = [{"days": x["days"], "df": x["df"]} for x in dfs]
                    df.insert(0, {"days": 0.00, "df": 1.00})
                else:
                    df = [{"times": x["time"], "df": x["df"]} for x in dfs]
                    df.insert(0, {"times": 0.00, "df": 1.00})
                results.append(df)

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
