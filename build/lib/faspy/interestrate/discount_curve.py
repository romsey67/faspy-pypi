#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:49:32 2020

@author: RMS671214
"""
from .rmp_curves import *
from .fas_ircls import STRates, LTRates, Rates
from .conventions import *
from .rmp_dates import tenor_to_maturity as ttm, generate_dates as gen_dates
from numpy import datetime64 as dt64
import math
from collections import deque

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


def flat_curve(start_date, end_date, rate, rate_basis="Money Market",
               day_count="Actual/365", bus_day="No Adjustment",
               tenors=None, return_type="time"):
    if tenors is None:
        mytenors = [x for x in std_tenors if x != "12M"]
    else:
        mytenors = list(tenors)
    sdate = dt64(start_date, "D")
    edate = dt64(end_date, "D")
    time = day_cf("Actual/365", sdate, edate)
    dcf = day_cf(day_count, sdate, edate)
    if rate_basis == "Money Market":
        df = _mmr2df(float(rate), float(dcf))
    elif rate_basis == "Discount Rate":
        df = _dr2df(float(rate), float(dcf))
    else:
        return None

    time = day_cf("Actual/365", sdate, edate)
    crate = -math.log(df)/time
    dates = [ttm(sdate, tenor, day_count, bus_day) for tenor in mytenors]

    times = [{"times": float(day_cf("Actual/365", sdate, date)) } for date in dates]
    timelen = len(times)
    ret_array=[]
    if return_type == "time":
        dfs = [{"df": math.exp(-crate * time["times"])} for time in times]
        for i in range(timelen):
            time = times[i]
            df = dfs[i]
            time.update(df)
            ret_array.append(time)

    elif return_type == "days":
        days = [{"days": (date-sdate).astype("float") } for date in dates]
        dfs = [{"df": math.exp(-crate * time["times"])} for time in times]
        for i in range(timelen):
            day = days[i]
            df = dfs[i]
            day.update(df)
            ret_array.append(day)
    return ret_array


def discount_factor_from_ytm(value_date, maturity, day_count, frequency,
                             business_day, ytm):
    # dates is a deque
    dates = gen_dates(value_date, maturity, issueDate=value_date,
                      frequency=frequency, business_day=business_day,
                      method=date_gen_method[1])

    start_dates = deque(dates)
    end_dates = deque(dates)
    start_dates.pop()
    end_dates.popleft()
    noofcpns = len(start_dates)
    data = deque()
    for no in range(noofcpns):
        structure = {"start_date": start_dates[no], "end_date": end_dates[no]}
        data.append(structure)

    df_curve =  discount_factor_from_ytm(value_date, list(data), day_count,
                                        frequency, business_day, ytm)
    return df_curve


def discount_factor_from_ytm_using_structures(value_date, date_structure,
                                              day_count,
                                              frequency, business_day, ytm):

    df = 1
    maturity = date_structure[-1]["end_date"]

    data = [{"start_date": x["start_date"], "end_date": x["end_date"]}
            for x in date_structure if value_date < dt64(x["end_date"])]

    for datum in data:
        if dt64(datum["start_date"]) >= value_date:
            datum["dcf"] = day_cf(day_count, datum["start_date"], datum["end_date"],
                                  bondmat_date=maturity,
                                  next_coupon_date=datum["end_date"],
                                  business_day=business_day,
                                  Frequency=12/frequencies[frequency])
            datum["time"] = day_cf("Actual/365", value_date, datum["end_date"],
                                  bondmat_date=maturity,
                                  next_coupon_date=datum["end_date"],
                                  business_day=business_day,
                                  Frequency=12/frequencies[frequency])

        elif datum["start_date"] < value_date:
            datum["dcf"] = day_cf(day_count, value_date, datum["end_date"],
                                  bondmat_date=maturity,
                                  next_coupon_date=datum["end_date"],
                                  business_day=business_day,
                                  Frequency=12/frequencies[frequency])
            datum["time"] = day_cf("Actual/365", value_date, datum["end_date"],
                                  bondmat_date=maturity,
                                  next_coupon_date=datum["end_date"],
                                  business_day=business_day,
                                  Frequency=12/frequencies[frequency])

        datum["period_df"] = 1/(1 + ytm * datum["dcf"] / 100)
        df = datum["df"] = df * datum["period_df"]

    disfac = list(map(lambda datum: {"times": datum["time"], "df": datum["df"]},
                      data))
    disfac.insert(0,{"times":0, "df": 1})
    return disfac


@numba.njit('float64(float64, float64)')
def _df2mmr(df, dcf):
    rate = ((1 / df) - 1) / (0.01 * dcf)
    return rate


@numba.njit('float64(float64, float64)')
def _mmr2df(rate, dcf):
    df = 1 / (1 + rate * 0.01 * dcf)
    return df


@numba.njit('float64(float64, float64)')
def _df2dr(df, dcf):
    rate = (1 - df) / (0.01 * df)
    return rate


@numba.njit('float64(float64, float64)')
def _dr2df(rate, dcf):
    df = 1 - rate * 0.01 * dcf
    return df
