#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 09:57:01 2020

@author: RMS671214
"""

from numpy import datetime64 as dt64
from .rmp_dates import generate_dates as gen_dates, \
    frequencies as fre,  day_count_factor as day_cf
from .conventions import start_basis, frequencies
from .rmp_curves import interpolation, calc_shortrate_from_df, \
calc_df_from_shortrate, calc_fwd_df, shift_curve, discount_factor_from_zspread
from .discount_curve import discount_factor_from_ytm as ytm_df, \
    discount_factor_from_ytm_using_structures as ytm_df_struct, \
    discount_factor_gen as df_gen
from collections import deque
import numpy as np
import math


def fixbond(bond):

    # create structure for the bond
    structure = fixbond_structures(bond)
    risks = fixbond_value(bond["value_date"], structure, bond["ytm"],
                          bond["day_count"], bond["frequency"],
                          bond["business_day"])

    return {"structure": structure, "risks": risks}


def fixbond_value(value_date, structures, yld, day_count, frequency,
                  business_day="No Adjustment"):
    const = 0.01
    try:
        ytm = float(yld)
        ytm1 = ytm + const
        ytm2 = ytm1 + const
    except Exception:
        return None
    newstructures = [dict(x) for x in structures]
    maturity = newstructures[-1]["end_date"]

    # Curves to be used in calculation of duration, convexity and pvbp01
    df_curve = ytm_df_struct(value_date, newstructures, day_count, frequency,
                             business_day, ytm)
    x_axis = [x["times"] for x in df_curve]
    y_axis = [x["df"] for x in df_curve]
    ifunc = interpolation(x_axis, y_axis, float(1/366), is_function=True)

    df_curve1 = ytm_df_struct(value_date, newstructures, day_count, frequency,
                              business_day, ytm1)
    x_axis1 = [x["times"] for x in df_curve1]
    y_axis1 = [x["df"] for x in df_curve1]
    ifunc1 = interpolation(x_axis1, y_axis1, float(1/366), is_function=True)

    df_curve2 = ytm_df_struct(value_date, newstructures, day_count, frequency,
                              business_day, ytm2)
    x_axis2 = [x["times"] for x in df_curve2]
    y_axis2 = [x["df"] for x in df_curve2]
    ifunc2 = interpolation(x_axis2, y_axis2, float(1/366), is_function=True)

    # covert data into list
    dates = [{"start_date": x["start_date"], "end_date": x["end_date"]}
             for x in structures if value_date < dt64(x["end_date"])]
    cfs = np.asarray([x["cash_flow"] for x in structures
                      if value_date < dt64(x["end_date"])])

    times = np.asarray([day_cf("Actual/365", value_date, x["end_date"])
                        for x in dates])
    # interpolating discount factors
    dfs = [ifunc(x) for x in times]
    # calculating the present values
    pvs = cfs * dfs
    weighted_pvs = pvs * times
    value = np.sum(pvs)
    weighted_value = np.sum(weighted_pvs)

    # interpolating 2nd discount factors
    dfs1 = [ifunc1(x) for x in times]
    # calculating the 2nd present values
    pvs1 = cfs * dfs1
    value1 = np.sum(pvs1)

    # interpolating 3rd discount factors
    dfs2 = [ifunc2(x) for x in times]
    # calculating the 3rd present values
    pvs2 = cfs * dfs2
    value2 = np.sum(pvs2)

    mac_dur = weighted_value/value
    compound = 12 / frequencies[frequency]
    mod_dur = mac_dur / pow((1 + ytm / (compound * 100)), compound)
    pvbp01 = value1 - value

    # 1st derivative at ytm using forward differential
    der0 = (value1 - value) / (ytm1 - ytm)
    # 1st derivative at ytm1 using forward differential
    der1 = (value2 - value1) / (ytm2 - ytm1)
    # 2nd derivative at ytm using forward differential
    conv = (der1 - der0) / (ytm2-ytm1)

    return {"macaulay_duration": mac_dur, "modified_duration": mod_dur,
            "pvbp01": pvbp01, "convexity": conv, "value": value}


def fixbond_structures(bond):
    """
    Generate bonds coupon structures.
            Parameters:
                bonds: a dictionary with the following keys - value_date,
                business_day, issue_date, value_date, maturity,
                frequency, day_count, date_generation, face_value, coupon,
                ytm and type.
            Returns:
                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """
    results = list(_fixbond_gen_structure(bond))

    return results


def _fixbond_gen_structure(bond):

    dates = _dates_gen_structure(bond)
    dates = deque(dates)

    noofcpns = len(dates)
    for no in range(noofcpns):
        structure = {}
        structure["cpn_dcf"] = day_cf(bond["day_count"],
                                      dates[no]["start_date"],
                                      dates[no]["end_date"],
                                      bondmat_date=bond["maturity"],
                                      next_coupon_date=dates[no]["end_date"],
                                      business_day=bond["business_day"],
                                      Frequency=12/frequencies[bond["frequency"]])
        structure["coupon"] = bond["coupon"]
        structure["face_value"] = bond["face_value"]
        structure["coupon_interest"] = (structure["cpn_dcf"] *
                                        structure["coupon"] *
                                        structure["face_value"] / 100)
        if no != noofcpns - 1:
            structure["fv_flow"] = 0
        else:
            structure["fv_flow"] = structure["face_value"]
        structure["cash_flow"] = (structure["coupon_interest"] +
                                  structure["fv_flow"])

        dates[no].update(structure) # merging the dictionary

    if (bond.get("ytm") and bond["value_date"]):
        fixbond_value(bond["value_date"], dates, bond["ytm"], bond["day_count"],
                      bond["frequency"])

    return dates


def date_structures(bonds):
    """
    Generate coupon structures.

            Parameters:
                bonds: a dictionary with the following keys - value_date,
                business_day, issue_date, value_date, maturity,
                frequency, day_count, date_generation, face_value, coupon,
                ytm and type.

            Returns:
                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """
    if isinstance(bonds, dict):
        results = _dates_gen_structure(bonds)

    elif isinstance(bonds, list):
        results = []
        for bond in bonds:
            if isinstance(bond, dict):
                results.append(_dates_gen_structure(bond))
    return results


def _dates_gen_structure(bond):
    """
    Generate date structures. Function can be used for all coupon bearing
    products with bullet principal repayment

            Parameters:
                bonds: a dictionary with the following keys - value_date,
                business_day, issue_date, value_date, maturity,
                frequency, day_count, date_generation.

            Returns:
                a dictionary or an array of dictionaries with the following
                keys - "start_date" and "end_date"
    """
    bus_day = None
    if bond['business_day'] == 'NULL':
        bus_day = 'No Adjustment'
    else:
        bus_day = bond['business_day']


    start_date = bond.get('issue_date')
    if start_date is not None:
        start_date = dt64(start_date, 'D')
    value_date = bond.get('value_date')

    if value_date is not None:
        value_date = dt64(value_date, 'D')
    if start_date is not None and value_date is not None:
        use_date = start_date
    elif start_date is not None:
        use_date = start_date
    elif value_date is not None:
        use_date = value_date
    else:
        raise Exception('Both value_date and issue_date do not have any value')
    # dates is a deque
    dates = gen_dates(use_date, bond['maturity'], issueDate=start_date,
                      frequency=fre[bond['frequency']],
                      business_day=bus_day, method=bond['date_generation'])

    start_dates = deque(dates)
    end_dates = deque(dates)
    start_dates.pop()
    end_dates.popleft()
    noofcpns = len(start_dates)
    structures = deque()
    for no in range(noofcpns):
        structure = {"start_date": start_dates[no], "end_date": end_dates[no]}
        structures.append(structure)

    newstructure = list(map(lambda sdate, edate: {"start_date": sdate,
                                                  "end_date": edate},
                            start_dates, end_dates))

    return newstructure


def floatbond(bond, ref_curve, val_curve=None):
    """
    Generate bonds coupon structures.
            Parameters:
                bonds: a dictionary with the following keys - value_date,
                business_day, issue_date, value_date, maturity,
                frequency, day_count, date_generation, face_value, coupon,
                margin.
            Returns:
                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """
    # create structure for the bond
    structure = list(_floatbond_gen_structure(bond))

    risks = floatbond_value(bond["value_date"], structure, bond["spread"],
                            bond["day_count"], bond["frequency"],
                            ref_curve=ref_curve, val_curve=val_curve)

    return {"structure": structure, "risks": risks}


def _floatbond_gen_structure(bond, holidays=[]):

    # generate the face value and coupons
    dates = _dates_gen_structure(bond)
    dates = deque(dates)

    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)
    for date in dates:
        offset = -start_basis[bond["fixing_basis"]]
        date["fixing_date"] = np.busday_offset(date["start_date"], offset,
                                               roll='backward', busdaycal=bdc)
        date["face_value"] = bond["face_value"]
        date["margin"] = bond["margin"]
        date["cpn_dcf"] = day_cf(bond["day_count"],
                                 date["start_date"],
                                 date["end_date"],
                                 bondmat_date=bond["maturity"],
                                 next_coupon_date=date["end_date"],
                                 business_day=bond["business_day"],
                                 Frequency=bond["frequency"])
        # fixing date is a forward date
        if date["fixing_date"] > bond["value_date"]:
            date["is_fixed"] = False
            date["accrued"] = 0

        # fixing date is for the current coupon period
        elif (date["fixing_date"] <= bond["value_date"] and
              date["end_date"] > bond["value_date"]):
            date["is_fixed"] = True
            date["coupon"] = bond["current_coupon"]
            date["coupon_interest"] = (date["cpn_dcf"] *
                                       date["coupon"] *
                                       date["face_value"] / 100)
            if bond["value_date"] > date["start_date"]:
                acc_dcf = day_cf(bond["day_count"],
                                 date["start_date"],
                                 bond["value_date"],
                                 bondmat_date=bond["maturity"],
                                 next_coupon_date=date["end_date"],
                                 business_day=bond["business_day"],
                                 Frequency=bond["frequency"])
                date["accrued"] = acc_dcf / date["cpn_dcf"] * date["coupon_interest"]
            else:
                date["accrued"] = 0

        else:
            date["is_fixed"] = True
            date["accrued"] = 0
        date["fv_flow"] = 0
    dates[-1]["fv_flow"] = bond["face_value"]

    return dates


def floatbond_value(value_date, structures, spread, day_count, frequency,
                    ref_curve, val_curve=None):
    if val_curve is None:
        risks = _float_val_use_zspread(value_date, structures, spread,
                                       day_count, frequency, ref_curve)
    else:
        risks = _float_val_use_curve(value_date, structures, day_count,
                                     frequency, ref_curve, val_curve)

    return risks


def _float_val_use_zspread(value_date, structures, zspread, day_count,
                           frequency, ref_curve):

    # copy only active structure period
    data = [dict(x) for x in structures if value_date < x["end_date"]]
    # covert data into list
    times = np.asarray([day_cf("Actual/365", value_date, x["end_date"])
                        for x in data])

    # Original Curves
    ref_df_curve = df_gen(ref_curve, return_type="times")
    val_df_curve = discount_factor_from_zspread(value_date, structures,
                                                day_count, frequency,
                                                ref_df_curve, zspread)

    # Shift the original curves by 1 bp
    ref_curve1 = dict(ref_curve)
    ref_curve1["rates"] = shift_curve(ref_curve1["rates"])
    ref_df_curve1 = df_gen(ref_curve1, return_type="times")
    val_df_curve1 = discount_factor_from_zspread(value_date, structures,
                                                 day_count, frequency,
                                                 ref_df_curve1, zspread)

    # Shift the curves by another 1bp
    ref_curve2 = dict(ref_curve1)
    ref_curve2["rates"] = shift_curve(ref_curve2["rates"])
    ref_df_curve2 = df_gen(ref_curve2, return_type="times")
    val_df_curve2 = discount_factor_from_zspread(value_date, structures,
                                                 day_count, frequency,
                                                 ref_df_curve2, zspread)

    pvs0 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve, val_df_curve)
    value = np.sum(pvs0)

    pvs1 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve1, val_df_curve1)
    value1 = np.sum(pvs1)

    pvs2 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve2, val_df_curve2)
    value2 = np.sum(pvs2)
    pvbp01 = value1 - value
    duration = pvbp01 / 0.01
    duration2 = (value2 - value1) / 0.01
    convexity = (duration2 - duration) / 0.01

    return {"duration": duration, "pvbp01": pvbp01, "convexity": convexity,
            "value": value}


def _float_val_use_curve(value_date, structures, day_count, frequency,
                         ref_curve, val_curve):

    # copy only active structure period
    data = [dict(x) for x in structures if value_date < x["end_date"]]
    # covert data into list
    times = np.asarray([day_cf("Actual/365", value_date, x["end_date"])
                        for x in data])

    # Original Curves
    ref_df_curve = df_gen(ref_curve, return_type="times")
    val_df_curve = df_gen(val_curve, return_type="times")

    # Shift the original curves by 1 bp
    ref_curve1 = dict(ref_curve)
    ref_curve1["rates"] = shift_curve(ref_curve1["rates"])
    ref_df_curve1 = df_gen(ref_curve1, return_type="times")

    val_curve1 = dict(val_curve)
    val_curve1["rates"] = shift_curve(val_curve1["rates"])
    val_df_curve1 = df_gen(val_curve1, return_type="times")

    # Shift the curves by another 1bp
    ref_curve2 = dict(ref_curve1)
    ref_curve2["rates"] = shift_curve(ref_curve2["rates"])
    ref_df_curve2 = df_gen(ref_curve2, return_type="times")

    val_curve2 = dict(val_curve1)
    val_curve2["rates"] = shift_curve(val_curve2["rates"])
    val_df_curve2 = df_gen(val_curve2, return_type="times")

    pvs0 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve, val_df_curve)
    value = np.sum(pvs0)

    pvs1 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve1, val_df_curve1)
    value1 = np.sum(pvs1)

    pvs2 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve2, val_df_curve2)
    value2 = np.sum(pvs2)
    pvbp01 = value1 - value
    duration = pvbp01 / 0.01
    duration2 = (value2 - value1) / 0.01
    convexity = (duration2 - duration) / 0.01

    return {"duration": duration, "pvbp01": pvbp01, "convexity": convexity,
            "value": value}


def _float_recalc_loop(value_date, data, day_count, frequency,
                       times,  df_curve, val_df_curve):

    x_axis = [x["times"] for x in df_curve]
    y_axis = [x["df"] for x in df_curve]
    ifunc = interpolation(x_axis, y_axis, float(1/366), is_function=True)

    x_val = [x["times"] for x in val_df_curve]
    y_val = [x["df"] for x in val_df_curve]
    ifunc_val = interpolation(x_val, y_val, float(1/366), is_function=True)
    val_dis_curve = [ifunc_val(x) for x in times]

    for datum in data:
        end_time = day_cf("Actual/365", value_date, datum["end_date"])

        if not datum["is_fixed"]:  # coupon has not been fixed
            start_time = day_cf("Actual/365", value_date,
                                datum["start_date"])

            fwd_df = calc_fwd_df(start_time, end_time, ifunc=ifunc)
            fwd_rate = calc_shortrate_from_df(datum["start_date"],
                                              datum["end_date"],
                                              fwd_df, day_count)
            datum["coupon"] = fwd_rate + datum["margin"]
            datum["coupon_interest"] = (datum["coupon"] * datum["cpn_dcf"] *
                                        datum["face_value"] / 100)

        if datum.get("coupon_interest"):
            datum["cash_flow"] = datum["coupon_interest"] + datum["fv_flow"]
        else:
            datum["cash_flow"] = datum["fv_flow"]

    cash_flows = np.asarray([datum["cash_flow"] for datum in data])
    pvs = cash_flows * val_dis_curve

    return pvs




def fixleg(data, val_curve=None):
    structures = _fixleg_gen_structure(data)
    if not data["principal_exchange"]:
        for structure in structures:
            structure["cash_flow"] = structure["coupon_interest"]

    if val_curve is not None:
        risks = fixleg_value(data["value_date"], structures, data["day_count"],
                             data["frequency"], val_curve)
    else:
        risks = None

    return {"structure": structures, "risks": risks}


def _fixleg_gen_structure(data):

    dates = _dates_gen_structure(data)
    dates = deque(dates)

    noofcpns = len(dates)
    for no in range(noofcpns):
        structure = {}
        structure["cpn_dcf"] = day_cf(data["day_count"],
                                      dates[no]["start_date"],
                                      dates[no]["end_date"],
                                      bondmat_date=data["maturity"],
                                      next_coupon_date=dates[no]["end_date"],
                                      business_day=data["business_day"],
                                      Frequency=12/frequencies[data["frequency"]])
        structure["coupon"] = data["coupon"]
        structure["face_value"] = data["face_value"]
        structure["coupon_interest"] = (structure["cpn_dcf"] *
                                        structure["coupon"] *
                                        structure["face_value"] / 100)
        if no != noofcpns - 1:
            structure["fv_flow"] = 0
        else:
            structure["fv_flow"] = structure["face_value"]
        structure["cash_flow"] = (structure["coupon_interest"] +
                                  structure["fv_flow"])

        dates[no].update(structure) # merging the dictionary


    return dates


def fixleg_value(value_date, structures, day_count, frequency, val_curve):

    # copy only active structure period
    data = [dict(x) for x in structures if value_date < x["end_date"]]
    # covert data into list
    times = np.asarray([day_cf("Actual/365", value_date, x["end_date"])
                        for x in data])

    # Original Curves
    val_df_curve = df_gen(val_curve, return_type="times")

    # Shift the original curves by 1 bp
    val_curve1 = dict(val_curve)
    val_curve1["rates"] = shift_curve(val_curve1["rates"])
    val_df_curve1 = df_gen(val_curve1, return_type="times")

    # Shift the curves by another 1bp
    val_curve2 = dict(val_curve1)
    val_curve2["rates"] = shift_curve(val_curve2["rates"])
    val_df_curve2 = df_gen(val_curve2, return_type="times")

    pvs0 = _fixleg_recalc_loop(value_date, data, day_count, frequency,
                               times, val_df_curve)
    value = np.sum(pvs0)

    pvs1 = _fixleg_recalc_loop(value_date, data, day_count, frequency,
                               times, val_df_curve1)
    value1 = np.sum(pvs1)

    pvs2 = _fixleg_recalc_loop(value_date, data, day_count, frequency,
                               times, val_df_curve2)
    value2 = np.sum(pvs2)
    pvbp01 = value1 - value
    duration = pvbp01 / 0.01
    duration2 = (value2 - value1) / 0.01
    convexity = (duration2 - duration) / 0.01

    return {"duration": duration, "pvbp01": pvbp01, "convexity": convexity,
            "value": value}


def _fixleg_recalc_loop(value_date, data, day_count, frequency,
                        times, val_df_curve):

    x_val = [x["times"] for x in val_df_curve]
    y_val = [x["df"] for x in val_df_curve]
    ifunc_val = interpolation(x_val, y_val, float(1/366), is_function=True)
    val_dis_curve = [ifunc_val(x) for x in times]

    cash_flows = np.asarray([datum["cash_flow"] for datum in data])
    pvs = cash_flows * val_dis_curve

    return pvs


def floatleg(data, ref_curve, val_curve=None):
    structures = _floatleg_gen_structure(data)
    if not data["principal_exchange"]:
        for structure in structures:
            structure["fv_flow"] = 0

    if val_curve is not None:
        risks = _floatleg_value(data["value_date"], structures,
                                data["day_count"], data["frequency"],
                                ref_curve, val_curve)
    else:
        risks = None

    return {"structure": structures, "risks": risks}


def _floatleg_gen_structure(data, holidays=[]):

    # generate the face value and coupons
    dates = _dates_gen_structure(data)
    dates = deque(dates)

    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)
    for date in dates:
        offset = -start_basis[data["fixing_basis"]]
        date["fixing_date"] = np.busday_offset(date["start_date"], offset,
                                               roll='backward', busdaycal=bdc)
        date["face_value"] = data["face_value"]
        date["margin"] = data["margin"]
        date["cpn_dcf"] = day_cf(data["day_count"],
                                 date["start_date"],
                                 date["end_date"],
                                 bondmat_date=data["maturity"],
                                 next_coupon_date=date["end_date"],
                                 business_day=data["business_day"],
                                 Frequency=data["frequency"])
        # fixing date is a forward date
        if date["fixing_date"] > data["value_date"]:
            date["is_fixed"] = False
            date["accrued"] = 0

        # fixing date is for the current coupon period
        elif (date["fixing_date"] <= data["value_date"] and
              date["end_date"] > data["value_date"]):
            date["is_fixed"] = True
            date["coupon"] = data["current_coupon"]
            date["coupon_interest"] = (date["cpn_dcf"] *
                                       date["coupon"] *
                                       date["face_value"] / 100)
            if data["value_date"] > date["start_date"]:
                acc_dcf = day_cf(data["day_count"],
                                 date["start_date"],
                                 data["value_date"],
                                 bondmat_date=data["maturity"],
                                 next_coupon_date=date["end_date"],
                                 business_day=data["business_day"],
                                 Frequency=data["frequency"])
                date["accrued"] = (acc_dcf / date["cpn_dcf"] *
                                   date["coupon_interest"])
            else:
                date["accrued"] = 0

        else:
            date["is_fixed"] = True
            date["accrued"] = 0
        date["fv_flow"] = 0
    dates[-1]["fv_flow"] = data["face_value"]

    return dates


def _floatleg_value(value_date, structures, day_count, frequency,
                    ref_curve, val_curve):

    # copy only active structure period
    data = [dict(x) for x in structures if value_date < x["end_date"]]
    # covert data into list
    times = np.asarray([day_cf("Actual/365", value_date, x["end_date"])
                        for x in data])

    # Original Curves
    ref_df_curve = df_gen(ref_curve, return_type="times")
    val_df_curve = df_gen(val_curve, return_type="times")

    # Shift the original curves by 1 bp
    ref_curve1 = dict(ref_curve)
    ref_curve1["rates"] = shift_curve(ref_curve1["rates"])
    ref_df_curve1 = df_gen(ref_curve1, return_type="times")

    val_curve1 = dict(val_curve)
    val_curve1["rates"] = shift_curve(val_curve1["rates"])
    val_df_curve1 = df_gen(val_curve1, return_type="times")

    # Shift the curves by another 1bp
    ref_curve2 = dict(ref_curve1)
    ref_curve2["rates"] = shift_curve(ref_curve2["rates"])
    ref_df_curve2 = df_gen(ref_curve2, return_type="times")

    val_curve2 = dict(val_curve1)
    val_curve2["rates"] = shift_curve(val_curve2["rates"])
    val_df_curve2 = df_gen(val_curve2, return_type="times")

    pvs0 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve, val_df_curve)
    value = np.sum(pvs0)

    pvs1 = _float_recalc_loop(value_date, data, day_count, frequency,
                              times, ref_df_curve1, val_df_curve1)
    value1 = np.sum(pvs1)

    pvs2 = _floatleg_recalc_loop(value_date, data, day_count, frequency,
                                 times, ref_df_curve2, val_df_curve2)
    value2 = np.sum(pvs2)
    pvbp01 = value1 - value
    duration = pvbp01 / 0.01
    duration2 = (value2 - value1) / 0.01
    convexity = (duration2 - duration) / 0.01

    return {"duration": duration, "pvbp01": pvbp01, "convexity": convexity,
            "value": value}


def _floatleg_recalc_loop(value_date, data, day_count, frequency,
                          times,  df_curve, val_df_curve):

    x_axis = [x["times"] for x in df_curve]
    y_axis = [x["df"] for x in df_curve]
    ifunc = interpolation(x_axis, y_axis, float(1/366), is_function=True)

    x_val = [x["times"] for x in val_df_curve]
    y_val = [x["df"] for x in val_df_curve]
    ifunc_val = interpolation(x_val, y_val, float(1/366), is_function=True)
    val_dis_curve = [ifunc_val(x) for x in times]

    for datum in data:
        end_time = day_cf("Actual/365", value_date, datum["end_date"])

        if not datum["is_fixed"]:  # coupon has not been fixed
            start_time = day_cf("Actual/365", value_date,
                                datum["start_date"])

            fwd_df = calc_fwd_df(start_time, end_time, ifunc=ifunc)
            fwd_rate = calc_shortrate_from_df(datum["start_date"],
                                              datum["end_date"],
                                              fwd_df, day_count)
            datum["coupon"] = fwd_rate + datum["margin"]
            datum["coupon_interest"] = (datum["coupon"] * datum["cpn_dcf"] *
                                        datum["face_value"] / 100)

        if datum.get("coupon_interest"):
            datum["cash_flow"] = datum["coupon_interest"] + datum["fv_flow"]
        else:
            datum["cash_flow"] = datum["fv_flow"]

    cash_flows = np.asarray([datum["cash_flow"] for datum in data])
    pvs = cash_flows * val_dis_curve

    return pvs


def swaps(fxleg, flleg, fxlval_curve=None, flref_curve=None,
          fllval_curve=None, fxrate=1.00):

    flleg = floatleg(flleg, ref_curve=flref_curve, val_curve=fllval_curve)
    fxleg = fixleg(fxleg, val_curve=fxlval_curve)
    fxleg_risks = fxleg["risks"]
    flleg_risks = flleg["risks"]

    swapvalue = fxleg_risks["value"] * fxrate - flleg_risks["value"]
    swappvbp01 = fxleg_risks["pvbp01"] * fxrate - flleg_risks["pvbp01"]

    return {"value": swapvalue, "pvbp01": swappvbp01}

def loans_structures(loans):
    """
    Generate loan schedule.

            Parameters:
                loans: a dictionary with the following keys - value_date,
                business_day, start_date, value_date, maturity,
                frequency, day_count, date_generation, face_value, rate,
                rate_type and rate_compounding.

            Returns:
                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """

    if isinstance(loans, dict):
        results = _loan_gen_structure(loans)

    elif isinstance(loans, list):
        results = []
        for loan in loans:
            if isinstance(loan, dict):
                results.append(_fixbond_gen_structure(loan))

    return results


def _loan_gen_structure(loan):
    value_date = loan.get("value_date")
    structures = _dates_gen_structure(loan)
    if loan["rate_type"] == "fixed":
        for structure in structures:
            pass
    elif loan["rate_type"] == float:
        pass

    return structures


def calc_customfix_structures(structures, day_count, frequency, business_day):
    """
    Calculate all the values of the custome structure.

            Parameters:
                value_date: numpy.datetime64.
                structures: a list of dictionaries with the following keys: start_date, end_date, face_value, coupon, fv_flow.
                day_count: str
                frequeny: str
            Returns:
                a dictionary with the following keys - value_date, start_date, end_date, face_value, coupon, coupon_interest, cash_flow and cpn_dcf.
    """

    for structure in structures:
        structure["cpn_dcf"] = day_cf(day_count,
                                      structure["start_date"],
                                      structure["end_date"],
                                      bondmat_date=structure["end_date"],
                                      next_coupon_date=structure["end_date"],
                                      business_day=business_day,
                                      Frequency=frequency)

        structure["coupon_interest"] = (structure["cpn_dcf"] *
                                        structure["coupon"] *
                                        structure["face_value"] / 100)
        if structure.get("fv_flow"):
            structure["cash_flow"] = (structure["coupon_interest"] +
                                            structure["fv_flow"])
        else:
            structure["fv_flow"] = 0.00
            structure["cash_flow"] = (structure["coupon_interest"])

    return structures


def value_customfix_structures(value_date, structures, day_count, frequency,
                               dis_curve):
    """
    Calculate all the values of the custome structure.

            Parameters:
                value_date: numpy.datetime64.
                structures: a list of dictionaries with the following keys:
                    start_date, end_date, face_value, coupon and fv_flow.
                day_count: str
                frequeny: str
                dis_curve: list of ditionaries with the following keys: time
                and df
            Returns:
                a dictionary with the following keys - value_date, start_date,
                end_date, face_value, coupon, coupon_interest, cash_flow and
                cpn_dcf.
    """
    newstructures = [dict(x) for x in structures]

    if isinstance(dis_curve, list):
        xaxis = [x["times"] for x in dis_curve]
        yaxis = [x["df"] for x in dis_curve]
        ifunc = interpolation(xaxis, yaxis, 1, model='chip', is_function=True)

    for structure in newstructures:
        temp_structure = {}
        temp_structure["cpn_dcf"] = day_cf(day_count,
                                      structure["start_date"],
                                      structure["end_date"],
                                      bondmat_date=structure["end_date"],
                                      next_coupon_date=structure["end_date"],
                                      Frequency=frequency)

        temp_structure["coupon_interest"] = (temp_structure["cpn_dcf"] *
                                        structure["coupon"] *
                                        structure["face_value"] / 100)
        if structure.get("fv_flow"):
            temp_structure["cash_flow"] = (temp_structure["coupon_interest"] +
                                            structure["fv_flow"])
        else:
            structure["fv_flow"] = 0.00
            temp_structure["cash_flow"] = (temp_structure["coupon_interest"])
        temp_structure["time"] = float(day_cf(day_count,
                                      value_date,
                                      structure["end_date"]))
        structure.update(temp_structure) # merging the dictionary

    # Discount curve was provided
    if isinstance(dis_curve, list):
        for structure in newstructures:
            if structure["time"] > 0:
                structure["df"] = float(ifunc(structure["time"]))
                structure["pv"] = structure["df"] * structure["cash_flow"]
            else:
                structure["df"] = 0
                structure["pv"] = 0

    # yield to maturity was provided instead
    elif isinstance(dis_curve, float):
        newstructures = fixbond_value(value_date, newstructures, dis_curve,
                                      day_count, frequency)

    return newstructures


def _calc_customfix_risks(value_date, structures, day_count, frequency,
                               dis_curve):
    data = list(structures)
    mylist = []
    for datum in data:
        mydatum = dict(datum)
        mylist.append

    if isinstance(dis_curve, list):
        pass
    elif isinstance(dis_curve, float):
        ytm = dis_curve
    # Calculate the PV weight
    maturity = mylist[-1]["end_date"]
    for datum in data:
            datum["time"] = day_cf("Actual/365",
                                   value_date,
                                   datum["end_date"],
                                   bondmat_date=maturity,
                                   next_coupon_date=datum["end_date"],
                                   business_day="No Adjustment",
                                   Frequency=frequency)
            datum["cf_weight"] = datum["pv"] * datum['time']
            datum["period_df1"] = 1 / (1 + datum["ytm_dcf"] *
                                       (bond["ytm"] + 0.01) / 100)

    mac_dur = 0
    val = 0
    for datum in data:
        if datum["end date"] > bond["value_date"]:
            mac_dur += datum["cf_weight"]
            val += datum["pv"]
    mac_dur = mac_dur / val
    mod_dur = mac_dur(1 + bond["ytm"] /
                      (12 / frequencies[bond["frequency"]]))


def create_structures_from_dates(dates, coupons, face_values, fv_flows):
    """
    Create structures from date.

            Parameters:
                dates: list of dict with the following keys - start_date and
                end_date
                coupons: list of float.
                face_values: list of float
                fv_flows: list of float.

            Returns:
                list of dictionaries with the following keys - start_date, end_date, coupon, face_value, fv_flow
    """

    structures = list(map( lambda date, coupon, face_value, fv_flow: {"start_date": date["start_date"], "end_date": date["end_date"], "coupon": coupon, "face_value": face_value, "fv_flow": fv_flow}, dates, coupons, face_values, fv_flows))

    return structures


def mm_matamount(start, maturity, rate, rate_basis, day_count, principal=10_000_000):
    dcf = day_cf(day_count, start, maturity)
    if rate_basis == 'Simple':
        matamount = principal * (1 + rate * dcf/100)
    elif rate_basis == 'Continuous':
        matamount = principal * math.exp(rate * dcf)

    return matamount
    
