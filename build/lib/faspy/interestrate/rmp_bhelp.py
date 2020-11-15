# !/usr/bin/env python3.
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 00:08:48 2020

@author: ROSLI MOHD SANI
"""
from faspy.interestrate.rmp_dates import generate_datesv102 as gen_dates, \
    frequencies as fre
from faspy.interestrate.rmp_dates import day_count_factor as day_cf

from numpy import datetime64 as dt64, busday_offset as bus_off
import numpy as np
from faspy.interestrate.conventions import start_basis

import functools
from faspy.interestrate.rmp_solvers import solver_bond_derivatives as sbd, bond_risk
from numba import njit, jit


bondtype = ('Discount Bond', 'Convertible Discount Bond',
            'Callable Discount Bond', 'Callable Convertible Discount Bond',
            'Bullet Bond', 'Callable Bullet Bond',
            'Fixed Rate Bond', 'Fixed Rate ABS', 'Fixed Rate MBS',
            'Bond with Warrants', 'Convertible Bond',
            'Exchangeable Bond', 'Callable Bond', 'Callable ABS',
            'Callable MBS', 'Callable Convertible Bond',
            'Callable Exchangeable Bond', 'Amortising Bond',
            'Fixed Rate Amortising ABS Bond', 'Callable Amortising Bond',
            'Stepping FRB', 'Convertible Stepping Bond',
            'Exchangeable Stepping Bond', 'Callable Stepping Bond',
            'Callable Convertible Stepping Bond', 'Stepping Amortising Bond',
            'Callable Stepping Amortising Bond',
            'Floating Rate Note', 'Floating Rate ABS',
            'Floating Rate MBS', 'Floating Convertible Bond',
            'Floating Exchangeable Bond', 'Floating Callable Bond',
            'Floating Callable ABS', 'Floating Callable MBS',
            'Floating Callable Convertible Bond',
            'Floating Callable Exchangeable', 'Bond',
            'Floating Amortising Note', 'Stepping FRN',
            'Fixed Rate Bond', 'Discount Bond', 'Floating Rate Note',
            'Bond with Secondary Notes', 'Stepping Bonds with Secondary Notes',
            'Callable Convertible Discount Bond',
            'Callable Bond with Secondary Notes',
            'Callable Discount Bond', 'Callable Bond', 'Callable Bullet Bond',
            'Convertible Discount Bond', 'Convertible Bond',
            'Callable Convertible Stepping Bond',
            'Callable Stepping Bonds with Secondary Notes',
            'Callable Stepping Amortising Bonds with Secondary Notes',
            'Callable Stepping Bond', 'Amortising Bonds with Secondary Notes',
            'Bullet Bond', 'Stepping FRB', 'Callable Stepping Amortising Bond',
            'Stepping Amortising Bonds with Secondary Notes',
            'Floating Callable Bond', 'Convertible Stepping Bond',
            'Exchangeable Stepping Bond', 'Exchangeable Bond',
            'Callable Convertible Bond', 'Callable Exchangeable Bond',
            'Amortising Bond',
            'Callable Amortising Bond', 'Floating Amortising Note')
bondtype = set(bondtype)


# In[1]:

def construct_fixbond(bond):
    bus_day = None
    if bond['business_day'] == 'NULL':
        bus_day = 'No Adjustment'
    else:
        bus_day = bond['business_day']

    if bond['day_count'] == 'NULL':
        day_count = 'Actual/Actual ICMA'
    else:
        day_count = bond['day_count']

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

    dates = gen_dates(use_date, bond['maturity'], issueDate=start_date,
                      frequency=fre[bond['frequency']],
                      business_day=bus_day, method=bond['date_generation'])

    #print(dates)
    # *********************************************************************************
    #
    #               CALCULATING THE DATA FOR FULL STRUCTURES
    #
    # *********************************************************************************
    fs = list(dates)
    flen = len(fs)
    for i in range(flen):
        fs_sgl = fs[i]
        fs_sgl['fv_flow'] = 0.00
        fs_sgl['fv'] = bond['face_value']
        if i > 0:
            fs_prev = fs[i-1]
            fs_sgl['cpn_dcf'] = day_cf(day_count, fs_prev['date'],
                                       fs_sgl['date'],
                                       bondmat_date=fs[-1]['date'],
                                       next_coupon_date=fs_sgl['date'],
                                       Frequency=12 /
                                       fre[bond['frequency']])
            fs_sgl['coupon'] = bond['coupon']
            fs_sgl['time'] = day_cf('Actual/365', use_date,
                                   fs_sgl['date'],
                                   bondmat_date=fs[-1]['date'],
                                   next_coupon_date=fs_sgl['date'],
                                   Frequency=12 /
                                   fre[bond['frequency']])
            fs_sgl['cf'] = (fs_sgl['coupon'] * fs_sgl['cpn_dcf'] * 0.01
                            * fs_sgl['fv'])
        else:
            fs_sgl['cpn_dcf'] = 0.00
            fs_sgl['coupon'] = 0.00
            fs_sgl['fv'] = 0.00
            fs_sgl['cf'] = 0.00

    fs[-1]['fv_flow'] = bond['face_value']
    bond['full_structures'] = fs

    # *********************************************************************************
    #
    #               CALCULATING THE DATA FOR ACTIVE STRUCTURES
    #
    # *********************************************************************************
    if value_date is None:
        bond['active_structures'] = {}

    elif value_date is not None and start_date is None:
        bond['active_structures'] = bond['full_structures']

    elif value_date < start_date:
        bond['active_structures'] = bond['full_structures']

    else:
        booldate = list(map(lambda x: x['date'] > value_date, fs))
        index = booldate.index(True)  # exclude paid coupon period/start_date
        acs = list(fs[index - 1:])
        acs[0]['dcf'] = 0
        acs[1]['dcf'] = day_cf(day_count, bond['value_date'], acs[1]['date'],
                               bondmat_date=acs[-1]['date'],
                               next_coupon_date=acs[1]['date'],
                               Frequency=12 / fre[bond['frequency']])
        alen = len(acs)
        for i in range(alen):
            acs[i]['time'] = day_cf('Actual/365', bond['value_date'],
                                    acs[i]['date'],
                                    bondmat_date=acs[-1]['date'],
                                    next_coupon_date=acs[i]['date'])

        bond['active_structures'] = acs


# In[]:
def construct_discountbond(bond):

    value_date = bond.get('value_date')
    issue_date = bond.get('issue_date')
    use_date = None
    if value_date is not None and issue_date is not None:
        if issue_date >= value_date:
            use_date = issue_date
        else:
            use_date = value_date
    elif value_date is not None:
        use_date = value_date
    elif issue_date is not None:
        use_date = issue_date

    dcfs = []
    times = []
    cfs = []
    cpns = []
    dates2 = []
    if use_date is not None:
        dcf = day_cf('Actual/365', use_date, bond['maturity'],
                     bondmat_date=bond['maturity'])
        dcfs.append(dcf)
        times.append(dcf)

    cfs.append(bond['position'])
    dates2.append(bond['maturity'])

    # cash_flows = {}
    cash_flows = {'dates': dates2, 'dcfs': dcfs, 'cfs': cfs, 'times': times}
    bond['cash_flows'] = cash_flows
    bond['coupon_structure'] = cpns
    # return bond


# %%
def construct_frn(bond, holidays=[]):
    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)

    bus_day = None
    if bond['business_day'] == 'NULL':
        bus_day = 'No Adjustment'
    else:
        bus_day = bond['business_day']

    if bond['day_count'] == 'NULL':
        day_count = 'Actual/365'
    else:
        day_count = bond['day_count']

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

    dates = gen_dates(use_date, bond['maturity'], issueDate=start_date,
                      frequency=fre[bond['frequency']],
                      business_day=bus_day, method=bond['date_generation'])

    #print(dates)
    # *********************************************************************************
    #
    #               CALCULATING THE DATA FOR FULL STRUCTURES
    #
    # *********************************************************************************
    fs = list(dates)

    flen = len(fs)
    for i in range(flen):
        row = fs[i]
        row['cf'] = 0
        row['time'] = day_cf('Actual/365', use_date, row['date'],
                               bondmat_date=fs[-1]['date'],
                               next_coupon_date=row['date'])
        if i == 0:
            row['dcf'] = 0
            row['y_dcf'] = 0
            row['margin'] = 0
            row['fv'] = 0
            row['fv_flow'] = 0
            row['fixing_date'] = None

        else:
            row['dcf'] = day_cf(day_count, fs[i-1]['date'], row['date'],
                                  bondmat_date=fs[-1]['date'],
                                  next_coupon_date=row['date'],
                                  Frequency=12 / fre[bond['frequency']])
            row['y_dcf'] = row['dcf']
            row['margin'] = bond['margin']
            row['fv'] = bond['face_value']
            row['fv_flow'] = 0
            if i == flen - 1:
                row['fv_flow'] = bond['face_value']
            row['fixing_date'] = bus_off(fs[i-1]['date'],
                                           -start_basis[bond['fixing_basis']],
                                           roll='backward', busdaycal=bdc)

    bond['full_structures'] = fs

    # *********************************************************************************
    #
    #               CALCULATING THE DATA FOR ACTIVE STRUCTURES
    #
    # *********************************************************************************
    if value_date is None:
        bond['active_structures'] = {}

    elif value_date is not None and start_date is None:
        bond['active_structures'] = bond['full_structures']

    elif value_date < start_date:
        bond['active_structures'] = bond['full_structures']

    else:
        # print(value_date,start_date)
        booldate = list(map(lambda x: x['date'] > value_date, fs))
        index = booldate.index(True)  # exclude paid coupon period/start_date
        acs = list(fs[index - 1:])
        acs[0]['y_dcf'] = 0
        acs[1]['y_dcf'] = day_cf(day_count, bond['value_date'], acs[1]['date'],
                               bondmat_date=acs[-1]['date'],
                               next_coupon_date=acs[1]['date'],
                               Frequency=12 / fre[bond['frequency']])
        alen = len(acs)
        for i in range(alen):
            acs[i]['time'] = day_cf('Actual/365', bond['value_date'],
                                    acs[i]['date'],
                                    bondmat_date=acs[-1]['date'],
                                    next_coupon_date=acs[i]['date'])

        bond['active_structures'] = acs

    return bond


# In[]:
def construct_bond(bond):

    if bond['type'] not in bondtype:
        raise Exception(bond['type'] + ' is not in the bondtype list')

    if bond['type'] == 'Fixed Rate Bond':
        construct_fixbond(bond)

    elif bond['type'] == 'FRN':
        pass

    elif bond['type'] == 'Discount Bond':
        construct_discountbond(bond)

    elif bond['type'] == 'Stepping Bond':
        pass


def fixbond_price(value_date, bond, yld):

    try:
        val_date = dt64(value_date)

    except Exception:
        raise Exception('Value date is not a date')

    try:
        acs = bond['active_structures']
    except Exception:
        raise Exception('full_structures does not have the required key(s)')

    try:
        freq = bond['frequency']
    except Exception:
        raise Exception('Error getting frequency of the floating leg')

    if bond['value_date'] is None:
        raise Exception('Value date for the bond is not set')

    bond['last_coupon_date'] = acs[0]['date']
    bond['accrued_interest'] = 0
    next_cpn = acs[1]
    # Calculation for fixed coupon
    if bond['value_date'] >= acs[0]['date']:
        coupon = bond.get('coupon')
        if coupon is not None:
            #print(next_cpn)
            next_cpn['coupon'] = coupon
            next_cpn['cf'] = coupon * next_cpn['dcf'] * next_cpn['fv'] * 0.01
            y_dcf = day_cf(bond['day_count'], val_date, next_cpn['date'],
                         bondmat_date=acs[-1]['date'],
                         next_coupon_date=next_cpn['date'],
                         Frequency=12 / fre[bond['frequency']])
            next_cpn['y_dcf'] = y_dcf
            df = 1 / (1 + yld * y_dcf * 0.01)
            next_cpn['pv'] = (next_cpn['cf'] + next_cpn['fv_flow']) * df
            next_cpn['df'] = df
            next_cpn['p_df'] = df

            # Calculation for PVBP01
            df01 = 1 / (1 + (yld + 0.01) * y_dcf * 0.01)
            next_cpn['pv01'] = (next_cpn['cf'] + next_cpn['fv_flow']) * df01
            next_cpn['df01'] = df01
            next_cpn['p_df01'] = df01
            #print(next_cpn)
        accrued_dcf = day_cf(bond['day_count'], bond['last_coupon_date'],
                             bond['value_date'], bondmat_date=acs[-1]['date'],
                             next_coupon_date=next_cpn['date'],
                             Frequency=12 / fre[bond['frequency']])
        # print(bond.get('current_coupon'), next_cpn)
        bond['accrued'] = accrued_dcf
        bond['accrued_interest'] = (accrued_dcf * next_cpn['coupon'] * 0.01 *
                                    next_cpn['fv'])

        lacs = len(acs)
        for i in range(2, lacs, 1):
            acs_sgl = acs[i]
            lastrow = acs[i-1]
            acs_sgl['y_dcf'] = acs_sgl['cpn_dcf']
            acs_sgl['p_df'] = 1 / (1 + yld * acs_sgl['y_dcf'] * 0.01)
            acs_sgl['p_df01'] = 1 / (1 + (yld + 0.01) * acs_sgl['y_dcf']
                                     * 0.01)
            acs_sgl['df'] = acs_sgl['p_df'] * lastrow['df']
            acs_sgl['df01'] = acs_sgl['p_df01'] * lastrow['df01']
            acs_sgl['pv'] = ((acs_sgl['cf'] + acs_sgl['fv_flow'])
                             * acs_sgl['df'])
            acs_sgl['pv01'] = ((acs_sgl['cf'] + acs_sgl['fv_flow'])
                               * acs_sgl['df01'])

        sum_pv = 0
        sum_pv01 = 0
        for i in range(1, lacs, 1):
            sgl_acs = acs[i]
            sum_pv += sgl_acs['pv']
            sum_pv01 += sgl_acs['pv01']
        bond['pvbp01'] = sum_pv01 - sum_pv
        bond['proceed'] = sum_pv
        bond['price per 100'] = ((bond['proceed'] - bond['accrued_interest'])
                                 / bond['face_value'] * 100)

        macD = 0
        for i in range(1, lacs, 1):
            acs_sgl = acs[i]
            macD = (macD + acs_sgl['time'] *
                    (acs_sgl['cf'] + acs_sgl['fv_flow']) *
                    acs_sgl['df'])

        macD = macD / bond['proceed']
        bond['modified duration'] = -macD / (1 + (yld * 0.01) / (12 / fre[freq]))
        bond['macaulay duration'] = macD
        #print('acs===>', acs)
        derivatives = bond_risk(acs, yld)
        bond['duration'] = derivatives['duration']
        bond['convexity'] = derivatives['convexity']
        bond['active_structures'] = acs
    else:
        # *******************************************************************
        # else clause process the bond as a forward issuance bond
        # disfac (discount factor) has to be provided to properly price the
        # bond) yld is taken to be the fwd yield of the bond
        # ********************************************************************
        raise Exception('Forward bond issuance is yet to be done')

# %%
def frn_price(bond, disc_curve):

    try:
        val_date = dt64(bond.get('value_date'))

    except Exception:
        raise Exception('Value date is not a date')

    try:
        acs = bond['active_structures']

    except Exception:
        raise Exception('full_structures does not have the required key(s)')

    bond['last_coupon_date'] = acs[0]['date']
    bond['accrued_interest'] = 0
    next_cpn = acs[1]
    interp = disc_curve.interpolate(1, is_function=True)
    # Calculation for fixed coupon
    if bond['value_date'] >= acs[0]['date']:
        coupon = next_cpn['coupon']
        if coupon is not None:
            gen_time = next_cpn['time']
            next_cpn['cf'] = coupon * next_cpn['dcf'] * next_cpn['fv'] * 0.01
            dcf = next_cpn['dcf']
            df = interp(gen_time)
            next_cpn['ref_rate'] = (1 / df - 1) / (dcf * 0.01)
            next_cpn['ref_with_spread'] = (next_cpn['ref_rate']
                                           + bond['spread'])
            df = 1 / (1 + next_cpn['ref_with_spread'] * dcf * 0.01)
            next_cpn['pv'] = (next_cpn['cf'] + next_cpn['fv_flow']) * df
            next_cpn['df'] = df
            # Calculation for PVBP01
            next_cpn['coupon01'] = coupon
            next_cpn['cf01'] = next_cpn['cf']
            df = interp(gen_time)
            next_cpn['ref_rate01'] = (1 / df - 1) / (dcf * 0.01) + 0.01
            next_cpn['ref01_with_spread'] = (next_cpn['ref_rate01']
                                             + bond['spread'])
            #print(next_cpn['ref_rate01'], bond['spread'])
            df01 = 1 / (1 + next_cpn['ref01_with_spread'] * dcf * 0.01)
            next_cpn['pv01'] = (next_cpn['cf01'] + next_cpn['fv_flow']) * df01
            next_cpn['df01'] = df01

    accrued_dcf = day_cf(bond['day_count'], bond['last_coupon_date'],
                         bond['value_date'], bondmat_date=acs[-1]['date'],
                         next_coupon_date=next_cpn['date'],
                         Frequency=12 / fre[bond['frequency']])
    #print(bond.get('current_coupon'), next_cpn)
    bond['accrued'] = accrued_dcf
    bond['accrued_interest'] = (accrued_dcf * next_cpn['coupon'] * 0.01 *
                                next_cpn['fv'])


    lacs = len(acs)
    # Calculation for unfixed coupon
    for i in range(2, lacs, 1):
        prev_acs = acs[i-1]
        sgl_acs = acs[i]
        #sgl_acs['fixed'] = False
        time0 = prev_acs['time']
        time1 = sgl_acs['time']
        df0 = disc_curve.interpolate(time0)
        df1 = disc_curve.interpolate(time1)
        dcf = sgl_acs['dcf']
        fwd_df = df1/df0
        sgl_acs['ref_rate'] = (1 / fwd_df - 1) / (dcf * 0.01)
        sgl_acs['coupon'] = (sgl_acs['ref_rate'] + sgl_acs['margin'])
        sgl_acs['cf'] = sgl_acs['coupon'] * dcf * 0.01 * sgl_acs['fv']
        sgl_acs['ref_with_spread'] = (sgl_acs['ref_rate'] + bond['spread'])
        sgl_acs['fwd_df'] = 1 / (1 + sgl_acs['ref_with_spread'] * dcf * 0.01)
        sgl_acs['df'] = sgl_acs['fwd_df'] * prev_acs['df']
        sgl_acs['pv'] = (sgl_acs['cf'] + sgl_acs['fv_flow']) * sgl_acs['df']
        # Calculation for PVBP01
        sgl_acs['ref_rate01'] = (1 / fwd_df - 1) / (dcf * 0.01) + 0.01
        sgl_acs['coupon01'] = (sgl_acs['ref_rate01'] + sgl_acs['margin'])
        sgl_acs['cf01'] = sgl_acs['coupon01'] * dcf * 0.01 * sgl_acs['fv']
        sgl_acs['ref01_with_spread'] = (sgl_acs['ref_rate01'] +
                                        bond['spread'])
        sgl_acs['fwd_df01'] = 1 / (1 + sgl_acs['ref01_with_spread'] *
                                   dcf * 0.01)
        sgl_acs['df01'] = sgl_acs['fwd_df01'] * prev_acs['df01']
        sgl_acs['pv01'] = ((sgl_acs['cf01'] + sgl_acs['fv_flow']) *
                           sgl_acs['df01'])

    sum_pv = 0
    sum_pv01 = 0
    for i in range(1, lacs, 1):
        sgl_acs = acs[i]
        sum_pv += sgl_acs['pv']
        sum_pv01 += sgl_acs['pv01']
    bond['pvbp01'] = sum_pv01 - sum_pv
    bond['proceed'] = sum_pv
    bond['price per 100'] = ((bond['proceed'] - bond['accrued_interest'])
                             / bond['face_value'] * 100)
    # derivatives = sbd(bond['face_value'], active_cfs, active_dcfs,
    # #                 active_fvs_flow, yld)
    # bond['duration'] = derivatives['duration']
    # bond['convexity'] = derivatives['convexity']



# %%
def bond_price(value_date, bond, freq=1):

    try:
        myyield = bond.get('ytm')
    except Exception:
        raise ValueError('ytm is required for valueation')

    bistype = bond.get('instrument_subtype')
    btype = bond.get('type')
    type_li = [bistype, btype]

    if bistype not in bondtype and btype not in bondtype:
        raise Exception(bistype + ' or ' + btype +
                        ' are not in the bondtype list')

    if 'Fixed Rate Bond' in type_li:
        fixbond_price(value_date, bond, myyield)

    elif 'Discount Bond' in type_li:
        cf = bond['cash_flows']
        bond['proceed'] = cf['cfs'][0] * (1 /
                                          (1 + (myyield * 0.01) /
                                           freq))**(cf['times'][0] * freq)
        value2 = cf['cfs'][0] * (1 /
                                 (1 + ((myyield + 0.01) * 0.01) /
                                  freq)) ** (cf['times'][0] * freq)
        bond['pvbp01'] = value2 - bond['proceed']
        bond['accrued interest'] = 0
        bond['price per 100'] = bond['proceed']/cf['cfs'][0] * 100
        bond['risk_stats'] = {}
        bond['risk_stats']['macaulay duration'] = cf['times'][0]
        bond['risk_stats']['modified duration'] = (cf['times'][0] /
                                                   (1 + (myyield * 0.01) /
                                                    (freq)))
        derivatives = sbd(bond['position'], cf['cfs'], cf['times'], myyield)
        bond['risk_stats']['duration'] = derivatives['duration']
        bond['risk_stats']['convexity'] = derivatives['convexity']

    else:
        raise Exception(btype + ' and ' + bistype + ' are not supported yet')


@jit
def _calc_df(active_cfs, active_dcfs, accrued_dcf, yld):
    dfs = []
    dfs01 = []
    totalcf = len(active_cfs)
    for i in range(totalcf):
        df = 1
        df01 = 1
        for ind in range(i+1):
            if ind == 0:
                df = df * 1 / (1 + yld *
                               (active_dcfs[ind] - accrued_dcf) * 0.01)
                df01 = df01 * 1 / (1 + (yld + 0.01) *
                                   (active_dcfs[ind] - accrued_dcf) * 0.01)
            else:
                df = df * 1 / (1 + yld * active_dcfs[ind] * 0.01)
                df01 = df01 * 1 / (1 + (yld + 0.01) * active_dcfs[ind] * 0.01)
        # print(df)
        dfs.append(df)
        dfs01.append(df01)

    return (dfs, dfs01)


@njit
def _calc_risk(active_cfs, dfs, dfs01, active_fvs_flow):
    # present value of each coupon flow
    dfs_len = len(dfs)
    active_cpns_values = []
    active_cpns_values01 = []
    active_pflow_values = []
    active_pflow_values01 = []
    for i in range(dfs_len):
        df = dfs[i]
        df01 = dfs01[i]
        cf = active_cfs[i]
        flow = active_fvs_flow[i]
        active_cpns_values.append(float(cf * df))
        active_pflow_values.append(float(flow * df))
        active_cpns_values01.append(float(cf * df01))
        active_pflow_values01.append(float(flow * df01))

    # net present value of all the principal flows
    pflow_value = 0.00
    pflow_len = len(active_pflow_values)
    for i in range(pflow_len):
        pflow_value = pflow_value + float(active_pflow_values[i])

    # summing the coupon pv to get the value of the bond
    cpn_values = 0.00
    cpn_len = len(active_cpns_values)
    for i in range(cpn_len):
        cpn_values = cpn_values + float(active_cpns_values[i])

    bond_values = cpn_values + pflow_value

    # net present value of all the principal flows
    pflow_value01 = 0.00
    pflow_len = len(active_pflow_values01)
    for i in range(pflow_len):
        pflow_value01 = pflow_value01 + float(active_pflow_values01[i])

    # summing the coupon pv to get the value of the bond
    cpn_values01 = 0.00
    cpn_len = len(active_cpns_values01)
    for i in range(cpn_len):
        cpn_values01 = cpn_values01 + float(active_cpns_values01[i])

    bond_values01 = cpn_values01 + pflow_value01
    pvbp01 = float(bond_values01 - bond_values)

    return (bond_values, pvbp01, active_cpns_values, active_pflow_values)
