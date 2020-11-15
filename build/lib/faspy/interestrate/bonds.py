
import numpy as np
from numpy import datetime64 as dt64
from sympy import Symbol, diff, solve
from faspy.interestrate.rmp_dates import generate_dates as gen_dates
from faspy.interestrate.rmp_dates import frequencies as fre
from faspy.interestrate.rmp_dates import day_count_factor as day_cf, convs
import functools
from faspy.interestrate.rmp_solvers import solver_bond_derivatives as sbd
from numba import njit


def init_bond_cfs(face_value, value_date, maturity, convention,
                  freq, cpns=[], issue_date=None, business_day='No Adjustment',
                  date_method='Backward from maturity date',
                  holidays=[], yld=0):

    dates = gen_dates(value_date, maturity, issueDate=issue_date,
                      frequency=fre[freq], business_day=business_day,
                      method=date_method, holidays=holidays)
    return dates


def bond_cfs_values(face_value, value_date, convention, freq='Semi-Annual',
                    cpn_dates=[], cpns=[], business_day='No Adjustment',
                    holidays=[], yld=0):
    if len(cpn_dates) <= len(cpns):
        return None
    elif len(cpn_dates) > len(cpns) + 1:
        return None

    results = {}
    noofdates = len(cpn_dates)
    cpn_dcfs = []
    yld_dcfs = []
    fwd_dfs = []
    fwd_dfs01 = []
    dfs = []
    dfs01 = []
    cpn_interests = []
    times = []
    accrued_interest = 0
    maturity = cpn_dates[-1]
    for i in range(1, noofdates, 1):
        dcf = None
        dcf_unaccr = None
        time = None
        start_date = cpn_dates[i-1]
        end_date = cpn_dates[i]
        # day count factor for coupon interest
        dcf = day_cf(convention, start_date, end_date,
                     bondmat_date=cpn_dates[-1],
                     next_coupon_date=end_date,
                     business_day=business_day,
                     Frequency=12 / fre[freq])
        time = day_cf('Actual/365', value_date, end_date,
                      bondmat_date=cpn_dates[-1],
                      next_coupon_date=end_date,
                      business_day=business_day,
                      Frequency=12 / fre[freq])

        cpn_interest = cpns[i-1] * 0.01 * dcf * face_value
        # calculate day count factor for discount factor
        if start_date < value_date and end_date > value_date:
            dcf_accr = day_cf(convention, start_date, value_date,
                              bondmat_date=maturity,
                              next_coupon_date=end_date,
                              business_day=business_day,
                              Frequency=12 / fre[freq])
            dcf_unaccr = dcf - dcf_accr
            accrued_interest = dcf_accr * face_value * cpns[i-1] * 0.01

        elif end_date <= value_date:
            dcf_unaccr = 0
        else:  # start_date >= value_date:
            dcf_unaccr = dcf

        fwd_df = 1 / (1 + yld * 0.01 * dcf_unaccr)
        fwd_df01 = 1/(1 + (yld + 0.01) * 0.01 * dcf_unaccr)

        times.append(time)
        fwd_dfs.append(fwd_df)
        fwd_dfs01.append(fwd_df01)
        cpn_dcfs.append(dcf)
        yld_dcfs.append(dcf_unaccr)
        cpn_interests.append(cpn_interest)

    results['cpn_interest'] = cpn_interests
    results['fwd_dfs'] = fwd_dfs
    results['cpn_dcfs'] = cpn_dcfs
    results['yld_dcfs'] = yld_dcfs

    fwd_dfs_len = len(fwd_dfs)
    # calculating df from fwd df
    df = None
    df01 = None
    for i in range(fwd_dfs_len):
        # df= None
        if i == 0:
            df = fwd_dfs[i]
            df01 = fwd_dfs01[i]
        else:
            for k in range(i+1):
                if k == 0:
                    df = fwd_dfs[k]
                    df01 = fwd_dfs01[k]
                else:
                    df = df * fwd_dfs[k]
                    df01 = df01 * fwd_dfs01[k]
        dfs.append(df)
        dfs01.append(df01)

    results['dfs'] = dfs
    results['dfs01'] = dfs01
    # should use reduce
    pv = 0
    pv01 = 0
    dfs_len = len(dfs)
    # calculate the bond value
    for i in range(dfs_len):
        if yld_dcfs[i] != 0:
            pv = pv + cpn_interests[i]*dfs[i]
            pv01 = pv01 + cpn_interests[i]*dfs01[i]

    pv = pv + face_value * dfs[dfs_len - 1]
    pv01 = pv01 + face_value * dfs01[dfs_len - 1]
    results['pv'] = pv
    results['pv01'] = pv01

    # calculate the macaulay duration
    macD = 0
    for i in range(dfs_len):
        if yld_dcfs[i] != 0:
            macD = macD + times[i]*cpn_interests[i]*dfs[i]
    macD = ((macD + times[dfs_len - 1] * face_value *
             dfs[dfs_len - 1]) / pv)

    results['modified duration'] = macD/(1 + yld / 12 / fre[freq])
    results['macaulay duration'] = macD
    results['accrued'] = accrued_interest

    derivatives = solver_bond_duration(face_value, cpns, results['cpn_dcfs'],
                                       results['yld_dcfs'], yld)
    results['duration'] = derivatives['duration']
    results['convexity'] = derivatives['convexity']

    return results


# slow solver but working
def solver_bond_yield(fv, pv, cpns=[], cpn_dcfs=[], yld_dcfs=[]):
    yld = Symbol('yld')
    df = Symbol('df')
    expr = Symbol('expr')
    expr = None
    dcfs_len = len(yld_dcfs)
    for i in range(dcfs_len):
        ci = None
        df = None
        if yld_dcfs[i] != 0:
            ci = cpns[i] * cpn_dcfs[i] * 0.01 * fv

            if df is None:
                df = (1 + yld * yld_dcfs[i] * 0.01)

            for k in range(i):
                if yld_dcfs[k] != 0:
                    df = df * (1 + yld * yld_dcfs[k] * 0.01)

            if expr is None:
                expr = ci/df
            else:
                expr = expr + ci/df

            if i == dcfs_len - 1:
                expr = expr + fv/df - (pv)

    solved_yield = solve(expr, yld)
    new_yields = []
    for i in range(len(solved_yield)):
        try:
            float(solved_yield[i])
        except Exception:
            continue
        new_yields.append(solved_yield[i])

    new_yields = list(filter(lambda x: x > 0, new_yields))
    return new_yields


def solver_bond_duration(fv, cpns=[], cpn_dcfs=[], yld_dcfs=[],
                         my_yield=5):
    yld = Symbol('yld')
    df = Symbol('df')
    expr = Symbol('expr')
    expr = None
    dcfs_len = len(yld_dcfs)
    for i in range(dcfs_len):
        ci = None
        df = None
        if yld_dcfs[i] != 0:
            ci = cpns[i] * cpn_dcfs[i] * 0.01 * fv

            if df is None:
                df = (1 + yld * yld_dcfs[i] * 0.01)

            for k in range(i):
                if yld_dcfs[k] != 0:
                    df = df * (1 + yld * yld_dcfs[k] * 0.01)

            if expr is None:
                expr = ci/df
            else:
                expr = expr + ci/df

            if i == dcfs_len - 1:
                expr = expr + fv/df

    expr_duration = diff(expr, yld)
    duration = expr_duration.evalf(subs={yld: my_yield})
    # print('Duration >>>>',duration)
    expr_convex = diff(expr_duration, yld)
    convexity = expr_convex.evalf(subs={yld: my_yield})
    # print('Convexity >>>>>>',convexity)
    return {'duration': duration, 'convexity': convexity}



def create_bond_dates(bond):
    """
    Create list of coupon dates.

    Parameters
    ----------
    bond : dict
        A dictionary containing the following keys - 'day_count', 'postion',
        'frequency', and 'coupon'.

    Returns
    -------
    dates : list of numpy.datetime64
        list of coupon dates

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

    dates = gen_dates(use_date, bond['maturity'], issueDate=start_date,
                      frequency=fre[bond['frequency']],
                      business_day=bus_day, method=bond['date_generation'])

    return dates


def create_bond_full_structures(bond, dates):
    """
    Create full structure from coupon dates.

    Parameters
    ----------
    bond : dict
        A dictionary containing the following keys - 'day_count', 'postion',
        'frequency', and 'coupon'.
    dates : array of numpy.datetime64
        full list of coupon dates.

    Returns
    -------
    dict
        a dictionary with the following keys: 'dates', 'dcfs', 'cfs',
            'times', 'coupons', 'fvs','fvs_flow'.

    """
    day_count = bond.get('day_count')
    if day_count == 'NULL' or day_count is None:
        day_count = 'Actual/Actual'

    face_value = bond.get('position')
    if face_value is None:
        face_value = 100.00
    else:
        try:
            face_value = float(face_value)
        except Exception:
            face_value = 100.00

    full_dates = list(dates)

    dates1 = list(dates)
    dates2 = list(dates)
    dates1.pop(-1)
    dates2.pop(0)
    full_dcfs = list(map(lambda x1, x2: float(day_cf(day_count, x1, x2,
                                               bondmat_date=dates2[-1],
                                               next_coupon_date=x2,
                                               Frequency=12 /
                                               fre[bond['frequency']])),
                         dates1, dates2))
    full_dcfs.insert(0, 0.0)
    full_cpns = list(map(lambda x: float(bond['coupon']), full_dates))
    full_cpns[0] = 0.0
    full_fvs = list(map(lambda x: float(face_value), full_dates))
    full_fvs[0] = 0.0
    full_fvs_flow = list(map(lambda x: 0, full_dates))
    full_fvs_flow[-1] = float(face_value)
    # print(face_value,full_cpns,full_dcfs)
    full_cfs = list(map(lambda cpn, dcf: float(cpn * dcf * 0.01 * face_value),
                        full_cpns, full_dcfs))
    full_times = list(map(lambda x: day_cf('Actual/365', dates[0], x,
                                           bondmat_date=dates2[-1],
                                           next_coupon_date=x), full_dates))
    return {'dates': full_dates, 'dcfs': full_dcfs, 'cfs': full_cfs,
            'times': full_times, 'coupons': full_cpns, 'fvs': full_fvs,
            'fvs_flow': full_fvs_flow}


def create_bond_act_structures(value_date, full_structures,
                               day_count='Actual/Actual', freq='Semi-Annual'):
    """
    Create active structure from full coupon structure.
    An active structure consist of coupon information after value date.

    Parameters
    ----------
    value_date : numpy datetime64
        value date of the bond.
    full_structures : dict
        full_structures must have the following keys:
            'dates', 'dcfs', 'cfs', 'times', 'coupons', 'fvs' and 'fvs_flow'.
    day_count: str
        day count convention and must be in the list of accepted values.
    freq: str
        coupon frequency and must be in the list of accepted values.

    Returns
    -------
    dict
        containing the following information/keys:
            'accrued' -  accrued
            'accrued_interest' - accrued_interest,
            'last_coupon_date' - last_coupon_date,
            'coupons' - an array of coupon rate (%) for each coupon date,
            'coupon_flow' -  array of coupon interest,
            'dates' - coupon dates after value date,
            'fvs' - face value applicable for the coupon date,
            'fvs_flow' -  principal flows after value date if any,
            'times' - time in fraction of a year from value date,
            'dcfs' - day count factor for the coupon period.

    """
    if day_count not in convs:
        # print(day_count, convs)
        raise Exception('day count error in create_bond_active_structures')
    dates = full_structures.get('dates')
    # print('full_structures==>>', full_structures)
    full_dates = full_structures['dates']
    full_coupons = full_structures['coupons']
    full_dcfs = full_structures['dcfs']
    full_cfs = full_structures['cfs']
    # full_times = full_structures['times']
    full_fvs = full_structures['fvs']
    full_fvs_flow = full_structures['fvs_flow']
    # pflow = full_structures['fvs_flow']
    date_bool = list(map(lambda x: value_date < x, dates))
    index = date_bool.index(True)

    active_dates = list(full_dates[index:])
    # print(active_dates)
    active_coupons = list(full_coupons[index:])
    active_dcfs = list(full_dcfs[index:])
    active_cfs = list(full_cfs[index:])
    active_fvs = list(full_fvs[index:])
    active_fvs_flow = list(full_fvs_flow[index:])
    accrued = 0
    last_coupon_date = None
    accrued_interest = 0
    if index >= 0:
        next_date = dt64(full_dates[index])
        last_coupon_date = dt64(full_dates[index - 1])
        accrued_dcf = float(day_cf(day_count, last_coupon_date, value_date,
                             bondmat_date=full_dates[-1],
                             next_coupon_date=next_date,
                             Frequency=12/fre[freq]))
        if value_date < next_date:
            accrued = accrued_dcf
            accrued_interest = float((accrued_dcf * full_coupons[index]
                                * full_fvs[index] * 0.01))
        else:
            last_coupon_date = next_date
            accrued = 0
            accrued_interest = 0
    # recalculating time from value date
    active_times = list(map(lambda x: day_cf('Actual/365', value_date, x,
                                             bondmat_date=active_dates[-1],
                                             next_coupon_date=x),
                            active_dates))

    return {'accrued': accrued, 'accrued_interest': accrued_interest,
            'last_coupon_date': last_coupon_date, 'coupons': active_coupons,
            'coupon_flow': active_cfs,
            'dates': active_dates, 'fvs': active_fvs,
            'fvs_flow': active_fvs_flow, 'times': active_times,
            'dcfs': active_dcfs}

#@njit(parallel=True)
def value_active_structures(value_date, yld, face_value, structure,
                            day_count='Actual/Actual', freq='Semi-Annual'):
    """
    Value an active cash flow structure of a bond.

    Parameters
    ----------
    value_date : numpy.datetime64
        Valuation date
    yld : float
        Yield to maturity of the cash flows
    face_value : float
        Face value of the bond
    structure : dict
        The active cash flow structure. The structure should have keys
         and data structure similar to those return by
         create_bond_act_structures.
    day_count : str, optional
        Day count convention. The default is 'Actual/Actual'.
    freq : TYPE, optional
        Coupon frequecy. The default is 'Semi-Annual'.


    Returns
    -------
    risk_stats : dict
        Keys to the dictionary are 'proceed', 'pvbp01', 'price_per_100',
        'modified duration', 'macaulay duration', 'duration', 'convexity',
        and 'active_structures'.

    """
    try:
        val_date = dt64(value_date)
    except Exception:
        raise Exception('Value date in value_active_structures is not a date')

    active_dates = list(structure.get('dates'))
    active_coupons = list(structure.get('coupons'))
    active_dcfs = list(structure.get('dcfs'))
    active_cfs = list(structure.get('coupon_flow'))
    active_fvs = list(structure.get('fvs'))
    active_fvs_flow = list(structure.get('fvs_flow'))
    last_coupon_date = structure.get('last_coupon_date')
    dfs = []
    dfs01 = []

    totalcf = len(active_cfs)
    accrued_dcf = day_cf(day_count, last_coupon_date, val_date,
                         bondmat_date=active_dates[-1],
                         next_coupon_date=active_dates[0],
                         Frequency=12/fre[freq])
    #dfs, dfs01 = _calc_df(active_cfs, active_dcfs, accrued_dcf, yld, )
    # Calculating the discount factor
    # print(valid_cfs)
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

    # present value of each coupon flow
    active_cpns_values = list(map(lambda cf, df: cf * df, active_cfs, dfs))
    # present value of each principal flow
    active_pflow_values = list(map(lambda df, flow:  flow * df, dfs,
                                   active_fvs_flow))
    # net present value of all the principal flows
    pflow_value = functools.reduce(lambda a, b: a + b, active_pflow_values)
    # summing the coupon pv to get the value of the bond
    cpn_values = functools.reduce(lambda a, b: a + b, active_cpns_values)
    bond_values = cpn_values + pflow_value

    # Values after one bp shift
    active_cpns_values01 = list(map(lambda cf, df: cf * df, active_cfs, dfs01))

    # present value of each principal flow
    active_pflow_values01 = list(map(lambda df, flow:  flow * df, dfs01,
                                     active_fvs_flow))

    # net present value of all the principal flows
    pflow_value01 = functools.reduce(lambda a, b: a + b, active_pflow_values01)

    # summing the coupon pv to get the value of the bond
    cpn_values01 = functools.reduce(lambda a, b: a + b, active_cpns_values01)
    bond_values01 = cpn_values01 + pflow_value01
    risk_stats = {}
    risk_stats['proceed'] = bond_values
    # present value of 1 bp
    risk_stats['pvbp01'] = bond_values01 - bond_values
    risk_stats['price_per_100'] = ((risk_stats['proceed'] -
                                    structure['accrued_interest']) /
                                   face_value * 100)
    # recalculating time from value date
    active_times = list(map(lambda x: day_cf('Actual/365', val_date, x,
                                             bondmat_date=active_dates[-1],
                                             next_coupon_date=x),
                            active_dates))
    macD = 0
    for i in range(len(active_times)):
        macD = (macD + active_times[i] *
                (active_cfs[i] + active_fvs_flow[i]) *
                dfs[i])

    macD = macD/risk_stats['proceed']

    risk_stats['modified duration'] = macD / (1 + (yld * 0.01)
                                              / (12 / fre[freq]))
    risk_stats['macaulay duration'] = macD

    derivatives = sbd(face_value, active_cfs, active_dcfs,
                      active_fvs_flow, yld)
    risk_stats['duration'] = derivatives['duration']
    risk_stats['convexity'] = derivatives['convexity']
    risk_stats['active_structures'] = {'coupons': active_coupons,
                                       'coupon_flow': active_cfs,
                                       'coupon_values': active_cpns_values,
                                       'dates': active_dates,
                                       'fvs': active_fvs,
                                       'fvs_flow': active_fvs_flow,
                                       'fvs_flow_value': active_pflow_values,
                                       'times': active_times, 'dfs': dfs,
                                       'dcfs': active_dcfs}
    return risk_stats


def value_active_structures_jit(value_date, yld, face_value, structure,
                                day_count='Actual/Actual', freq='Semi-Annual'):
    """
    Value an active cash flow structure of a bond.

    Parameters
    ----------
    value_date : numpy.datetime64
        Valuation date
    yld : float
        Yield to maturity of the cash flows
    face_value : float
        Face value of the bond
    structure : dict
        The active cash flow structure. The structure should have keys
         and data structure similar to those return by
         create_bond_act_structures.
    day_count : str, optional
        Day count convention. The default is 'Actual/Actual'.
    freq : TYPE, optional
        Coupon frequecy. The default is 'Semi-Annual'.


    Returns
    -------
    risk_stats : dict
        Keys to the dictionary are 'proceed', 'pvbp01', 'price_per_100',
        'modified duration', 'macaulay duration', 'duration', 'convexity',
        and 'active_structures'.

    """
    try:
        val_date = dt64(value_date)
    except Exception:
        raise Exception('Value date in value_active_structures is not a date')

    active_dates = np.array(structure.get('dates'))
    active_coupons = np.array(structure.get('coupons'))
    active_dcfs = np.array(structure.get('dcfs'))
    active_cfs = np.array(structure.get('coupon_flow'))
    active_fvs = np.array(structure.get('fvs'))
    active_fvs_flow = np.array(structure.get('fvs_flow'))
    last_coupon_date = structure.get('last_coupon_date')

    # totalcf = len(active_cfs)
    accrued_dcf = day_cf(day_count, last_coupon_date, val_date,
                         bondmat_date=active_dates[-1],
                         next_coupon_date=active_dates[0],
                         Frequency=12/fre[freq])
    dfs, dfs01 = _calc_df(np.array(active_cfs), np.array(active_dcfs),
                          accrued_dcf, yld)

    # print(_calc_risk(active_cfs, dfs, dfs01, active_fvs_flow))
    proceed, pvbp01, active_cpns_values, active_pflow_values = _calc_risk(
        np.array(active_cfs), np.array(dfs), np.array(dfs01),
        np.array(active_fvs_flow))

    risk_stats = {}
    risk_stats['proceed'] = proceed
    # present value of 1 bp
    risk_stats['pvbp01'] = pvbp01
    risk_stats['price_per_100'] = ((risk_stats['proceed'] -
                                    structure['accrued_interest']) /
                                   face_value * 100)
    # recalculating time from value date
    active_times = list(map(lambda x: day_cf('Actual/365', val_date, x,
                                             bondmat_date=active_dates[-1],
                                             next_coupon_date=x),
                            active_dates))
    macD = 0
    for i in range(len(active_times)):
        macD = (macD + active_times[i] *
                (active_cfs[i] + active_fvs_flow[i]) *
                dfs[i])

    macD = macD / risk_stats['proceed']

    risk_stats['modified duration'] = macD / (1 + (yld * 0.01)
                                              / (12 / fre[freq]))
    risk_stats['macaulay duration'] = macD

    derivatives = sbd(face_value, active_cfs, active_dcfs,
                      active_fvs_flow, yld)
    risk_stats['duration'] = derivatives['duration']
    risk_stats['convexity'] = derivatives['convexity']
    risk_stats['active_structures'] = {'coupons': active_coupons,
                                       'coupon_flow': active_cfs,
                                       'coupon_values': active_cpns_values,
                                       'dates': active_dates,
                                       'fvs': active_fvs,
                                       'fvs_flow': active_fvs_flow,
                                       'fvs_flow_value': active_pflow_values,
                                       'times': active_times, 'dfs': dfs,
                                       'dcfs': active_dcfs}
    return risk_stats


@njit
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
