import numpy as np
from numpy import datetime64 as dt64
from .rmp_dates import tenor_to_maturity as ttm, \
day_count_factor as day_cf
from .rmp_dates import generate_dates as gen_dates
from scipy import interpolate
import sympy as sy
import numba
import math
import time
from collections import deque
from .conventions import frequencies



def generate_st_df(value_date, curve, day_count, business_day,
                   rate_basis='Simple', holidays=[]):
    new_curve = []
    for k in curve:
        tenor = {}
        tenor['tenor'] = k
        tenor['rate'] = curve[k]
        tenor['date'] = ttm(value_date, k, convention=day_count,
                                business_day=business_day,
                                holidays=holidays)
        tenor['dcf'] = day_cf(day_count, value_date, tenor['date'])
        tenor['time'] = day_cf('Actual/365', value_date, tenor['date'])
        tenor["days"] = (tenor["date"] - value_date).astype("int")
        if rate_basis == 'Simple':
            tenor['df'] = _mmr2df(tenor['rate'], tenor['dcf'])
        else:
            tenor['df'] = _dr2df(tenor['rate'], tenor['dcf'])
        new_curve.append(tenor)

    return new_curve


def generate_st_df_bymaturity(value_date, maturity, rate, convention,
                              business_day, rate_basis='Simple',
                              holidays=[]):
    result = {}
    dcf = day_cf(convention, value_date, maturity)

    time = day_cf('Actual/365', value_date, maturity)
    #print('dcf: ',dcf)
    df = None
    if rate_basis == 'Money Market':
        df = 1/(1 + rate * 0.01 * dcf)
    elif rate_basis == 'Discount Rate':
        df = 1 - rate * 0.01 * dcf
    #print('dfs in generate_st_df: ', dfs)

    result['dates'] = maturity
    result['dcfs'] = dcf
    result['dfs'] = df
    result['times'] = time

    return result


# In[10]:
def generate_fulldf(value_date, st_curve, st_daycount, st_business_day,
                    st_rate_basis, lt_curve, lt_daycount,
                    lt_business_day, frequency=6,
                    method="Forward from issue date",
                    holidays=[]):
    par_time = deque()
    par_rate = deque()
    df_time = deque()
    df_df = deque()
    st_ = generate_st_df(value_date, st_curve, st_daycount,
                         st_business_day, rate_basis=st_rate_basis)

    len_st = len(st_)
    st2compound = convert_shortrate_to_compounding
    for i in range(len_st):
        lst_ = st_[i]
        lst_['par_rate'] = st2compound(lst_['rate'], value_date,
                                       lst_['date'],
                                       frequency=12,
                                       compound_busday=lt_business_day,
                                       compound_dc=lt_daycount,
                                       compound_frequency=frequency)
        par_time.append(lst_['time'])
        par_rate.append(lst_['par_rate'])
        df_time.append(lst_['time'])
        df_df.append(lst_['df'])

    lt_ = deque()

    for k in lt_curve:
        tenor = {}
        tenor['tenor'] = k
        tenor['rate'] = lt_curve[k]
        tenor['date'] = ttm(value_date, k, convention=lt_daycount,
                                business_day=lt_business_day,
                                holidays=holidays)
        tenor['time'] = day_cf('Actual/365', value_date, tenor['date'])
        tenor["days"] = (tenor["date"] - value_date).astype("int")
        time_bool = list(map(lambda ptime: tenor['time'] >= ptime, par_time))

        try:
            index = time_bool.index(False)
            if tenor['time'] not in par_time:
                par_time.insert(index, tenor['time'])
                par_rate.insert(index, tenor['rate'])
        except ValueError:
            # print(tenor['time'] not in par_time, tenor['time'])
            if tenor['time'] not in par_time:
                par_time.append(tenor['time'])
                par_rate.append(tenor['rate'])
        lt_.append(tenor)
        # print(tenor)

    # print(lt_)
    lt_last_maturity = lt_[-1]['date']

    lt_alldates = gen_dates(value_date, lt_last_maturity, issueDate=value_date,
                            frequency=frequency, business_day=lt_business_day,
                            method=method, holidays=[])
    len_dates = len(lt_alldates)
    lt_all = deque()
    interp_par = interpolation(par_time, par_rate, 1, model='chip',
                               is_function=True)
    for i in range(len_dates):
        lt_sgl = {}
        thedate = lt_alldates[i]
        if thedate != value_date:
            lt_sgl['date'] = thedate
            lt_sgl["days"] = (thedate - value_date).astype("int")
            lt_sgl['time'] = day_cf('Actual/365', value_date, thedate)
            if i != len_dates - 1:
                lt_sgl['dcf'] = day_cf(lt_daycount, lt_alldates[i-1],
                                       thedate, bondmat_date=lt_alldates[-1],
                                       next_coupon_date=lt_alldates[i+1])
            else:
                lt_sgl['dcf'] = day_cf(lt_daycount, lt_alldates[i-1],
                                       thedate, bondmat_date=lt_alldates[-1],
                                       next_coupon_date=lt_alldates[-1])
            lt_sgl['rate'] = float(interp_par(lt_sgl['time']))
            lt_all.append(lt_sgl)

    # *********************************************************************
    #   BOOTSTRAPPING STARTS HERE
    # *********************************************************************
    len_all = len(lt_all)
    lt_all = list(lt_all)
    for i in range(len_all):
        lt_sgl = lt_all[i]
        ctime = lt_sgl['time']
        if ctime <= max(df_time):
            idf = interpolation(df_time, df_df, ctime)
            lt_sgl['df'] = idf
        else:

            dcfs = np.array([x['dcf'] * x['df'] for x in lt_all[:i]])
            new_dcfs = dcfs * (lt_sgl['rate'] * 0.01)
            values = np.sum(new_dcfs)
            cdf = (1 - values)/(1 + lt_sgl['rate'] * 0.01 * lt_sgl['dcf'])
            lt_sgl['df'] = cdf
            df_time.append(lt_sgl['time'])
            df_df.append(cdf)

    # merging the discount factors
    new_df = [*st_, *lt_all]
    def sorter(items):
        return items['time']

    new_df = sorted(new_df, key=sorter)
    lnew_df = len(new_df)
    dfs = deque()
    dfs.append(new_df[0])
    for i in range(1, lnew_df, 1):
        sgl_newdf = new_df[i]
        if sgl_newdf['time'] != new_df[i-1]['time']:
            dfs.append(sgl_newdf)

    #new_df = list(set(new_df))
    # return lt_all
    return dfs


def interpolation(x_axis, y_axis, x_value, model='chip', method=None,
                  is_function=False):
    methods = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
               'previous', 'next']
    # print(type(x_value),type(min(x_axis)))
    models = ['akima', 'chip', 'interp1d', 'cubicspline', 'krogh',
              'barycentric']
    # print('x_axis => ',x_axis)

    result = None
    f = None
    if model is None or model == 'interp1d':
        if method in methods:
            f = interpolate.interp1d(x_axis, y_axis, kind=method)
        else:
            f = interpolate.interp1d(x_axis, y_axis, kind='cubic')

    elif model in models:
        if model == 'akima':
            f = interpolate.Akima1DInterpolator(x_axis, y_axis)

        elif model == 'chip':
            f = interpolate.PchipInterpolator(x_axis, y_axis)

        elif model == 'cubicspline':
            f = interpolate.CubicSpline(x_axis, y_axis)

        elif model == 'krogh':
            f = interpolate.KroghInterpolator(x_axis, y_axis)

        elif model == 'barycentric':
            f = interpolate.BarycentricInterpolator(x_axis, y_axis)
    else:
        f = interpolate.PchipInterpolator(x_axis, y_axis)

    if is_function == True:
        return f
    else:
        if not isinstance(x_value, list):
            # if x_value <min(x_axis) or x_value >max(x_axis):
            #    raise Exception('interpolation error: value requested is outside of range')
            #    return result
            try:
                result = float(f(x_value))
            except:
                return result

        else:
            result = list(map(lambda x: float(f(x)),x_value))

        return result


def convert_shortrate_to_compounding(rate, start, end, frequency=12,
                                     convention='Actual/365',
                                     compound_busday='No Adjustment',
                                     rate_basis='Money Market',
                                     compound_dc='Actual/365',
                                     compound_frequency=12,
                                     method="Forward from issue date"):
    result = None
    fre = 12/frequency
    compound_fre = 12/compound_frequency

    dcf = day_cf(convention, start, end, Frequency=fre)
    df = None
    if rate_basis == 'Money Market':
        df = 1/(1 + rate * 0.01 * dcf)
    else:
        df = 1-rate * 0.01 * dcf
    two_year_date = ttm(start, '2Y',
                                          convention=compound_dc,
                                          business_day=compound_busday,
                                          holidays=[])
    li_cpn_dcf = []

    dates = gen_dates(start, two_year_date, issueDate=start,
                               frequency=compound_frequency,
                               business_day=compound_busday,
                               method=method, holidays=[])

    no_of_cpn = len(dates)

    for cpn_no in range(no_of_cpn-1):
        prev_cpn = dates[cpn_no]
        next_cpn = dates[cpn_no +1]
        if end <= prev_cpn:
            break
        elif next_cpn >= end:
            cur_date = end
        else:
            cur_date = next_cpn

        cmp_dcf = day_cf(compound_dc,prev_cpn, cur_date,
                               bondmat_date = dates[-1], next_coupon_date = next_cpn,
                               business_day = compound_busday,
                               Frequency = compound_fre)

        #print('start=>',start,'end=>',end,'cpn_no==>',cpn_no,'prev_cpn=>',prev_cpn,'cur_date==>', cur_date, 'nextcpn==>', next_cpn)
        #print('cmp_dcf==>>', cmp_dcf)
        li_cpn_dcf.append(cmp_dcf)

    #print('dates==>',dates,'enddate===>',end,'startdate===>',start,
    #      'df:', df,'li:',li_cpn_dcf)
    result = solver_rate_from_compounded_df(df,li_cpn_dcf)
    #print(type(result),len(result),result)
    result = float(result[0])

    #print(start,end,result)

    return result


def solver_rate_from_compounded_df(dis_factor,daycount_factors):
    rate = sy.Symbol('rate')
    fv = sy.Symbol('fv')
    #df = sy.Symbol('df')
    dcfs = daycount_factors
    #df = dis_factor

    for d in range(len(dcfs)):
        if d == 0:

            fv = (1 + (rate* 0.01*dcfs[d]) )
        else:
            fv = fv * (1 + (rate* 0.01*dcfs[d]) )

    fv = dis_factor * fv - 1
    #solved_rates = sy.solve(fv,rate)
    solved_rates = sy.solveset(fv, rate, domain=sy.S.Reals)
    new_rates = []
    #print('solved_rates==> ',solved_rates, type(solved_rates),list(solved_rates))
    solved_rates = list(solved_rates)
    for i in range(len(solved_rates)):
        try:
            float(solved_rates[i])
        except:
            continue
        new_rates.append(solved_rates[i])

    new_rates = list(filter(lambda x: x>0, new_rates))
    return new_rates


def calc_fwd_df(start, end, time_axis=None, df_axis=None, ifunc=None):

    if ifunc is not None:
        df1 = ifunc(start)
        df2 = ifunc(end)
        return df2 / df1

    elif (time_axis is not None and df_axis is not None):
        interp = interpolation(time_axis, df_axis, start, is_function=True)
        df1 = interp(start)
        df2 = interp(end)
        return df2 / df1
    else:
        return None


def calc_shortrate_from_df(startdate, enddate, df, day_count,
                           rate_basis="Money Market"):
    sdate = dt64(startdate, "D")
    edate = dt64(enddate, "D")
    dcf = day_cf(day_count, sdate, edate, bondmat_date=edate,
                 next_coupon_date=edate)
    if rate_basis == "Discount Rate":
        return _df2dr(df, dcf)
    else:
        return _df2mmr(float(df), float(dcf))
    return None


def calc_df_from_shortrate(startdate, enddate, rate, day_count,
                           rate_basis="Money Market"):
    sdate = dt64(startdate, "D")
    edate = dt64(enddate, "D")
    dcf = day_cf(day_count, sdate, edate, next_coupon_date=edate)
    if rate_basis == "Discount Rate":
        return _dr2df(rate, dcf)
    else:
        return _mmr2df(rate, dcf)
    return None


def continuous_rate(value_date, rate, tenor, day_count="Actual/365",
                    rate_basis="Money Market"):
    maturity = ttm(value_date, tenor, convention=day_count)
    time = day_cf("Actual/365", value_date, maturity)
    dcf = day_cf(day_count, value_date, maturity)
    if rate_basis == "Money_Market":
        df = _mmr2df(float(rate), float(dcf))
    elif rate_basis == "Discount Rate":
        df = _dr2df(float(rate), float(dcf))

    crate = -math.log(df)/time
    return crate


def shift_curve(curve, bp=0.01):
    
    rates = dict(curve)
    for key in rates:
        rates[key] += bp
    return rates


def discount_factor_from_zspread(value_date, date_structure, day_count,
                                 frequency, df_curve, zspread):
    data = [{"start_date": x["start_date"], "end_date": x["end_date"]}
            for x in date_structure if value_date < x["end_date"]]

    x_axis = [x["times"] for x in df_curve]
    y_axis = [x["df"] for x in df_curve]
    ifunc = interpolation(x_axis, y_axis, 1/366, is_function=True)
    maturity = data[-1]["end_date"]
    zdfs = [{"times": 0, "df": 1}]
    for datum in data:
        if value_date > datum["start_date"]:
            time1 = day_cf("Actual/365", value_date, datum["end_date"],
                           bondmat_date=maturity,
                           next_coupon_date=datum["end_date"])
            df = ifunc(time1)
            fwd_rate = calc_shortrate_from_df(value_date,
                                              datum["end_date"],
                                              df, day_count)
            zfwd_rate = fwd_rate + zspread
            zfwd_df = calc_df_from_shortrate(value_date, datum["end_date"],
                                             zfwd_rate, day_count)
            zdfs.append({"times": time1, "zfwd_df": zfwd_df})

        else:
            start_time = day_cf("Actual/365", value_date,
                                datum["start_date"],
                                bondmat_date=maturity,
                                next_coupon_date=datum["end_date"])
            end_time = day_cf("Actual/365", value_date, datum["end_date"],
                              bondmat_date=maturity,
                              next_coupon_date=datum["end_date"])
            fwd_df = calc_fwd_df(start_time, end_time, ifunc=ifunc)
            fwd_rate = calc_shortrate_from_df(datum["start_date"],
                                              datum["end_date"],
                                              fwd_df, day_count)

            zfwd_rate = fwd_rate + zspread
            zfwd_df = calc_df_from_shortrate(datum["start_date"],
                                             datum["end_date"],
                                             zfwd_rate, day_count)
            zdfs.append({"times": end_time, "zfwd_df": zfwd_df})

    data_len = len(zdfs)

    for i in range(1, data_len, 1):
        
        datum = zdfs[i]
        pdatum = zdfs[i-1]
        datum["df"] = pdatum["df"] * datum["zfwd_df"]

    zdf_curve = [{"times": x["times"], "df": x["df"]} for x in zdfs]
    

    return zdf_curve


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
