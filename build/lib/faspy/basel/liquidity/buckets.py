import numba
from numpy import datetime64 as dt64
from .rmp_dates import day_count_factor as day_cf
from faspy.interestrate.rmp_curves import interpolation
import math


def bucketing(data, buckets, value_date, dfcurve=None):
    """
    Generate bonds coupon structures.

            Parameters:
                data: a list containing dictionary with the following keys -
                "id", "date" and "amount"
                buckets: a list containg dictionaries with the following keys -
                "name", "from" and "to". The value for each key is number of days


            Returns:
                a dictionary with the following keys - date, dcf, time,
                days, df, rate
    """
    vdate = dt64(value_date, "D")
    mydata = list(data)
    calc_data = [{"days": _datediff(vdate, dt64(x["date"], "D")).astype("int"),
                  "tenor": day_cf("Actual/365", vdate, dt64(x["date"], "D"))}
                 for x in mydata]
    if dfcurve is not None:
        df_xaxis = [x["times"] for x in dfcurve]
        df_yaxis = [x["df"] for x in dfcurve]
        i_func = interpolation(df_xaxis, df_yaxis, 1, is_function=True)
        df_array = [{"df": float(i_func(x["tenor"]))} for x in calc_data]


    dtlen = len(mydata)
    newdata = []
    for i in range(dtlen):
        old = dict(mydata[i])
        new = dict(calc_data[i])
        old.update(new)

        if dfcurve is not None:
            df_dic = df_array[i]
            old.update(df_dic)
            old["pv"] = old["df"] * old["amount"]
        newdata.append(old)

    mdur = [{"mod_duration": _mod_duration(float(x["tenor"]),
                                           float(x["df"]))}
            for x in newdata]
    bucket = [{"bucket": assign_bucket(x["tenor"], buckets)} for x in newdata]

    for i in range(dtlen):
        old = newdata[i]
        new = mdur[i]
        buc = bucket[i]
        old.update(new)
        old.update(buc)

    sortedbucket = sorted(buckets, key=sorter)
    buclen = len(sortedbucket)
    bucket_list = {}
    # return sortedbucket
    for i in range(buclen):
        bucket_list[sortedbucket[i]["name"]] = 0

    for datum in newdata:
        bucket_list[datum["bucket"]] += datum["pv"]

    return (newdata, bucket_list)


def sorter(item):
    return item["from"]


def assign_bucket(data, bucket):

    for buck in bucket:
        if data > buck["from"] and data <= buck["to"]:
            return buck["name"]

    return None


@numba.njit('float64(float64, float64)')
def _mod_duration(time, df):
    crate = -math.log(df)/time
    df01 = math.exp(-time * (crate + 0.0001))
    return (df01-df) / 0.0001


@numba.njit('float64(float64, float64)')
def _cal_crate(time, df):
    crate = -math.log(df)/time
    return float(crate)


@numba.njit
def _datediff(date1, date2):
    return date2 - date1
