import numpy as np
from numpy import datetime64 as dt64, timedelta64 as td64
import re
from datetime import datetime
import numba
from .conventions import *
from collections import deque


_a365 = ('Actual/365', 'Act/365', 'A/365', 'English')
_a365f = ('Actual/365 Fixed', 'Actual/365F', 'Act/365 Fixed', 'A/365 Fixed',
          'A/365F')
_a365l = ('Actual/365L', 'Act/365L', 'A/365L')
_a364 = ('Actual/364', 'Act/364', 'A/364')
_a360 = ('Actual/360', 'Act/360', 'A/360', 'French')
_a30_360 = ('30/360 Bond Basis', '30/360 US', '30E/360', '30E/360 ISDA')
_o30_360bb = ('30/360 Bond Basis', '30A/360')
_o30_360us = ('30/360', '30/360 US')
_a30e_360 = ('30E/360', '30/360 ICMA', '30S/360', 'Eurobond basis (ISDA 2006)',
             'Special German')
_a30e_360i = ('30E/360 ISDA', 'Eurobond basis (ISDA 2000)', 'German')
_A_A_icma = ('Actual/Actual', 'Actual/Actual ICMA', 'Act/Act ICMA', 'ISMA-99',
             'Act/Act ISMA')
_A_A_isda = ('Actual/Actual ISDA', 'Actual/365', 'Act/365')
# ,'Actual/Actual','Act/Act')
_A_A_afb = ('Actual/Actual AFB')


@numba.njit('float64(int32, int32, int32, int32, int32, int32)')
def _dcf_360(d1, m1, y1, d2, m2, y2):
    return (360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1))/360


def day_count_factor(convention, prev_date, current_date, bondmat_date=None,
                     next_coupon_date=None, business_day=None, Frequency=1):

    if convention in _a365:
        return _dcf_a365(prev_date, current_date)

    elif convention in _a360:
        return _dcf_a360(prev_date, current_date)

    elif convention in _a365f:

        return _dcf_a365f(dt64(prev_date, 'D'), dt64(current_date, 'D'))

    elif convention in _a364:
        return _dcf_a364(prev_date, current_date)

    elif convention in _a365l:
        return _dcf_a365l(prev_date, current_date, Frequency=Frequency)

    elif (convention in _a30_360 or convention in _o30_360bb or
          convention in _o30_360us or convention in _a30e_360 or
          convention in _a30e_360i):
        year1 = np.datetime64(prev_date, 'Y')
        month1 = np.datetime64(prev_date, 'M')
        d1 = np.int32(prev_date - np.datetime64(month1)) + 1
        m1 = np.int32(month1 - np.datetime64(year1, 'M')) + 1
        y1 = np.int64(year1)

        year2 = np.datetime64(current_date, 'Y')
        month2 = np.datetime64(current_date, 'M')
        d2 = np.int32(current_date - np.datetime64(month2)) + 1
        m2 = np.int32(month2 - np.datetime64(year2, 'M')) + 1
        y2 = np.int64(year2)

        if convention in _o30_360bb:
            d1 = min(d1, 30)
            if d1 >= 30:
                d2 = min(d2, 30)
            return _dcf_360(d1, m1, y1, d2, m2, y2)

        elif convention in _o30_360us:
            str_date1 = str(prev_date)
            str_date2 = str(current_date)
            arr_date1 = str_date1.split('-')
            arr_date2 = str_date2.split('-')

            date1m = np.datetime64(prev_date, 'M')
            days = np.int32(prev_date - np.datetime64(date1m, 'D'))
            addmonth = date1m + 1
            nextmonth = np.datetime64(addmonth, 'D') + days
            daysinamonth1 = np.int32(nextmonth - prev_date)

            date2m = np.datetime64(current_date, 'M')
            days2 = np.int32(current_date - np.datetime64(date2m, 'D'))
            addmonth2 = date2m + 1
            nextmonth2 = np.datetime64(addmonth2, 'D') + days2
            daysinamonth2 = np.int32(nextmonth2 - current_date)

            if (np.int32(arr_date1[1]) == 2 and np.int32(arr_date2[1]) == 2
                and days == daysinamonth1 and days2 == daysinamonth2 and
                business_day == 'EOM'):
                d2 = 30
            if (np.int32(arr_date1[1]) == 2 and days == daysinamonth1 and
                business_day == 'End of Month'):
                d1 = 30
            if d2 == 31 and d1 >= 30:
                d2 = 30
            if d1 == 31:
                d1 = 30

            return _dcf_360(d1, m1, y1, d2, m2, y2)

        elif convention in _a30e_360:
            if d1 == 31:
                d1 = 30
            if d2 == 31:
                d2 = 30
            value = _dcf_360(d1, m1, y1, d2, m2, y2)
            return value

        elif convention in _a30e_360i:
            daysinamonth = daysinthemonth(prev_date)
            monthend = np.datetime64(month1, 'D') + daysinamonth - 1
            if monthend == prev_date:
                d1 = 30

            daysinamonth2 = daysinthemonth(current_date)
            monthend2 = np.datetime64(month2, 'D') + daysinamonth2 - 1
            if bondmat_date != monthend2 and m2 != 2:
                d2 = 30

            return _dcf_360(d1, m1, y1, d2, m2, y2)

    elif convention in _A_A_icma:
        if next_coupon_date is None:
            return None
        else:
            accrued_days = _datediff(prev_date, current_date).astype('int')
            coupon_days = _datediff(prev_date, next_coupon_date).astype('int')
            dcf = accrued_days/(Frequency * coupon_days)

            return dcf

    elif convention in _A_A_isda:
        return _dcf_a365(prev_date, current_date, business_day=business_day,
                         Frequency=Frequency)

    elif convention in _A_A_afb:
        return "In development"


def _dcf_a365_4ayear(start, end, business_day=None, Frequency=1):

    year1 = np.datetime64(start, 'Y')
    year2 = np.datetime64(end, 'Y')
    if year2 == year1:
        daysinayear = daysintheyear(start)
        days = np.int64(end - start)
        return days/daysinayear
    else:
        daysinayear1 = daysintheyear(start)
        end1 = np.datetime64(year1 + 1, 'D')
        days1 = np.int64(end1 - start)

        daysinayear2 = daysintheyear(end)
        start2 = np.datetime64(year2, 'D')
        days2 = np.int64(end - start2)
        return days1/daysinayear1 + days2/daysinayear2


def _dcf_a365(prev_date, current_date, business_day=None, Frequency=1):
    year1 = np.datetime64(prev_date, 'Y')
    year2 = np.datetime64(current_date, 'Y')
    years = np.int64(year2 - year1)

    if years == 0:
        return _dcf_a365_4ayear(prev_date, current_date,
                                business_day=business_day,
                                Frequency=Frequency)
    else:
        start_m = np.datetime64(prev_date, 'M')
        start_day = np.int64(np.datetime64(prev_date, 'D') -
                             np.datetime64(start_m, 'D'))

        years_enddate_m = np.datetime64(prev_date, 'M') + years * 12
        years_enddate = np.datetime64(years_enddate_m, 'D') + start_day

        extra_days = np.int64(current_date - years_enddate)

        if years == 1:
            if extra_days == 0:
                return years
            elif extra_days < 0:
                return _dcf_a365_4ayear(prev_date, current_date,
                                        business_day=business_day,
                                        Frequency=Frequency)
            else:
                yearfrac = _dcf_a365_4ayear(years_enddate, current_date,
                                            business_day=business_day,
                                            Frequency=Frequency)
                return years + yearfrac
        else:
            if extra_days == 0:
                return years
            else:
                yearfrac = _dcf_a365_4ayear(years_enddate, current_date,
                                            business_day=business_day,
                                            Frequency=Frequency)
                return years + yearfrac


def _dcf_a365f(prev_date, current_date):

    datediff = _datediff(prev_date, current_date)
    days = datediff.astype('int')
    return days/365


def _dcf_a364(prev_date, current_date, business_day=None):
    days = np.int32(current_date - prev_date)
    return days/364


# Cant use this => @numba.njit('float64(numba.datetime, numba.datetime)')

def _dcf_a360(prev_date, current_date):
    datediff = _datediff(prev_date, current_date)
    days = datediff.astype('int')
    return days/360


def _dcf_a365l(prev_date, current_date, business_day=None, Frequency=1):

    year1 = np.datetime64(prev_date, 'Y')
    year2 = np.datetime64(current_date, 'Y')
    days = np.int32(current_date - prev_date)
    # check for leap year
    if Frequency == 1:
        start = np.datetime64(year1, 'D')
        end = np.datetime64(year1 + 1, 'D')
        daysinayear = np.int32(end - start)
        if daysinayear == 365:
            return days/365
        else:
            leapdate = np.datetime64(year1, 'M') + 1
            leapdate = np.datetime64(leapdate, 'D') + 28
            if prev_date < leapdate and current_date > leapdate:
                return days/366
            else:
                return days/365
    else:
        start2 = np.datetime64(year2, 'D')
        end2 = np.datetime64(year2 + 1, 'D')
        daysinayear2 = np.int32(end2 - start2)

    return days/daysinayear2


def daysinthemonth(date1):
    start = np.datetime64(date1, 'M')
    end = start + 1
    datediff = _datediff(np.datetime64(start, 'D'), np.datetime64(end, 'D'))
    return datediff.astype(int)


def daysintheyear(date1):
    start = np.datetime64(date1, 'Y')
    end = start + 1
    return np.int32(np.datetime64(end, 'D') - np.datetime64(start, 'D'))


def npdate_to_datetime(npdate):
    strdate = str(npdate)
    dtdate = datetime.strptime(strdate, "%Y-%m-%d").date()
    return dtdate


def movedatebyyear(your_date, no_of_year=1):
    start_m = np.datetime64(your_date, 'M')
    end_m = start_m + 12 * no_of_year
    days = np.datetime64(your_date, 'D') - np.datetime64(start_m, 'D')
    enddate = np.datetime64(end_m, 'D') + days
    return enddate


def movedatebymonth(your_date, no_of_month=1, business_day='No Adjustment',
                    holidays=[]):

    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)
    start_m = dt64(your_date, 'M')
    day = dt64(your_date, 'D') - dt64(start_m, 'D')
    end_m = start_m + no_of_month
    next_m = end_m + 1
    days_in_end_m = dt64(next_m, 'D') - dt64(end_m, 'D')

    if days_in_end_m > day:
        enddate = dt64(end_m, 'D') + day
    else:
        enddate = dt64(end_m, 'D') + days_in_end_m - 1

    if business_day == 'No Adjustment' or business_day is None:
        return enddate

    elif business_day == 'Following':
        return np.busday_offset(enddate, 0, roll='forward', busdaycal=bdc)

    elif business_day == 'Preceeding':
        return np.busday_offset(enddate, 0, roll='backward', busdaycal=bdc)

    elif business_day == 'Modified Following':
        new_date = np.busday_offset(enddate, 0, roll='forward', busdaycal=bdc)
        if dt64(new_date, 'M') > end_m:
            return np.busday_offset(enddate, 0, roll='backward', busdaycal=bdc)
        else:
            return new_date


def generate_dates(currentDate, maturity, issueDate=None, frequency=6,
                   business_day='No Adjustment',
                   method='Forward from issue date', holidays=[],
                   cpnrate=None, ytm=None):

    dates = deque()

    if method == 'Backward from maturity date':
        cpnDate = maturity
        dates.append(cpnDate)
        counter = 0
        if issueDate is None:
            while currentDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-counter * frequency,
                                          business_day=business_day,
                                          holidays=holidays)
                dates.insert(0, cpnDate)

        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")
            while issueDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if issueDate >= cpnDate:
                    dates.insert(0, issueDate)
                else:
                    dates.insert(0, cpnDate)

    elif method == 'Forward from issue date':
        if issueDate is None:
            cpnDate = currentDate
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(currentDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                dates.append(cpnDate)
        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")

            cpnDate = issueDate
            dates.append(cpnDate)
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(issueDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if maturity < cpnDate:
                    dates.append(maturity)
                else:
                    dates.append(cpnDate)

    return list(dates)



def generate_datesv102(currentDate, maturity, issueDate=None, frequency=6,
                   business_day='No Adjustment',
                   method='Forward from issue date', holidays=[],
                   cpnrate=None, ytm=None):

    dates = []

    if method == 'Backward from maturity date':
        cpnDate = maturity
        dates.append({'date': cpnDate})
        counter = 0
        if issueDate is None:
            while currentDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-counter * frequency,
                                          business_day=business_day,
                                          holidays=holidays)
                dates.insert(0, {'date': cpnDate})

        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")
            while issueDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if issueDate >= cpnDate:
                    dates.insert(0, {'date': issueDate})
                else:
                    dates.insert(0, {'date': cpnDate})

    elif method == 'Forward from issue date':
        if issueDate is None:
            cpnDate = currentDate
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(currentDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                dates.append({'date': cpnDate})
        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")
            cpnDate = issueDate
            dates.append({'date': cpnDate})
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(issueDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if maturity < cpnDate:
                    dates.append({'date': maturity})
                else:
                    dates.append({'date': cpnDate})

    return dates

def generate_dates_v101(currentDate, maturity, issueDate=None, frequency=6,
                   business_day='No Adjustment',
                   method='Forward from issue date', holidays=[],
                   cpnrate=None, ytm=None):

    dates = []

    if method == 'Backward from maturity date':

        cpnDate = maturity
        dates.append(cpnDate)
        counter = 0
        if issueDate is None:
            while currentDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-counter * frequency,
                                          business_day=business_day,
                                          holidays=holidays)
                dates.insert(0, cpnDate)

        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")
            while issueDate < cpnDate:
                counter += 1
                cpnDate = movedatebymonth(maturity,
                                          no_of_month=-int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if issueDate >= cpnDate:
                    dates.insert(0, issueDate)
                else:
                    dates.insert(0, cpnDate)

    elif method == 'Forward from issue date':
        if issueDate is None:
            cpnDate = currentDate
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(currentDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                dates.append(cpnDate)
        else:
            if type(issueDate) != np.datetime64:
                raise Exception("issueDate is not a date!")
            cpnDate = issueDate
            dates.append(cpnDate)
            counter = 0
            while maturity > cpnDate:
                counter += 1
                cpnDate = movedatebymonth(issueDate,
                                          no_of_month=int(counter * frequency),
                                          business_day=business_day,
                                          holidays=holidays)
                if maturity < cpnDate:
                    dates.append(maturity)
                else:
                    dates.append(cpnDate)

    return dates



def tenor_to_maturity(start_date, tenors, convention='Actual/365 Fixed',
                      business_day='No Adjustment', holidays=[]):

    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)

    if isinstance(tenors, list):
        result = []
        for tenor in tenors:
            mat_date = None

            if 'w' in tenor or 'W' in tenor:
                tenor_len = [float(number) for number in re.findall(r'-?\d+\.?\d*',
                                                                    tenor)]
                if len(tenor_len) == 0:
                    mat_date = None
                    result.append(mat_date)
                else:
                    number = tenor_len[0]
                    number = np.int32(number)
                    mat_date = np.datetime64(start_date, 'D') + 7 * number
                    mat_date = np.busday_offset(mat_date, 0, roll='forward',
                                                busdaycal=bdc)
                    result.append(mat_date)

            elif 'm' in tenor or 'M' in tenor:
                mat_date = mat_tenor_by_month(start_date, tenor,
                                              convention=convention,
                                              business_day=business_day,
                                              holidays=holidays)
                result.append(mat_date)

            elif 'y' in tenor or 'Y' in tenor:
                mat_date = mat_tenor_by_month(start_date, tenor,
                                              convention=convention,
                                              business_day=business_day,
                                              holidays=holidays,
                                              multiplier=12)
                result.append(mat_date)
            elif tenor == 'O/N':
                mat_date = np.busday_offset(start_date, 1, roll='forward',
                                            busdaycal=bdc)
                result.append(mat_date)

        return result

    else:
        result = tenor_to_maturity(start_date, [tenors], convention=convention,
                                   business_day=business_day,
                                   holidays=holidays)
        return result[0]


def mat_tenor_by_month(start_date, tenor, convention=None, business_day=None,
                       holidays=[], multiplier=1):
    tenor_len = [float(number) for number in re.findall(r'-?\d+\.?\d*', tenor)]
    if len(tenor_len) == 0:
        return None
    else:
        number = tenor_len[0]

        try:
            number = np.int32(number) * multiplier
        except:
            return None


        new_date = movedatebymonth(start_date, no_of_month=number,
                                   business_day=business_day,
                                   holidays=holidays)
        return new_date


@numba.njit
def _datediff(date1, date2):
    return date2 - date1


def adjustdates(months, dates, business_day, holidays=[]):

    bdc = np.busdaycalendar(weekmask='1111100', holidays=holidays)
    datalen = len(dates)
    newdates = []

    for i in range(datalen):
        month = months[i]
        dmonth = dt64(dates[i], 'M')
        if dmonth == month:
            if business_day == 'Following':
                dates[i] = np.busday_offset(dates[i], 0, roll='forward',
                                           busdaycal=bdc)

            elif business_day == 'Preceeding':
                dates[i] = np.busday_offset(dates[i], 0, roll='backward',
                                           busdaycal=bdc)

            elif business_day == 'Modified Following':
                new_date = np.busday_offset(dates[i], 0, roll='forward',
                                            busdaycal=bdc)
                if dt64(new_date, 'M') > month:
                    dates[i] = np.busday_offset(dates[i], 0, roll='backward',
                                            busdaycal=bdc)
        elif dmonth > month:
            #print(month,dmonth)
            nextmonth = month + td64(1, 'M')
            nextmonthdate = dt64(nextmonth, 'D')
            days = _datediff(dt64(month, 'D'), nextmonthdate)
            dates[i] = dt64(month, 'D') + days - 1
            if business_day == 'Following':
                dates[i] = np.busday_offset(dates[i], 0, roll='forward',
                                           busdaycal=bdc)

            elif business_day == 'Preceeding':
                dates[i] = np.busday_offset(dates[i], 0, roll='backward',
                                           busdaycal=bdc)

            elif business_day == 'Modified Following':
                new_date = np.busday_offset(dates[i], 0, roll='forward',
                                            busdaycal=bdc)
                if dt64(new_date, 'M') > month:
                    dates[i] = np.busday_offset(dates[i], 0, roll='backward',
                                            busdaycal=bdc)
