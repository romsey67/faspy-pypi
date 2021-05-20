#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:39:29 2020

@author: RMS671214
"""

from .rmp_dates import  day_count_factor as day_cf
from .rmp_dates import tenor_to_maturity as ttm
from .rmp_curves import generate_st_df as gen_stdf
from numpy import datetime64 as dt64, busdaycalendar as cal
from numpy import busday_offset as bus_off
from .rmp_curves import generate_fulldf as gen_fdf, \
    interpolation
from numba import njit
from .conventions import *

class Data:
    """
    A class to store single data by means of a single level/depth dictionary

    Attributes:
    ===========
    data: dict
        Single level dictionary. There can be various number of keys.
        Values can be any object required subject to the limitation that
        the class is meant for single level dictionary

    Methods:
    ========
    single_data:
        return the value for the key provided.
    """

    def __init__(self, data=None):
        # data should be a dictionary like {'O/N': True, 'T/N':False}
        # the dict may hold value for rate, boolean etc depending on the intended
        # use of the class
        self.data = data

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, newdata):
        if isinstance(newdata, dict):
            self.__data = newdata
        else:
            self.__data = {}


    def get(self, data_name):
        """
        Parameters
        ----------
        data_name : str
            key for the data attributes

        Returns
        -------
        obj
            Returns the object held by the key in the atributes. If the key
            is not found, it will return None.

        """
        return self.__data.get(data_name)




class HistoricalData:
    """

    A collection class for Data. It hold multiple values/records of Data object
    Attributes:
    ===========
    historical: dict
        A dictionary containing Data object

    Methods:
        addnewdata:
            method to add new data to the class

        view_data_as_dict:
            return multilevel dictionary of historical attribute

    """
    def __init__(self):
        # historical is a dictionary of Data
        self.historical = {}

    @property
    def historical(self):

        return self.__historical

    @historical.setter
    def historical(self,historical):

        if historical is None or isinstance(historical,dict):
            self.__historical = historical
        else:
            raise Exception ('historical data is not in a dictionary')

    def addnewdata(self,newdata):

        #historical is a dictionary of Data
        if isinstance(newdata,dict):
            for key in newdata:
                self.__historical[key] = newdata[key]
        else:
            raise Exception('HistoricalRates error: data given for addnewdata is not a dictionary')


    def view_data_as_dict(self):
        """
        converts historical attribute to full-fledge dictionary
        Returns
        -------
        result : dict
            The whole content of historical attributes as dict object

        """
        result = {}
        for key in self.__historical:
            tenor_data = self.__historical[key]
            if isinstance(tenor_data, Data):
                result[key] = tenor_data.data
            else:
                result[key] = self.__historical[key]

        return result


# %%

class STRate():
    _basis = ('Simple', 'Discount Rate')

    def __init__(self, start_date, maturity, rate,
                 day_count='Actual/365 Fixed', rate_basis='Money Market',
                 business_day='No Adjustment'):
        self.start_date = start_date
        self.maturity = maturity
        self.day_count = day_count
        self.rate_basis = rate_basis
        self.business_day = business_day
        self.rate = rate


    def discount_factor(self):
        if self.day_count not in day_counts:
            return None
        elif self.rate_basis not in self._basis:
            return None
        elif self.business_day not in business_days:
            return None

        d_count = dc_fac(self.day_count, self.start_date, self.maturity)
        result = None
        if self.rate_basis == self._basis[0]:
            result = 1/(1 + self.rate * d_count * 0.01)
        elif self.rate_basis == self._basis[1]:
            result = 1 - self.rate * d_count * 0.01


        return result

# %%
class STRt(Data):
    __basis = ('Simple', 'Discount Rate', 'Continuous')
    __termsorder = {'O/N': 1, 'T/N': 2, 'S/N': 2.5, '1W': 3, '2W': 4, '3W': 5, '1M': 6,
                    '2M': 7, '3M': 8, '4M': 9, '5M': 10, '6M': 11, '9M': 15,
                    '12M': 20}
    __termsrates = {'O/N': None, 'T/N': None, 'S/N': None, '1W': None, '2W': None,
                    '3W': None, '1M': None, '2M': None, '3M': None,
                    '4M': None, '5M': None, '6M': None, '9M': None,
                    '12M': None}
    __spotbasis = ('Same Day', 'Tom', 'Spot')

    def __init__(self):
        self.value_date = None
        self.rate_basis = 'Simple'
        self.day_count = None
        self.business_day = None
        self.data = {}
        self.__ratedict = dict(self.__termsrates)

    @property
    def current_date(self):
        return self.__current_date

    @current_date.setter
    def current_date(self, current_date):
        if current_date is not None:
            try:
                self.__current_date = dt64(current_date)
            except:
                raise ValueError('current date is not a date')
        else:
            self.__current_date = current_date

    @property
    def value_date(self):
        return self.__value_date

    @value_date.setter
    def value_date(self, value_date):
        if value_date is not None:
            try:
                self.__value_date = dt64(value_date)
            except:
                raise ValueError('start_date is not a date')
        else:
            self.__value_date = value_date

    @property
    def rate_basis(self):
        return self.__rate_basis

    @rate_basis.setter
    def rate_basis(self, rate_basis):

        if rate_basis is not None:
            if rate_basis in self.__basis:
                self.__rate_basis = rate_basis
            else:
                raise ValueError('rate_basis is not in the list of acceptable values')
        else:
            self.__rate_basis = 'Money Market'

    @property
    def start_basis(self):
        return self.__start_basis

    @start_basis.setter
    def start_basis(self, start_basis):

        if start_basis is not None:
            if start_basis in self.__spotbasis:
                self.__start_basis = start_basis
            else:
                raise ValueError('start_basis is not in the list of acceptable values')
        else:
            self.__start_basis = 'Same Day'

    @property
    def day_count(self):
        return self.__day_count

    @day_count.setter
    def day_count(self, day_count):
        if day_count is not None:
            if day_count in day_counts:
                self.__day_count = day_count
            else:
                raise ValueError('day_count is not in the list of acceptable values')
        else:
            self.__day_count = 'Actual/365'

    @property
    def business_day(self):
        return self.__business_day

    @business_day.setter
    def business_day(self, business_day):
        if business_day is not None:
            if business_day in business_days:
                self.__business_day = business_day
            else:
                raise ValueError('business_day is not in the list of acceptable values')
        else:
            self._business_day = 'Modified Following'

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, newdata):
        if isinstance(newdata, dict):
            tenors = []
            therates = []
            for k in newdata:
                if newdata.get(k) is not None:
                    self.__ratedict[k] = newdata[k]
                    tenors.append(k)
                    therates.append(newdata[k])
            self.__data = {'tenors': tenors, 'rates': therates}

        elif isinstance(newdata, Data):
            self.__data = newdata.data
        else:
            self.__data = {}

    def calcdf(self, holidays=[]):
        setting = {'cdate': self.__current_date,
                   'start_basis': self.start_basis,
                   'day_count': self.__day_count,
                   'bus_day': self.__business_day,
                   'rate_basis': self.__rate_basis}
        tenors =[]
        therates =[]
        for k in self.__ratedict:
            tenors.append(k)
            therates.append(self.__ratedict[k])
        result = _calc_shorttenor(tenors, therates,
                                  setting, holidays=holidays)
        return result


# %%
class STRates(Data):
    __basis = ('Simple', 'Discount Rate', 'Continuous')
    __termsrates = {'1W': None, '2W': None,
                    '3W': None, '1M': None, '2M': None, '3M': None,
                    '4M': None, '5M': None, '6M': None, '9M': None,
                    '12M': None}

    def __init__(self):
        self.start_date = None
        self.rate_basis = 'Simple'
        self.day_count = None
        self.business_day = None
        self.data = {}
        self.rates = dict(self.__termsrates)


    @property
    def start_date(self):
        return self.__start_date

    @start_date.setter
    def start_date(self, start_date):
        if start_date is not None:
            try:
                self.__start_date = dt64(start_date)
            except:
                raise ValueError('start_date is not a date')
        else:
            self.__start_date = start_date

    @property
    def rate_basis(self):
        return self.__rate_basis

    @rate_basis.setter
    def rate_basis(self, rate_basis):

        if rate_basis is not None:
            if rate_basis in self.__basis:
                self.__rate_basis = rate_basis
            else:
                raise ValueError('rate_basis is not in the list of acceptable values')
        else:
            self.__rate_basis = 'Simple'

    @property
    def day_count(self):
        return self.__day_count

    @day_count.setter
    def day_count(self, day_count):
        if day_count is not None:
            if day_count in day_counts:
                self.__day_count = day_count
            else:
                raise ValueError('day_count is not in the list of acceptable values')
        else:
            self.__day_count = 'Actual/365'

    @property
    def business_day(self):
        return self.__business_day

    @business_day.setter
    def business_day(self, business_day):
        if business_day is not None:
            if business_day in business_days:
                self.__business_day = business_day
            else:
                raise ValueError('business_day is not in the list of'
                                 + ' acceptable values')
        else:
            self._business_day = 'Modified Following'

    @property
    def rates(self):
        return self.__rates

    @rates.setter
    def rates(self, newdata):
        self.__rates = dict(self.__termsrates)
        if isinstance(newdata, dict):
            tenors = []
            therates = []
            for k in self.__termsrates:
                if newdata.get(k) is not None:
                    try:
                        float(newdata[k])
                        tenors.append(k)
                        therates.append(newdata[k])
                        self.__rates[k] = newdata[k]
                    except:
                        pass

        elif isinstance(newdata, Data):
            self.__rates = newdata.data

    def getdata(self):
        strates = {}
        for k in self.__rates:
            if self.__rates.get(k) is not None:
                strates[k] = self.__rates[k]

        return strates


    def calc_df(self):
        strates = self.getdata()
        result = gen_stdf(self.start_date, strates, self.day_count,
                          self.business_day, rate_basis=self.rate_basis)
        return result



# %%
class LTRates(Data):

    __termsrates = {'1Y': None, '2Y': None,
                    '3Y': None, '4Y': None, '5Y': None, '6Y': None,
                    '7Y': None, '10Y': None, '15Y': None, '20Y': None,
                    '30Y': None}

    def __init__(self):
        self.start_date = None
        self.frequency = 'Semi-Annual'
        self.day_count = None
        self.business_day = None
        self.__rates = dict(self.__termsrates)


    @property
    def start_date(self):
        return self.__start_date

    @start_date.setter
    def start_date(self,start_date):
        if start_date is not None:
            try:
                self.__start_date = dt64(start_date)
            except:
                raise ValueError('start_date is not a date')
        else:
            self.__start_date = start_date


    @property
    def frequency(self):
        return self.__frequency

    @frequency.setter
    def frequency(self,frequency):
        if frequency is not None:
            if frequency in frequencies:
                self.__frequency = frequency
            else:
                raise ValueError('frequency is not in the list of acceptable values')
        else:
            self.__frequency = 'Semi-Annual'


    @property
    def day_count(self):
        return self.__day_count

    @day_count.setter
    def day_count(self,day_count):
        if day_count is not None:
            if day_count in day_counts:
                self.__day_count = day_count
            else:
                raise ValueError('day_count is not in the list of acceptable values')
        else:
            self.__day_count = 'Actual/365'



    @property
    def business_day(self):
        return self.__business_day

    @business_day.setter
    def business_day(self,business_day):
        if business_day != None:
            if business_day in business_days:
                self.__business_day = business_day
            else:
                raise ValueError('business_day is not in the list of acceptable values')
        else:
            self._business_day = 'Modified Following'


    @property
    def rates(self):
        return self.__rates

    @rates.setter
    def rates(self, newdata):
        self.__rates = dict(self.__termsrates)
        if isinstance(newdata, dict):
            tenors = []
            therates = []
            for k in self.__termsrates:
                if newdata.get(k) is not None:
                    try:
                        float(newdata[k])
                        tenors.append(k)
                        therates.append(newdata[k])
                        self.__rates[k] = newdata[k]
                    except:
                        pass


        elif isinstance(newdata, Data):
            self.__rates = newdata.data


    def getdata(self):
        strates = {}
        for k in self.__rates:
            if self.__rates.get(k) is not None:
                strates[k] = self.__rates[k]

        return strates


#%%
class Rates:
    def __init__(self):
        self.strates = STRates()
        self.ltrates = LTRates()
        self.__value_date = None
        self.date_gen_method = date_gen_method[0]
        self.holidays = []
        self.__df = DFactor()
        self.__par = {}

    @property
    def value_date(self):
        return self.__value_date

    @value_date.setter
    def value_date(self, value_date):
        if value_date is not None:
            try:
                self.__value_date = dt64(value_date, "D")
            except:
                raise ValueError('value_date is not a date')
        else:
            self.__value_date = value_date

    @property
    def df(self):
        return self.__df

    @property
    def par(self):
        return self.__par


    def calcdf(self):
        st = self.strates
        lt = self.ltrates

        st_curve = st.getdata()
        lt_curve = lt.getdata()

        all_df = gen_fdf(self.value_date, st_curve, st.day_count,
                          st.business_day, st.rate_basis, lt_curve,
                          lt.day_count, lt.business_day,
                          frequencies[lt.frequency],
                          method=self.date_gen_method,
                          holidays=self.holidays)

        self.df.data = list(all_df)
        # print(list(all_df))
        #self.__par = all_df['par_rates']


    def interpolate(self, x_value, model='chip', method='slinear'):
        disfac = self.df.data
        x_axis = [x['time'] for x in disfac]
        #x_axis = self.par['times']
        y_axis = [y['df'] for y in disfac]
        #y_axis =self.par['rates']
        return interpolation(x_axis, y_axis, x_value, model=model,
                             method=method)


#%%
class DFactor(Data):
    def __init__(self):
        self.times = None
        self.day_counts = None
        self.rates = None
        self.dfs = None
        self.dates = None
        self.tenors = None
        self.data = []

    @property
    def day_counts(self):
        return self.__day_counts

    @day_counts.setter
    def day_counts(self,day_counts):
        self.__day_counts = day_counts

    @property
    def times(self):
        return self.__times

    @times.setter
    def times(self,times):
        self.__times = times


    @property
    def dfs(self):
        return self.__dfs

    @dfs.setter
    def dfs(self,dfs):
        self.__dfs = dfs

    @property
    def rates(self):
        return self.__rates

    @rates.setter
    def rates(self, rates):
        self.__rates = rates



    @property
    def dates(self):
        return self.__dates

    @dates.setter
    def dates(self, dates):
        self.__dates = dates


    @property
    def tenors(self):
        return self.__tenors

    @tenors.setter
    def tenors(self, tenors):
        self.__tenors = tenors



    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, newdata):
        if isinstance(newdata, list):
            self.__data = newdata
        elif isinstance(newdata, Data):
            self.__data = newdata.data
        else:
            self.__data = []

        #self._assign_attrs()



    def _assign_attrs(self):

        self.__times = self.__data.get('times')
        self.__dates = self.__data.get('dates')
        self.__rates = self.__data.get('rates')
        self.__dfs = self.__data.get('dfs')
        self.__tenors = self.__data.get('tenors')
        self.__day_counts = self.__data.get('dcfs')


    def interpolate(self, x_value, x_axis='times', model='chip',
                    method='slinear',is_function=False):
        if x_axis == 'times':
            myx = list(map(lambda x: x['time'], self.__data))
        elif x_axis == 'day_count':
            myx = list(map(lambda x: x['dcf'], self.__data))

        y_axis = list(map(lambda x: x['df'], self.__data))
        return interpolation(myx, y_axis, x_value, model=model, method=method, is_function=is_function)


class SwapPoint():
    __termsorder = {'O/N': 1, 'T/N': 2, 'S/N': 2.5, '1W': 3, '2W': 4, '3W': 5, '1M': 6,
                    '2M': 7, '3M': 8, '4M': 9, '5M': 10, '6M': 11, '9M': 15,
                    '12M': 20}
    __spotbasis = ('Same Day', 'Tom', 'Spot')
    __termsrates = {'O/N': None, 'T/N': None, 'S/N': None, '1W': None, '2W': None,
                    '3W': None, '1M': None, '2M': None, '3M': None,
                    '4M': None, '5M': None, '6M': None, '9M': None,
                    '12M': None}
    def __init__(self):
        self.multiplier = 10000
        self.__ratedict = dict(self.__termsrates)

    @property
    def fxrate(self):
        return self.__fxrate

    @fxrate.setter
    def fxrate(self, fxrate):
        try:
            self.__fxrate = float(fxrate)
        except:
            raise Exception('fxrate must be a number')

    @property
    def current_date(self):
        return self.__current_date

    @current_date.setter
    def current_date(self, current_date):
        if current_date is not None:
            try:
                self.__current_date = dt64(current_date)
            except:
                raise ValueError('current date is not a date')
        else:
            self.__current_date = current_date

    @property
    def money(self):
        return self.__money

    @money.setter
    def money(self, money):
        if isinstance(money, STRt):
            self.__money = money
        else:
            self.__money = None

    @property
    def base(self):
        return self.__base

    @base.setter
    def base(self, base):
        if isinstance(base, STRt):
            self.__base = base
        else:
            self.__base = None

    @property
    def multiplier(self):
        return self.__multiplier

    @multiplier.setter
    def multiplier(self, multiplier):
        try:
            self.__multiplier = float(multiplier)
        except:
            raise Exception('multiplier is not a number')


    @property
    def swappoint(self, swappoint):
        return self.__swappoint

    @swappoint.setter
    def swappoint(self, swappoint):
        if isinstance(swappoint, dict):
            tenors = []
            therates = []
            for k in swappoint:
                if swappoint.get(k) is not None and self.__termsorder.get(k) is not None:
                    self.__ratedict[k] = swappoint[k]
                    tenors.append(k)
                    therates.append(swappoint[k])
            self.__data = {'tenors': tenors, 'rates': therates}

        elif isinstance(swappoint, Data):
            self.__data = swappoint.data
        else:
            self.__data = {}

    @property
    def start_basis(self):
        return self.__start_basis

    @start_basis.setter
    def start_basis(self, start_basis):

        if start_basis is not None:
            if start_basis in self.__spotbasis:
                self.__start_basis = start_basis
            else:
                raise ValueError('start_basis is not in the list of acceptable values')
        else:
            self.__start_basis = 'Spot'

    @property
    def business_day(self):
        return self.__business_day

    @business_day.setter
    def business_day(self, business_day):
        if business_day is not None:
            if business_day in business_days:
                self.__business_day = business_day
            else:
                raise ValueError('business_day is not in the list of acceptable values')
        else:
            self._business_day = 'Modified Following'


    def calc_impliedswap(self, holidays = []):
        mclass = self.__money
        bclass = self.__base

        moneydf = mclass.calcdf()
        basedf = bclass.calcdf()

        mcurve = _getcurve_dict(moneydf)
        bcurve = _getcurve_dict(basedf)


        tenors =[]
        therates =[]
        for k in self.__ratedict:
            tenors.append(k)
            therates.append(self.__ratedict[k])

        tomdate, spotdate, spotnext = _calc_ondate(self.current_date, holidays=holidays)
        swap_arr = []
        for k in self.__ratedict:
            swap_dict = {}
            swap_dict['tenor'] = k
            swap_dict['point'] = self.__ratedict[k]

            if swap_dict['tenor'] == 'O/N':
                swap_dict['start'] = self.__current_date
                swap_dict['maturity'] = tomdate
                swap_dict['stime'] = 0
            elif swap_dict['tenor'] == 'T/N':
                swap_dict['start'] = tomdate
                swap_dict['maturity'] = spotdate
                swap_dict['stime'] = float(day_cf('Actual/365', self.__current_date, tomdate))
            elif swap_dict['tenor'] == 'S/N':
                swap_dict['start'] = spotdate
                swap_dict['maturity'] = spotnext
                swap_dict['stime'] = float(day_cf('Actual/365', self.__current_date, spotdate))
            else:
                if self.__start_basis == 'Same Day':
                    swap_dict['start'] = self.__current_date
                    swap_dict['stime'] = 0
                elif self.__start_basis == 'Tom':
                    swap_dict['start'] = tomdate
                    swap_dict['stime'] = float(day_cf('Actual/365', self.__current_date, tomdate))
                elif self.__start_basis == 'Spot':
                    swap_dict['start'] = spotdate
                    swap_dict['stime'] = float(day_cf('Actual/365', self.__current_date, spotdate))
                swap_dict['maturity'] = ttm(swap_dict['start'], swap_dict['tenor'], 'Actual/365', self.__business_day)

            stime = float(day_cf('Actual/365', self.__current_date, swap_dict['start']))
            time = float(day_cf('Actual/365', self.__current_date, swap_dict['maturity']))

            swap_dict['time'] = float(day_cf('Actual/365', self.__current_date, swap_dict['maturity']))
            mon_df = interpolation(mcurve['time'], mcurve['df'], swap_dict['time'])/interpolation(mcurve['time'], mcurve['df'], swap_dict['stime'])
            mon_dcf = day_cf(mclass.day_count, swap_dict['start'], swap_dict['maturity'] )
            swap_dict['mrate'] = (1/mon_df - 1) / (mon_dcf * 0.01)

            bas_df = interpolation(bcurve['time'], bcurve['df'], swap_dict['time'])/interpolation(bcurve['time'], bcurve['df'], swap_dict['stime'])
            bas_dcf = day_cf(bclass.day_count, swap_dict['start'], swap_dict['maturity'] )
            swap_dict['brate'] = (1/bas_df - 1) / (bas_dcf * 0.01)

            try:
                spt = float(swap_dict['point'])
                swap_dict['imrate']= _calc_impliedmoney(self.__fxrate, spt/ self.__multiplier, bas_df, mon_dcf)
                swap_dict['ibrate']= _calc_impliedbase(self.__fxrate, spt/ self.__multiplier, mon_df, bas_dcf)
            except:
                swap_dict['imrate'] = None
                swap_dict['ibrate']= None
            forward = self.__fxrate * bas_df/mon_df
            point = (forward - self.__fxrate) * self.__multiplier
            swap_dict['impoint'] = point



            swap_arr.append(swap_dict)

        return swap_arr


@njit('float64(float64, float64, float64, float64)')
def _calc_impliedmoney(fxrate, swappoint, basedf, moneydcf):
    # fwdrate = fxrate * basedf/mon_df
    # mondf = (fxrate * basedf)/fwdrate
    # 1 / (1 + monrate * moneydcf * 0.01) = (fxrate * basedf)/fwdrate
    # (1 + monrate * moneydcf * 0.01) = fwdrate/(fxrate * basedf)
    fwdrate = fxrate + swappoint
    return (fwdrate/(fxrate * basedf) - 1)/ ( moneydcf * 0.01)


@njit('float64(float64, float64, float64, float64)')
def _calc_impliedbase(fxrate, swappoint, moneydf, basedcf):
    # fwdrate = fxrate * basedf / mon_df
    # basedf = fwdrate * mondf / fxrate
    # 1 / (1 + baserate * basedcf * 0.01) = fwdrate * mondf / fxrate
    # (1 + baserate * basedcf * 0.01) = fxrate / (fwdrate * mondf)
    fwdrate = fxrate + swappoint
    return (fxrate / (fwdrate * moneydf) - 1 ) / (basedcf * 0.01)


def _calc_impliedswap(mcurve, bcurve, starttime, endtime, fxrate):
    basedf = list(bcurve['df'])
    basedf.insert(0,1)
    basetime = list(bcurve['time'])
    basetime.insert(0,0)

    moneydf = list(mcurve['df'])
    moneydf.insert(0,1)
    moneytime = list(mcurve['time'])
    moneytime.insert(0,0)

    points = []
    if isinstance(starttime,list):
        base_startdf = interpolation(basetime, basedf, starttime, model='chip', method=None )
        base_enddf = interpolation(basetime, basedf, endtime, model='chip', method=None )
        df_base = list(map(lambda start, end: end/start, base_startdf, base_enddf))

        money_startdf = interpolation(moneytime, moneydf, starttime, model='chip', method=None )
        money_enddf = interpolation(moneytime, money_startdf, endtime, model='chip', method=None )
        df_money = list(map(lambda start, end: end/start, money_startdf, money_enddf))

        forwards = list(map(lambda bdf,mdf: fxrate * bdf / mdf, df_base, df_money))
        points = list(map(lambda fwd: fwd - fxrate, forwards))
    return points

def _getcurve_dict(data_array):
    """
    An internal function to convert array of dictionary (eg [{'time': value0, 'df_tocdate': value1},...]) to dictionary containing array (eg {'time': [], 'df': []}). Used to converts arrays of discount factors to dictionary of discount factors. Source must have the following keys:
    i. time
    ii. df_tocdate
    Returns:
    ========
    dict: {'time':[], 'df': []}
    """
    datalen = len(data_array)
    time = []
    df =[]
    for i in range(datalen):
        datum = data_array[i]
        time.append(datum['time'])
        df.append(datum['df_tocdate'])
    time.insert(0,0)
    df.insert(0,1)

    return {'time': time, 'df': df}


def _calc_maturities(start_basis, day_count, bus_day, dates, tenors):
    if start_basis == 'Same Day':
        return ttm(dates['current'], tenors, day_count,
                         bus_day)
    elif start_basis == 'Tom':
        return ttm(dates['Tom'], tenors, day_count,
                         bus_day)
    elif start_basis == 'Spot':
        return ttm(dates['Spot'], tenors, day_count,
                         bus_day)


def _calc_ondate(current_date, holidays =[]):
    bdc = cal(weekmask='1111100', holidays=holidays)
    tomdate = bus_off(current_date, 1, roll='forward', busdaycal=bdc)
    spotdate = bus_off(tomdate, 1, roll='forward', busdaycal=bdc)
    spotnext = bus_off(spotdate, 1, roll='forward', busdaycal=bdc)
    return (tomdate,spotdate,spotnext)

# %%
def _calc_shorttenor(tenors, rates, setting, holidays=[]):

    v_date = dt64(setting['cdate'])
    tenors1 = tenors[3:]
    tenors2 = tenors[:3]
    tomdate, spotdate, spotnext = _calc_ondate(v_date, holidays=holidays)
    dates = {'current': v_date, 'Tom': tomdate, 'Spot': spotdate}
    data_len = len(tenors)
    data = []
    rates1 = rates[3:]
    rates2 = rates[:3]

    ratebasis = 0
    if setting['rate_basis'] == 'Discount Rate':
        ratebasis = 1

    data_len = len(tenors2)
    for i in range(data_len):
        datum = {}
        datum['tenor'] = tenors2[i]
        datum['rate'] = rates2[i]
        if datum['tenor'] == 'O/N':
            datum['start'] = v_date
            datum['maturity'] = tomdate
        elif datum['tenor'] == 'T/N':
            datum['start'] = tomdate
            datum['maturity'] = spotdate
        elif datum['tenor'] == 'S/N':
            datum['start'] = spotdate
            datum['maturity'] = spotnext

        datum['stime'] = float(day_cf('Actual/365', v_date,
                                    datum['start']))
        datum['dcf'] = float(day_cf(setting['day_count'], datum['start'],
                                    datum['maturity']))
        if datum['rate'] is not None:
            datum['df'] = _strate_to_df(datum['rate'], datum['dcf'], ratebasis)
        else:
            datum['df'] = None
        datum['time'] = day_cf('Actual/365', v_date, datum['maturity'])
        data.append(datum)

    maturities = _calc_maturities(setting['start_basis'], setting['day_count'], setting['bus_day'], dates, tenors1)
    start_date = None
    if setting['start_basis'] == 'Same Day':
        start_date = v_date
    elif setting['start_basis'] == 'Tom':
        start_date = tomdate
    elif setting['start_basis'] == 'Spot':
        start_date = spotdate

    stime = day_cf('Actual/365', v_date, start_date)
    data_len = len(tenors1)
    for i in range(data_len):
        datum = {}
        datum['tenor'] = tenors1[i]
        datum['rate'] = rates1[i]
        datum['maturity'] = maturities[i]
        datum['start'] = start_date
        datum['stime'] = stime
        datum['dcf'] = float(day_cf(setting['day_count'], datum['start'],
                                    datum['maturity']))
        if datum['rate'] is not None:
            datum['df'] = _strate_to_df(datum['rate'], datum['dcf'], ratebasis)
        else:
            datum['df'] = None
        datum['time'] = float(day_cf('Actual/365', datum['start'],
                                    datum['maturity']))
        data.append(datum)

    # Calculate df to current date
    data_len = len(data)
    if setting['start_basis'] == 'Same Day':
        for i in range(data_len):
            datum = data[i]
            if datum['tenor'] == 'O/N':
                if datum.get('df') is not None:
                    datum['df_tocdate'] = datum['df']
                else:
                    datum['df_tocdate'] = None
            elif datum['tenor'] == 'T/N' and datum.get('df') is not None:
                on = dict(data[i-1])
                if on.get('df') is not None:
                    datum['df_tocdate'] = on['df'] * datum['df']
                else:
                    on['rate'] = datum['rate']
                    on['df'] = 1 / (1 + on['rate'] * 0.01 * on['dcf'])
                    datum['df_tocdate'] = on['df'] * datum['df']

            elif datum['tenor'] == 'S/N' and datum.get('df') is not None:
                tn = dict(data[i-1])
                if tn.get('df_tocdate') is not None:
                    datum['df_tocdate'] = tn['df_tocdate'] * datum['df']
                else:
                    on = dict(data[i-2])
                    on['rate'] = datum['rate']
                    on['df'] = 1 / (1 + on['rate'] * 0.01 * on['dcf'])
                    on['df_tocdate'] = on['df']
                    tn['rate'] = datum['rate']
                    tn['df'] = 1 / (1 + tn['rate'] * 0.01 * tn['dcf'])
                    tn['df_tocdate'] = on['df'] * tn['df']
                    datum['df_tocdate'] = tn['df'] * datum['df']
            else:
                datum['df_tocdate'] = datum['df']

    elif setting['start_basis'] == 'Spot':
        _adjustdf_fromspot(data)

    elif setting['start_basis'] == 'Tom':
        df_totom = None
        ondf = None
        ontime = None

        for i in range(0, 1, 1):
            datum = data[i]

            if datum['tenor'] == 'O/N':
                ondf = datum['df']
                ontime = datum['time']
                if ondf is not None:
                    datum['df_tocdate'] = ondf

        tomtime = ontime
        if ondf is None:
            df_totom = ondf
            # get discount factor to interpolate
            x_axis = [0]
            y_axis = [1]
            for n in range(1, data_len, 1):
                datum1 = data[n]
                if datum1['df'] is not None:
                    x_axis.append(datum1['time'])
                    y_axis.append(datum1['df'])

            if ondf is None:
                x_value = tomtime
                df_totom = interpolation(x_axis, y_axis, x_value)


        for i in range(1, data_len, 1):
            datum = data[i]
            if datum['df'] is not None:
                datum['df_tocdate'] = datum['df'] * df_totom

    return data



def _adjustdf_fromspot(data):
    df_tospot = None
    ondf = None
    tndf = None
    ontime = None
    tntime = None
    data_len = len(data)
    for i in range(0, 2, 1):
        datum = data[i]
        if datum['tenor'] == 'O/N':
            ondf = datum['df']
            if ondf is not None:
                datum['df_tocdate'] = ondf
            ontime = datum['time']
        if datum['tenor'] == 'T/N':
            tndf = datum['df']
            tntime = datum['time']
            if ondf is not None and tndf is not None:
                datum['df_tocdate'] = ondf * tndf


    spottime = tntime
    if ondf is not None and tndf is not None:
        df_tospot = ondf * tndf
    else:
        # get discount factor to interpolate
        x_axis = [0]
        y_axis = [1]
        for n in range(2, data_len, 1):
            if datum1['df'] is not None:
                x_axis.append(datum1['time'])
                y_axis.append(datum1['df'])

        if ondf is None and tndf is None:
            x_value = spottime
            df_tospot = interpolation(x_axis, y_axis, x_value)

        elif ondf is None:
            x_value = ontime
            ondf = interpolation(x_axis, y_axis, x_value)
            df_tospot = ondf * tndf

        elif tndf is None:
            x_value = tntime - ontime
            tndf = interpolation(x_axis, y_axis, x_value)
            df_tospot = ondf * tndf

    for i in range(2, data_len, 1):
        datum = data[i]

        if datum['df'] is not None:
            datum['df_tocdate'] = datum['df'] * df_tospot
            if datum['tenor'] != 'S/N':
                datum['time'] = datum['time'] + datum['stime']


@njit('float64(float64, float64, int16)')
def _strate_to_df(rate, dcf, ratebasis):
    # Discount Rate
    if ratebasis == 1:
        return (1 - rate * 0.01 * dcf)
    # Money Market
    else:
        return (1 / (1 + rate * 0.01 * dcf))
