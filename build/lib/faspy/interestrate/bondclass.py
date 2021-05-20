#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:29:52 2020

@author: RMS671214
"""
from numpy import datetime64 as dt64
from .conventions import *
from .rmp_dates import day_count_factor as day_cf
from .rmp_bhelp import *



class Bond:
    def __init__(self):
        self.__risk_stats = {}

    @property
    def issue_date(self):
        return self.__issue_date

    @issue_date.setter
    def issue_date(self, issue_date):
        if issue_date is None:
            self.__issue_date = None
        else:
            try:
                self.__issue_date = dt64(issue_date)
            except:
                raise ValueError('Bond class error: issue_date is not a date')

    @property
    def last_coupon_date(self):
        return self.__last_coupon_date

    @last_coupon_date.setter
    def last_coupon_date(self, last_coupon_date):
        if last_coupon_date is None:
            self.__last_coupon_date = None
        else:
            try:
                self.__last_coupon_date = dt64(last_coupon_date)
            except:
                raise ValueError('Bond class error: last_coupon_date is '
                                 + 'not a date')

    @property
    def value_date(self):
        return self.__value_date

    @value_date.setter
    def value_date(self, value_date):
        if value_date is None:
            self.__value_date = None
        else:
            try:
                self.__value_date = dt64(value_date)

            except:
                raise ValueError('Bond class error: ' + value_date + ' passed'
                                 + ' for value date is not a date')


    @property
    def maturity(self):
        return self.__maturity

    @maturity.setter
    def maturity(self, maturity):
        if maturity is None:
            self.__maturity = maturity
        else:
            try:
                self.__maturity = dt64(maturity)
            except:
                raise ValueError('Bond class error: ' + maturity + ' passed'
                                 + ' for date argument is not a date')

    @property
    def day_count(self):
        return self.__day_count

    @day_count.setter
    def day_count(self, day_count):
        if day_count is None or day_count in day_counts:
            self.__day_count = day_count
        else:
            raise ValueError('Bond class error: ' + day_count + ' is not in'
                             + ' the list of conventions')


    @property
    def business_day(self):
        return self.__business_day

    @business_day.setter
    def business_day(self, business_day):
        if business_day is None or business_day in business_days:
            self.__business_day = business_day
        else:
            raise ValueError('Bond class error: ' + business_day + ' is not in'
                             + ' the list of conventions')


    @property
    def date_gen(self):
        return self.__date_gen

    @date_gen.setter
    def date_gen(self, date_gen):
        if date_gen is None or date_gen in date_gen_method:
            self.__date_gen = date_gen
        else:
            raise ValueError('Bond class error: ' + date_gen + ' is not'
                             + ' in the list of methods')


    @property
    def frequency(self):
        return self.__frequency

    @frequency.setter
    def frequency(self, frequency):
        if frequency is None or frequency in frequencies:
            self.__frequency = frequency
        else:
            raise ValueError('Bond class error: ' + frequency + ' is not in'
                             + ' the list of frequencies')


class FixBond(Bond):

    def __init__(self):
        super().__init__()
        self.__risk_stats = {}

    @property
    def coupon(self):
        return self.__coupon

    @coupon.setter
    def coupon(self, coupon):
        if coupon is None:
            self.__coupon = coupon
        else:
            try:
                self.__coupon = float(coupon)
            except:
                raise ValueError('Bond class error: ' + coupon + ' is not'
                                 + ' a number')


    @property
    def accrued_interest(self):
        return self.__accrued_interest

    @property
    def ytm(self):
        return self.__ytm

    @ytm.setter
    def ytm(self, ytm):
        if ytm is None:
            self.__ytm = ytm
        else:
            try:
                self.__ytm = float(ytm)
            except:
                raise ValueError('Bond class error: ' + ytm + ' is not'
                                 + ' a number')

    @property
    def priceper100(self):
        return self.__priceper100

    @priceper100.setter
    def priceper100(self, priceper100):
        if priceper100 is None:
            self.__priceper100 = priceper100
        else:
            try:
                self.__priceper100 = float(priceper100)
            except:
                raise ValueError('Bond class error: ' + priceper100
                                 + ' is not a number')

    @property
    def face_value(self):
        return self.__face_value

    @face_value.setter
    def face_value(self, face_value):
        if face_value is None:
            self.__face_value = face_value
        else:
            try:
                self.__face_value = float(face_value)
            except:
                raise ValueError('Bond class error: ' + face_value
                                 + ' is not a valid number for position')

    @property
    def proceed(self):
        return self.__proceed

    @property
    def risk_stats(self):
        return self.__risk_stats

    @risk_stats.setter
    def risk_stats(self, risk_stats):
        if risk_stats is None:
            self.__risk_stats = {}
        else:
            if isinstance(risk_stats, dict):
                self.__risk_stats = risk_stats
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a dictionary.')

    @property
    def full_structures(self):
        return self.__full_structures

    @full_structures.setter
    def full_structures(self, full_structures):
        if full_structures is None:
            self.__full_structures = []
        else:
            if isinstance(full_structures, list):
                self.__full_structures = full_structures
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a dictionary.')

    @property
    def active_structures(self):
        return self.__active_structures

    @active_structures.setter
    def active_structures(self, active_structures):
        if active_structures is None:
            self.__structures = []
        else:
            if isinstance(active_structures, list):
                self.__active_structures = active_structures
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a Structures class.')

    def __get_bonddict(self):
        bonddict = {}
        bonddict['issue_date'] = self.issue_date
        bonddict['maturity'] = self.maturity
        bonddict['value_date'] = self.value_date
        bonddict['date_generation'] = self.date_gen
        bonddict['frequency'] = self.frequency
        bonddict['day_count'] = self.day_count
        bonddict['business_day'] = self.business_day
        bonddict['coupon'] = self.coupon
        bonddict['face_value'] = self.face_value
        return bonddict


    def construct_bond(self):
        bonddict = self.__get_bonddict()
        try:
            construct_fixbond(bonddict)
            # print(bonddict)

            self.full_structures = bonddict['full_structures']
            self.active_structures = bonddict['active_structures']
        except ValueError as error:
            print(str(error))


    def calculate(self, reconstruct=False):
        bonddict = self.__get_bonddict()
        bonddict['full_structures'] = self.full_structures
        if reconstruct is True:
            self.construct_bond()

        if self.__ytm is None:
            raise ValueError('ytm value is not assigned')

        try:
            self.recalc_act_structures()
            bonddict['active_structures'] = self.active_structures
            fixbond_price(self.value_date, bonddict, self.ytm)
        except ValueError as myerror:
            print(str(myerror))

        self.__accrued_interest = bonddict['accrued_interest']
        self.__risk_stats['pvbp01'] = bonddict['pvbp01']
        self.__risk_stats['modified_duration'] = bonddict['modified duration']
        self.__risk_stats['macaulay_duration'] = bonddict['macaulay duration']
        self.__risk_stats['duration'] = bonddict['duration']
        self.__risk_stats['convexity'] = bonddict['convexity']
        self.__active_structures = bonddict['active_structures']
        self.__priceper100 = bonddict['price per 100']
        self.__proceed = bonddict['proceed']
        self.last_coupon_date = bonddict['last_coupon_date']


    def recalc_act_structures(self):
        if self.value_date is None:
            raise Exception("Fixed Bond need value date for constructing" +
                            " active structures")

        # an error when issue date is none is not handle
        elif self.value_date < self.issue_date:
            acs = list(self.full_structures)

        else:
            fs = self.full_structures
            booldate = list(map(lambda x: x['date'] > self.value_date, fs))
            index = booldate.index(True)  # exclude paid coupon period
            acs = list(fs[index - 1:])

            acs[0]['dcf'] = 0
            acs[1]['dcf'] = day_cf(self.day_count, self.value_date,
                                   acs[1]['date'],
                                   bondmat_date=acs[-1]['date'],
                                   next_coupon_date=acs[1]['date'],
                                   Frequency=12 /
                                   frequencies[self.frequency])
            alen = len(acs)

            for i in range(alen):
                acs[i]['time'] = day_cf('Actual/365', self.value_date,
                                        acs[i]['date'],
                                        bondmat_date=acs[-1]['date'],
                                        next_coupon_date=acs[i]['date'])
            self.active_structures = acs


class FRN(Bond):

    def __init__(self):
        super().__init__()
        self.__risk_stats = {}
        self.__full_structures = {}
        self.average_period = 1
        self.fixing_basis = 'Same Day'


    @property
    def margin(self):
        return self.__margin

    @margin.setter
    def margin(self, margin):
        if margin is None:
            self.__margin = margin
        else:
            try:
                self.__margin = float(margin)
            except:
                self.__margin = None
                raise ValueError('FRN class error: ' + margin + ' is not'
                                 + ' a number')

    @property
    def spread(self):
        return self.__spread

    @spread.setter
    def spread(self, spread):
        if spread is None:
            self.__spread = spread
        else:
            try:
                self.__spread = float(spread)
            except:
                raise ValueError('FRN class error: ' + spread + ' is not'
                                 + ' a number')

    @property
    def fixing_basis(self):
        return self.__fixing_basis

    @fixing_basis.setter
    def fixing_basis(self, fixing_basis):
        if fixing_basis is None:
            self.__fixing_basis = 'Same Day'
        elif fixing_basis in start_basis:
            self.__fixing_basis = fixing_basis
        else:
            self.__fixing_basis = 'Same Day'

    @property
    def average_period(self):
        return self.__average_period

    @average_period.setter
    def average_period(self, average_period):
        if average_period is None:
            self.__average_period = average_period
        else:
            try:
                self.__average_period = int(average_period)
            except:
                raise ValueError('FRN class error: ' + average_period +
                                 ' is not a number')

    @property
    def current_coupon(self):
        return self.__current_coupon

    @current_coupon.setter
    def current_coupon(self, current_coupon):
        if current_coupon is None:
            self.__current_coupon = current_coupon
        else:
            try:
                self.__current_coupon = float(current_coupon)
            except:
                raise ValueError('FRN class error: ' + current_coupon +
                                 ' is not a number')

    @property
    def accrued_interest(self):
        return self.__accrued_interest

    @property
    def priceper100(self):
        return self.__priceper100

    @priceper100.setter
    def priceper100(self, priceper100):
        if priceper100 is None:
            self.__priceper100 = priceper100
        else:
            try:
                self.__priceper100 = float(priceper100)
            except:
                raise ValueError('FRN class error: ' + priceper100
                                 + ' is not a number')

    @property
    def face_value(self):
        return self.__face_value

    @face_value.setter
    def face_value(self, face_value):
        if face_value is None:
            self.__face_value = face_value
        else:
            try:
                self.__face_value = float(face_value)
            except:
                raise ValueError('Bond class error: ' + face_value
                                 + ' is not a valid number for face value')

    @property
    def proceed(self):
        return self.__proceed

    @property
    def risk_stats(self):
        return self.__risk_stats

    @risk_stats.setter
    def risk_stats(self, risk_stats):
        if risk_stats is None:
            self.__risk_stats = {}
        else:
            if isinstance(risk_stats, dict):
                self.__risk_stats = risk_stats
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a dictionary.')

    @property
    def full_structures(self):
        return self.__full_structures

    @full_structures.setter
    def full_structures(self, full_structures):
        if full_structures is None:
            self.__full_structures = []
        else:
            if isinstance(full_structures, list):
                self.__full_structures = full_structures
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a dictionary.')

    @property
    def active_structures(self):
        return self.__active_structures

    @active_structures.setter
    def active_structures(self, active_structures):
        if active_structures is None:
            self.__active_structures = []
        else:
            if isinstance(active_structures, list):
                self.__active_structures = active_structures
            else:
                raise ValueError('Bond class error: argument passed is not'
                                 + ' a Structures class.')

    def __get_bonddict(self):
        bonddict = {}
        bonddict['issue_date'] = self.issue_date
        bonddict['maturity'] = self.maturity
        bonddict['value_date'] = self.value_date
        bonddict['date_generation'] = self.date_gen
        bonddict['frequency'] = self.frequency
        bonddict['day_count'] = self.day_count
        bonddict['business_day'] = self.business_day
        bonddict['face_value'] = self.face_value
        bonddict['margin'] = self.margin
        bonddict['spread'] = self.spread
        bonddict['current_coupon'] = self.current_coupon
        bonddict['fixing_basis'] = self.fixing_basis
        bonddict['average_period'] = self.average_period
        return bonddict


    def construct_bond(self):
        bonddict = self.__get_bonddict()
        try:
            construct_frn(bonddict)
            self.__full_structures = bonddict['full_structures']
            self.__active_structures = bonddict['active_structures']
        except ValueError as error:
            print(str(error))


    def recalc_act_structures(self):
        if self.value_date is None:
            raise Exception("FRN need value date for constructing" +
                            " active structures")

        elif self.value_date < self.issue_date:
            acs = list(self.full_structures)

        else:
            fs = self.full_structures
            booldate = list(map(lambda x: x['date'] > self.value_date, fs))
            index = booldate.index(True)  # exclude paid coupon period
            acs = list(fs[index - 1:])

            acs[0]['y_dcf'] = 0
            acs[0]['fixed'] = None
            acs[1]['y_dcf'] = day_cf(self.day_count, self.value_date,
                                   acs[1]['date'],
                                   bondmat_date=acs[-1]['date'],
                                   next_coupon_date=acs[1]['date'],
                                   Frequency=12 /
                                   frequencies[self.frequency])
            alen = len(acs)

            for i in range(alen):
                row = acs[i]
                row['time'] = day_cf('Actual/365', self.value_date,
                                        row['date'],
                                        bondmat_date=acs[-1]['date'],
                                        next_coupon_date=row['date'])
                if row['fixing_date'] is not None:
                    if ((row['fixing_date'] <= self.value_date) and (
                            self.value_date < row['date'])):
                        row['fixed'] = True
                        row['coupon'] = self.current_coupon
                    else:
                        row['fixed'] = False

            self.active_structures = acs


    def calculate(self, disc_curve, reconstruct=False):
        """


        Parameters
        ----------
        disc_curve : DFactor
            DFactor class containg information on the discount factors.
        reconstruct : boolean, optional
            Constructing the bonds coupon structure if set to true.
            The default is False.

        Returns
        -------
        None

        """
        bonddict = self.__get_bonddict()
        if reconstruct is True:
            self.construct_bond()
        try:
            self.recalc_act_structures()
            bonddict['active_structures'] = self.active_structures
            frn_price(bonddict, disc_curve)
        except ValueError as myerror:
            print(str(myerror))

        self.__accrued_interest = bonddict['accrued_interest']
        self.__risk_stats['pvbp01'] = bonddict['pvbp01']
        #self.__risk_stats['modified_duration'] = bonddict['modified duration']
        #self.__risk_stats['macaulay_duration'] = bonddict['macaulay duration']
        #self.__risk_stats['duration'] = bonddict['duration']
        #self.__risk_stats['convexity'] = bonddict['convexity']
        #self.__active_structures = bonddict['active_structures']
        self.__priceper100 = bonddict['price per 100']
        self.__proceed = bonddict['proceed']
