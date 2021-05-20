#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 14:08:14 2020

@author: RMS671214
"""

class Level:

    def  __init__(self, description, factor):
        self.__factor = factor
        self.__description = description
        self.__child = {}
        # self.__parent = None
        self.__amount = 0.00

    @property
    def factor(self):
        return self.__factor

    @factor.setter
    def factor(self, factor):
        self.__factor = factor

    @property
    def description(self):
        return self.__description

    @description.setter
    def description(self, description):
        self.__description = description

    @property
    def child(self):
        return self.__child

    @property
    def amount(self):
        return self.__amount

    @amount.setter
    def amount(self, amount):
        self.__amount = amount

    def value(self):
        val = 0
        level_size = len(self.__child)
        if level_size == 0 and self.__factor is not None:
            val += self.__amount * self.__factor

        else:
            for key in self.__child:
                val += self.__child[key].value()
                if isinstance(self.__factor, float):
                    val += self.__amount * self.__factor
        return val


    def get_array(self):

        level_size = len(self.__child)
        if level_size == 0:
            mydict = {}
            mydict["description"] = self.__description
            mydict["amount"] = self.__amount
            mydict["factor"] = self.__factor
            mydict["value"] = self.value()
            return mydict

        else:
            myarr = {"description": self.__description}
            for key in self.__child:
                myarr[key] = self.__child[key].get_array()
            return myarr
