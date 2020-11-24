.. faspy_docs documentation master file, created by
   sphinx-quickstart on Tue Nov 17 19:34:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============
Discount Curve
==============
The submodule is named discount_curve as can be imported in the following manner:

::
  import faspy.interestrate.discount_curve as dc

or:

::
  from faspy.interestrate import discount_curve as dc

There are several functions available in the submodule. The primary function in the module is *discount_factor_gen()*.

.. py:function:: discount_factor_gen(rates, return_type="time")

   Generate discount factors from the provided rates

   :param rates: a dictionary or an array of dictionaries with the following
                keys - value_date, st_busday, st_ratebasis, st_daycount, lt_busday, lt_frequency, lt_daycount and rates. 'rates' key in the dictionary is a dictionary of interest rate in percentage having the following keys - 1W,  2W, 3W, 1M, 2M, 3M, 4M, 5M, 6M, 9M, 12M, 1Y, 2Y, 3Y, 4Y, 5Y, 6Y, 7Y, 10Y, 15Y, 20Y, 30Y. Non listed rates' key will be ignored and ommission of several keys are acceptable.

   :param return_type: a string of either "time" or "days". Default is "time"

   :rtype: An array of dictionary with keys "times" and "df" (for
            return_type="time" or "days" and "df" for return_type="days". If rates is an array, the function will return an array of arrays.
