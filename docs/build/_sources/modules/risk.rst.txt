.. faspy_docs documentation master file, created by
   sphinx-quickstart on Tue Nov 17 19:34:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================
Risk
=================
The initial purpose of the module is to calculate value at risk. In order to do so, various matrices of statistics measure have to be calculated  powered by numpy. Relevant functions are in a submodule named *var_utils*.
::
    import faspy.risk.var_utils as fvar

or:
::
    from faspy.risk import var_utils as fvar

The matrices calculated are:

* historical returns
* volatilities and correlation
* asset/cash flows weigthing

Historical Returns
-------------------

**Log Returns From Prices**

Calculate the log returns of prices. For interest rate product, please use discount factors as prices.

.. py:function:: returns_from_prices(prices, period=1)

  Calculate the log returns of prices

  :param prices: *numpy.array()*.

  :param period: *int*.

  :rtype: *numpy.array()*. An array of log returns

**Returns From Rates**

Calculate the returns using rates. If this method is used, it is recommended that the continuous rate is used instead of the normally quoted rate.

.. py:function:: returns_from_rates(prices, period=1)

  Calculate the rate difference between the comparison period

  :param prices: *numpy.array()*.

  :param period: *int*.

  :rtype: *numpy.array()*. An array of log returns


Volatilies and Correlations
---------------------------

**Volatilities**

Calculate volatilities (standard deviation) of returnS for the purpose of calculating value at risk.

.. py:function:: volatilities(data)

  Calculate the volatilities of data. Variables are the columns (*n*)

  :param data: *numpy.array()*. *m* x *n* array

  :rtype: *numpy.array()*. 1 x *n* array

**Correlations**

Calculates correlations of data. For the purpose of calculating value at risk the data are historical returns.

.. py:function:: correlation(data)

  Calculate the correlation of data.

  :param data: *numpy.array()*. *m* x *n* array

  :rtype: *numpy.array()*. n x *n* array


Weighting and Mapping
---------------------
As part of the process to calculate value at risk, cash flows of interest rate products have to be mapped to proper vertices. Once they are properly mapped, they will be weighted by the confidence interval, foreign exchange rate (in case of the foreign currencies cash flows) and basis point sensitivity (for interest rate products)

**Mapping Interest Rate Cash Flows**


.. py:function:: map_cf_to_var_vertices(cashflows, disfac, var_vertices, corr, vola)

  Map the cash flows provided vertices weighted by correlation, volatility and discount factor.

  :param cashflows: *list[{"times":, "pv":}]*.

  :param disfac: *list[{"times":, "pv":}]*.

  :param var_vertices: *list*. 1 x *n* array. The float value refers to time as fraction of a year or days.

  :param corr: *list*. *n* x *n* array of correlation.

  :param vola: *list*. *1* x *n* array of volatility.

  :rtype: *list*. 1 x *n* array of mapped cash flows.


**Weighting the Cash Flows**

  The weighting is required as part of the calculation of value at risk which use the following formula:

.. math::

  \begin{align}
    V&=\sqrt{WRW^T} \\
    \text{where } V &= \text{value at risk} \\
    W &= \text{ weighted cashflows/assets} \\
      &= \begin{bmatrix} w_1, w_2, ... w_n \end{bmatrix} \\
      w_i &= \text{ asset}_i \times z \times \sigma_{asset_i} \times fx \\
      z &= \text{ inverse cumulative distribution for the confidence interval }\\
      fx &= \text{ foreign exchange rate  as per unit of var currency} \\
      \sigma_{asset_i} &= \text{ volatility of} asset_i { returns} \\
      R &= \text{ correlation matrix} \\
      W_T &= \text{ transpose of matrix } W
  \end{align}

.. py:function:: var_asset_weightbyprice(assets, confidence, vol, fxrate=1)

  Weight the asset by confidence intervals, volatility and fx rate. Applicable for assets' volatilities measured by the return of their prices .

  :param assets: *ndarray*. 1 x n array of present values.
  :param confidence: *ndarray*. 1 x n array of confidence intervals.
  :param vol: *ndarray*. 1 x *n* array. 1 x n array of volatilities
  :param fxrate: *float*. Foreign exchange rate  as per unit of var currency. Default value is 1
  :rtype: *ndarray*. 1 x *n* array of weighted asset values

|

.. py:function:: var_asset_weightbyrate(assets, confidence, vol, bpvs, fxrate=1)

  Weight the asset by confidence intervals, volatility, present value of one basis point and fx rate. Applicable for assets' volatilities measured by the interest rate differentials over the observe period.

  :param assets: *ndarray*. 1 x n array of present values.
  :param confidence: *ndarray*. 1 x n array of confidence intervals.
  :param vol: *ndarray*. 1 x *n* array. 1 x n array of volatilities
  :param bpvs: *ndarray*. 1 x n array of present value of one basis point
  :param fxrate: *float*. Foreign exchange rate  as per unit of var currency. Default value is 1

  :rtype: *ndarray*. 1 x *n* array of weighted asset values
