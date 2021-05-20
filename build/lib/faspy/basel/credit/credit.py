#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 16:20:19 2020

@author: RMS671214
"""

from faspy.interestrate.rmp_stats import rmpcdf, rmpppf
import math
from numba import njit
from faspy.interestrate import rmp_curves as rcurve
from collections import deque

def expected_loss(prob_default, exposure_at_default, recovery_rate=None,
                  loss_given_default=None):
    """
    Calculate the expected loss.

    Parameters
    ----------
    prob_default : float
        Probability of default.
    exposure_at_default : float
        Exposure at default.
    recovery_rate : float
        Recovery rate. Defaulted at None. If None 'loss_given_default' must be provided
    loss_given_default : float
        Loss given default. Defaulted at None. If None, recovery rate must be provided.

    Returns
    -------
    float
        Expected Loss.

    """
    try:
        rr = float(recovery_rate)
        try:
            el = prob_default * (1 - rr) * exposure_at_default
            return el
        except Exception as exc:
            return str(exc)
    except:
        try:
            lgd = float(loss_given_default)
            try:
                el = prob_default * lgd * exposure_at_default
                return el
            except Exception as exc:
                return str(exc)
        except Exception as exc:
            return str(exc)


def mean_rr_with_betadist(alpha, beta):
    return alpha / (alpha + beta)


def basel_cap_req(loss_given_default, prob_default, corr, confidence=0.999,
                 mat_adj=1):
    d1 = (rmpppf(prob_default) * math.sqrt(1 / (1 - corr)) +
          rmpppf(confidence) * math.sqrt(corr / (1 - corr)))
    cond_exp_loss = loss_given_default * rmpcdf(d1)
    exp_loss = loss_given_default * prob_default
    cap_req = (cond_exp_loss - exp_loss) * mat_adj
    return cap_req


def maturity_adjustment(prob_default, maturity, fixed_maturity=2.5):
    b_pd = math.pow(0.11852 - 0.05478 * math.log(prob_default), 2)
    mat_adj = (1 + (maturity - fixed_maturity) + b_pd) / (1 - 1.5 * b_pd)
    return mat_adj


def rwa(exp_at_default, car, cap_req):
    return cap_req * exp_at_default / car


def credit_correlation(corr_lo_pd, corr_hi_pd, prob_default, kfactor,
                       size_adj=0):
    factor_weight = ((1 - math.exp(-kfactor * prob_default)) /
                     (1 - math.exp(-kfactor)))
    w1 = corr_hi_pd * factor_weight
    w2 = corr_lo_pd * factor_weight

    return w1 + w2 - size_adj


def size_adjustment(sales, min_sales=5, max_sales=50):
    if sales <= min_sales:
        return 0.04
    elif sales >= max_sales:
        return 0.00s
    else:
        return 0.04 * (1 - (sales - min_sales)/(max_sales - min_sales))


def calc_prob_default(data):
    """
    Calculate probability of default and its cumulative.

            Parameters:
                data: is an array of dictionary. The dictionary must have the
                following keys - tenor, spread, rec_rate.

            Returns:
                an array of dictionary with the following keys - tenor, spread
                prob_default and cum_prob_default
    """
    datas = list(data)
    for datum in data:
        datum["prob_default"] = 1 - math.exp((-datum["spread"]) /
                                             ((1 - datum["rec_rate"]) * 10000))
        datum["cum_prob_default"] = 1 - math.exp((-datum["spread"] *
                                                  datum["tenor"]) /
                                                 ((1 - datum["rec_rate"]) *
                                                  10000))

    return datas


def prob_default_interpolation(defaults, tenors):
    """
    Interpolation of default probability

    Parameters
    ----------
    defaults : [{"tenor":float, "prob_default": float, "cum_prob_default": float}]
        The probability default curve
    tenors : [float]
        Tenor stated in the same manner as "tenor" key of "defaults" e.g time as fraction of year or number of days.

    Returns
    -------
    [float]
        Default probability for each tenor.

    """
    x_axis = [x["tenor"] for x in defaults]
    x_axis.insert(0, 0)

    def_y_axis = [x["prob_default"] for x in defaults]
    def_y_axis.insert(0, 0)

    cum_def_y_axis = [x["cum_prob_default"] for x in defaults]
    cum_def_y_axis.insert(0, 0)

    def_func = rcurve.interpolation(x_axis, def_y_axis, 1, is_function=True)
    cum_def_func = rcurve.interpolation(x_axis, cum_def_y_axis, 1, is_function=True)

    results = deque()
    for tenor in tenors:
        result = {}
        result["tenor"] = tenor
        result["prob_default"] = def_func(tenor)
        result["cum_prob_default"] = cum_def_func(tenor)
        results.append(result)
    return list(results)
