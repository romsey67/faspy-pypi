#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 14:11:55 2020

@author: RMS671214
"""
from sympy import Symbol, diff, solve


def solver_bond_derivatives(fv, cfs=[], yld_dcfs=[], fvs=[], my_yield=5):
    yld = Symbol('yld')
    df = Symbol('df')
    expr = Symbol('expr')
    expr = None
    dcfs_len = len(yld_dcfs)
    for i in range(dcfs_len):
        ci = None
        df = None
        if yld_dcfs[i] != 0:
            ci = (cfs[i] + fvs[i]) / fv * 100

            if df is None:
                df = (1 + yld * yld_dcfs[i] * 0.01)

            for k in range(i):
                if yld_dcfs[k] != 0:
                    df = df * (1 + yld * yld_dcfs[k] * 0.01)

            if expr is None:
                expr = ci/df
            else:
                expr = expr + ci/df

    expr_duration = diff(expr, yld)
    duration = expr_duration.evalf(subs={yld: my_yield})
    # print(type(duration))
    # print('Duration >>>>',duration)
    expr_convex = diff(expr_duration, yld)
    convexity = expr_convex.evalf(subs={yld: my_yield})
    # print('Convexity >>>>>>',convexity)
    return {'duration': float(duration), 'convexity': float(convexity)}


def bond_risk(structure, my_yield=5):
    yld = Symbol('yld')
    df = Symbol('df')
    expr = Symbol('expr')
    expr = None
    dcfs_len = len(structure)
    fv = structure[-1]['fv']
    #print('structure ==>>', structure)
    for i in range(1, dcfs_len, 1):
        row = structure[i]
        #print(row)
        ci = None
        df = None
        #print('i======>>>>', i, row)
        if row['time'] > 0:
            ci = (row['cf'] + row['fv_flow']) / fv * 100

            if df is None:
                df = (1 + yld * row['y_dcf'] * 0.01)

            for k in range(1, i, 1):
                row2 = structure[k]
                if row2['y_dcf'] != 0:
                    df = df * (1 + yld * row2['y_dcf'] * 0.01)

            if expr is None:
                expr = ci/df
            else:
                expr = expr + ci/df

    expr_duration = diff(expr, yld)
    duration = expr_duration.evalf(subs={yld: my_yield})
    # print(type(duration))
    # print('Duration >>>>',duration)
    expr_convex = diff(expr_duration, yld)
    convexity = expr_convex.evalf(subs={yld: my_yield})
    # print('Convexity >>>>>>',convexity)
    return {'duration': float(duration), 'convexity': float(convexity)}
