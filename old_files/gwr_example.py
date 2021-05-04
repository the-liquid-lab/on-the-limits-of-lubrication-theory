# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:14:29 2021

@author: C6139404
"""
from functools import lru_cache

from mpmath import mp

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 16

from gwr_inversion import gwr

def fun(s: float):
    # unit step transform
    return 1 / s * mp.exp(-((s**2+mp.fraction(37, 100)*s+1)/(s**2 + s + mp.pi))**mp.fraction(1,2))

time = 10 ** np.linspace(-1, 2, 101)
a = 1.0

M = 32
inv1 = gwr(fun, time, M)
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.plot(time, inv1, 'o', c='C3', ms=10, mfc='w', label=r'Unit Step GWR')
ax.set(xscale='log', xlim=(1e-1, 1e2), ylim=(-1, 2))

