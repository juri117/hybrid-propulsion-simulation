'''
Created on 15.12.2015

@author: Juri Bieler
'''

import numpy as np
import csv
import math
import os
import itertools
import operator

'''
rechnet um: N [1/min] zu rad/s
'''
def N2omega(N):
    omega = 2 * np.pi * N / 60
    return omega