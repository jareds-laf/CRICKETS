from analysis_functions import get_tavg_kurtosis
from plotting_functions import plot_tavg_kurtosis
import matplotlib.pyplot as plt
# import pylab as plt
import blimpy
from blimpy import Waterfall
from blimpy import calcload
import os
import glob
import numpy as np
import time
from scipy.stats import norm, kurtosis
import scipy
import numpy.ma as ma
import pandas as pd

t0=time.time()
print(f'Start time: {time.strftime("%H:%M:%S", time.localtime())}')

# water0 = Waterfall(os.path.normpath(file_list[0]), max_load = ml_list[0])
water1 = Waterfall(os.path.normpath(file_list[1]), max_load = ml_list[1])
# water2 = Waterfall(os.path.normpath(file_list[2]), max_load = ml_list[2])
# water3 = Waterfall(os.path.normpath(file_list[3]), max_load = ml_list[3])
# water4 = Waterfall(os.path.normpath(file_list[4]), max_load = ml_list[4])
# water5 = Waterfall(os.path.normpath(file_list[5]), max_load = ml_list[5])
# water6 = Waterfall(os.path.normpath(file_list[6]), max_load = ml_list[6])

# water_list = [water0, water1, water2, water3, water4, water5, water6]
# water_list = [0, water1, water2, water3, water4, 0, 0]
# print(water_list)

t1 = time.time()
print(f'Elapsed time: {t1 - t0}')

get_tavg_kurtosis(water1)
plot_tavg_kurtosis(water1)
