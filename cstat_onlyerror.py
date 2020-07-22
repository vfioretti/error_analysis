"""
 cstat_error.py  -  pyXSPEC cstat error evaluation
 ---------------------------------------------------------------------------------
 Author: Mehdy Lefkir
 Modified by: V. Fioretti (INAF/OAS) valentina.fioretti@inaf.it
 ---------------------------------------------------------------------------------
 Dependencies:
 - python 2.7
 - numpy
 - scipy
 - matplotlib
 - PyPDF2
 - pyXSPEC running on Python 2.7
  ---------------------------------------------------------------------------------
 Parameters:
 - model = name of the model (the .xcm file must be <model>.xcm)
 - input_level = chi2 level to evaluate a confidence interval
 - statistic = 'cstat' or 'chi' (Statistic of the fit method)
 - blacklist = list of strings. Contains the list of the parameters's name to be frozen before the error computation.
 For instance if blacklist=['Sigma'] all parameters named 'Sigma' will be frozen and the error on those parameters won't be computed. Default is ['']
 - n_cores = float. Number of cores to set for the XSPEC parallel variable.
 - level : float. Chi2 level to evaluate a confidence interval. Default is chi2 = 2.706 for a 90% confidence interval.
 +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
 |   sigma   |  1.00s  |  1.28s  |  1.64s   |  1.96s  |  2.00s  |  2.58s  |  3.00s  |  3.29s  |  4.00s  |
 +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
 | conf_int  | 68.27%  | 80.00%  | 90.00%   | 95.00%  | 95.45%  | 99.00%  | 99.73%  | 99.90%  | 99.99%  |
 | level     | 1.000   | 1.642   | 2.706    | 3.841   | 4.000   | 6.635   | 9.000   | 10.828  | 16.000  |
 +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
 
 - plot_statistic = bool. Used to plot or not the statistic graphs.
 - interp_method : str. Method for the interpolation of the statistic. "linear" uses a linear interpolation and the brentq method
 for a the root finding. "spline" uses a spline interpolation and find the roots with a scipy method in the spline class.
 - selection = 'all' or list of numbers separated with space (e.g. "1 2 3"). Select the parameters used to get the errors with a list of parameter number.
 - [optional] = if selection is 'some',
 

 --------------------------------------------------------------------------------
 Usage example:
 > python cstat_error.py base10 2.706 
"""

import matplotlib.pyplot as plt
import numpy as np
from xspec import *
import sys, os
from scipy import interpolate  
import matplotlib.backends.backend_pdf
from matplotlib import rc
import shutil
import datetime


import pyXIFU as px

# Import the input parameters
arg_list = sys.argv
model = arg_list[1]
statistic = arg_list[2]
selection = arg_list[3]
blacklist = arg_list[4]
n_cores = np.int(arg_list[5])
input_level = np.float(arg_list[6])
plot_statistic = bool(arg_list[7])
interp_method = arg_list[8]


if selection == 'all':
    selection_input = selection

if selection == 'some':
    selection_input = list(arg_list[8].split(" "))


px.ml_get_errors(model+"_error",'cstat',selection = selection_input, blacklist = blacklist, n_cores=n_cores,level=input_level, plot_statistic = plot_statistic, interp_method = interp_method)
