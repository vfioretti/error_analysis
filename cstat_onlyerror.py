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

    +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
    |   sigma   |  1.00s  |  1.28s  |  1.64s   |  1.96s  |  2.00s  |  2.58s  |  3.00s  |  3.29s  |  4.00s  |
    +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
    | conf_int  | 68.27%  | 80.00%  | 90.00%   | 95.00%  | 95.45%  | 99.00%  | 99.73%  | 99.90%  | 99.99%  |
    | level     | 1.000   | 1.642   | 2.706    | 3.841   | 4.000   | 6.635   | 9.000   | 10.828  | 16.000  |
    +-----------+---------+---------+----------+---------+---------+---------+---------+---------+---------+
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
input_level = np.float(arg_list[2])

# simulation
#Xset.restore(model+".xcm")
#Xset.restore("simula.xcm")

#Xset.save(model+"_error.xcm",info="a")

px.ml_get_errors(model+"_error",'cstat',n_cores=8,level=input_level) 
