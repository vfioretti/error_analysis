"""
 cstat_simula.py  -  PyXSPEC simulation
 ---------------------------------------------------------------------------------
 Author: V. Fioretti (INAF/OAS) valentina.fioretti@inaf.it
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
 --------------------------------------------------------------------------------
 Usage example:
 > python cstat_simula.py base10 
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



# Import the input parameters
arg_list = sys.argv
model = arg_list[1]

# simulation
Xset.restore(model+".xcm")
Xset.restore("simula.xcm")

Xset.save(model+"_error.xcm",info="a") 
