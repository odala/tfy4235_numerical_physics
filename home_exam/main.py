# Diffusion
# from myIterativeMethods import explicit_euler, implicit_euler, ...

import sys  # for error massages
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time
import scipy.linalg
from scipy.special import erf

# --- Import my functions.
#from tdma import tdma