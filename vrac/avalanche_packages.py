import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from operator import mul
from functools import reduce
from scipy.stats import pearsonr
import pickle, random, itertools, sys
from math import log
from reliability.Fitters import Fit_Weibull_2P, Fit_Weibull_3P


