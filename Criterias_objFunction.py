#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:47:57 2022

@author: mikaelw
"""
import numpy as np
import matplotlib.pyplot as plt
from gr4j_model import gr4j



def rve(Qsim, Qobs):
    """Relative Volume Error function"""
    return sum(Qsim-Qobs)/sum(Qobs) * 100

def nse(Qsim, Qobs):
    """Nash–Sutcliffe model efficiency coefficient function"""
    return 1 - (sum((Qobs - Qsim)**2)) / (sum((Qobs-Qobs.mean())**2))

def kge(Qsim, Qobs):
    """Kling-Gupta efficiency function"""
    alpha = np.std(Qsim)/np.std(Qobs)
    beta = np.mean(Qsim) / np.mean(Qobs)
    r = np.corrcoef(Qsim,Qobs, rowvar = False)[0,1]  # rowvar = False because each column represents a variable, while the rows contain observations
    return 1 - np.sqrt((r-1)**2 + (alpha-1)**2 + (beta-1)**2)

def func_obj_calib(x0, Qobs_calib, inputX_calib, sign = -1): 
    """Compute the model with a specific set of parameters and use KGE function to evaluate the model"""
    X = [0]
    X.extend(x0)
    Qsim_calib = gr4j(X, inputX_calib)  
    return sign * kge(Qsim_calib, Qobs_calib)