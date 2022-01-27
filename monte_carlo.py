#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 10:16:08 2022

@author: mikaelw
"""
import numpy as np
import matplotlib.pyplot as plt
from gr4j_model import gr4j

def monte_carlo(X, inputX, Y):
    """A Monte Carlo simulation for the 2 most senstive parameters of the GR4J model"""
    N_sim = 500
    plt.figure(2)
    plt.show()
    for n in range(N_sim):
        X[2] = np.random.uniform(-10,5)  # Assuming a uniform law for errors in parameters
        X[3] = np.random.uniform(1,500)
        Qsim_MC = gr4j(X, inputX)  # Running the model
        t_plot_MC = np.arange(0, len(Qsim_MC), 1).reshape(len(Qsim_MC), 1)
        plt.plot(t_plot_MC, Qsim_MC, ".", markersize = 0.2)
        plt.xlabel("t (day)")
        plt.ylabel("Q (mm/day)")
    t = np.arange(0, len(Qsim_MC), 1).reshape(len(Qsim_MC), 1)
    plt.plot(t, Y, 'r', linewidth = 0.3) # Observed discharge
    plt.title("Monte Carlo simulations")

