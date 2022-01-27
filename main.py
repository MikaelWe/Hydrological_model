import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gr4j_model import gr4j
import Criterias_objFunction as cof
from sklearn.model_selection import train_test_split
from scipy.optimize import minimize      
from monte_carlo import monte_carlo



surface_watershed = 526 * (10 ** 3) ** 2;  # (m2)
# Data from 1984 to 1998 for Viroin watershed
data = pd.read_excel("/Users/mikaelw/PycharmProjects/GR4J/dataset/Observed time series 1984-1998modified.xlsx", 
                     sheet_name="data", skiprows=1 ).to_numpy()
inputX = data[:, 0:2]  # Precipitation and evap
Y = data[:, 2].reshape((-1, 1)) * ((3600*24*1000)/surface_watershed)  # Observed discharge (mm/day)

inputX_calib, inputX_valid, Y_calib, Y_valid = train_test_split(inputX, Y, test_size=0.2, shuffle=False)  # Splitting the data for calibration and validation periods
Qobs_calib = Y_calib

#  Initial parameters of the model :
x0 = [320, 2.42, 60, 1.3]  # X1 X2 X3 X4

bnds = ((1, 1500), (-10, 5), (1, 500), (0.5, 4))  # bounds for the 4 parameters


# Calibration :
res_calib = minimize(cof.func_obj_calib, x0, (Qobs_calib, inputX_calib), bounds = bnds, method = 'nelder-mead',  # Use of an optimizer algorithm to find the parameters that maximize KGE criteria
                options = {'xtol': 1e-7, 'disp': False})

# Validation : Split-sample test
X = [0]
x0 = res_calib.x  
X.extend(x0)  # X with the optimized parameters from the calibration period
Qsim_valid = gr4j(X, inputX_valid)
    # Evaluation functions for validation period
nse = cof.nse(Qsim_valid, Y_valid)
rve = cof.rve(Qsim_valid, Y_valid)
kge = cof.kge(Qsim_valid, Y_valid)
print("Validation :")
print("kge", kge, "nse", nse, "rve", rve)

# Plots
plt.figure(1)
t_plot = np.arange(0, len(Qsim_valid), 1).reshape(len(Qsim_valid), 1)  # t pour le plot
plt.plot(t_plot, Y_valid, label="Qobs validation")
plt.plot(t_plot, Qsim_valid, linewidth = 0.3, label="Qsim validation")
plt.xlabel("t (day)")
plt.ylabel("Q (mm/day)")
plt.title("Discharge according to days")
plt.legend()
plt.show()

# Monte Carlo Simulation 
monte_carlo(X, inputX, Y)  # An example of a monte carlo simulation