import numpy as np
from math import tanh

def gr4j(X, data):
    """GR4J hydrological model (rainfall-runoff model)
    
Perrin, C., C. Michel, et al. (2003). "Improvement of a parsimonious model for streamflow simulation." Journal of Hydrology 279(1-4): 275-289 
https://webgr.inrae.fr/modeles/journalier-gr4j-2/ """
    # input data
    P = data[:, 0].reshape(-1, 1)  # Precipitation (mm/day)  
    E = data[:, 1].reshape(-1, 1)  # Potentcial evapotranspiration (mm/day)
    
    start_date_Qsim_fullperiod = 20  # Starting at day x allows the construction of the unit hydrographs    
    
    # Initial conditions
    S = np.zeros((len(P), 1))  # Production store level (mm)
    R = np.zeros((len(P), 1))  # Rounting store level (mm)
    S[0] = 0.6 * X[1]  # S is calculated also during the first days of the simulation (allow ti-o construct the unit hydrographs)
    R[start_date_Qsim_fullperiod] = 0.7 * X[3] 
    
    # Initializations variables
    l = int(X[4]) + 1; m = int(2*X[4]) + 1
    UH1 = np.zeros((l+1, 1)); UH2 = np.zeros((m+1, 1))
    Pn = np.zeros((len(P), 1)); En = np.zeros((len(P), 1)); Ps = np.zeros((len(P), 1)); Es = np.zeros((len(P), 1)); Perc = np.zeros((len(P), 1))
    Pr = np.zeros((len(P), 1)); Q1 = np.zeros((len(P), 1)); Q9 = np.zeros((len(P), 1)); Qsim_fullperiod = np.zeros((len(P), 1));

    # Construction of the unit hydrographs
    for i in range(1, l + 1):
        UH1[i] = scurve_1(i, X, l) - scurve_1(i - 1, X, l)  # ordinate of HU1 
    for i in range(1, m + 1):
        UH2[i] = scurve_2(i, X, m) - scurve_2(i - 1, X, m)  # ordinate of HU2
        
    # Running the model
    for t in range(1, len(P)):  
        if P[t] > E[t]:
            Pn[t] = P[t] - E[t]
            En[t] = 0
            Ps[t] = (X[1] * (1-(S[t-1]/X[1])**2)*tanh(Pn[t]/X[1]))/(1+S[t-1]/X[1] * tanh(Pn[t]/X[1]))
        else:
            Pn[t] = 0
            En[t] = E[t] - P[t]
            Es[t] = (S[t-1]*(2-S[t-1]/X[1])*tanh(En[t]/X[1]))/(1+(1-S[t-1]/X[1])*tanh(En[t]/X[1]))
        S[t] = S[t-1] - Es[t] + Ps[t]
        Perc[t] = S[t] * (1 - (1 + (4/9 * S[t]/X[1])**4)**(-1/4))
        S[t] = S[t] - Perc[t]
        Pr[t] = Perc[t] + Pn[t] - Ps[t]

        if t > start_date_Qsim_fullperiod:  
            sumq9 = 0  
            for k in range(1, l+1): 
                sumq9 = sumq9 + UH1[k] * Pr[t-k+1]
            Q9[t] = 0.9 * sumq9
            sumq1 = 0
            for k in range(1, m+1): 
                sumq1 = sumq1 + UH2[k] * Pr[t - k + 1]
            Q1[t] = 0.1 * sumq1

            F = X[2] * (R[t-1]/X[3])**(7/2)  # groundwater exchange term
            R[t] = max(0, R[t-1] + Q9[t] + F)
            Qr = R[t] * ( 1 - (1+(R[t]/X[3])**4)**(-1/4))
            R[t] = R[t] - Qr
            Qd = max(0, Q1[t] + F)
            Qsim_fullperiod[t] = Qr + Qd
    return Qsim_fullperiod

def scurve_1(t, X, l):
    SH1 = np.zeros((l+1, 1))
    if t < X[4]:
        SH1[t] = (t/X[4])**(5/2)
    else:
        SH1[t] = 1
    return SH1[t]

def scurve_2(t, X, m):
    SH2 = np.zeros((m+1, 1))
    if t < X[4]:
        SH2[t] = (t/X[4])**(5/2)
    elif X[4]<=t and t<(2*X[4]):
        SH2[t] = 1 - 0.5 * (2 - t/X[4])**(5/2)
    else:
        SH2[t] = 1
    return SH2[t]









