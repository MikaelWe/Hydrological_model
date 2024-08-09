import numpy as np
from math import *
import matplotlib.pyplot as plt
#paramètres rentrés par l'utilisateur
tOn = 7
tRecouvrement = 0
tTransition = 3
stepTime = 0.1
nbRubanLed = 3
#Calcul temps séquence et discrétisation du tableau de la sequence
tSeq = nbRubanLed * (tOn  + tRecouvrement )
t = int(tSeq / stepTime)
seqRef = np.zeros(t)
#remplissage de la transition montee
nbStepTransition = int(tTransition / stepTime)
stepTransitionDebut = 0
stepDebutTOnMax = nbStepTransition + 1
for step in range(stepTransitionDebut,stepDebutTOnMax,1):
    seqRef[step] = (sin(- pi / 2 + step * pi / nbStepTransition) + 1) * 255 / 2
#remplissage durant pendant niveau 255
steptOnMax = int((tOn - 2 * tTransition) / stepTime) + stepDebutTOnMax
for step in range(stepDebutTOnMax,steptOnMax ,1):
    seqRef[step] = 255
#remplissage de la transition descente
stepFinSeq = steptOnMax + nbStepTransition
k = 0
for step in range(steptOnMax,stepFinSeq + 1,1):
    seqRef[step] = (sin(pi / 2 + k * pi / nbStepTransition) + 1) * 255 / 2
    k = k + 1
plt.plot(seqRef,"-*")
plt.show()



#seq = np.zeros((nbRubanLed,t))



