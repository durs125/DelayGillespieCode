# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import Module3 as M3
import numpy as np

alpha = 300
beta = .1
R0 = 1
C0 = 10
yr = 80
mu = 5
cv = .5
production = M3.Gill.Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, C0, 2], 0, [mu, mu * cv])
dilution = M3.Gill.Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1,[0])  # np array used to allow expandability to multi species
enzymatic_degradation = M3.Gill.Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])


tempReactionList = np.array([production, enzymatic_degradation, dilution])
tempStopTime = 400
tempInitialState = np.array([160], dtype=int)
tempRunCount = 1
tempBurnInTime = 200
tempSampleRate = 10

M3.run_pipeline(tempReactionList, tempStopTime, tempInitialState, tempRunCount, tempBurnInTime, tempSampleRate)





print ("ALL DONE")
