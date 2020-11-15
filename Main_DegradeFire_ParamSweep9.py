#!/usr/bin/env python
#yr

import multiprocessing as mp
safeProcessors = max(1, mp.cpu_count() - 1)
pool2 = mp.Pool(safeProcessors)
# Degrade and fire
from  Functions_Gillespie import *
from  Classes_Gillespie import *
#Import after you start the multiprocess system
import pandas as pd
import os
import math
import time
from pathlib import Path
import numpy as np

mean_range =np.array([7.5, 0])
cv_range = np.linspace(0, .5, 16)
alpha = 300
beta = .1
R0 = 1
C0 = 10
yr =80
par_range = np.linspace(5, 20, 2)  # Co
param = 'Co' #THE parameter looped over

#This loops over THE parameter in question to first name the files according to which parameter value they have, an then loops over 
for par in par_range:
    
    C0 = par #This shoud be the parameter
    path1 = 'PostProcessing/Simulations/{}{}'.format(param,par)
    Path(path1).mkdir(parents=True, exist_ok=True)
    
#This generates and organizes meta data
    pd.DataFrame([mean_range]).to_csv(path1 + '/0metadata.csv', header=False, index=False)
    pd.DataFrame([cv_range]).to_csv(path1 + '/0metadata.csv', mode='a', header=False, index=False)
    row_of_file_names = []
    for mu in mean_range:
        if mu >0:
            for cv in cv_range:
                file_name = path1+ '/mean=' + str(mu) + '_CV=' + str(cv)  + '.csv'
                row_of_file_names.append(file_name)
            pathFile = path1 + '/1metadata.csv'
            pd.DataFrame([row_of_file_names]).to_csv(pathFile,  mode='a', header=False, index=False)
            row_of_file_names = [] # This includes the directory containing the files 
   #This ends the metadata creation
    dilution = Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])#np array used to allow expandability to multi species
    enzymatic_degradation = Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])
    #Loop over delay parameters only, mean and coeficient of varience (CV)
    for mu2 in mean_range:
        if mu2 >0:
            pool2.starmap(gillespie_sim, [(mu2, cv2,alpha,beta,R0 ,C0,yr,param,par,dilution,enzymatic_degradation) for cv2 in cv_range])

pool2.close()
pool2.join()
