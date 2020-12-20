#!/usr/bin/env python

import Functions_DegradeFire4a as Gill
import Mod3Gillespie as Mod
import multiprocessing as mp

safeProcessors = 1# max(1, int(mp.cpu_count() * .5))
pool2 = mp.Pool(safeProcessors)
# Degrade and fire
from Functions_DegradeFire4a import *
from Classes_DegradeFire4 import *

import pandas as pd
import os
import math
import time
from pathlib import Path
import numpy as np

from numpy import random

from Functions_DegradeFire4a import gillespie_sim

# main

def Initialize_Classes(initialize_Gillespie):
    [alpha, beta, yr, r0, c0, mu, cv] = initialize_Gillespie  # remove runcount, remove indexing, get rid of extra shit in main.py

    dilution = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])

    enzymatic_degradation = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, r0, 1], 1, [0])
    print("degrade")
    production = Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, c0, 2], 0, [mu, mu * cv])

    reactionList = np.array([production, enzymatic_degradation, dilution])
    return reactionList

mean_range = np.linspace(5, 10, 16)
cv_range = np.linspace(0, .5, 16)
alpha = 300
beta = .1
R0 = 1
C0 = 10
yr = 80
par_range1 = np.linspace(150, 600, 2)  # alpha
par_range2 = np.linspace(.05, .2, 2)  # beta
param1 = 'alpha1'
param2 = 'beta'

# have list of parameters that vary and values, use double for loop to loop over these parameters Talk to SPC

path1 = "PostProcessing/Simulations/params_{}_{}/".format(param1, param2)
Path(path1).mkdir(parents=True, exist_ok=True)
# Add headers to the output files now before anything else is written to them.

'''path1 = 'PostProcessing/Simulations/'
temp_signalcsv = path1 + 'temp_signal.csv'
temp_peakscsv = path1 + "temp_peaks.csv"
statscsv = path1 + "stats.csv"


tempOutFileSignal = open(temp_signalcsv, 'a')
tempOutFilePeaks = open(temp_peakscsv, 'a')
OutFileStats = open(statscsv, 'a')

header_peaks = ["alpha" ,"beta, yr, r0, c0, mu, cv, initialState, run_number, peaks"]
header_stats = ["alpha",  "beta, yr, r0, c0, mu, cv, initialState, run_number, 
Mean Amplitude, CV Amp, Mean Period, CV Per"]

np.savetxt("temp_peaks.csv", header_peaks)
np.savetxt("stats.csv", header_stats) #, delimiter=","
tempOutFilePeaks.close()
OutFileStats.close()
'''

'''def gillespie_sim(mu, cv, alpha, beta, R0, C0, yr,param,par,dilution,enzymatic_degradation):
    #initialize(0)
    init_Protein = (alpha - yr) * (
                mu - C0 * (math.sqrt(alpha / yr) - 1) / yr)  # calculate the avg peak to initialize at a peak
    path1 = 'PostProcessing/Simulations/'+ param + str(par)
    #dilution = Reaction(np.array([-1, 0], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])
    #enzymatic_degradation = Reaction(np.array([-1, 0], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])
    production = Reaction(np.array([1, 0], dtype=int), 0, 1, [alpha, C0, 2], 0, [mu, mu * cv])
    time_series = gillespie(np.array([production, enzymatic_degradation, dilution]), 1,
                                np.array([init_Protein, 0], dtype=int))
    #file_name =  path1+ '/mean=' + str(mu) + '_CV=' + str(cv) + "?"+str(alpha) +'.csv'
    file_name =  path1+ '/mean=' + str(mu) + '_CV=' + str(cv) + '.csv'
    #badWords = open(file_name, 'w')  # make sure there is an empty file there
    #badWords.close()
    pd.DataFrame(time_series).to_csv(file_name, header=False, index=False)
    #time.sleep(140)
    pass'''
Initialize_Gillespie = []
iterations = 32
for par1 in par_range1:
    alpha = par1

    # path2 = 'Simulations/'+ param + str(par)

    for par2 in par_range2:
        beta = par2

        # path2 = 'Simulations/'+ param + str(par)

        pd.DataFrame([mean_range]).to_csv(path1 + '0metadata.csv', header=False, index=False)
        pd.DataFrame([cv_range]).to_csv(path1 + '0metadata.csv', mode='a', header=False, index=False)

        row_of_file_names = []

        for mu in mean_range:
            for cv in cv_range:
                for run in range(iterations):  # list all parameters here, they will be made into a list # [alpha,
                    # the following initial state is only for the degrade and fire model
                    initialState = int((alpha - yr) * (mu - C0 * (math.sqrt(alpha / yr) - 1) / yr) * (.1 / (.1 + cv)))
                    initialState = np.array([initialState], dtype=int)

                    dilution = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])
                    enzymatic_degradation = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])
                    production = Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, C0, 2], 0, [mu, mu * cv])
                    vector_Initialize_Gillespie = [[alpha, beta, yr, R0, C0, mu, cv, initialState]]
                    # order of entry vecor
                    Initialize_Gillespie = Initialize_Gillespie + vector_Initialize_Gillespie
                file_name = 'mean=' + str(mu) + '_CV=' + str(cv) + '.csv'
                row_of_file_names.append(file_name)
            paths = path1 + '1metadata.csv'
            pd.DataFrame([row_of_file_names]).to_csv(paths, mode='a', header=False, index=False)
            # pd.DataFrame([row_of_file_names]).to_csv('PostProcessing/Simulations/yr/1metadata.csv', mode='a',
            # header=False, index=False)
            row_of_file_names = []

        # for yr in yr_range:
        ''' 
        for mu in mean_range:
            run_pipeline(initialize_Gillespie, reactionList, stopTime, runCount, burnInTime, sampleRate, path1)

            # gillespie_sim(mu, cv, alpha, beta, R0, C0, yr,param,par,dilution,enzymatic_degradation)
            pool2.starmap(gillespie_sim, [(mu, cv, alpha, beta, R0, C0, yr, param1, par1, dilution,
                                           enzymatic_degradation) for cv in cv_range])
'''
try:
    time_chunk = 4
    stopTime = 40
    burnInTime = 4
    sampleRate = .1
    runCount = 1
    print(Initialize_Gillespie[2])
    runVector = [stopTime, runCount, burnInTime, sampleRate]

   # Mod.run_pipeline_splits(Initialize_Gillespie[2], runVector, path1)
    pool2.starmap(Mod.run_pipeline_single2,
                      [(Gillespie_vector, runVector, path1) for Gillespie_vector in Initialize_Gillespie])
finally:
    pool2.close()
    pool2.join()

