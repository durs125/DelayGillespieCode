# correct the for loop it loops over the wrong opject. It has the wrong thing as a variable name

#!/usr/bin/env python
#         gamma_r40.0
import multiprocessing as mp
safeProcessors = max(1, int(mp.cpu_count()/2 ))
pool2 = mp.Pool(safeProcessors)
import PostProcessing_Functions8Short as Fun
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os
import imp
imp.reload(Fun)

#directory2 = ["/scratch/Infor/__pycache__/PostProcessing/Simulations/Co20.0/","/scratch/Infor/__pycache__/PostProcessing/Simulations/Co5.0/", "/scratch/Infor/__pycache__/PostProcessing/Simulations/Ro0.5/","/scratch/Infor/__pycache__/PostProcessing/Simulations/Ro2.0/"]
directory1 = ["/home/david/BioResearch/python/PostProcessing/Simulations/beta0.2/","/home/david/BioResearch/python/PostProcessing/Simulations/Ro2.0/"]
#directory2 = "/scratch/Infor/__pycache__/"
#freqs = np.linespace(1,3,10)
for directory in directory1:

    file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))
    heat_map_axes = [16]
    heat_map_axes0 = [1]
    heat_map_matrices = np.zeros([4, file_names.shape[0], file_names.shape[1]])


    directory2 = "scratch/Infor/"
    file_names = np.array(pd.read_csv(directory+'1metadata.csv', header=None))


    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(file_names[mean_axis,nameRange],  600,  10, 0) for nameRange in range(file_names.shape[1])])

    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])

    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Once_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Once_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Once_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Once_amplitude_CV.csv', mode = "a")

    for mean_axis in range(file_names.shape[0]):

        pool2.starmap(Fun.cleanStatsHeatMap,[(file_names,  600,  10, 0)  for file_names in range(file_names.shape[1])])

    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'UniPass_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'UniPass_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'UniPass_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'UniPass_amplitude_CV.csv', mode = "a")

    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(file_names , 600,10, 2) for file_names in range(file_names.shape[1])])


    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Triang_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Triang_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Triang_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Triang_amplitude_CV.csv', mode = "a")
    for mean_axis in range(file_names.shape[0]):
        pool2.starmap(Fun.cleanStatsHeatMap, [(file_names, 600,10, 3) for file_names in range(file_names.shape[1])])

    mean_period = pd.DataFrame(heat_map_matrices[0, :, :])
    mean_amplitude = pd.DataFrame(heat_map_matrices[1, :, :])
    period_CV = pd.DataFrame(heat_map_matrices[2, :, :])
    amplitude_CV = pd.DataFrame(heat_map_matrices[3, :, :])

    mean_period.to_csv(directory + 'Binom_mean_period.csv', mode = "a")
    mean_amplitude.to_csv(directory + 'Binom_mean_amplitude.csv', mode = "a")
    period_CV.to_csv(directory + 'Binom_period_CV.csv', mode = "a")
    amplitude_CV.to_csv(directory + 'Binom_amplitude_CV.csv', mode = "a")
pool2.close()
pool2.join()
