import Functions_DegradeFire4a as Gill
import PostProcessing_Functions as Post
import numpy as np
import Classes_DegradeFire4 as Classy

'''if 1 ==0:
    mean_range =np.array([7.5,0]) # np.linspace(5, 10, 16)
    cv_range = np.linspace(0, .5, 16)
    alpha = 300
    beta = .1
    R0 = 1
    C0 = 10
    yr =80
    mean_range =np.array([7.5,0]) # np.linspace(5, 10, 16)
    cv_range = np.linspace(0, .5, 16)
    reactionList = (mu2, cv2,alpha,beta,R0 ,C0,yr,param,par,dilution,enzymatic_degradation)
    initialState = (alpha - yr) * ( mu - C0 * (math.sqrt(alpha / yr) - 1) / yr)
    stopTime = 40'''

import matplotlib.pyplot as plt


# add ,filterParamiters)


def run_pipeline_splits(initialize_Gillespie,  runVector, path1):
    # Prepare the reactions classes.
    [alpha, beta, yr, r0, c0, mu, cv, initialState, run_count ] = initialize_Gillespie [0:9]
    initialState = int (initialState)
    [length_run, split_Size, burnInTime, sampleRate] = runVector
    [stopTime, runCount, burnInTime, sampleRate] = runVector
    # print(run_count)
    dilution = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])

    enzymatic_degradation = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, r0, 1], 1, [0])
    print("degrade")
    production = Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, c0, 2], 0, [mu, mu * cv])

    reactionList = np.array([production, enzymatic_degradation, dilution])
    temp_signalcsv = path1 + 'temp_signal.csv'
    temp_peakscsv = path1 + "temp_peaks.csv"
    statscsv = path1 + "stats.csv"

    temp_signalcsv = path1 + '{}temp_signal.csv'.format( initialize_Gillespie )
    temp_peakscsv = path1 + "{}temp_peaks.csv".format( initialize_Gillespie )
    statscsv = path1 + "{}stats.csv".format( initialize_Gillespie )

    tempOutFileSignal = open(temp_signalcsv, 'a')
    tempOutFilePeaks = open(temp_peakscsv, 'a')
    OutFileStats = open(statscsv, 'a')

    header = "alpha, beta, yr, r0, c0, mu, cv, initialState, peaks"
    np.savetxt("temp_peaks.csv", header, delimiter=",")

    tempSignal = Gill.gillespie(reactionList, split_Size, initialState)
    np.savetxt("temp_signal.csv", tempSignal, delimiter=",")
    burned_in_signal = Post.slice_up_time_series(tempSignal, burnInTime,
                                                 length_run)  # do this as list and then iterate over list
    for x in burned_in_signal:
        uniform_sampling = Post.uniformly_sample(burned_in_signal, sampleRate)  # this needs to be passed by reference.
        peaks = Post.detect_peaks(uniform_sampling)
        # Fix1, need to make it such that the things are appended to a something taht can be saved to file

    # to_peak_file = initialize_Gillespie.extend(peaks)  # DMI  FIX!, issue: the peaks are in the form of time and height. Not good for 1 column format
    # recomended fix in format of making the peak data have a different set of deliminators than the other data
    # other fix, just put the parameters near the peak amp and time, along with a slice number
    np.savetxt("temp_peaks.csv", to_peak_file, delimiter=",")

    to_stat_file = Post.run_statistics(peaks)  # send the stats here
    np.savetxt("stats.csv", to_stat_file, delimiter=",")

    tempOutFileSignal.close()
    tempOutFilePeaks.close()
    OutFileStats.close()

    plt.plot(tempSignal[:, 0], tempSignal[:, 1])
    plt.show()

    return


def run_pipeline_single(initialize_Gillespie,  runVector, path1):
    print(initialize_Gillespie)
    ''' Make 87-95 it's own function. To be written in main.py 
        -SPC
        ''' 
    [alpha, beta, yr, r0, c0, mu, cv, initialState, run_count] = initialize_Gillespie [0:9]  # remove runcount, remove indexing, get rid of extra shit in main.py
    
    dilution = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])

    enzymatic_degradation = Classy.Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, r0, 1], 1, [0])
    print("degrade")
    production = Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, c0, 2], 0, [mu, mu * cv])

    reactionList = np.array([production, enzymatic_degradation, dilution])

    initialState = np.array([initialState], dtype=int)  #Do we need to recreate initialState? Probably put this in main.py
    [stopTime, runCount, burnInTime, sampleRate] = runVector # remove runcount
    
    ''' Rename temp_* to just * and make filenames "human readable" in lines 102,103
        -SPC'''
    temp_signalcsv = path1 + '{}temp_signal.csv'.format( initialize_Gillespie )
    temp_peakscsv = path1 + "{}temp_peaks.csv".format( initialize_Gillespie )
    ''' move line 105 to main.py and put all stats in a single file. -SPC '''
    statscsv = path1 + "{}stats.csv".format( initialize_Gillespie )

    tempOutFileSignal = open(temp_signalcsv, 'a')  # refactor temp and open with 'w' instead of 'a'
    tempOutFilePeaks = open(temp_peakscsv, 'a')  # refactor temp and open with 'w' instead of 'a'
    OutFileStats = open(statscsv, 'a')  # VERY IMPORTANT TO KEEP AS 'a'

    print("Gillespie start")

    tempSignal = Gill.gillespie(reactionList, stopTime, initialState)  # refactor tempSignal to Signal
    np.savetxt(temp_signalcsv, initialize_Gillespie, delimiter=",")  # maybe make this a header
    np.savetxt(temp_signalcsv, tempSignal, delimiter=",")
    print(tempSignal)
    burned_in_signal = Post.burn_in_time_series(tempSignal, burnInTime)  # this needs to be passed by reference.
    uniform_sampling = Post.uniformly_sample(burned_in_signal, sampleRate)  # this needs to be passed by reference.
    ''' Break up time series here '''
    peaks = Post.detect_peaks(uniform_sampling)
    to_peak_file = initialize_Gillespie.extend(peaks)  # DMI  FIX!, issue: the peaks are in the form of time and height.

    # recommended fix in format of making the peak data have a different set of deliminators than the other data
    # other fix, just put the parameters near the peak amp and time, along with a slice number
    #np.savetxt("temp_peaks.csv", to_peak_file, delimiter=",")  # might want to try json, might work
    # know np.save(filename.npy,list) can save list to file. Cannot append multiple files
    to_stat_file = initialize_Gillespie
    to_stat_file[len(to_stat_file):] = Post.run_statistics(peaks)
    to_stat_file = Post.run_statistics(peaks)  # send the stats here
    print(peaks)
    np.savetxt(statscsv, to_stat_file, delimiter=",")

    tempOutFileSignal.close()
    tempOutFilePeaks.close()
    OutFileStats.close()
# plotting
    ''' """
    plt.plot(tempSignal[:, 0], tempSignal[:, 1])
    plt.show()
    '''
    return

