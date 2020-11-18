import Functions_Gillespie as Gill
import PostProcessing_Functions as Post
import numpy as np
import matplotlib.pyplot as plt

# add ,filterParamiters)
def run_pipeline(reactionList, stopTime, initialState, runCount, burnInTime, sampleRate):
    tempOutFileSignal = open('temp_signal.csv', 'a')
    tempOutFilePeaks = open("temp_peaks.csv", 'a')

    for i in range(runCount):
        tempSignal = Gill.gillespie(reactionList, stopTime, initialState)
        np.savetxt("temp_signal.csv",tempSignal,delimiter=",")
        Post.burn_in_time_series(tempSignal, burnInTime)
        Post.uniformly_sample(tempSignal, sampleRate)
        peaks = Post.detect_peaks(tempSignal)
        np.savetxt("temp_peaks.csv",peaks, delimiter=",")

    tempOutFileSignal.close()
    tempOutFilePeaks.close()

    plt.plot(tempSignal[:,0],tempSignal[:,1])
    plt.show()


    return

# init_Protein = (alpha - yr) * ( mu - C0 * (math.sqrt(alpha / yr) - 1) / yr)   # calculate the avg peak to initialize at a peak
# time_series = gillespie(np.array([production, enzymatic_degradation, dilution]), timeRun,np.array([init_Protein], dtype=int))
