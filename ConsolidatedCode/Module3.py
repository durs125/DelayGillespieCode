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
        burntIn = Post.burn_in_time_series(tempSignal, burnInTime)
        uSampled = Post.uniformly_sample(burntIn, sampleRate)
        peaks = Post.detect_peaks(uSampled)
        np.savetxt("temp_peaks.csv",peaks, delimiter=",")

    tempOutFileSignal.close()
    tempOutFilePeaks.close()

    plt.plot(tempSignal[:,0],tempSignal[:,1])
    plt.show()

    plt.plot(uSampled[:,0],uSampled[:,1])
    plt.scatter(peaks[:,0],peaks[:,1])
    plt.show()

    return

