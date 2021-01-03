#!/usr/bin/env python
import numpy as np


def burn_in_time_series(signal, burn_in_time):
    temp_signal = signal
    temp_signal[:, 0] = signal[:, 0] - burn_in_time

    new_start_time = np.where(temp_signal[:, 0] > 0)[0][0]
    new_start_state = np.zeros(temp_signal.shape[1])
    new_start_state[1:] = temp_signal[new_start_time - 1, 1:]
    temp_signal[new_start_time - 1, :] = new_start_state
    burned_in_signal = temp_signal[(new_start_time - 1):, :]
    return burned_in_signal


def uniformly_sample(signal, rate=0, number_of_samples=0):
    n = min(np.shape(signal))
    end = signal[np.shape(signal)[0] - 1, 0]
    if rate > 0:
        number_of_samples = int(end * rate)
    elif rate < 0 or number_of_samples < 1:
        raise ValueError("No samples specified or no sampling rate specified")
    uniform_sampling = np.float32(np.zeros([number_of_samples, n]))
    uniform_timestamps = np.linspace(0, end, number_of_samples)
    uniform_sampling[:, 0] = uniform_timestamps
    counter = 0
    for index in range(number_of_samples):
        while counter < max(np.shape(signal)):
            if signal[counter, 0] > uniform_timestamps[index]:
                uniform_sampling[index, 1:n] = signal[counter - 1, 1:n]
                break
            counter += 1
    return uniform_sampling


def is_max_in_window(signal, length_of_signal, window_size):
    def window_checker(index):
        if index <= window_size or index >= length_of_signal - window_size:
            return np.float32(1 + np.random.uniform())
        elif signal[index, 1] < signal[index - window_size, 1] \
                or signal[index, 1] < signal[index + window_size, 1]:
            return np.float32(1 + np.random.uniform())
        else:
            return np.float32(0)

    return window_checker


def compute_optimal_time_window(signal):
    """ first half of the peak detection algorithm """
    n = max(np.shape(signal))
    rows = int(np.ceil(n / 2) - 1)
    lms = np.zeros((rows, n), dtype="float16")
    for x in range(0, rows):
        lms[x, :] = np.array(list(map(is_max_in_window(signal, n, x + 1), range(n))))
    row_sum = np.sum(lms, axis=1)
    gamma = np.where(row_sum == np.amin(row_sum))
    rescaled_lms = np.vsplit(lms, gamma[0] + 1)[0]
    return rescaled_lms


def detect_peaks(signal):
    column_sd = np.std(compute_optimal_time_window(signal), axis=0)
    peaks_index = np.where(column_sd == 0)
    peaks = signal[peaks_index, :]
    peaks = peaks[0, :, :]
    return peaks


def chop_peaks(signal, filename, chop_size=2000):
    with open(filename, 'w') as record_peaks:
        for index in range(int(max(np.shape(signal))/chop_size)):
            start = index * (chop_size + 1)
            peak_data = detect_peaks(signal[start:(start + chop_size), :])
            delta_time = np.diff(peak_data[:, 0])
            full_peak_data = np.zeros([np.shape(peak_data)[0] - 1, 3])
            full_peak_data[:, 0] = delta_time
            full_peak_data[:, 1:] = peak_data[1:, :]
            np.savetxt(record_peaks, full_peak_data, delimiter=',')
