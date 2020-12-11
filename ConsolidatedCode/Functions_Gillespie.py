import Classes_DegradeFire4 as Classy
import numpy as np
import pandas as pd
import math
from pathlib import Path
from numpy import random

''' main Distributed Delay Stochastic Simulation Algorithm 
 http://localhost:8888/edit/BioResearch/python/Functions.py#   
 NOTE: There is no checking for negative values in this version.'''


def gillespie(reactions_list, stop_time, initial_state_vector):
    [state_vector, current_time, service_queue, time_series] = initialize(initial_state_vector)#changed
    subIterations = 220

#computationally, appending one row to a time series 100k times is expensive. It is much cheaper to append 100 rows to the time series 1k times
#As such, a subseries has
# the equation for how many computations are required =
#0.5*(data_rows/subIterations+1)*data_rows+(subIterations+2)*(subIterations+1)*data_rows/(2*data_rows) This is optimized for about 50000 rows, actual use varies but for the estimated values this results in only .7% time used
    time_series = pd.DataFrame([[current_time, state_vector]], columns=["time", "state"])
    iter = 0
    while current_time < stop_time:
        time_subseries = pd.DataFrame( columns=["time", "state"])
        if iter ==0:
            time_subseries =pd.DataFrame([[current_time, state_vector]], columns=["time", "state"])
        iter = 0
        while current_time < stop_time and iter<subIterations:

            cumulative_propensities = calculate_propensities(state_vector, reactions_list)
            next_event_time = draw_next_event_time(current_time, cumulative_propensities)
            if reaction_will_complete(service_queue, next_event_time):
                [state_vector, current_time] = trigger_next_reaction(service_queue, state_vector)
                time_subseries = update_time_series(time_subseries, current_time, state_vector)

                continue
            current_time = next_event_time
            next_reaction = choose_reaction(cumulative_propensities, reactions_list)
            processing_time = next_reaction.distribution()
            if processing_time == 0: # DMI
                state_vector = state_vector + next_reaction.change_vec
                time_subseries = update_time_series(time_subseries, current_time, state_vector)
            else:
                add_reaction(service_queue, current_time + processing_time, next_reaction)
            if(len(time_subseries)>0):
                iter = len(time_subseries)
        print(type(time_series))
        print(type(time_subseries))
        print(time_subseries)
        if len(time_series)<3:
            time_series = time_subseries
        else:
            time_series = time_series.append(time_subseries, ignore_index=True)# append the subseries to the series
    print("time_series")
    print(time_series)
    return dataframe_to_numpyarray(time_series)



def initialize(initial_state_vector):
    state_vector = initial_state_vector
    current_time = 0
    service_queue = []
    time_series = pd.DataFrame([[current_time, state_vector]], columns=['time', 'state'])
    return [state_vector, current_time, service_queue, time_series]


''' calculate_propensities creates an array with the cumulative sum of the propensity functions. '''


def calculate_propensities(x, reactions_list):
    propensities = np.zeros(np.shape(reactions_list))
    for index in range(np.size(reactions_list)):
        propensities[index] = reactions_list[index].propensity(x)
    return np.cumsum(propensities)


def reaction_will_complete(queue, next_event_time):
    if len(queue) > 0:
        if next_event_time > queue[0].comp_time:
            return True
    return False


def draw_next_event_time(current_time, cumulative_propensities):
    # print(cumulative_propensities[0])# debug
    return current_time + np.random.exponential(scale=(1 / cumulative_propensities[-1]))


''' choose_reaction rolls a biased die to determine which reaction will take place or be scheduled next. '''


def choose_reaction(cumulative_propensities, reactions_list):
    u = np.float32(np.random.uniform())
    next_reaction_index = min(
        np.where(np.float32(cumulative_propensities) > np.float32(cumulative_propensities[-1]) * u)[0])
    return reactions_list[next_reaction_index]


''' add_reaction, while not a pure function, does what it is supposed to,
    inserts into the queue a new delayed reaction sorted by completion time. '''


def add_reaction(queue, schedule_time, next_reaction):
    reaction = Classy.ScheduleChange(schedule_time, next_reaction.change_vec)
    if len(queue) == 0:
        return queue.append(reaction)
    else:
        for k in range(len(queue)):
            if reaction.comp_time < queue[k].comp_time:
                return queue.insert(k, reaction)
    return queue.append(reaction)


''' trigger_next_reaction has the side effect of removing the first entry of the queue it was passed. '''


def trigger_next_reaction(queue, state_vector):
    next_reaction = queue.pop(0)
    state_vector = state_vector + next_reaction.change_vec
    current_time = next_reaction.comp_time
    return [state_vector, current_time]# made data frame


def update_time_series(time_subseries, current_time, state_vector):
    '''print("time fresh")
    print(time_subseries)
    print("time Append")
    print(time_subseries.append(pd.DataFrame([[current_time, state_vector]])))'''
    return time_subseries.append (pd.DataFrame([[current_time, state_vector]],columns = ["time",  "state"]))
    # DMI columns=['time', 'state']) , ignore_index=True


# we are not sure if it is memory efficient to hae this function or if it is better to reove the content of
# this function and move to where the function is called

''' dataframe_to_numpyarray allows us to use the more efficient DataFrame class to record time series
    and then convert that object back into a usable numpy array. '''


def dataframe_to_numpyarray(framed_data):
    timestamps = np.array(framed_data[['time']])
    states = framed_data[['state']]
    arrayed_data = np.zeros([max(np.shape(timestamps)), np.shape(states.iloc[0, 0])[0] + 1])
    arrayed_data[:, 0] = timestamps.transpose()
    for index in range(max(np.shape(timestamps))):
        arrayed_data[index, 1:] = states.iloc[index, 0]
    return arrayed_data


def gillespie_sim(mu, cv, alpha, beta, R0, C0, yr, param, par, dilution, enzymatic_degradation):
    # model parameters
    init_Protein = (alpha - yr) * (mu - C0 * (math.sqrt(alpha / yr) - 1) / yr) * (
                .1 / (.1 + cv))  # calculate the avg peak to initialize at a peak

    dilution = Reaction(np.array([-1], dtype=int), 0, 0, [0, beta, 1, 0], 1, [0])
    enzymatic_degradation = Reaction(np.array([-1], dtype=int), 0, 0, [0, yr, R0, 1], 1, [0])
    production = Classy.Reaction(np.array([1], dtype=int), 0, 1, [alpha, C0, 2], 0, [mu, mu * cv])
    timeRun = 4
    # Naming files and paths
    path1 = 'PostProcessing/Simulations/{}{}'.format(param, par)
    file_name = '{}/mean={}_CV={}.csv'.format(path1, mu, cv)

    # Gillespie
    time_series = gillespie(np.array([production, enzymatic_degradation, dilution]), timeRun,
                            np.array([init_Protein], dtype=int))

    pd.DataFrame(time_series).to_csv(file_name, header=False, index=False)
    return file_name
