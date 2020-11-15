from numpy import random
import numpy as np
#from numba import jit, jitclass, types, typed

'''spec = [
    ('value', int32),               # a simple scalar field
    ('array', float32[:]),          # an array field
]

@jitclass(spec)
class Bag(object):
    def __init__(self, value):
        self.value = value
        self.array = np.zeros(value, dtype=np.float32)'''

# DO NOT EDIT THIS SECTION '''

'''spec = [
    ('l', types.ListType(types.pyobject)),             
    ('l', types.ListType(types.pyobject)),  
    ('value', types.int32),         
]'''

class Reaction:
    propensities_list = ['mobius_propensity', 'decreasing_hill_propensity',
                         'increasing_hill_propensity', 'mobius_sum_propensity',
                         'dual_feedback_decreasing_hill_propensity', 'dual_feedback_increasing_hill_propensity']
    distributions_list = ['gamma_distribution', 'trivial_distribution', 'bernoulli_distribution']
    system_size = 1




#!jit error
    def __init__(self, state_change_vector, parts_of_state_vector,
                 propensity_id, propensity_params,
                 distribution_id, distribution_params):
        self.change_vec = state_change_vector
        self.parts_of_vec = parts_of_state_vector
        if type(propensity_id) == str:
            self.prop_id = Reaction.propensities_list.index(propensity_id)
        else:
            self.prop_id = propensity_id
        self.prop_par = propensity_params
        if type(distribution_id) == str:
            self.dist_id = Reaction.distributions_list.index(distribution_id)
        else:
            self.dist_id = distribution_id
        self.dist_par = distribution_params

    ''' DO NOT EDIT THIS SECTION '''
#!jit error
    def propensity(self, x):
        return getattr(self, Reaction.propensities_list[self.prop_id])(x * Reaction.system_size) / Reaction.system_size
#!jit error
    def distribution(self):
        return getattr(self, Reaction.distributions_list[self.dist_id])()

    ''' Propensities start here.
    After adding a propensity function to the list of definitions, 
    append the name to props_list. '''
    #!jit error
    def mobius_propensity(self, state_vector):
        """ For a constant function f(x) = c, assign the vector [c,0,1,0]
            For a linear map f(x) = a * x,    assign the vector [0,a,1,0] """
        x = state_vector[self.parts_of_vec]
        return (self.prop_par[0] + self.prop_par[1] * x) / (self.prop_par[2] + self.prop_par[3] * x)
    #!jit error
    def decreasing_hill_propensity(self, state_vector):
        x = state_vector[self.parts_of_vec]
        scale = self.prop_par[0]
        threshold = self.prop_par[1]
        exponent = self.prop_par[2]
        return scale * (threshold / (x + threshold)) ** exponent
    #!jit error
    def increasing_hill_propensity(self, state_vector):
        x = state_vector[self.parts_of_vec]
        scale = self.prop_par[0]
        threshold = self.prop_par[1]
        exponent = self.prop_par[2]
        return scale * (x / (x + threshold)) ** exponent
    #@jit
    def mobius_sum_propensity(self, state_vector):
        total_species = np.sum(state_vector)
        x = state_vector[self.parts_of_vec]
        return (self.prop_par[0] + self.prop_par[1] * x) / (self.prop_par[2] + self.prop_par[3] * total_species)
    #@jit
    def dual_feedback_decreasing_hill_propensity(self, state_vector):
        scale = self.prop_par[0]
        threshold0 = self.prop_par[1]
        threshold1 = self.prop_par[2]
        factor = self.prop_par[3]
        exponent = self.prop_par[4]
        return scale * (threshold0 / (threshold0 + state_vector[self.parts_of_vec])) ** exponent * \
               (threshold1 / factor + state_vector[(self.parts_of_vec + 1) % 2]) / \
               (threshold1 + state_vector[(self.parts_of_vec + 1) % 2])
    #@jit
    def dual_feedback_increasing_hill_propensity(self, state_vector):
        scale = self.prop_par[0]
        threshold0 = self.prop_par[1]
        threshold1 = self.prop_par[2]
        factor = self.prop_par[3]
        exponent = self.prop_par[4]
        return scale * (threshold0 / (threshold0 + state_vector[(self.parts_of_vec + 1) % 2])) ** exponent * \
               (threshold1 / factor + state_vector[self.parts_of_vec]) / \
               (threshold1 + state_vector[self.parts_of_vec])
    #@jit
    def heavyside_propensity(self, state_vector):
        scale = self.prop_par[0]
        threshold = self.prop_par[1]
        

    ''' Distributions start here.
    After adding a distribution to the list of definitions,
    append the name to distr_list. '''
   #!jit error
    def gamma_distribution(self):
        if self.dist_par[1] == 0:
            return self.trivial_distribution()
        mean = self.dist_par[0]
        stdev = self.dist_par[1]
        return random.gamma(shape=(mean / stdev) ** 2, scale=stdev ** 2 / mean)
    #!jit error
    def trivial_distribution(self):
        mean = self.dist_par[0]
        return mean
    #@jit #errored
    def bernoulli_distribution(self):
        mean = self.dist_par[0]
        stdev = self.dist_par[1]
        return mean - 1 + 2 * stdev * random.randint(2)


class ScheduleChange:
    #jit error
    def __init__(self, completion_time, change_vector):
        self.comp_time = completion_time
        self.change_vec = change_vector
