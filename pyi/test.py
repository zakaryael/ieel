#!/opt/homebrew/anaconda3/bin/python
import os
import numpy as np 

Ns = 200

def mkwdir(where='./output_data'):
    counter = 0
    while os.path.isdir(where + '/wdir' + str(counter)): #improvement: use os.path.join()
        counter += 1
    path = where + '/wdir' + str(counter)
    os.makedirs(where + '/wdir' + str(counter))
    return path
    

def swim(wdir, iteration, action, nswim=10, nout=10):
    os.system('../runs/bin/swim -dir ' + wdir + ' -it ' + str(iteration) + ' -a ' + str(action) + ' -nsw ' + str(nswim) + ' -nout ' + str(nout))



def get_center(state):
    return state[0].mean(axis=1)

def get_orientation(state):
    return [state[0][0,0] - get_center(state)[0], state[0][1,0] - get_center(state)[1]]
    
def compute_observation(state, u0=0.1):
    center = get_center(state)
    orientation = get_orientation(state)[0] > 0
    wind = (fluid[0,0] > -u0) + (fluid[0,0] > u0)
    return wind + 3 * orientation

def get_center(state):
    return state[0].mean(axis=1)

def get_orientation(state):
    return [state[0][0,0] - get_center(state)[0], state[0][1,0] - get_center(state)[1]]

def compute_observation(state, u0=0.1):
    center = get_center(state)
    orientation = get_orientation(state)[0] > 0
    wind = (state[-1][0,0] > -u0) + (state[-1][0,0] > u0)
    return wind + 3 * orientation

def det_to_sto_policy(pi, number_of_actions = 7):
    number_of_states = len(pi)
    Pi = np.zeros((number_of_states, number_of_actions))
    for s in range(number_of_states):
        Pi [s, pi[s]] = 1 
    return Pi

def select_action(state, Pi):
    #dependent on number of actions
    prob = Pi[state]
    return np.random.choice(np.arange(0,7), p=prob)

def get_state(iteration):
    with open(wdir + "/fiber" + str(iteration) + ".ff", 'r') as file:
        data = np.fromfile(file, dtype=float)
        position = data[0:2*(Ns+1)]
        #velocity = data[2*(Ns+1): 4*(Ns+1)]
        fluid = data[4*(Ns+1):]
        position = position.reshape(Ns+1,2).T
        fluid = fluid.reshape(Ns+1, 2).T
        return (position, position, fluid)


wdir = mkwdir()
pi = [0, 0, 0, 0, 6, 6]
obs = 4
nout = 10
for i in range(10000):
    Pi = det_to_sto_policy(pi)
    action = select_action(obs, Pi)
    if i % nout == 0: print(f"action: {action}, state: {obs}")
    swim(wdir, i, action)
    state = get_state(i+1)
    obs = compute_observation(state)