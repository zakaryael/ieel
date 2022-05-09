import os
import numpy as np
import subprocess

def mkwdir(folder_name='wdir', where='../output_data/'):
    counter = 0
    while os.path.isdir(where + folder_name + str(counter)): #improvement: use os.path.join()
        counter += 1
    path = where + folder_name + str(counter)
    os.makedirs(where + folder_name + str(counter))
    return path
    
def swim(wdir, iteration, action, nswim=10, nout=10):
    subprocess.Popen('../runs/bin/swim -dir ' + wdir + ' -it ' + str(iteration) + ' -a ' + str(action) + ' -nsw ' + str(nswim) + ' -nout ' + str(nout), shell=True).wait()

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

def select_action(state, Pi, na=7):
    #dependent on number of actions
    prob = Pi[state]
    return np.random.choice(np.arange(0,na), p=prob)

def get_state(iteration, wdir):
    Ns = 200
    with open(wdir + "/fiber" + str(iteration) + ".ff", 'r') as file:
        data = np.fromfile(file, dtype=float)
        position = data[0:2*(Ns+1)]
        #velocity = data[2*(Ns+1): 4*(Ns+1)]
        fluid = data[4*(Ns+1):]
        position = position.reshape(Ns+1,2).T
        fluid = fluid.reshape(Ns+1, 2).T
        return (position, position, fluid)

def write_command(exec,wdir, **kwargs):
    file = open(wdir+"/arguments.txt", "w")
    str_dictionary = repr(kwargs)
    file.write(str_dictionary + "\n")
    file.close()
    command = exec
    for key, value in kwargs.items():
        command = command + " --" + key + " " + str(value)
    return command

def save_csv(array, outdir):
    """ saves a numpy array  to a csv file"""
    np.savetxt(outdir + ".csv", array, delimiter=',')
    print("array has been saved to" +  outdir + ".csv")