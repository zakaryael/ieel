#loading libraries
#! 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import netCDF4 as nc
from scipy.stats import norm
import scipy
from tqdm import tqdm

def arrangeQ(Q):
    '''arranges the values of a misread Q matrix (from row major to col major?)'''
    number_of_states, number_of_actions = Q.shape
    return Q.flatten().reshape(number_of_actions, number_of_states).T

def load_data(dir, nout=400, T=600, dt=5e-4, progress=False, number_of_actions=8, number_of_states=6):
    '''load simulation data from nc files '''
    numberOfFiles = int(T / (nout * dt))
    Qs = np.zeros((numberOfFiles, number_of_states, number_of_actions))
    Vel = np.zeros(numberOfFiles)
    Pos = np.zeros(numberOfFiles)
    States = np.zeros(numberOfFiles)
    Actions = np.zeros(numberOfFiles)
    number_of_time_steps = int(T / dt)
    for i in tqdm(range(0, number_of_time_steps, nout)):
        ds = nc.Dataset(dir+"/fiber"+str(i)+".nc")
        Q = ds['Qtable'][:].flatten().reshape(number_of_actions, number_of_states).T
        Qs[i // nout] = Q
        Pos[i // nout] = np.array(ds['Pos'])[0,-1]
        Vel[i // nout] = np.array(ds['Vel'])[0,0]
        States[i // nout] = int(ds["state"][:])
        Actions[i // nout] = int(ds["action"][:])
        ds.close()
        d = {"state":States, "action":Actions, "pos":Pos, "vel":Vel}
    return pd.DataFrame(data=d), Qs

def get_greedy_policy(Q):
    '''extract the greedy policy from a matrix Q'''
    number_of_states = Q.shape[0]
    number_of_actions = Q.shape[1]
    pi = 10 * np.ones(number_of_states)
    for i in range(number_of_states):
        pi[i] = np.argmax(Q[i])
    return pi

def moving_average(x, w):
    '''computes the moving average of a vector x over a sliding winow of size w'''
    return np.convolve(x, np.ones(w), 'valid') / w

#statistical analysis functions
def rewards_dist(current_state, action, data, plot=True):
    """computes and plots the distribution of conditional reward on the current state and action"""
    Actions = data["action"]
    States = data["state"]
    Pos = data["pos"]
    dist = []
    rewards = np.diff(Pos)
    N = Pos.shape[0]
    for i in range(N-1):
        if States[i] == current_state and Actions[i] == action:
            dist.append([rewards[i], i])
    if dist == []: dist.append([10, N])
    dist = np.array(dist).T
    if plot: plt.hist(dist[0], bins=1000, density=True)
    return dist

def rewards_dist_next(current_state, action, next_state, data):
    dist_next = []
    dist = rewards_dist(current_state, action, data)
    for i in dist[1]:
        if next_state == data["state"][i+1]:
            dist_next.append(dist[:, i])
    return np.array(dist_next).T

#action-state occupancy rate
def occupancy_matrix(data, freq=True):
    counts = data[["state", "action"]].value_counts(normalize=freq)
    counts_sorted = counts.sort_index()
    occupancy_matrix = np.array(counts_sorted).reshape((12,8))
    occupancy_matrix_wm = np.vstack((occupancy_matrix, occupancy_matrix.sum(axis=0)))
    occupancy_matrix_wmm = np.hstack((occupancy_matrix_wm, np.array([occupancy_matrix_wm.sum(axis=1)]).T))
    return occupancy_matrix_wmm

def transition_dist(state, action, data):
    print("hello!")
    dist = []
    for i in range(data.shape[0] - 1):
        if data.state[i] == state and data.action[i] == action:
            dist.append([i, data.state[i+1]])
    return np.array(dist).T

def transition_matrices(data, number_of_actions, number_of_states):
    P = np.zeros((number_of_states, number_of_actions, number_of_states)) #initialize transitions array
    for s in range(number_of_states):
        for a in range(number_of_actions):
            d = pd.DataFrame(transition_dist(s, a, data)[1])
            freq = d.value_counts(normalize=True)
            for state in range(number_of_states):
                if state not in freq.index:
                    freq.loc[state] = 0
            P[s, a, :] = freq.sort_index()
    return P

def Qevol(Qs):
    number_of_actions = Qs.shape[2]
    number_of_states = Qs.shape[1]
    plt.figure(figsize=(25,16))
    for s in range(number_of_states):
        for a in range(number_of_actions):
            #plt.subplot(number_of_actions,number_of_states, s * number_of_actions + a + 1)
            plt.subplot(number_of_states // 3 ,3, s+1)
            plt.plot(Qs[:, s, a], label="action " + str(a))
            plt.title("state {}".format(s))
            plt.legend()

def get_reward(state, data):
    return data[data.state==state].reward

def saveQ(Q, outputname='Q'):
    """ saves a numpy array Q to a csv file"""
    np.savetxt(outputname + ".csv", Q, fmt='%.8f', delimiter=',')
    print("Q matrix has been saved to ./" + outputname + ".csv")

def reduce_data(data):
    '''fuse buckled and non buckled states'''
    buckledr = [3, 4, 5]
    buckledl = [9, 10, 11]
    data_reduced = data.copy()
    for s in buckledr:
        data_reduced['state'] = data_reduced['state'].replace([s], s-3)
    for s in [6, 7, 8]:
            data_reduced['state'] = data_reduced['state'].replace([s], s-3)
    for s in buckledl:
        data_reduced['state'] = data_reduced['state'].replace([s], s-6)
    return data_reduced

def Qlearning(data, Qinit, lrate0, gamma):
    n = data.shape[0]
    number_of_states, number_of_actions = Qinit.shape
    Qmats = np.zeros((n, number_of_states, number_of_actions))
    Qmats[0] = Qinit
    number_of_visits = np.zeros((number_of_states, number_of_actions))
    if 'reward' not in data.columns:
        data['reward'] = np.diff(data.pos, append=0)
    for i in range(1,n):
        di = data.iloc[i-1]
        state = int(di.state)
        action = int(di.action)
        number_of_visits[state, action] += 1
        lrate = lrate0 #/ np.log(number_of_visits[state, action]+1)
        next_state = int(data.iloc[i].state)
        Qmats[i] = Qmats[i-1]
        Qmats[i, state, action] = (1 - lrate0) * Qmats[i, state, action] +   + lrate * (di.reward + gamma *  Qmats[i, next_state].max())
    return Qmats

def load_from_binary_file(dir, number_of_states=6, number_of_actions=8):
    file = open(dir, 'r')
    data = np.fromfile(file, dtype=float)
    number_of_entries = int(data.shape[0] / (6 + number_of_actions * number_of_states))
    data = data.reshape(number_of_entries, (6 + number_of_actions * number_of_states))
    Qs = data[:, 6:]
    header = ['time', 'step', 'position', 'reward', 'state', 'action']
    df = pd.DataFrame(data[:, :6], columns=header)
    return df.convert_dtypes(), Qs

def MC_rewards_matrix(df, number_of_actions, number_of_states):
    R = np.zeros((number_of_states, number_of_actions))
    for s in range(number_of_states):
        for a in range(number_of_actions):
            R[s,a] = np.mean(rewards_dist(s, a, df, plot=False)[0])
    return R

def value_iteration(gamma, n_iter, number_of_states, transition_matrix):
    ## Value iteration algorithm
    V = np.zeros((n_iter, number_of_states))
    for k in range(n_iter-1):
        for state in range(number_of_states):
            V[k+1, state] = np.max(R[state] + gamma * P[state].dot(V[k]))
    return V

def get_Q_from_v(v, R, P):
    return R + gamma * P[s].dot(v)

def mc_evol(y, level=0.99):
    """computes monte carlo evolution with confidence level bounds"""
    n = y.shape[0]
    nvec = np.arange(1, n+1)
    delta = np.cumsum(y) / nvec
    Var = np.zeros(n)
    Var[1:] =  (np.cumsum(y**2)[1:] - nvec[1:] * (delta[1:] ** 2)) / nvec[:n-1]; Var[0] = 0
    q = norm.interval(level)[1]
    i1 = delta - q * np.sqrt (Var / nvec)
    i2 = delta + q * np.sqrt (Var / nvec)
    #return {"delta": delta, "i1": i1, "i2": i2}
    return delta, i1, i2, nvec
    
def mc_evol_plot(evol, level, label=""):
    """plots mc evolution with confidence bounds """
    estimation, lower_bound, upper_bound, xx = evol
    plt.plot(xx, estimation, label=label)
    #plt.fill_between(xx, lower_bound, upper_bound, alpha=.2, label=str(100*level)+'% level confidence interval')
    plt.legend()

#compute transition matrices
    
def transition_matrices(data, number_of_actions, number_of_states):
    """Computes the transition matrices"""
    P = np.zeros((number_of_states, number_of_actions, number_of_states)) #initialize transitions array
    for s in range(number_of_states):
        for a in range(number_of_actions):
            d = data.next_state[(data.current_state == s) & (data.action == a+1)]
            freq = d.value_counts(normalize=True)
            for state in range(number_of_states):
                if state not in freq.index:
                    freq.loc[state] = 0
            P[s, a, :] = freq.sort_index()
    return P

def mixed_policy(pi, number_of_actions = 8, sim=False):
    pi = pi.astype(int)
    number_of_states = len(pi)
    Pi = np.zeros((number_of_states, number_of_actions))
    for s in range(number_of_states):
        Pi [s, pi[s] - 1 * sim] = 1 
    return Pi

#simulate learning
def get_next_state(current_state, action, data):
    return np.random.choice(data.next_state[(data.current_state == current_state) & (data.action == action)]) 

def get_reward(current_state, action, next_state, data):
   return np.random.choice(data.reward[(data.current_state == current_state) & (data.action == action) & (data.next_state == next_state)])

def select_action(state, Pi):
    prob = Pi[state]
    return np.random.choice(np.arange(1,8), p=prob)

def update_Q(Q, current_state, action, next_state, reward, gamma, alpha):
    action -= 1
    Q[current_state, action] = (1 - alpha) * Q[current_state, action] + alpha * (reward + gamma * Q[next_state].max())
    return Q

def update_policy(Q, pi, state, epsilon):
    a = Q[state].argmax()
    pi[state] = epsilon / 6
    pi[state, a] = 1-epsilon
    return pi

def sim_ql(Qinit, Pinit, gamma, alpha0, epsilon0, T, learn, data, decreasing_lrate=True):
    print("Starting QL simulation: \n initial policy: {} \n initial Q matrix: {} \n gamma = {}, epsilon0 = {}, alpha0 = {}, learn = {}".format(np.round(Pinit, 2), np.round(Qinit, 2), gamma, epsilon0, alpha0, learn ))
    pi = Pinit
    Q = Qinit
    alpha = alpha0
    epsilon = epsilon0
    number_of_states, number_of_actions = Q.shape
    number_of_visits = np.zeros((number_of_states, number_of_actions))
    Qmats = np.zeros((T+1, number_of_states, number_of_actions))
    Qmats[0] = Q
    current_state = data.current_state.iloc[0] #always start from the same initial state (does it matter?)
    list = []
    for t in range(T):
        action = select_action(current_state, pi)
        next_state = get_next_state(current_state, action, data)
        reward = get_reward(current_state, action, next_state, data)
        number_of_visits[current_state, action-1] += 1
        if decreasing_lrate: alpha = alpha0 / np.sqrt(number_of_visits[current_state, action-1])
        Q = update_Q(Q, current_state, action, next_state, reward, gamma, alpha)
        epsilon = epsilon0 #/ number_of_visits[current_state].sum()
        if learn: pi = update_policy(Q, pi, current_state, epsilon)
        Qmats[t+1] = Q
        list.append([t, current_state, action, next_state, reward, number_of_visits[current_state, action-1]])
        current_state = next_state
    return pd.DataFrame(list, columns=['time', 'current_state', 'action', 'next_state', 'reward', 'state-action number of visits']), Qmats
