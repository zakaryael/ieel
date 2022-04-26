#!/opt/homebrew/anaconda3/bin/python
from utils import *

wdir = mkwdir()
print(wdir)
pi = [0, 0, 0, 0, 6, 6]
Pi = det_to_sto_policy(pi)
obs = 4 #initial state
nout = 10
for i in range(11):
    action = select_action(obs, Pi)
    if i % nout == 0: print(f"action: {action}, state: {obs}")
    swim(wdir, i, action, nout=1)
    state = get_state(i+1, wdir)
    obs = compute_observation(state)