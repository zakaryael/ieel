
from utils import *

#parameters:
phys_params = {'length': 1.0, 'zeta': 2500.0, 'EI': 1.0,
    'penalisation': 400, 'timestep': 0.001, 'Ns': 200, 'wavenumber': 2, 'frequency': 20, 'alpha': 1,
    'Velocity': 0.5, 'speed': 0.1}

pg_params = {'learning_param1':5e-3, 'learning_param2':2.5e-2}

common_params = {
    'make_output': 0, 'step_out':10, 'step_learning':10, 'step_save_learning':10
}

#make a working directory
folder_name = 'pg'
wdir = mkwdir(folder_name)

#create theta0.csv where it should be created
pi0 = [3, 3, 3, 3, 6, 6]
Pi0 = det_to_sto_policy(pi0)
theta0 = 1 * Pi0
save_csv(theta0, wdir + "/theta0")
save_csv(theta0, '../input_data/theta0')

#save the various parameters
command = write_command('../runs/bin/policy_gradient', wdir, time=1e6, outdir=wdir, **{**phys_params, **pg_params, **common_params})

#run the command
os.system(command)

