from utils import *

#parameters:
phys_params = {'length': 1.0, 'zeta': 2500.0, 'EI': 1.0,
    'penalisation': 400, 'timestep': 0.001, 'Ns': 200, 'wavenumber': 2, 'frequency': 20, 'alpha': 1,
    'Velocity': 0.5, 'speed': 0.1}

specific_params = {'learnrate':5e-3, 'learning': 0, 'angle': -45}

common_params = {
    'make_output': 1, 'step_out':10, 'step_learning':10, 'step_save_learning':10
}

#make a working directory
wdir = mkwdir('cellflow')

pi0 = [3, 3, 3, 3, 6, 6]
Pi0 = det_to_sto_policy(pi0)
save_csv(Pi0, '../input_data/Pi')
save_csv(pi0, wdir + "/policy")
#save the various parameters
#angle = np.random.rand() * 180
#print(angle)
command = write_command('../runs/bin/learn2D_cellflow', wdir, outdir=wdir, time=10000, **{**phys_params, **specific_params, **common_params})

#run the command
os.system(command)