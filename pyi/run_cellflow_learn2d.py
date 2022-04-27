
from utils import *

#parameters:
length = 1.0
zeta = 2.5e3
EI = 1.0
beta = 400
timestep = 1e-3
Ns = 200
wavenumber = 2
frequency = 20
alpha = 1
Velocity = 0.5
speed = 0.1
#
learning = 0
beta0 = 0.01
alpha0 = 0.1

#make a working directory
wdir = mkwdir()

#create theta0.csv where it should be created
pi0 = [0, 0, 0, 6, 6, 6]
Pi0 = det_to_sto_policy(pi0)
theta0 = 1 * Pi0
save_csv(Pi0, '../input_data/Pi')

#save the various parameters
command = write_command('../runs/bin/learn2D_cellflow', wdir, outdir=wdir,\
            time=1e4,length=length, zeta=zeta, EI=EI, penalisation=beta, timestep=timestep,\
            Ns=Ns, wavenumber=wavenumber, frequency=frequency, alpha=alpha,\
            Velocity=Velocity, learning=learning,\
            step_learning=10, step_save_learning=10)

#run the command
os.system(command)

