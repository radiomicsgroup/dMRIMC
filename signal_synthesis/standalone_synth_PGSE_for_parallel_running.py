import numpy as np
import glob as gb
import pickle as pk
import time
import os
from tqdm import tqdm
from numba import njit
import argparse
import subprocess
from collections import deque

@njit
def calculate_phase_accrual(Gtime_row, Xpos, Ypos, gammar, dT, Nmolecules):
    """Calculate phase accrual for all molecules
    
    Args:
        Gtime_row: 1D array, the gradient waveform (already extracted from Gtime[0,:])
        Xpos, Ypos: 2D arrays (Nmolecules, Nsteps)
    """
    phix = np.zeros(Nmolecules)
    phiy = np.zeros(Nmolecules)
    
    for nn in range(Nmolecules):
        phix[nn] = -gammar * dT * np.sum(Gtime_row * Xpos[nn, :])
        phiy[nn] = -gammar * dT * np.sum(Gtime_row * Ypos[nn, :])
    
    return phix, phiy

def load_trajectories_2D(Nspins, Nsteps, Tdur, file):
    # Create array of time
    T = np.linspace(0,Tdur,Nsteps) # in seconds 
    # Load trajectories
    # The bin files of the trajectories need to be imported as float32
    # or you get garbage
    data = np.array(np.fromfile(file, dtype="float32"))
    # The simulator returns the numbers in mm so we make them m
    data = 0.001*data
    # Trajectories will be stored in a matrix with dims:
    # Nspins x Nsteps
    # Remember the way the data are stored from the simulator!!

    ###### x-component
    x = data[0:np.copy(data).size:3]
    X = np.zeros((Nspins, Nsteps))
    deltaX = np.zeros((Nspins, Nsteps))
    # For every time step, create a row containing all the positions of
    # spin j
    for j in range(Nsteps):
        X[:,j] = x[j:x.size:Nsteps]
    # Do the same for the Δx for every spin
    for k in range(Nspins):
        deltaX[k,:] = X[k,:] - X[k,0]
    
    ###### y-component
    y = data[1:np.copy(data).size:3]
    Y = np.zeros((Nspins,Nsteps))
    deltaY = np.zeros((Nspins,Nsteps))
    for tt in range(0,Nsteps):
        Y[:,tt] = y[tt:y.size:Nsteps]
    for nn in range(0,Nspins):
        deltaY[nn,:] = Y[nn,:] - Y[nn,0]

    return deltaX, deltaY, T


def construct_signal_PGSE(inlist):
    """
    Simulate an MRI signal from a set of random walks given the files containing
    the trajectories of the spins and information about the MRI sequence
    """
    snr = inlist[0]             # Signal-to-noise ratio (SNR = S(0) / sigma)
    b_value = inlist[1]*1e6     # b-value in sec/m^2
    delta1 = inlist[2]*1e-3     # Gradient duration 1 in sec
    delta2 = inlist[3]*1e-3     # Gradient duration 2 in sec
    Delta12 = inlist[4]*1e-3    # Gradient separation1,2 in sec

    # Get info for loading trajectories
    Nspinfile = inlist[5]   # number of spins per files
    Nsteps = inlist[6]      # number of time steps
    Tdur = inlist[7]       # total simulation duration in sec
    file = inlist[8]       # File with MRI trajectories

    # Load files
    Xpos, Ypos, T = load_trajectories_2D(Nspinfile,Nsteps,Tdur,file)
    Nmolecules = Xpos.shape[0]      # Total number of water molecules
    dT = T[1] - T[0]                # duration of time step in sec
    # Calculate gradient strength
    gammar = 267.522187e6           # Gyromagnetic ratio in 1/(sec T)
    Nmeas = b_value.size            # Number of measurements
    Gval = np.zeros(Nmeas)
    # Places where δ > Δ are invalid
    if np.any(delta1 > Delta12):
        invalid_indexes = np.where(delta1 > Delta12)
    
    for mm in range(Nmeas):
        if(b_value[mm]>0.0):
            pasta = (gammar**2) * (delta1[mm]**2) * (Delta12[mm] - (delta1[mm]/3))
            Gval[mm] = np.sqrt( b_value[mm]/ pasta)
    # This does not get rid of the issue, it is here for emphasis
    if np.any(delta1 > Delta12):
        Gval[invalid_indexes] = np.nan
    # Construct the MRI signal
    mri_sig = np.zeros((b_value.shape))
    Ggs = []
    for mm in range(Nmeas):
        # b_value = 0 image
        # T2* weighted image, no diffusion attenuation of the signal
        if b_value[mm] == 0.0:
            mri_sig[mm] = 1.0
        
        # Diffusion-weighted image
        # b != 0

        else:
            # Generate gradient waveform
            # This is the from the graph of the sequence, the peaks of the
            # magnetic fields during the times that the coils are activating
            Gtime = np.zeros((1,T.size))
            Gtime[0,T<=delta1[mm]] = Gval[mm]
            Gtime[0,T>Delta12[mm]] = -Gval[mm]
            Gtime[0,T>Delta12[mm]+delta2[mm]] = 0.0
            # Make sure they add up to zero
            more_neg = len(np.argwhere(Gtime < 0)) > len(np.argwhere(Gtime > 0))
            more_pos = len(np.argwhere(Gtime > 0)) > len(np.argwhere(Gtime < 0))
            eq = len(np.argwhere(Gtime > 0)) == len(np.argwhere(Gtime < 0))
            if not eq:
                if more_neg:
                    while more_neg:
                        # find last positive entry and make the next one the same
                        indx = np.argwhere(Gtime[0] > 0)[-1]
                        Gtime[0][indx+1] = Gtime[0][indx][0]
                        more_neg = len(np.argwhere(Gtime < 0)) > len(np.argwhere(Gtime > 0))
                elif more_pos:
                    while more_pos:
                        # same procedure for negative entries
                        indx = np.argwhere(Gtime[0] < 0)[-1]
                        Gtime[0][indx+1] = Gtime[0][indx][0]
                        more_pos = len(np.argwhere(Gtime > 0)) > len(np.argwhere(Gtime < 0))

            Ggs.append(Gtime)
            # Allocate variables to store phase accruals for 2 orthogonal gradients
            # phix = np.zeros(Nmolecules)
            # phiy = np.zeros(Nmolecules)

            # Calculate phase accruals and magnetic moments for gradient directions along x, y, z
            # for nn in range(0,Nmolecules):
            #     # phase acrual for nn-th spin when gradient direction is [1 0 0]
            #     phix[nn] = -1.0*gammar*dT*np.sum(Gtime[0,:]*Xpos[nn,:])
            #     # phase acrual for nn-th spin when gradient direction is [0 1 0]
            #     phiy[nn] = -1.0*gammar*dT*np.sum(Gtime[0,:]*Ypos[nn,:])
            # Numba version
            phix, phiy = calculate_phase_accrual(Gtime[0,:], Xpos, Ypos, gammar, dT, Nmolecules)
            # Get complex magnetic moments for all spins
            momx = np.exp(-1j*phix)
            momy = np.exp(-1j*phiy)

            # Average over spin ensemble
            # There are a few cases of NaN values in the deltaX/Y in some
            # substrates so we go for nanmean
            sx = np.nanmean(momx)
            sy = np.nanmean(momy)
            
            # Store synthetic magnitude MRI signal
            mri_sig[mm] = (np.abs(sx) + np.abs(sy))/2.0
    
    # Now we are ok
    if np.any(delta1 > Delta12):
        mri_sig[invalid_indexes] = np.nan
    # Add noise
    mri_sig_noisy = np.sqrt(( mri_sig + (1.0/snr)*np.random.randn(mri_sig.size))**2  +  ((1.0/snr)*np.random.randn(mri_sig.size) )**2)

    return [mri_sig,mri_sig_noisy,Gval,Ggs]


def synthesis_misc_PGSE(Nspins, Ntime, Tdur, substrate, SEQUENCE_NAME):
    # Simulation parameters
    Nsteps = Ntime + 1  # number of simulation interations
    # Random walks directory
    trajdir = f"./simulations/{substrate}"
    # Get all trajectory files
    traj_files = sorted(gb.glob(f"{trajdir}/random_walks/*_0.traj"))
    # Sequence parameters
    snr = 20 # we wont keep the noisy signals anyway
    sequence = SEQUENCE_NAME
    outstring = f"./synthesized_signals/{sequence}/{substrate}"
    if not os.path.exists(outstring):
        os.makedirs(outstring)

    sequence_params_path = f"./sequence_parameters/{sequence}"
    bvalseq = np.loadtxt(f"{sequence_params_path}/custom.bval", delimiter=' ')            # List of b-value in sec/mm2
    delta1seq = np.loadtxt(f"{sequence_params_path}/custom.gdur1", delimiter=' ')         # List of delta1 in msec
    delta2seq = np.loadtxt(f"{sequence_params_path}/custom.gdur2", delimiter=' ')         # List of delta2 in msec
    Delta12seq = np.loadtxt(f"{sequence_params_path}/custom.gsep12", delimiter=' ')       # List of Delta1,2 in msec
    ########## START OF PROCESSING ##########
    print(f"Starting signal synthesis")
    print(f"Signals to construct: {len(traj_files)}")
    start = time.time()
    # Log file for errors
    errors = 0
    f = open(f"{outstring}/errors.txt",'at')
    for j in (pbar :=tqdm(traj_files,total=len(traj_files))):
        # print(j)
        traj_file = j.replace(trajdir, '')
        p = traj_file.replace(f"{trajdir}/",'')
        p = p.replace("_0.traj",'')
        p = p.replace("/random_walks/",'')
        pbar.set_description(f"Processing trajectories for file {p}")
        # print(p)
        # Parameter list
        processinglist = [snr,bvalseq,delta1seq,delta2seq,Delta12seq,Nspins,Nsteps,Tdur,j]
        # Compute MRI signals
        try:
            fitresults = construct_signal_PGSE(processinglist)
            outlist = [fitresults,p]
        except Exception as e:
            print(f"ERROR: {e}")
            errors += 1
            f.write(f"{p}")
        # Save results
        # print(f'{outstring}/{p}.signals.bin')
        hfile = open(f'{outstring}/{p}.signals.bin','wb')
        pk.dump(outlist,hfile,protocol=5) 
        hfile.close()
    f.close()
    end = time.time()
    print(f"Done")
    print(f"Total time: {(end-start):.2f}s")
    print(f"Total errors: {errors}")


def synth_run_bash_commands_in_parallel_with_timeout(commands, n_parallel, t, retry=False):
    """
    Run a list of bash commands in parallel with maximum number of processes
    and kill the process if it takes more than t seconds.
    """
    # Log files for errors and timeouts
    error_file = open("PARAL_SYNTH_errors.txt", 'wt')
    timeout_file = open("PARAL_SYNTH_timed_out.txt", 'wt')

    waiting = deque([(command, 1) for command in commands])
    running = deque()
    cnt_errs = 0
    cnt_timeouts = 0
    cnt_completed = 0
    while len(waiting) > 0 or len(running) > 0:
        print(f"Running: {len(running)}, Waiting: {len(waiting)}, Completed: {cnt_completed}, Errors: {cnt_errs}, Timeouts: {cnt_timeouts}", end='\r')
        while len(waiting) > 0 and len(running) < n_parallel:
            command, tries = waiting.popleft()
            try:
                print(command)
                process = subprocess.Popen(command, stdout=subprocess.DEVNULL)
                running.append((process, command, tries, time.time()))
                print(f"Working on {command[9]}")
            except OSError:
                print(f"{command[9]} failed to start")
                error_file.write(f"{command}\n")
                cnt_errs += 1

        for _ in range(len(running)):
            process, command, tries, start_time = running.popleft()
            ret = process.poll()
            if ret is None:
                if time.time() - start_time > t:
                    process.kill()
                    timeout_file.write(f"{command}\n")
                    print(f"{command} timed out")
                    cnt_timeouts += 1
                else:
                    running.append((process, command, tries, start_time))
            elif ret != 0:
                print(f"{command[9]} returned an error")
                error_file.write(f"{command}\n")
                cnt_errs += 1
            else:
                exec_time = np.round((time.time() - start_time)/60, decimals=4)
                print(f"{command[9]} finished successfully after {exec_time} minutes")
                cnt_completed += 1
        time.sleep(0.5)

    print('All tasks done')
    error_file.close()
    timeout_file.close()

parser = argparse.ArgumentParser(description="Synthesize PGSE signal from trajectory file")

# Required arguments
parser.add_argument("--Nspins", type=int, required=True, help="Number of spins")
parser.add_argument("--Ntime", type=int, required=True, help="Number of time steps")
parser.add_argument("--Tdur", type=float, required=True, help="Duration of the simulation")
parser.add_argument("--substrate", type=str, required=True, help="Substrate type")
parser.add_argument("--SEQUENCE_NAME", type=str, required=True, help="Name of the sequence")

args = parser.parse_args()

# Example usage of the arguments
print(f"Running synthesis_misc_PGSE with:")
print(f"Nspins: {args.Nspins}")
print(f"Ntime: {args.Ntime}")
print(f"Tdur: {args.Tdur}")
print(f"Substrate: {args.substrate}")
print(f"SEQUENCE_NAME: {args.SEQUENCE_NAME}")

synthesis_misc_PGSE(Nspins=int(args.Nspins), Ntime=int(args.Ntime), Tdur=float(args.Tdur), substrate=str(args.substrate), SEQUENCE_NAME=str(args.SEQUENCE_NAME))