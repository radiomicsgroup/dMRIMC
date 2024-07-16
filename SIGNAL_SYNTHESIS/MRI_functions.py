import numpy as np
import glob as gb
import pickle as pk
import time
import os
from tqdm import tqdm


def load_trajectories_2D(Nspins, Nsteps, Tdur, file):
    # Create array of time
    T = np.linspace(0, Tdur, Nsteps)  # in seconds
    # Load trajectories
    # The bin files of the trajectories need to be imported as float32
    # or you get garbage
    data = np.array(np.fromfile(file, dtype="float32"))
    # The simulator returns the numbers in mm so we make them m
    data = 0.001 * data
    # Trajectories will be stored in a matrix with dims:
    # Nspins x Nsteps

    # x-component
    x = data[0:np.copy(data).size:3]
    X = np.zeros((Nspins, Nsteps))
    deltaX = np.zeros((Nspins, Nsteps))
    # For every time step, create a row containing all the positions of
    # spin j
    for j in range(Nsteps):
        X[:, j] = x[j:x.size:Nsteps]
    # Do the same for the Δx for every spin
    for k in range(Nspins):
        deltaX[k, :] = X[k, :] - X[k, 0]

    # y-component
    y = data[1:np.copy(data).size:3]
    Y = np.zeros((Nspins, Nsteps))
    deltaY = np.zeros((Nspins, Nsteps))
    for tt in range(0, Nsteps):
        Y[:, tt] = y[tt:y.size:Nsteps]
    for nn in range(0, Nspins):
        deltaY[nn, :] = Y[nn, :] - Y[nn, 0]

    return deltaX, deltaY, T


def construct_signal_PGSE(inlist):
    """
    Simulate an MRI signal from a set of random walks given the files containing
    the trajectories of the spins and information about the MRI sequence

    Returns: MRI signal at infinite SNR; MRI signal at finite SNR,
    the value of the gradient, the timings of the gradient

    Note: In the current pipeline we are polluting the signal with noise
    later on but the snr variable is left here for compatibility reasons
    """
    snr = inlist[0]             # Signal-to-noise ratio (SNR = S(0) / sigma)
    b_value = inlist[1] * 1e6     # b-value in sec/m^2
    delta1 = inlist[2] * 1e-3     # Gradient duration 1 in sec
    delta2 = inlist[3] * 1e-3     # Gradient duration 2 in sec
    Delta12 = inlist[4] * 1e-3    # Gradient separation1,2 in sec

    # Get info for loading trajectories
    Nspinfile = inlist[5]   # number of spins per files
    Nsteps = inlist[6]      # number of time steps
    Tdur = inlist[7]       # total simulation duration in sec
    file = inlist[8]       # File with MRI trajectories

    # Load files
    Xpos, Ypos, T = load_trajectories_2D(Nspinfile, Nsteps, Tdur, file)
    Nmolecules = Xpos.shape[0]      # Total number of water molecules
    dT = T[1] - T[0]                # duration of time step in sec

    # Calculate gradient strength
    gammar = 267.522187e6           # Gyromagnetic ratio in 1/(sec T)
    Nmeas = b_value.size            # Number of measurements
    Gval = np.zeros(Nmeas)
    for mm in range(Nmeas):
        if(b_value[mm] > 0.0):
            pasta = (gammar**2) * (delta1[mm]**2) * \
                (Delta12[mm] - (delta1[mm] / 3))
            Gval[mm] = np.sqrt(b_value[mm] / pasta)
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
            Gtime = np.zeros((1, T.size))
            Gtime[0, T <= delta1[mm]] = Gval[mm]
            Gtime[0, T > Delta12[mm]] = -Gval[mm]
            Gtime[0, T > Delta12[mm] + delta2[mm]] = 0.0

            # Make sure they add up to zero
            more_neg = len(np.argwhere(Gtime < 0)) > len(
                np.argwhere(Gtime > 0))
            more_pos = len(np.argwhere(Gtime > 0)) > len(
                np.argwhere(Gtime < 0))
            eq = len(np.argwhere(Gtime > 0)) == len(np.argwhere(Gtime < 0))
            if not eq:
                if more_neg:
                    while more_neg:
                        # find last positive entry and make the next one the same
                        indx = np.argwhere(Gtime[0] > 0)[-1]
                        Gtime[0][indx + 1] = Gtime[0][indx][0]
                        more_neg = len(np.argwhere(Gtime < 0)) > len(
                            np.argwhere(Gtime > 0))
                elif more_pos:
                    while more_pos:
                        # same procedure for negative entries
                        indx = np.argwhere(Gtime[0] < 0)[-1]
                        Gtime[0][indx + 1] = Gtime[0][indx][0]
                        more_pos = len(np.argwhere(Gtime > 0)) > len(
                            np.argwhere(Gtime < 0))

            Ggs.append(Gtime)
            # Allocate variables to store phase accruals for 2 orthogonal gradients
            phix = np.zeros(Nmolecules)
            phiy = np.zeros(Nmolecules)

            # Calculate phase accruals and magnetic moments for gradient directions along x, y, z
            for nn in range(0, Nmolecules):
                # phase acrual for nn-th spin when gradient direction is [1 0 0]
                phix[nn] = -1.0 * gammar * dT * \
                    np.sum(Gtime[0, :] * Xpos[nn, :])
                # phase acrual for nn-th spin when gradient direction is [0 1 0]
                phiy[nn] = -1.0 * gammar * dT * \
                    np.sum(Gtime[0, :] * Ypos[nn, :])
            # Get complex magnetic moments for all spins
            momx = np.exp(-1j * phix)
            momy = np.exp(-1j * phiy)

            # Average over spin ensemble
            # There are a few cases of NaN values in the deltaX/Y in some
            # substrates so we go for nanmean
            sx = np.nanmean(momx)
            sy = np.nanmean(momy)

            # Store synthetic magnitude MRI signal
            mri_sig[mm] = (np.abs(sx) + np.abs(sy)) / 2.0
    # Add noise
    mri_sig_noisy = np.sqrt((mri_sig + (1.0 / snr) * np.random.randn(mri_sig.size))
                            ** 2 + ((1.0 / snr) * np.random.randn(mri_sig.size))**2)

    return [mri_sig, mri_sig_noisy, Gval, Ggs]


def synthesis_misc_PGSE(Nspins, Ntime, Tdur, substrate, substrate_type, snr, sequence):

    # Simulation parameters
    Nsteps = Ntime + 1  # number of simulation interations
    # Random walks directory
    if substrate_type == "cells":
        trajdir = f"../SIMULATIONS/{substrate}_CELLS"
    elif substrate_type == "EXTRA":
        trajdir = f"../SIMULATIONS/{substrate}_EXTRA"
    # Get all trajectory files
    traj_files = sorted(gb.glob(f"{trajdir}/random_walks/*_0.traj"))
    # Sequence parameters
    outstring = f"{sequence}/{sequence}_{substrate}_{substrate_type}"
    if not os.path.exists(outstring):
        os.makedirs(outstring)

    sequence_params_path = f"./sequence_parameters/{sequence}"
    
    # List of b-value in sec/mm2
    bvalseq = np.loadtxt(f"{sequence_params_path}/custom.bval", delimiter=' ')
    # List of delta1 in msec
    delta1seq = np.loadtxt(
        f"{sequence_params_path}/custom.gdur1", delimiter=' ')
    # List of delta2 in msec
    delta2seq = np.loadtxt(
        f"{sequence_params_path}/custom.gdur2", delimiter=' ')
    # List of Delta1,2 in msec
    Delta12seq = np.loadtxt(
        f"{sequence_params_path}/custom.gsep12", delimiter=' ')
    
    ########## START OF SYNTHESIS ##########
    print(f"Starting signal synthesis")
    print(f"Signals to construct: {len(traj_files)}")
    start = time.time()
    # Log file for errors
    errors = 0
    f = open(f"{outstring}/errors.txt", 'at')
    for j in (pbar := tqdm(traj_files, total=len(traj_files))):
        traj_file = j.replace(trajdir, '')
        pbar.set_description(f"Processing trajectories for file {traj_file}")
        p = traj_file.replace(f"{trajdir}/", '')
        p = p.replace("_0.traj", '')
        p = p.replace("/random_walks", '')
        processinglist = [snr, bvalseq, delta1seq,
                          delta2seq, Delta12seq, Nspins, Nsteps, Tdur, j]
        
        # Compute MRI signals
        try:
            fitresults = construct_signal_PGSE(processinglist)
            outlist = [fitresults, p]
            # Save results
            hfile = open(f'{outstring}/{p}.signals.bin', 'wb')
            pk.dump(outlist, hfile, protocol=5)
            hfile.close()
        except:
            print("ERROR")
            errors += 1
            f.write(f"{p}\n")
    f.close()
    end = time.time()
    print(f"Done")
    print(f"Total time: {(end-start):.2f}s")
    print(f"Total errors: {errors}")
