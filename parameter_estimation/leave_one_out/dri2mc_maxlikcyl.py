import argparse, os, sys
import nibabel as nib
import numpy as np
from scipy.optimize import minimize
import multiprocessing as mp
import time



def _get_roots():
	# First ten roots of J'1(x) = 0, where J'1(x) is the derivative of the Bessel function of first kind and order 1
	rootvals = np.array([1.84118378, 5.33144277, 8.53631637, 11.7060049, 14.86358863, 18.01552786, 21.16436986, 24.31132686, 27.45705057, 30.60192297])	
	return rootvals



def _adccylperpGPA(mripar,tpar):

	# Get input stuff
	b_smm2 = mripar[0,:]           # b-value in s/mm2
	gd_ms_array = mripar[1,:]      # gradient duration delta in ms
	gs_ms_array = mripar[2,:]      # gradient separation Delta in ms
	b_msum2 = b_smm2/1000.0        # b-value in ms/um2	
	L_um = tpar[0]                 # Cylinder diameter L in um
	R_um = L_um/2.0                # Cylinder radius R = L/2.0 in um
	D0_um2ms = tpar[1]             # Intrinsic cylinder diffusivity D0 in um2/ms
	nmeas = mripar.shape[1]        # Number of measurements 

	# Get intra-axonal apparent diffusion coefficient ADCintra in um2/ms (for large axons)
	jroots = _get_roots()                # Roots of J'1(x) = 0, where J'1(x) is the derivative of the Bessel function of first kind and order 1
	Nseries = jroots.size                # Number of terms to be kept in the truncated series
	ADCi_um2ms_array = np.zeros(nmeas)   # Allocate output array for ADCin
	for nn in range(0,nmeas):
		wfactor_um2xms = 0.0         # Weighting factor in um2*ms2
		gd_ms = gd_ms_array[nn]      # Gradient duration delta in ms of current measurement
		gs_ms = gs_ms_array[nn]      # Gradient separation Delta in ms of current measurement
		bv = b_msum2[nn]             # b-value in ms/um2 of current measurement
		if( (gd_ms==0) or (gs_ms==0) or (bv==0) ):
			ADCi_um2ms_array[nn] = np.nan
		else:
			for mm in range(0,Nseries):
				a_um = jroots[mm]/R_um      # alpham expressed in 1/um
				num_value =  2.0*D0_um2ms*a_um*a_um*gd_ms - 2.0 + 2.0*np.exp(-D0_um2ms*a_um*a_um*gd_ms) + 2.0*np.exp(-D0_um2ms*a_um*a_um*gs_ms) - np.exp(-D0_um2ms*a_um*a_um*(gs_ms - gd_ms)) - np.exp(-D0_um2ms*a_um*a_um*(gs_ms + gd_ms)) 
				den_value = a_um*a_um*a_um*a_um*a_um*a_um*(R_um*R_um*a_um*a_um - 1.0)
				wfactor_um2xms = wfactor_um2xms + num_value/den_value		
			ADCi_um2ms_array[nn] = ( 2.0/( D0_um2ms*D0_um2ms*gd_ms*gd_ms*( gs_ms - gd_ms/3.0 )  ) )*wfactor_um2xms     # ADC in um2/ms
		
	# Return intra-axonal perpendicular ADC in um2/ms
	# Same as in Van Geldern et al, "Evaluation of restricted diffusion in cylinders. Phosphocreatine in rabbit leg muscle"
	# J Magn Reson B 1994, 103(3):255-60, doi: 10.1006/jmrb.1994.1038
	return ADCi_um2ms_array



def _sigscylperpGPA(mripar,tpar):

	# Get input stuff
	b_smm2 = mripar[0,:]           # b-value in s/mm2
	gd_ms_array = mripar[1,:]      # gradient duration delta in ms
	gs_ms_array = mripar[2,:]      # gradient separation Delta in ms
	b_msum2 = b_smm2/1000.0        # b-value in ms/um2
		
	# Get ADCintra in um2/ms
	adcval = _adccylperpGPA(mripar,tpar)   # tpar = [L,D0] with L -> cell diameter L in um, D0 -> intrinsic cytosol diffusivity D0 in um2/ms
	
	# Get MRI signals
	sig = np.exp( - b_msum2*adcval )
	
	# Make sure signal = 1.0 when no diffusion gradients have been turned on
	sig[b_msum2==0] = 1.0
	sig[gs_ms_array==0] = 1.0
	sig[gd_ms_array==0] = 1.0
	
	# Return signals
	return sig




def _sigDinDexTD(acqpar,rweigh,tispar):

	# Get sequence parameters
	mybvals = acqpar[0,:]           # bvalues in s/mm2
	ldelta = acqpar[2,:]            # Gradient separation in ms

	# Get tissue parameters
	lval = tispar[0]                # cell size in um2/ms
	d0val = tispar[1]               # intrinsic cell diffusivity in um2/ms
	fval = tispar[2]                # intra-cellular signal fraction
	adcexinf = tispar[3]            # extra-cellular diffusivity in um2/ms
	Betaval = tispar[4]             # time-dependent extra-cellular coefficient in um2
	
	# Synthesise MRI signal
	dex = adcexinf + Betaval/ldelta
	mrisig = fval*_sigscylperpGPA(acqpar,np.array([lval,d0val])) +  ( 1 - fval )*np.exp(-(mybvals/1000.0)*dex)
	mrisig[mybvals==0] = 1.0
	
	# Add relaxation weighting if echo times had been varied during the acquisition
	if rweigh is not None:
		mrisig = mrisig*rweigh
	
	# Return MRI signals
	return mrisig
	
	
	
def _sigDinDex(acqpar,rweigh,tispar):

	# Get sequence parameters
	mybvals = acqpar[0,:]           # bvalues in s/mm2

	# Get tissue parameters
	lval = tispar[0]                # cell size in um2/ms
	d0val = tispar[1]               # intrinsic cell diffusivity in um2/ms
	fval = tispar[2]                # intra-cellular signal fraction
	adcexinf = tispar[3]            # extra-cellular diffusivity in um2/ms
	
	# Synthesise MRI signal
	dex = adcexinf
	mrisig = fval*_sigscylperpGPA(acqpar,np.array([lval,d0val])) +  ( 1 - fval )*np.exp(-(mybvals/1000.0)*dex)
	mrisig[mybvals==0] = 1.0
	
	# Add relaxation weighting if echo times had been varied during the acquisition
	if rweigh is not None:
		mrisig = mrisig*rweigh
	
	# Return MRI signals
	return mrisig
	
	

def _sigDin(acqpar,rweigh,tispar):

	# Get sequence parameters
	mybvals = acqpar[0,:]           # bvalues in s/mm2

	# Get tissue parameters
	lval = tispar[0]                # cell size in um2/ms
	d0val = tispar[1]               # intrinsic cell diffusivity in um2/ms
	fval = tispar[2]                # intra-cellular signal fraction
	
	# Synthesise MRI signal
	mrisig = fval*_sigscylperpGPA(acqpar,np.array([lval,d0val]))
	mrisig[mybvals==0] = 1.0
	
	# Add relaxation weighting if echo times had been varied during the acquisition
	if rweigh is not None:
		mrisig = mrisig*rweigh
	
	# Return MRI signals
	return mrisig
	
	

def _sig(acqpar,rweigh,tispar,modelid):

	# Get the right MRI signal prediction according to the specified model
	if(modelid=='DinDexTD'):
		mrs = _sigDinDexTD(acqpar,rweigh,tispar)
	elif(modelid=='DinDex'):
		mrs = _sigDinDex(acqpar,rweigh,tispar)
	elif(modelid=='Din'):
		mrs = _sigDin(acqpar,rweigh,tispar)
	else:
		raise RuntimeError('ERROR: signal model {} is unknown.'.format(modelid))
	
	# Return predicted MRI signals	
	return mrs


	
def _fobj(tisp,meas,acqp,rw,nf,modelstr):

	# Synthesise MRI signal
	measpred = _sig(acqp,rw,tisp,modelstr)
	
	# Calculate objective function value (MSE with offset noise floor model)
	fout = np.nansum( ( meas - np.sqrt(measpred**2 + nf**2) )**2 )

	# Return objective function
	return fout



def _procslice(inlist):

	# Input information for current MRI slice
	tevals = inlist[0]
	tevals_intra = inlist[1]
	bvals_intra = inlist[2]
	meas_b0 = inlist[3]
	navg = inlist[4]
	dict_pars = inlist[5]
	dict_sigs = inlist[6]
	k_data = inlist[7]
	m_data = inlist[8]
	m_data_intra = inlist[9]
	f_data = inlist[10]
	t2v_exist = inlist[11]
	t2t_data = inlist[12]
	t2v_data = inlist[13]
	sgm_exist = inlist[14] 
	sgm_data = inlist[15]
	slidx = inlist[16]
	dofit = inlist[17]
	gsep_intra = inlist[18]
	gdur_intra = inlist[19]
	tdiff_max = inlist[20]
	tdiff_min = inlist[21]
	prngmin = inlist[22]
	prngmax = inlist[23]
	sigmodel = inlist[24]
	optalgo = inlist[25]
		
	# Allocate output maps
	Nii = m_data.shape[0]
	Njj = m_data.shape[1]
	ntissue = prngmin.size   # Number of tissue parameters to estimate depending on the model
	pars_slice = np.zeros((Nii,Njj,ntissue))
	Exitmap_slice = np.zeros((Nii,Njj))      # Exit code will be 0 in the background, 1 where fitting was successful, -1 where fitting failed
	Fobj_slice = np.zeros((Nii,Njj))
	LogL_slice = np.zeros((Nii,Njj))
	BIC_slice = np.zeros((Nii,Njj))
	AIC_slice = np.zeros((Nii,Njj))
	
	# Loop through voxels
	for ii in range(0,Nii):
		for jj in range(0,Njj):
				
			## Process voxels within binary mask (if provided)
			if(k_data[ii,jj]==1):
				
				# Initialise exit code to 1 (success)
				Exflag = 1
				
				# Get actual MRI measurements
				mvox = m_data[ii,jj,:]                # MRI measurements
				mvox_intra = m_data_intra[ii,jj,:]    # MRI measurements with negligible vascular component
				fv = f_data[ii,jj]                    # Get vascular fraction
				
				# Estimate proton density s0
				if(t2v_exist):
					t2v = t2v_data[ii,jj]
					t2t = t2t_data[ii,jj]
					s0 = np.mean(mvox[meas_b0]/( fv*np.exp(-tevals[meas_b0]/t2v) + (1.0-fv)*np.exp(-tevals[meas_b0]/t2t) ))
				else:
					s0 = np.mean(mvox[meas_b0])
										
				
				# Normalise measurements
				mvox_intra_norm = mvox_intra/s0
				if(t2v_exist):
					mvox_intra_norm_relax = mvox_intra_norm/( (1.0-fv)*np.exp(-tevals_intra/t2t) )
					mvox_intra_norm_relax[mvox_intra_norm_relax<=0.0] = 1e-3
				else:
					mvox_intra_norm_relax = np.copy(mvox_intra_norm)
					mvox_intra_norm_relax[mvox_intra_norm_relax<=0.0] = 1e-3
				
				# Estimate tissue ADCc, to filter out some implausible candidate cytosolic diffusivities D0
				Nmeas = bvals_intra.size          # Number of measurements
				Bcol = np.reshape(bvals_intra/1000.0,(Nmeas,1))  # b-values as column array expressed in ms/um2
				Bmat = np.concatenate((1.0 + 0.0*Bcol,-1.0*Bcol),axis=1)    # Sequence parameter matrix for matrix multiplication
				Acol = np.reshape(np.log(mvox_intra_norm_relax),(Nmeas,1))  # Log-signals arranged as a column for matrix multiplication
				params = np.matmul( np.linalg.pinv( Bmat ) , Acol )         # Linear parameter estimation
				ADCtissue = params[1]                                        # Tissue ADC, to be used as lower bound for Din estimation
				           				
				# Make sure intrinsic intra-cellular diffusivity is at least as large as tissue ADC
				if( (sigmodel=='DinDexTD') or (sigmodel=='DinDex') or (sigmodel=='Din') ):   # For these models, D0 is in position 1
					unphys_idx =  dict_pars[:,1] < ADCtissue        # Indices of unphysical values
					dict_pars_phys = dict_pars[~unphys_idx,:]      # Get rid on unphysical D0 (lower than tissue ADC)
					dict_sigs_phys = dict_sigs[~unphys_idx,:]   # Get rid of corresponding synthetic MRI signals
					Nmicro_phys = dict_pars_phys.shape[0]      # Final number of potential candidate solutions
					mvox_intra_norm_mat = np.tile(mvox_intra_norm,(Nmicro_phys,1))    # Matrix of actual MRI measurements
				else:
					raise RuntimeError('ERROR: signal model {} is unknown.'.format(sigmodel))
				
				# Add relaxation weighting to the signal dictionary, if required
				if(t2v_exist):
					dict_sigs_phys.shape[0]
					relax_weight = (1.0-fv)*np.exp(-tevals_intra/t2t)
					relax_weight_mat = np.tile(relax_weight,(dict_sigs_phys.shape[0],1))
					dict_sigs_phys = dict_sigs_phys*relax_weight_mat
				else:
					relax_weight = None
						
				# Get noise level if it was provided
				if(sgm_exist):
				
					# Get noise standard deviation
					sgm_vox_raw = sgm_data[ii,jj]
					
					# Normalise noise standard deviation
					sgm_vox = sgm_vox_raw/s0
					
					
				## Grid search: fitting based on a discrete dictionary of signals	
				if(sgm_exist):   # Check whether noise floor modelling is required
					if( ~np.isnan(sgm_vox_raw) and ~np.isinf(sgm_vox_raw) ):
						nfloor = sgm_vox*np.sqrt(0.5*np.pi)*np.sqrt(float(navg))
					else:
						nfloor = 0.0		
					mse_array = np.nansum( ( mvox_intra_norm_mat - np.sqrt(dict_sigs_phys**2 + nfloor**2) )**2 , axis=1)
					Nmri = mvox_intra_norm_mat.shape[1]    # Number of measurements for log-likelihood computation
				else:
					nfloor = 0.0
					mse_array = np.nansum( (mvox_intra_norm_mat - dict_sigs_phys)**2 , axis=1)
						
				# Estimate tissue parameters via dictionary fitting
				try:
					# Get the position within the dictionary of the signal that most resembles our measurements
					min_idx = np.argmin(mse_array)
						
					# Get microstructural parameters corresponding to the selected signal 
					parest = dict_pars_phys[min_idx,:]
					
					# Get corresponding value of the objective function
					acpmat = np.zeros((3,bvals_intra.size))
					acpmat[0,:] = bvals_intra
					acpmat[1,:] = gdur_intra
					acpmat[2,:] = gsep_intra
					Fobj_grid = _fobj(parest,mvox_intra_norm,acpmat,relax_weight,nfloor,sigmodel)
						
					# Also perform objective function minimisation if non-linear fitting is required
					if(dofit):
						
						# Prepare data for function minimisation
						param_bound = []
						for yy in range(0,ntissue):
							param_bound.append((prngmin[yy],prngmax[yy]))
						param_init = np.copy(parest)
						
						
						# Minimise objective function
						modelfit = minimize(_fobj, param_init, method=optalgo, args=tuple([mvox_intra_norm,acpmat,relax_weight,nfloor,sigmodel]), bounds=param_bound)
						
						# Store estimated tissue parameters and value of the minimised objective function
						del parest
						parest = modelfit.x
						Fobj_new = modelfit.fun
						if(Fobj_new>Fobj_grid):
							Exflag = -1
						Fobj = 1.0*Fobj_new	
						
					else:
						Fobj = 1.0*Fobj_grid
						
					# Compute log-likelihood, BIC and AIC if an estimate for the noise standard deviation was provided
					if(sgm_exist):
						Nukn = float(ntissue)
						LogL = (-0.5/(sgm_vox*sgm_vox))*Fobj - 0.5*Nmri*np.log( np.sqrt(2*np.pi*sgm_vox*sgm_vox) )   # Log-likelihood
						BIC = -2.0*LogL + Nukn*np.log(Nmri)      # Bayesian information criterion
						AIC = -2.0*LogL + 2.0*Nukn               # Akaike information criterion
					else:
						LogL = np.nan
						BIC = np.nan
						AIC = np.nan

					# Check if values make sense - check if any NaNs were obtained in the output parameter maps
					if( np.isnan(np.sum(parest))  ):
						parest = np.nan*parest
						Fobj = np.nan
						LogL = np.nan
						BIC = np.nan
						AIC = np.nan
						Exflag = -1
							
				except:
					parest = np.nan*np.zeros(ntissue)
					Fobj = np.nan
					LogL = np.nan
					BIC = np.nan
					AIC = np.nan
					Exflag = -1
					

				# Store microstructural parameters for output
				for yy in range(0,ntissue):
					pars_slice[ii,jj,yy] = parest[yy]
				Exitmap_slice[ii,jj] = Exflag
				Fobj_slice[ii,jj] = Fobj
				LogL_slice[ii,jj] = LogL
				BIC_slice[ii,jj] = BIC
				AIC_slice[ii,jj] = AIC

	# Prepare output list and return it
	outlist = [pars_slice,Exitmap_slice,Fobj_slice,LogL_slice,BIC_slice,AIC_slice,slidx]
	return outlist
	




def run(mrifile, mriseq, output, maskfile=None, fvfile=None, T2Vfile=None, T2Tfile=None, noisefile=None, navg=1, bth=0.0, bmin=5.0, Nword=10, pmin=None, pmax=None, nthread=1, slicedim=2, nlinfit=True, nlinalgo='trust-constr', modstr='DinDex'):
	''' This tool computes 2- or 3-compartment model parameters from sets of diffusion MRI mesurements and pre-computed 
	    vascular tissue properties, using analytical signal models and maximum likelihood estimation.
	    
	    Third-party dependencies: nibabel, numpy, scipy. 
	    Developed and validated with versions: nibabel 3.2.1, numpy 1.21.5, scipy 1.7.3.
	    
	    Author: Francesco Grussu, Vall d Hebron Institute of Oncology (VHIO). 
	    Email: <francegrussu@gmail.com> <fgrussu@vhio.net>.
	    
	    USAGE
	    run(mrifile, mriseq, output, maskfile=None, fvfile=None, T2Vfile=None, T2Tfile=None, noisefile=None, navg=1, ...
	        bth=100.0, bmin=5.0, Nword=10, pmin=None, pmax=None, nthread=1, slicedim=2, ... 
	        nlinfit=True, nlinalgo='trust-constr', modstr='DinDex')
	    
	    * mrifile:    path of a 4D NIFTI file storing M diffusion MRI measurements acquired at multiple b-values 
	                  and/or multiple diffusion times and, potentially, multiple echo times TE
	                  
	    * mriseq:     path of a text file storing information on b-values, diffusion times and, optionally, echo times TE 
	                  corresponding to each volume of s_file. The acquisition must be a standard pulsed gradient spin echo
	                  (PGSE, also known as single linear diffusion encoding).
	                  This file must contain a space-separated array of 3 x M elements if echo
	                  time is not varied during the acquisition, or 4 x M is varied during the acquisition. 
	                  If the echo time is varied, the text file will have 3 lines (first line: b-values in s/mm2; 
	                  second line: gradient duration delta in ms; third line: gradient separation Delta in ms).
	                  If the echo time is varied, the text file will have 4 lines (first line: b-values in s/mm2; 
	                  second line: gradient duration delta in ms; third line: gradient separation Delta in ms; 
	                  fourth line: echo time TE in ms)
	                  
	    * output:     root file name of output files; output NIFTIs will be stored as double-precision floating point images 
	                  (FLOAT64), and the file names will end in *_Lum.nii (cell size L in um), 
	                  *_D0um2ms-1.nii (intrinsic intra-cellular diffusivity D0 in um2/ms), 
	                  *_Dexinfum2ms-1.nii (extra-cellular apparent diffusion coefficient parameter Dexinf in um2/ms), 
	                  *_Betaum2.nii (extra-cellular apparent diffusion coefficient parameter Beta -- note that the 
	                  extra-cellular apparent diffusion coefficient Dex is written as Dex = Dexinf + Beta/t, where t is 
	                  the gradient separation of measurements kept after b-value thresholding (see input parameter bth), 
	                  *_fin.nii (intra-cellular tissue fin fraction), 
	                  *_finx1-fv.nii (intra-cellular voxel Fin fraction), 
	                  *_fv.nii (vascular voxel fraction Fv - 0 everywhere if no vascular fraction map was provided with 
	                  input parameter fvfile, otherwise just a copy of the input provided with fvfile),
	                  *_cellsmm-2.nii (2D histology-like cellularity C in cells/mm2), 
	                  *_cellsmm-3.nii (cellularity C in cells/mm3), 
	                  *_exit.nii (voxel-wise exit code; -1: warning, failure in non-linear fitting; 0: background; 
	                  1 successful parameter estimation). 
	                  
	                  If a noise map was provided with the noisefile input parameter, additional output NIFTI files 
	                  storing quality of fit metrics are stored, i.e.: *_logL.nii (log-likelihood), 
	                  *_BIC.nii (Bayesian Information Criterion), and *_AIC.nii (Akaike Information Criterion). 
	                  
	                  If vascular/tissue T2 maps were provided with the T2Vfile and T2Tfile input parameters, copies of 
	                  those will also be stored as *_t2v.nii (vascular T2) and *_t2t.nii (tissue T2) so that all maps 
	                  are delivered together. 
	                  
	                  The number of parametric maps outputted depends on the model specified with input parameter modstr
	                  (see below). These will be: 
	                  L, D0, fin, Fin, Fv for model "Din"; 
	                  L, D0, fin, Dexinf, Fin, Fv for model "DinDex"; 
	                  L, D0, fin, Dexinf, Beta, Fin, Fv for model "DinDexTD"
	                  
	    * maskfile:   3D mask in NIFTI format (computation will be performed only in voxels where mask = 1)

	    * fvfile:     path of a 3D NIFTI file storing the vascular signal fraction map F, so that the tissue fraction 
	                  is equal to 1 - F (e.g. in liver MRI, F can be reasonably approximated using a IVIM model that 
	                  includes two water pools: vascular water and tissue water). If no file is provided, it is assumed
	                  that F = 0 throught the field-of-view and the modelling is effectively based on two compartments
	                  (2-compartment intra-cellular + extra-cellular model) 
	    
	    * T2Vfile:    3D vascular T2 map in NIFTI format, expressed in ms. If provided, measurements will be corrected 
	                  for T2-weighting before comparing MRI measurements and synthetic signals. 
	                  Compulsory if input T2Tfile is passed
	                  
	    * T2Tfile:    3D tissue T2 map in NIFTI format, expressed in ms, which is used as a proxy for intra-cellular T2 
	                  (this code hypothesis that intra-cellular and extra-cellular, extra-vascular T2 are approximately 
	                  the same). If provided, measurements will be corrected for T2-weighting before comparing MRI 
	                  measurements and synthetic signals. Compulsory if input T2Vfile is passed. If provided, the input text
	                  file scheme_file must contain 2 lines (with the second line expressing echo times TE in ms)
	                  
	    * noisefile:  3D noise standard deviation map in NIFTI format. If provided, the signal level will be compared to the
	                  expected Rician noise floor. Voxels where the mean intra-cellular signal overlaps considerably with the 
	                  noise floor will be flagged as -4 in the output exit code map (*_exit.nii)
	                  
	    * navg:       number of signal averages used for acquisition, which is used to get a more accurate estimate of the
	                  noise floor given the estimate of the noise standard deviation (default: 1; ignored if noisefile = None)
	                  Note that in some vendors, this parameter is referred to as number of excitations (NEX)
	    
	    * bth:        b-value threshold in s/mm2, so that b-values smaller than bth will not be used for model fitting 
	                 (default: 0 s/mm2, i.e., use all b-values). To fit the model on signals where vascular IVIM 
	                 contributions have been suppressed, set at least bth=100.0 (or higher)
	    
	    * bmin:       b-value in s/mm2 below which measurements will be used to estimate proton density (default: 5 s/mm2)
	                  
	    * Nword:      number of values to test for each tissue parameter grid search (default: 10)
	    
	    * pmin:       list or array storing the lower bounds for tissue parameters. These are: 
	                  L,D0,F for model "Din"; 
	                  L,D0,F,Dexinf for model "DinDex"; 
	                  L,D0,F,Dexinf,Beta for model "DinDexTD". 
	                  
	                  The symbols stand for: 
	                  L (cell size, in um), 
	                  D0 (intrinsic intra-cell diffusivity, um2/ms), 
	                  intra-cellular volume fraction F, 
	                  long-time extra-cellular apparent diffusion coefficient Dexinf (um2/ms), 
	                  extra-cellular apparent diffusion coefficient parameter Beta (um2). 
	                  Note that the extra-cellular apparent diffusion coeffficient is written as 
	                  Dex = Dexinf + Beta/t, where t is the diffusion time (~ gradient separation) 
	                  of measurements kept after b-value thresholding (see input parameter bth). 
	                  
	                  Default: 
	                  "8.0,0.8,0.0" for model "Din"; 
	                  "8.0,0.8,0.0,0.0" for model "DinDex"; 
	                  "8.0,0.8,0.0,0.0,0.0" for model "DinDexTD"
	                  
	                  For more information on the models, please look at input parameter modstr below.
	    
	    * pmax:       list or array storing the upper bounds for tissue parameters. These are: 
	                  L,D0,F for model "Din"; 
	                  L,D0,F,Dexinf for model "DinDex"; 
	                  L,D0,F,Dexinf,Beta for model "DinDexTD". 
	                  
	                  The symbols stand for: 
	                  L (cell size, in um), 
	                  D0 (intrinsic intra-cell diffusivity, um2/ms), 
	                  intra-cellular volume fraction F, 
	                  long-time extra-cellular apparent diffusion coefficient Dexinf (um2/ms), 
	                  extra-cellular apparent diffusion coefficient parameter Beta (um2). 
	                  Note that the extra-cellular apparent diffusion coeffficient is written as 
	                  Dex = Dexinf + Beta/t, where t is the diffusion time (~ gradient separation) 
	                  of measurements kept after b-value thresholding (see input parameter bth). 
	                  
	                  Default: 
	                  "40.0,3.0,1.0" for model "Din"; 
	                  "40.0,3.0,1.0,3.0" for model "DinDex"; 
	                  "40.0,3.0,1.0,3.0,10.0" for model "DinDexTD"
	                  
	                  For more information on the models, please look at input parameter modstr below.            
	    
	    * nthread:    number of threads to be used for computation (default: 1, single thread)
	    
	    * slicedim:   image dimension along which parallel computation will be exectued when nthread > 1
	                  (can be 0, 1, 2 -- default slicedim=2, implying parallel processing along 3rd 
	                  image dimension)
	                  
	    * nlinfit:    perform non-linear fitting via gradient descent (constrained objective function minimisation) if True, 
	                  perform only grid search if False (default: True)
	    
	    * nlinalgo:   algorithm to be used for constrained objective function minimisation in non-linear fitting (relevant if nlinfit = True). 
	                  Choose among: "Nelder-Mead", "L-BFGS-B", "TNC", "SLSQP", "Powell", and "trust-constr" 
	                  (default: "trust-constr" - see documentation of scipy.optimize.minimize for information on the optimisation algorithm)
	    
	    * modstr:     string specifying the signal model to fit. Choose among: 
	                  "Din" (extra-vascular signal dominated by intra-cellular diffusion), 
	                  "DinDex" (extra-vascular signal features both intra-cellular and extra-cellular contributions, without diffusion time 
	                           dependence in the extra-cellular ADC), 
	                  "DinDexTD" (extra-vascular signal features both intra-cellular and extra-cellular contributions, with diffusion time 
	                              dependence in the extra-cellular ADC). 
	                  
	                  Default: "DinDex". 
	                  
	                  Intra-cellular diffusion is modelled using the Gaussian Phase Approximation (GPA) formula for diffusion within sphere 
	                  of Van Geldern et al, "Evaluation of restricted diffusion in cylinders. Phosphocreatine in rabbit leg muscle"
	                  J Magn Reson B 1994, 103(3):255-60, doi: 10.1006/jmrb.1994.1038          
	    	    
	    Third-party dependencies: nibabel, numpy, scipy.
	    
	    Developed and validated with versions: 
	    - nibabel 3.2.1
	    - numpy 1.21.5
	    - scipy 1.7.3

	    Author: Francesco Grussu, Vall d'Hebron Institute of Oncology, November 2022
		    <fgrussu@vhio.net> <francegrussu@gmail.com>'''


	### Get time
	timeinitial = time.time()
	
	### Get number of threads
	nthread = int(nthread)
	ncpu = mp.cpu_count()
	if( (ncpu - 1) < nthread):
		nthread = ncpu - 1

	if(nthread<0):
		nthread = ncpu - 1
		print('WARNING: negative number of threads -- using {} instead'.format(nthread))

	## Get slice dimension for parallel processing
	if( not( (slicedim==0) or (slicedim==1) or (slicedim==2) ) ):
		slicedim = 2
		print('WARNING: invalid image dimension for parallel processing -- using dimension={} ({}-th dimension)'.format(slicedim,slicedim+1))		
	
	### Load MRI measurements and check for consistency
	print('')
	print('    ... loading MRI measurements')
	try:
		m_obj = nib.load(mrifile)
	except:
		print('')
		raise RuntimeError('ERROR: the 4D input file {} does not exist or is not in NIFTI format.'.format(mrifile))
	m_data = m_obj.get_fdata()
	m_size = m_data.shape
	m_size = np.array(m_size)
	m_header = m_obj.header
	m_affine = m_header.get_best_affine()
	m_dims = m_obj.shape
	if m_size.size!=4:
		print('')
		raise RuntimeError('ERROR: the 4D input file {} is not a 4D NIFTI.'.format(mrifile))				 

	### Load mask and check for consistency
	if (maskfile is not None):
		print('')
		print('    ... loading mask')
		try:
			k_obj = nib.load(maskfile)
		except:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} does not exist or is not in NIFTI format.'.format(maskfile))
		k_data = k_obj.get_fdata()
		k_size = k_data.shape
		k_size = np.array(k_size)
		k_header = k_obj.header
		k_affine = k_header.get_best_affine()
		k_dims = k_obj.shape
		if k_size.size!=3:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} is not a 3D NIFTI.'.format(maskfile))
		if ( (np.sum(k_affine==m_affine)!=16) or (k_dims[0]!=m_dims[0]) or (k_dims[1]!=m_dims[1]) or (k_dims[2]!=m_dims[2]) ):
			print('')
			raise RuntimeError('ERROR: the header geometry of {} and {} do not match.'.format(mrifile,maskfile))
		k_data[k_data>0] = 1
		k_data[k_data<=0] = 0	
	else:
		k_data = np.ones((m_size[0],m_size[1],m_size[2]))

	### Load vascular fraction and check for consistency
	if(fvfile is not None):
		print('')
		print('    ... loading vascular volume fraction map')
		try:
			f_obj = nib.load(fvfile)
		except:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} does not exist or is not in NIFTI format.'.format(fvfile))
		f_data = f_obj.get_fdata()
		f_size = f_data.shape
		f_size = np.array(f_size)
		f_header = f_obj.header
		f_affine = f_header.get_best_affine()
		f_dims = f_obj.shape
		if f_size.size!=3:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} is not a 3D NIFTI.'.format(fvfile))
		if ( (np.sum(f_affine==m_affine)!=16) or (f_dims[0]!=m_dims[0]) or (f_dims[1]!=m_dims[1]) or (f_dims[2]!=m_dims[2]) ):
			print('')
			raise RuntimeError('ERROR: the header geometry of {} and {} do not match.'.format(mrifile,fvfile))
	else:
		f_data = np.zeros((m_size[0],m_size[1],m_size[2]))
		f_size = f_data.shape
		
	### Load vascular and tissue T2, and check for consistency
	if ( (T2Vfile is not None) or (T2Tfile is not None)):
		print('')
		print('    ... loading vascular and tissue T2 maps')
		
		## Vascular T2 map
		try:
			t2v_obj = nib.load(T2Vfile)
		except:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} does not exist or is not in NIFTI format.'.format(T2Vfile))
		t2v_data = t2v_obj.get_fdata()
		t2v_size = t2v_data.shape
		t2v_size = np.array(t2v_size)
		t2v_header = t2v_obj.header
		t2v_affine = t2v_header.get_best_affine()
		t2v_dims = t2v_obj.shape
		if t2v_size.size!=3:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} is not a 3D NIFTI.'.format(T2Vfile))
		if ( (np.sum(t2v_affine==m_affine)!=16) or (t2v_dims[0]!=m_dims[0]) or (t2v_dims[1]!=m_dims[1]) or (t2v_dims[2]!=m_dims[2]) ):
			print('')
			raise RuntimeError('ERROR: the header geometry of {} and {} do not match.'.format(mrifile,T2Vfile))	
		
		## Tissue T2 map
		try:
			t2t_obj = nib.load(T2Tfile)
		except:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} does not exist or is not in NIFTI format.'.format(T2Tfile))
		t2t_data = t2t_obj.get_fdata()
		t2t_size = t2t_data.shape
		t2t_size = np.array(t2t_size)
		t2t_header = t2t_obj.header
		t2t_affine = t2t_header.get_best_affine()
		t2t_dims = t2t_obj.shape
		if t2t_size.size!=3:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} is not a 3D NIFTI.'.format(T2Tfile))
		if ( (np.sum(t2t_affine==m_affine)!=16) or (t2t_dims[0]!=m_dims[0]) or (t2t_dims[1]!=m_dims[1]) or (t2t_dims[2]!=m_dims[2]) ):
			print('')
			raise RuntimeError('ERROR: the header geometry of {} and {} do not match.'.format(mrifile,T2Tfile))	
		
		# Flags
		t2v_exist = True
		t2t_exist = True
	
	else:
		# Flags
		t2v_exist = False
		t2t_exist = False
		
		# Empty variables
		t2v_data = None
		t2t_data = None
	
	### Load noise standard deviation map
	if (noisefile is not None):
		print('')
		print('    ... loading noise map')
		try:
			sgm_obj = nib.load(noisefile)
		except:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} does not exist or is not in NIFTI format.'.format(noisefile))
		sgm_data = sgm_obj.get_fdata()
		sgm_size = sgm_data.shape
		sgm_size = np.array(sgm_size)
		sgm_header = sgm_obj.header
		sgm_affine = sgm_header.get_best_affine()
		sgm_dims = sgm_obj.shape
		if sgm_size.size!=3:
			print('')
			raise RuntimeError('ERROR: the 3D input file {} is not a 3D NIFTI.'.format(noisefile))
		if ( (np.sum(sgm_affine==m_affine)!=16) or (sgm_dims[0]!=m_dims[0]) or (sgm_dims[1]!=m_dims[1]) or (sgm_dims[2]!=m_dims[2]) ):
			print('')
			raise RuntimeError('ERROR: the header geometry of {} and {} do not match.'.format(mrifile,noisefile))
		
		# Flags
		sgm_exist = True
		sgm_exist = True
			
	else:
		# Flags
		sgm_exist = False
		sgm_exist = False
		
		# Empty variable
		sgm_data = None
	
	### Load MRI sequence parameters and check for consistency
	print('')
	print('    ... loading MRI sequence information')
	mripar = np.loadtxt(mriseq)
	if(t2v_exist):
		if(mripar.shape[0]!=4):
			print('')
			raise RuntimeError('ERROR: since T2 maps were provided, {} must have size 4 x M and list b-values, grad. dur., grad. sep., and TEs'.format(mriseq))
		if(m_size[3]!=mripar.shape[1]):
			print('')
			raise RuntimeError('ERROR: the number of measurements in {} and {} do not match'.format(mrifile,mriseq))
		bvals = mripar[0,:]
		gdur = mripar[1,:]
		gsep = mripar[2,:]
		tevals = mripar[3,:]
	else:
		if(mripar.shape[0]!=3):
			print('')
			raise RuntimeError('ERROR: since no T2 maps were provided, {} must have size 3 x M and list b-values, grad. dur., grad. sep.'.format(mriseq))
		if(m_size[3]!=mripar.shape[1]):
			print('')
			raise RuntimeError('ERROR: the number of measurements in {} and {} do not match'.format(mrifile,mriseq))
		bvals = mripar[0,:]
		gdur = mripar[1,:]
		gsep = mripar[2,:]
		tevals = None
		tevals_intra = None
	
	
		
	## Load tissue parameter lower bounds and assign default values if no bounds were provided
	if pmin is not None:
		pmin = np.array(pmin)
		pmin_size = pmin.size
		if(modstr=='DinDexTD'):
			if(pmin_size!=5):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmin = {} - for model {} it should be 5.'.format(pmin,modstr))
		elif(modstr=='DinDex'):
			if(pmin_size!=4):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmin = {} - for model {} it should be 4.'.format(pmin,modstr))
		elif(modstr=='Din'):
			if(pmin_size!=3):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmin = {} - for model {} it should be 3.'.format(pmin,modstr))
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))		
	else:
		if(modstr=='DinDexTD'):
			pmin = np.array([8.0,0.8,0.0,0.0,0.0])
		elif(modstr=='DinDex'):
			pmin = np.array([8.0,0.8,0.0,0.0])
		elif(modstr=='Din'):
			pmin = np.array([8.0,0.8,0.0])
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))
	
	
	## Load tissue parameter upper bounds and assign default values if no bounds were provided
	if pmax is not None:
		pmax = np.array(pmax)
		pmax_size = pmax.size
		if(modstr=='DinDexTD'):
			if(pmax_size!=5):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmax = {} - for model {} it should be 5.'.format(pmax,modstr))
		elif(modstr=='DinDex'):
			if(pmax_size!=4):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmax = {} - for model {} it should be 4.'.format(pmax,modstr))
		elif(modstr=='Din'):
			if(pmax_size!=3):
				raise RuntimeError('ERROR: wrong number of tissue parameters in lower bounds pmax = {} - for model {} it should be 3.'.format(pmax,modstr))
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))
	else:
		if(modstr=='DinDexTD'):
			pmax = np.array([40.0,3.0,1.0,3.0,10.0])
		elif(modstr=='DinDex'):
			pmax = np.array([40.0,3.0,1.0,3.0])
		elif(modstr=='Din'):
			pmax = np.array([40.0,3.0,1.0])
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))
	
	
	## Extract MRI measurement subset after thresholding b-values
	print('')
	print('    ... thresholding b-values')
	meas_keep = bvals>=bth
	if(np.sum(meas_keep)==0):
		print('')
		raise RuntimeError('ERROR: no MRI measurements are left after thresholding b-values. Try a lower b-value threshold.'.format(sigfile,mrifile))
	gsep_thr = gsep[meas_keep]
	gdur_thr = gdur[meas_keep]
	tdiff_min = np.min(gsep_thr[gsep_thr>0])   # Minimum diffusion time
	tdiff_max = np.max(gsep_thr[gsep_thr>0])   # Maximum diffusion time
	if(tdiff_max==tdiff_min):
		print('')
		raise RuntimeError('ERROR: your acquisition features a unique gradient separation value, equal to {} ms. The signal model cannot be fitted'.format(tdiff_max))	
		
	### Create dictionary of synthetic parameters and corresponding MRI signals
	print('')
	print('    ... creating dictionary of synthetic MRI signals')
	
	# Dictionary of tissue parameters
	if(modstr=='DinDexTD'):	
		L_list = np.linspace(pmin[0],pmax[0],Nword)
		Din_list = np.linspace(pmin[1],pmax[1],Nword)
		Fc_list = np.linspace(pmin[2],pmax[2],Nword)
		Dexref_list = np.linspace(pmin[3],pmax[3],Nword)
		Bex_list = np.linspace(pmin[4],pmax[4],Nword)
		L_array, Din_array, Fc_array, Dexref_array, Bex_array = np.meshgrid(L_list,Din_list,Fc_list,Dexref_list,Bex_list)
		L_array = L_array.flatten()
		Din_array = Din_array.flatten()
		Fc_array = Fc_array.flatten()
		Dexref_array = Dexref_array.flatten()
		Bex_array = Bex_array.flatten()
		Nmicro = L_array.size
		dict_pars = np.zeros((Nmicro,5))
		dict_pars[:,0] = L_array
		dict_pars[:,1] = Din_array
		dict_pars[:,2] = Fc_array
		dict_pars[:,3] = Dexref_array
		dict_pars[:,4] = Bex_array
		
	elif(modstr=='DinDex'):
		L_list = np.linspace(pmin[0],pmax[0],Nword)
		Din_list = np.linspace(pmin[1],pmax[1],Nword)
		Fc_list = np.linspace(pmin[2],pmax[2],Nword)
		Dexref_list = np.linspace(pmin[3],pmax[3],Nword)
		L_array, Din_array, Fc_array, Dexref_array = np.meshgrid(L_list,Din_list,Fc_list,Dexref_list)
		L_array = L_array.flatten()
		Din_array = Din_array.flatten()
		Fc_array = Fc_array.flatten()
		Dexref_array = Dexref_array.flatten()
		Nmicro = L_array.size
		dict_pars = np.zeros((Nmicro,4))
		dict_pars[:,0] = L_array
		dict_pars[:,1] = Din_array
		dict_pars[:,2] = Fc_array
		dict_pars[:,3] = Dexref_array
	
	elif(modstr=='Din'):
		L_list = np.linspace(pmin[0],pmax[0],Nword)
		Din_list = np.linspace(pmin[1],pmax[1],Nword)
		Fc_list = np.linspace(pmin[2],pmax[2],Nword)
		L_array, Din_array, Fc_array = np.meshgrid(L_list,Din_list,Fc_list)
		L_array = L_array.flatten()
		Din_array = Din_array.flatten()
		Fc_array = Fc_array.flatten()
		Nmicro = L_array.size
		dict_pars = np.zeros((Nmicro,3))
		dict_pars[:,0] = L_array
		dict_pars[:,1] = Din_array
		dict_pars[:,2] = Fc_array
	
	else:
		raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))
	
	# Dictionary of MRI signals
	npars = dict_pars.shape[1]	
	dict_sigs = np.zeros((Nmicro,mripar.shape[1]))	
	for qq in range(0,Nmicro):      # Loop through different microstructures in the dictionary
		# Synthesise signals without relaxation-weighting and store them
		synsigs = _sig(mripar,None,dict_pars[qq,:],modstr)
		dict_sigs[qq,:] = synsigs
	print('        ({} synthetic signals generated)'.format(Nmicro))


	# Remove measurements from the synthetic dictionary corresponding to b-values that are removed after thresholding
	m_data_intra = m_data[:,:,:,meas_keep]
	dict_sigs_intra = dict_sigs[:,meas_keep]
	if(t2v_exist):
		tevals_intra = tevals[meas_keep]
	bvals_intra = bvals[meas_keep]
	gsep_intra = gsep[meas_keep]
	gdur_intra = gdur[meas_keep]
	
	# Find out whether b-values where the vascular signal has been suppressed are good enough for model fitting (two non-zero b-values are needed)
	if(np.unique(bvals_intra).size==1):
		print('')
		raise RuntimeError('ERROR: you need at least two non-zero b-values where the vascular signal has been suppressed. The analysis cannot be run.')
	
	
	## Find indices of measurements to be used for proton density estimation
	meas_b0 = bvals<bmin
	
	### Allocate output parametric maps
	Tparmap = np.zeros((f_size[0],f_size[1],f_size[2],npars))     # Allocate output: parametric maps to be estimated
	Exitmap = np.zeros(f_size)    # Allocate output:  exit code map
	Fobjmap = np.zeros(f_size)    # Allolcate output: Fobj map
	LogLmap = np.zeros(f_size)    # Allolcate output: log-likelihood map
	BICmap = np.zeros(f_size)     # Allolcate output: Bayesian Information Criterion map
	AICmap = np.zeros(f_size)     # Allolcate output: Akaike Information Criterion map
	
	### Processing
	print('')
	print('    ... processing -- please wait')
	
	# Prepare information for current MRI slice
	inputlist = [] 
	for ww in range(0,m_size[slicedim]):

		if(slicedim==0):
		
			k_data_sl = k_data[ww,:,:]
			m_data_sl = m_data[ww,:,:,:]
			m_data_intra_sl = m_data_intra[ww,:,:,:]
			f_data_sl = f_data[ww,:,:]
			if(t2v_exist):
				t2t_data_sl = t2t_data[ww,:,:]
				t2v_data_sl = t2v_data[ww,:,:]
			else:
				t2t_data_sl = None
				t2v_data_sl = None
			if(sgm_exist):
				sgm_data_sl = sgm_data[ww,:,:]
			else:
				sgm_data_sl = None	
		
		elif(slicedim==1):
		
			k_data_sl = k_data[:,ww,:]
			m_data_sl = m_data[:,ww,:,:]
			m_data_intra_sl = m_data_intra[:,ww,:,:]
			f_data_sl = f_data[:,ww,:]
			if(t2v_exist):
				t2t_data_sl = t2t_data[:,ww,:]
				t2v_data_sl = t2v_data[:,ww,:]
			else:
				t2t_data_sl = None
				t2v_data_sl = None
			if(sgm_exist):
				sgm_data_sl = sgm_data[:,ww,:]
			else:
				sgm_data_sl = None			
		
		elif(slicedim==2):
		
			k_data_sl = k_data[:,:,ww]
			m_data_sl = m_data[:,:,ww,:]
			m_data_intra_sl = m_data_intra[:,:,ww,:]
			f_data_sl = f_data[:,:,ww]
			if(t2v_exist):
				t2t_data_sl = t2t_data[:,:,ww]
				t2v_data_sl = t2v_data[:,:,ww]
			else:
				t2t_data_sl = None
				t2v_data_sl = None
			if(sgm_exist):
				sgm_data_sl = sgm_data[:,:,ww]
			else:
				sgm_data_sl = None		
			
		else:
			raise RuntimeError('ERROR: invalid slice dimension slicedim = {}'.format(slicedim)) 

		
		sliceinfo = [tevals,tevals_intra,bvals_intra,meas_b0,navg,dict_pars,dict_sigs_intra,k_data_sl,m_data_sl,m_data_intra_sl,f_data_sl,t2v_exist,t2t_data_sl,t2v_data_sl,sgm_exist,sgm_data_sl,ww,nlinfit,gsep_intra,gdur_intra,tdiff_max,tdiff_min,pmin,pmax,modstr,nlinalgo]   
		inputlist.append(sliceinfo)	
	
	# Send slice to process in parallel if nthread>1	
	if(nthread>1):
		
		# Create the parallel pool and give jobs to the workers
		fitpool = mp.Pool(processes=nthread)  # Create parallel processes
		fitpool_pids_initial = [proc.pid for proc in fitpool._pool]  # Get PIDs
		fitresults = fitpool.map_async(_procslice,inputlist)          # Send off to parallel processing
			
		# Check whether results are ready
		while not fitresults.ready():
			fitpool_pids_new = [proc.pid for proc in fitpool._pool]  # Get PIDs again
			if fitpool_pids_new!=fitpool_pids_initial:
				print('')
				raise RuntimeError('ERROR: some processes died during parallel fitting.') # Error: some processes where killed, stop everything and throw an error 
			
		# Work done: get results
		fitlist = fitresults.get()

		# Collect fitting output and re-assemble MRI slices	## Recall: the output of _procslice() contains out = [pars_slice,Exitmap_slice,Fobj_slice,LogL_slice,BIC_slice,AIC_slice,slidx]	
		for ww in range(0, m_size[slicedim]):
			myresults = fitlist[ww]		
			allmaps_sl = myresults[0]	
			Exitmap_sl = myresults[1]
			Fobj_sl = myresults[2]
			LogL_sl = myresults[3]
			BIC_sl = myresults[4]
			AIC_sl = myresults[5]	
			myslice = myresults[6]		
			
			if(slicedim==0):
				for uu in range(0,npars):
					Tparmap[ww,:,:,uu] = allmaps_sl[:,:,uu]
				Exitmap[ww,:,:] = Exitmap_sl
				Fobjmap[ww,:,:] = Fobj_sl
				LogLmap[ww,:,:] = LogL_sl
				BICmap[ww,:,:] = BIC_sl
				AICmap[ww,:,:] = AIC_sl
				
			elif(slicedim==1):
				for uu in range(0,npars):
					Tparmap[:,ww,:,uu] = allmaps_sl[:,:,uu]
				Exitmap[:,ww,:] = Exitmap_sl
				Fobjmap[:,ww,:] = Fobj_sl
				LogLmap[:,ww,:] = LogL_sl
				BICmap[:,ww,:] = BIC_sl
				AICmap[:,ww,:] = AIC_sl			
		
			elif(slicedim==2):
				for uu in range(0,npars):
					Tparmap[:,:,ww,uu] = allmaps_sl[:,:,uu]
				Exitmap[:,:,ww] = Exitmap_sl
				Fobjmap[:,:,ww] = Fobj_sl
				LogLmap[:,:,ww] = LogL_sl
				BICmap[:,:,ww] = BIC_sl
				AICmap[:,:,ww] = AIC_sl
							
			else:
				raise RuntimeError('ERROR: invalid slice dimension slicedim = {}'.format(slicedim)) 		
			
		
	
	# Single CPU at work	
	else:
	
		for ww in range(0, m_size[slicedim]):
			myresults = _procslice(inputlist[ww]) 
			allmaps_sl = myresults[0]	
			Exitmap_sl = myresults[1]
			Fobj_sl = myresults[2]
			LogL_sl = myresults[3]
			BIC_sl = myresults[4]
			AIC_sl = myresults[5]	
			myslice = myresults[6]	
			
			if(slicedim==0):
				for uu in range(0,npars):
					Tparmap[ww,:,:,uu] = allmaps_sl[:,:,uu]
				Exitmap[ww,:,:] = Exitmap_sl
				Fobjmap[ww,:,:] = Fobj_sl
				LogLmap[ww,:,:] = LogL_sl
				BICmap[ww,:,:] = BIC_sl
				AICmap[ww,:,:] = AIC_sl
				
			elif(slicedim==1):
				for uu in range(0,npars):
					Tparmap[:,ww,:,uu] = allmaps_sl[:,:,uu]
				Exitmap[:,ww,:] = Exitmap_sl
				Fobjmap[:,ww,:] = Fobj_sl
				LogLmap[:,ww,:] = LogL_sl
				BICmap[:,ww,:] = BIC_sl
				AICmap[:,ww,:] = AIC_sl			
		
			elif(slicedim==2):
				for uu in range(0,npars):
					Tparmap[:,:,ww,uu] = allmaps_sl[:,:,uu]
				Exitmap[:,:,ww] = Exitmap_sl
				Fobjmap[:,:,ww] = Fobj_sl
				LogLmap[:,:,ww] = LogL_sl
				BICmap[:,:,ww] = BIC_sl
				AICmap[:,:,ww] = AIC_sl
							
			else:
				raise RuntimeError('ERROR: invalid slice dimension slicedim = {}'.format(slicedim)) 	
		
	
			
	### Save output NIFTIs				
	print('')
	print('    ... saving output files')
	buffer_header = m_obj.header
	buffer_header.set_data_dtype('float64')   # Make sure we save output files float64, even if input is not
	
	# Save exit code and objective function
	exit_obj = nib.Nifti1Image(Exitmap,m_obj.affine,buffer_header)
	nib.save(exit_obj, '{}_exit.nii'.format(output))
	
	fobj_obj = nib.Nifti1Image(Fobjmap,m_obj.affine,buffer_header)
	nib.save(fobj_obj, '{}_fobj.nii'.format(output))
	
	if(modstr=='DinDexTD'):
		map_strings = ['Lum.nii','D0um2ms-1','fin.nii','Dexinfum2ms-1','Betaum2.nii']
	elif(modstr=='DinDex'):
		map_strings = ['Lum.nii','D0um2ms-1','fin.nii','Dexinfum2ms-1']
	elif(modstr=='Din'):
		map_strings = ['Lum.nii','D0um2ms-1','fin.nii']
	else:
		raise RuntimeError('ERROR: signal model {} is unknown.'.format(modstr))
		
	# Save parametric maps
	for uu in range(0,npars):
	
		# Keep track of cell size
		if(map_strings[uu]=='Lum.nii'):
			Lmap = np.squeeze(Tparmap[:,:,:,uu])
		
		# Keep track of intra-cellular tissue fraction
		if(map_strings[uu]=='fin.nii'):
			Finmap = np.squeeze(Tparmap[:,:,:,uu])
	
		niiout_obj = nib.Nifti1Image(np.squeeze(Tparmap[:,:,:,uu]),m_obj.affine,buffer_header)
		nib.save(niiout_obj, '{}_{}'.format(output,map_strings[uu]))
	
	
	# Calculate intra-cellular voxel fraction and save it, jointly with the vascular voxel fraction map
	Finvoxmap = (1.0 - f_data[:,:,0])*Finmap
	
	finvox_obj = nib.Nifti1Image(Finvoxmap,m_obj.affine,buffer_header)
	nib.save(finvox_obj, '{}_finx1-fv.nii'.format(output))
	
	fv_obj = nib.Nifti1Image(f_data,m_obj.affine,buffer_header)
	nib.save(fv_obj, '{}_fv.nii'.format(output))
	
	# Calculate cellularity in cells/mm3 and save it
	cellmap = Finvoxmap/( (0.001*Lmap)*(0.001*Lmap)*(0.001*Lmap) )
	cellularity_obj = nib.Nifti1Image(cellmap,m_obj.affine,buffer_header)
	nib.save(cellularity_obj, '{}_cellsmm-3.nii'.format(output))

	# Calculate cellularity in cells/mm2 and save it
	cellmap2D = Finvoxmap/( (0.001*Lmap)*(0.001*Lmap) )
	cellularity2D_obj = nib.Nifti1Image(cellmap2D,m_obj.affine,buffer_header)
	nib.save(cellularity2D_obj, '{}_cellsmm-2.nii'.format(output))

		
	
	# If a noise map had been provided, we can get LogL, BIC and AIC from FObj and store them
	if(sgm_exist):
	
		logL_obj = nib.Nifti1Image(LogLmap,m_obj.affine,buffer_header)
		nib.save(logL_obj, '{}_logL.nii'.format(output))
		
		bic_obj = nib.Nifti1Image(BICmap,m_obj.affine,buffer_header)
		nib.save(bic_obj, '{}_BIC.nii'.format(output))
		
		aic_obj = nib.Nifti1Image(AICmap,m_obj.affine,buffer_header)
		nib.save(aic_obj, '{}_AIC.nii'.format(output))

	# Store copies of T2v and T2t maps if they were provided, so that all maps are nicely delivered together
	if(t2v_exist):
		t2v_obj = nib.Nifti1Image(t2v_data,m_obj.affine,buffer_header)
		nib.save(t2v_obj, '{}_t2v.nii'.format(output))
		
	if(t2t_exist):
		t2t_obj = nib.Nifti1Image(t2t_data,m_obj.affine,buffer_header)
		nib.save(t2t_obj, '{}_t2t.nii'.format(output))


	### Done
	timefinal = time.time()
	print('    ... done - it took {} sec'.format(timefinal - timeinitial))
	print('')




# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='This tool computes 2- or 3-compartment model parameters from sets of diffusion MRI mesurements and pre-computed vascular tissue properties, using analytical signal models and maximum likelihood estimation. Third-party dependencies: nibabel, numpy, scipy. Developed and validated with versions: nibabel 3.2.1, numpy 1.21.5, scipy 1.7.3. Author: Francesco Grussu, Vall d Hebron Institute of Oncology (VHIO). Email: <francegrussu@gmail.com> <fgrussu@vhio.net>.')
	parser.add_argument('s_file', help='path of a 4D NIFTI file storing M diffusion MRI measurements acquired at multiple b-values and/or multiple diffusion times and, potentially, multiple echo times TE')
	parser.add_argument('scheme_file', help='path of a text file storing information on b-values, diffusion times and, optionally, echo times TE corresponding to each volume of s_file. The acquisition must be a standard pulsed gradient spin echo (PGSE, also known as single linear diffusion encoding). This file must contain a space-separated array of 3 x M elements if echo time is not varied during the acquisition, or 4 x M is varied during the acquisition. If the echo time is varied, the text file will have 3 lines (first line: b-values in s/mm2; second line: gradient duration delta in ms; third line: gradient separation Delta in ms). If the echo time is varied, the text file will have 4 lines (first line: b-values in s/mm2; second line: gradient duration delta in ms; third line: gradient separation Delta in ms;  fourth line: echo time TE in ms).')
	parser.add_argument('out', help='root file name of output files; output NIFTIs will be stored as double-precision floating point images (FLOAT64), and the file names will end in *_Lum.nii (cell size L in um), *_D0um2ms-1.nii (intrinsic intra-cellular diffusivity D0 in um2/ms), *_Dexinfum2ms-1.nii (extra-cellular apparent diffusion coefficient parameter Dexinf in um2/ms), *_Betaum2.nii (extra-cellular apparent diffusion coefficient parameter Beta -- note that the extra-cellular apparent diffusion coefficient Dex is written as Dex = Dexinf + Beta/t, where t is the gradient separation of measurements kept after b-value thresholding (see option --bfil), *_fin.nii (intra-cellular tissue fin fraction), *_finx1-fv.nii (intra-cellular voxel Fin fraction), *_fv.nii (vascular voxel fraction Fv - 0 everywhere if no vascular fraction map was provided with option --fv, otherwise just a copy of the input provided with --fv), *_cellsmm-2.nii (2D histology-like cellularity C in cells/mm2), *_cellsmm-3.nii (cellularity C in cells/mm3), *_exit.nii (voxel-wise exit code; -1: warning, failure in non-linear fitting; 0: background; 1 successful parameter estimation). If a noise map was provided with the --noise option, additional output NIFTI files storing quality of fit metrics are stored, i.e.: *_logL.nii (log-likelihood), *_BIC.nii (Bayesian Information Criterion), and *_AIC.nii (Akaike Information Criterion). If vascular/tissue T2 maps were provided with the --t2v and --t2t options, copies of those will also be stored as *_t2v.nii (vascular T2) and *_t2t.nii (tissue T2) so that all maps are delivered together. The number of parametric maps outputted depends on the model specified with --modstr. These will be: L, D0, fin, Fin, Fv for model "Din"; L, D0, fin, Dexinf, Fin, Fv for model "DinDex"; L, D0, fin, Dexinf, Beta, Fin, Fv for model "DinDexTD"')
	parser.add_argument('--mask', metavar='<file>', help='3D mask in NIFTI format (computation will be performed only in voxels where mask = 1)')
	parser.add_argument('--fv', metavar='<file>', help='path of a 3D NIFTI file storing the vascular signal fraction map F, so that the tissue fraction is equal to 1 - F (e.g. in liver MRI, F can be reasonably approximated using an IVIM model that includes two water pools: vascular water and tissue water) -- if not specified, it is assumed that F = 0 everywhere and the modelling reduces to a 2-compartment representation of the tissue')
	parser.add_argument('--t2v', metavar='<file>', help='3D vascular T2 map in NIFTI format, expressed in ms. If provided, measurements will be corrected for T2-weighting before comparing MRI measurements and synthetic signals. Compulsory if --t2t input is provided.')
	parser.add_argument('--t2t', metavar='<file>', help='3D tissue T2 map in NIFTI format, expressed in ms, which is used as a proxy for intra-cellular T2 (this code hypothesis that intra-cellular and extra-cellular, extra-vascular T2 are approximately the same). If provided, measurements will be corrected for T2-weighting before comparing MRI measurements and synthetic signals. Compulsory if --t2v input is provided. If provided, the input text file scheme_file must contain 2 lines (with the second line expressing echo times TE in ms)')
	parser.add_argument('--noise', metavar='<file>', help='3D noise standard deviation map in NIFTI format. If provided, the signal level will be compared to the expected Rician noise floor. Voxels where the intra-cellular signal level overlaps considerably with the noise floor will be flagged as -4 in the output exit code map (*_exit.nii).')
	parser.add_argument('--savg', metavar='<num>', default='1', help='number of signal averages used for MRI data acquisition (default: 1). This parameter is used for the estimation of the noise floor (it is ignored if the option --noise is not used). Note that in some vendors, this parameter is also referred to as number of excitations (NEX).')
	parser.add_argument('--bfil', metavar='<num>', default='0.0', help='b-value threshold in s/mm2, so that b-values smaller than --bfil will not be used for model fitting (default: 0 s/mm2, i.e., use all b-values). To fit the model on signals where vascular IVIM contributions have been suppressed, set at least --bfil 100.0 (or higher).')
	parser.add_argument('--bmin', metavar='<num>', default='5.0', help='b-value in s/mm2 below which measurements will be used to estimate proton density (default: 5 s/mm2)')
	parser.add_argument('--nw', metavar='<num>', default='10', help='number of values to test for each unknown tissue parameter in the grid search (it must be an integer; default 10)')
	parser.add_argument('--pmin', metavar='<list>', help='comma-separaterd list storing the lower bounds for tissue parameters. These are: L,D0,F for model "Din"; L,D0,F,Dexinf for model "DinDex"; L,D0,F,Dexinf,Beta for model "DinDexTD". These stand for: L (cell size, in um), D0 (intrinsic intra-cell diffusivity, um2/ms), intra-cellular volume fraction F, long-time extra-cellular apparent diffusion coefficient Dexinf (um2/ms), extra-cellular apparent diffusion coefficient parameter Beta (um2). Note that the extra-cellular apparent diffusion coeffficient is written as Dex = Dexinf + Beta/t, where t is the diffusion time (~ gradient separation) of measurements kept after b-value thresholding (see option --bfil). Default: "8.0,0.8,0.0" for model "Din"; "8.0,0.8,0.0,0.0" for model "DinDex"; "8.0,0.8,0.0,0.0,0.0" for model "DinDexTD".')
	parser.add_argument('--pmax', metavar='<list>', help='comma-separaterd list storing the upper bounds for tissue parameters. These are: L,D0,F for model "Din"; L,D0,F,Dexinf for model "DinDex"; L,D0,F,Dexinf,Beta for model "DinDexTD". These stand for: L (cell size, in um), D0 (intrinsic intra-cell diffusivity, um2/ms), intra-cellular volume fraction F, long-time extra-cellular apparent diffusion coefficient Dexinf (um2/ms), extra-cellular apparent diffusion coefficient parameter Beta (um2). Note that the extra-cellular apparent diffusion coeffficient is written as Dex = Dexinf + Beta/t, where t is the diffusion time (~ gradient separation) of measurements kept after b-value thresholding (see option --bfil). Default: "40.0,3.0,1.0" for model "Din"; "40.0,3.0,1.0,3.0" for model "DinDex"; "40.0,3.0,1.0,3.0,10.0" for model "DinDexTD".')
	parser.add_argument('--ncpu', metavar='<num>', default='1', help='number of threads to be used for computation (default: 1, single thread)')
	parser.add_argument('--sldim', metavar='<num>', default='2', help='image dimension along which parallel computation will be exectued when nthread > 1 (it can be 0, 1, 2; default 2, implying parallel processing along the 3rd image dimension)')
	parser.add_argument('--nlfit', metavar='<num>', default='1', help='perform non-linear fitting via gradient descent (constrained objective function minimisation) if --nlfit 1, perform only grid search if --nlfit 0 (default: --nlfit 1)')
	parser.add_argument('--nlalgo', metavar='<string>', default='trust-constr', help='algorithm to be used for constrained objective function minimisation in non-linear fitting (relevant if --nlfit 1). Choose among: "Nelder-Mead", "L-BFGS-B", "TNC", "SLSQP", "Powell", and "trust-constr" (default: "trust-constr" - see documentation of scipy.optimize.minimize for information on the optimisation algorithm)')
	parser.add_argument('--modstr', metavar='<string>', default='DinDex', help='string specifying the signal model to fit. Choose among: Din (extra-vascular signal dominated by intra-cellular diffusion), DinDex (extra-vascular signal features both intra-cellular and extra-cellular contributions, without diffusion time dependence in the extra-cellular ADC), DinDexTD (extra-vascular signal features both intra-cellular and extra-cellular contributions, with diffusion time dependence in the extra-cellular ADC). Default: DinDex.                    Intra-cellular diffusion is modelled using the Gaussian Phase Approximation (GPA) formula for diffusion within sphere                     of Van Geldern et al, "Evaluation of restricted diffusion in cylinders. Phosphocreatine in rabbit leg muscle", J Magn Reson B 1994, 103(3):255-60, doi: 10.1006/jmrb.1994.1038')
	args = parser.parse_args()


	### Get input arguments
	sfile = args.s_file
	schfile = args.scheme_file
	outroot = args.out
	maskfile = args.mask
	ffile = args.fv
	t2vfile = args.t2v
	t2tfile = args.t2t
	sgmfile = args.noise
	nex = int(args.savg)
	bfilter = float(args.bfil)
	bmval = float(args.bmin)
	nword = int(args.nw)
	cpucnt = int(args.ncpu)
	slptr = int(args.sldim)
	nlflag = int(args.nlfit)
	nlopt = args.nlalgo
	mymodel = args.modstr
	if(nlflag==1):
		nlbool = True
	else:
		nlbool = False
	
	# Lower bounds
	if args.pmin is not None:
		pMin = (args.pmin).split(',')
		pMin = np.array( list(map( float, pMin )) )
	else:
		if(mymodel=='DinDexTD'):
			pMin = np.array([8.0,0.8,0.0,0.0,0.0])
		elif(mymodel=='DinDex'):
			pMin = np.array([8.0,0.8,0.0,0.0])
		elif(mymodel=='Din'):
			pMin = np.array([8.0,0.8,0.0])
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(mymodel))
	
	# Upper bounds
	if args.pmax is not None:
		pMax = (args.pmax).split(',')
		pMax = np.array( list(map( float, pMax )) )
	else:
		if(mymodel=='DinDexTD'):
			pMax = np.array([40.0,3.0,1.0,3.0,10.0])
		elif(mymodel=='DinDex'):
			pMax = np.array([40.0,3.0,1.0,3.0])
		elif(mymodel=='Din'):
			pMax = np.array([40.0,3.0,1.0])
		else:
			raise RuntimeError('ERROR: signal model {} is unknown.'.format(mymodel))
	
	
	### Print feedback
	print('')
	print('***********************************************************************')
	print('                             dri2mc_maxlik.py                          ')
	print('***********************************************************************')
	print('')
	print('** 4D NIFTI file with MRI measurements: {}'.format(sfile))
	print('** MRI sequence parameter text file: {}'.format(schfile))
	print('** b-value threshold (smaller b-values will not be included in the fitting): {} s/mm2'.format(bfilter))
	print('** b-value below which measurements will be treated as non-DW: {} s/mm2'.format(bmval))
	print('** Model to be fitted: {}'.format(mymodel))
	print('** Number of words for each tissue parameter grid search: {}'.format(nword))
	print('** Non-linear fitting via objective function minimisation required: {}'.format(nlbool))
	if(nlbool):
		print('** Non-linear fitting constrained minimisation algorithm: {}'.format(nlopt))
	print('** Number of threads for parallel slice processing: {}'.format(cpucnt))
	print('** Slice dimension for parallel processing: {}'.format(slptr))
	print('** Lower bound for tissue parameters: {}'.format(pMin))
	print('** Upper bound for tissue parameters: {}'.format(pMax))
	if( maskfile is not None ):
		print('** Optional binary mask file: {}'.format(maskfile))
	if( ffile is not None ):
		print('** Optional 3D NIFTI storing the vascular fraction map: {}'.format(ffile))
	if( t2vfile is not None ):
		print('** Optional 3D NIFTI file with vascular T2 map (ms): {}'.format(t2vfile))
	if( t2tfile is not None ):
		print('** Optional 3D NIFTI file with tissue T2 map (ms): {}'.format(t2tfile))
	if( sgmfile is not None ):
		print('** Optional 3D NIFTI file with noise standard deviation map: {}'.format(sgmfile))
		print('** Number of signal averages: {}'.format(nex))
	print('** Output root name: {}'.format(outroot))
		

	### Run computation
	run(sfile, schfile, outroot, maskfile=maskfile, fvfile=ffile, T2Vfile=t2vfile, T2Tfile=t2tfile, noisefile=sgmfile, navg=nex, bth=bfilter, bmin=bmval, Nword=nword, pmin=pMin, pmax=pMax, nthread=cpucnt, slicedim=slptr, nlinfit=nlbool, nlinalgo=nlopt, modstr=mymodel)
	
	### Job done
	sys.exit(0)

