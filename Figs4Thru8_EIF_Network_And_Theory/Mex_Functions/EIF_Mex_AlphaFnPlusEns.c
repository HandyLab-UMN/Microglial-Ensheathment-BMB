/*
**
** Written by Gregory Handy on 08/24/2021
** Adjusted and translated from the Julia code used in
** Litwin-Kumar et al., 2016 to mex
**
*/
#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"

/* returns a pseudorandom value drawn from the standard normal distribution */
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	int istim, T, NT, recstart; 
	int Npop, *Ncells, Ntot, m1, m2;
	double dt, *rates; 
	int *pinds;
	
	/* Synaptic parameters */
	double *p0, *J;
	double *tau_s;
	double *tau_m, *EL;
	double *vth, *vre, *tauref, *DeltaT, *VT;
	
	/* Connectivity matrix */
	int *wind, *wipost;
	double *wstr;
	//	, *wext, *rext;
	
	/* indexing variables */
	int tt, cc, kk, i, j, pp, jj;
	
   /******
    * Import variables from matlab
    * This is messy looking and is specific to mex.
    * Ignore if you're implementing this outside of mex.
    *******/
	T =  (int)mxGetScalar(prhs[0]);
	NT = (int)mxGetScalar(prhs[1]);
	recstart =  (int)mxGetScalar(prhs[2]);
	Npop = (int)mxGetScalar(prhs[3]);
	Ntot = (int)mxGetScalar(prhs[4]);
	
	/*
	  Needed for vectors of integers, which is important
	  since they will be used to access array locations
	  As a result, the mex call has been motified
	*/
	mxArray *temp[1];
	mxArray *lhs[1];
	temp[0] = (mxArray *) prhs[5];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	Ncells=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[5]);
	m2 = mxGetN(prhs[5]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("Ncells should be Npopx1 or 1xNpop");

	dt = mxGetScalar(prhs[6]);

	rates= mxGetPr(prhs[7]);
	m1 = mxGetM(prhs[7]);
	m2 = mxGetN(prhs[7]);
	if(m1!=Ntot && m2!=Ntot)
	    mexErrMsgTxt("rates should be 1xNtot or Ntotx1");

	/* load in synaptic parameters */
	p0 = mxGetPr(prhs[8]);
	m1 = mxGetM(prhs[8]);
	m2 = mxGetN(prhs[8]);
	if(m1!=Npop || m2!=Npop)
	    mexErrMsgTxt("p0 should be NpopxNpop");

	J = mxGetPr(prhs[9]);
	m1 = mxGetM(prhs[9]);
	m2 = mxGetN(prhs[9]);
	if(m1!=Npop || m2!=Npop)
	    mexErrMsgTxt("J should be NpopxNpop");

	/* Input 10 moved down past input 29 */

	tau_m= mxGetPr(prhs[11]);
	m1 = mxGetM(prhs[11]);
	m2 = mxGetN(prhs[11]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("tau_m should be 1xNpop or Npopx1");

	EL= mxGetPr(prhs[12]);
	m1 = mxGetM(prhs[12]);
	m2 = mxGetN(prhs[12]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("EL should be 1xNpop or Npopx1");

	vth= mxGetPr(prhs[13]);
	m1 = mxGetM(prhs[13]);
	m2 = mxGetN(prhs[13]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("vth should be 1xNpop or Npopx1");

	vre= mxGetPr(prhs[14]);
	m1 = mxGetM(prhs[14]);
	m2 = mxGetN(prhs[14]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("vre should be 1xNpop or Npopx1");

	tauref= mxGetPr(prhs[15]);
	m1 = mxGetM(prhs[15]);
	m2 = mxGetN(prhs[15]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("tauref should be 1xNpop or Npopx1");

	/*
	  Needed for vectors of integers, which is important
	  since they will be used to access array locations
	  As a result, the mex call has been motified
	*/
	temp[0] = (mxArray *) prhs[16];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	wind=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[16]);
	m2 = mxGetN(prhs[16]);
	if(m1!=(Ntot+1) && m2!=(Ntot+1))
	    mexErrMsgTxt("wind should be (Ntot+1)x1 or 1x(Ntot+1)");

	temp[0] = (mxArray *) prhs[17];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	wipost=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[17]);
	m2 = mxGetN(prhs[17]);

	wstr= mxGetPr(prhs[18]);
	m1 = mxGetM(prhs[18]);
	m2 = mxGetN(prhs[18]);

	temp[0] = (mxArray *) prhs[19];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	pinds=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[19]);
	m2 = mxGetN(prhs[19]);
	if(m1!=(Npop+1) && m2!=(Npop+1))
	    mexErrMsgTxt("pinds should be (Npop+1)x1 or 1x(Npop+1)");

	unsigned int iseed;
	iseed =  (int)mxGetScalar(prhs[20]);
	
	
	double *mu_vec, *sigma_vec;
	mu_vec = mxGetPr(prhs[21]);
	m1 = mxGetM(prhs[21]);
	m2 = mxGetN(prhs[21]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("mu_vec should be Npopx1 or 1xNpop");
	
	sigma_vec = mxGetPr(prhs[22]);
	m1 = mxGetM(prhs[22]);
	m2 = mxGetN(prhs[22]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("sigma_vec should be Npopx1 or 1xNpop");
	
	DeltaT = mxGetPr(prhs[23]);
	m1 = mxGetM(prhs[23]);
	m2 = mxGetN(prhs[23]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("DeltaT should be Npopx1 or 1xNpop");
	
	VT = mxGetPr(prhs[24]);
	m1 = mxGetM(prhs[24]);
	m2 = mxGetN(prhs[24]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("VT should be Npopx1 or 1xNpop");
	
	
	int *tauDelay;
	temp[0] = (mxArray *) prhs[25];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	tauDelay=(int*)mxGetData(lhs[0]);
	m1 = mxGetM(prhs[25]);
	m2 = mxGetN(prhs[25]);
	if(m1!=Npop && m2!=Npop)
	    mexErrMsgTxt("tauDelay should be Npopx1 or 1xNpop");
	
	int maxDelay;
	maxDelay = (int)mxGetScalar(prhs[26]);
	maxDelay = maxDelay+1;
	
	double shared_noise_mag;
	shared_noise_mag = mxGetScalar(prhs[27]);
	
	int maxns;
	maxns = (int)mxGetScalar(prhs[28]);
	
	int sValues;
	sValues=(int)mxGetScalar(prhs[29]);
	
	tau_s = mxGetPr(prhs[10]);
	m1 = mxGetM(prhs[10]);
	m2 = mxGetN(prhs[10]);
	if(m1!=Npop*sValues && m2!=Npop*sValues)
	    mexErrMsgTxt("tau_s should be sValues*Npopx1 or 1xNpop*sValues");
	
	int *ensh_realizations;
	temp[0] = (mxArray *) prhs[30];
	mexCallMATLAB(1,lhs,1,&temp[0],"int32");
	ensh_realizations = (int*)mxGetData(lhs[0]);
	m1 = mxGetN(prhs[30]);
	m2 = mxGetM(prhs[30]);
		
	/* Simple "srand()" seed: just use "time()" */
	/* iseed = (unsigned int)time(NULL); */
	// srand (iseed);
	
	/*	
	Code for mex printing
	mexPrintf("test \n");
	mexEvalString("drawnow;");
	*/
	
	/* NOT IMPORTED FROM MATLAB */
	/* state vectors */
	double *v, *lastSpike;
	int *whichpop;

	v = mxMalloc(Ntot * sizeof(double));
	lastSpike = mxMalloc(Ntot * sizeof(double));
	whichpop = mxMalloc(Ntot * sizeof(int));

	/* Store the synaptic decay parameters in three separate arrays*/
	double *tausyn_e,*tausyn_PV,*tausyn_SST;
	tausyn_e=mxMalloc(sValues*sizeof(double));
	tausyn_PV=mxMalloc(sValues*sizeof(double));
	tausyn_SST=mxMalloc(sValues*sizeof(double));
	for(kk=0; kk<sValues;kk++){
		tausyn_e[kk] = tau_s[kk];
		tausyn_PV[kk] = tau_s[kk+sValues];
		tausyn_SST[kk] = tau_s[kk+2*sValues];
	}

	// Allocate memory for the synaptic currents
	double *I_se, *I_sPV,*I_sSST,*y_se,*y_sPV,*y_sSST;
	I_se = mxMalloc(Ntot*sValues * sizeof(double));
	I_sPV = mxMalloc(Ntot*sValues * sizeof(double));
	I_sSST = mxMalloc(Ntot*sValues * sizeof(double));
	
	y_se = mxMalloc(Ntot*sValues * sizeof(double));
	y_sPV = mxMalloc(Ntot*sValues * sizeof(double));
	y_sSST = mxMalloc(Ntot*sValues * sizeof(double));

	int ns = 0; // counts the number of spikes	

	int *counts;
	counts = mxMalloc(Ntot * sizeof(int));
	
	int *spikeDelay;
	spikeDelay = mxMalloc(Ntot * maxDelay * sizeof(int));
	for(cc = 0; cc<Ntot;cc++){
		for(j = 0;j<maxDelay;j++){
			spikeDelay[j+cc*maxDelay] = 0;
		}
	}

	/* Allocate output vector */
	double *rates_official;
	plhs[0] = mxCreateDoubleMatrix(1, Ntot, mxREAL);
	rates_official=mxGetPr(plhs[0]);


	double *times;
	double *tinds;
	plhs[1] = mxCreateDoubleMatrix(1, maxns, mxREAL);
	times=mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(1, maxns, mxREAL);
	tinds=mxGetPr(plhs[2]);
	
	
	/* initial conditions */
	for(i = 0; i<Ntot; i++){
		v[i] = -65;
		lastSpike[i] = -100;

		/* save whichpop neuron i is */
		pp = 1;
		int initial_vT = 0;
		while(initial_vT == 0){
			if(i < pinds[pp]){
				whichpop[i] = pp-1;
				initial_vT = 1;
			}else{
				pp++;
			}

			if(pp > (Npop+1) ){
				mexErrMsgTxt("vT initialization failed");
			}
		}

		// These are the synaptic currents going into each neuron
		for(kk=0; kk<sValues;kk++){
			I_se[i*sValues+kk]=0;
			I_sPV[i*sValues+kk]=0;
			I_sSST[i*sValues+kk]=0;
			y_se[i*sValues+kk]=0;
			y_sPV[i*sValues+kk]=0;
			y_sSST[i*sValues+kk]=0;
		}
	}
	
	/* spike counts and times */
	ns = 0;
	for(i =0; i< maxns; i++){
		times[i] = 0;
		tinds[i] = 0;
	}
	for(cc = 0; cc < Ntot; cc++){
		counts[cc] = 0;
	}
		
	double sigma_hat[3];
	int pc;
	// for(pc = 0; pc<Npop; pc++){
	// 	// sigma_hat[pc] = var_vec[pc]*sqrt(tau_m[pc])/tau_s[pc]*sqrt(dt);
	// 	sigma_hat[pc] = var_vec[pc]*sqrt(tau_m[pc])*sqrt(dt);
	// }
	 	 
	int delayTracker = 0;	
	int delayLoc; 
		 
	/* simulation------------------------------------------------------------- */
	for(tt = 1; tt <= NT && ns<maxns; tt++){

		double t = dt*tt;
		
		double shared_noise = shared_noise_mag*randn(0,1)*sqrt(dt)/tau_m[pc];

		/* update synaptic/adaptation parameters */
		for(cc=0; cc<Ntot; cc++){
			int pc = whichpop[cc];

			// Uses Euler's method to approximate the exponential decay
			double I_rec = 0;
		  	for(kk=0; kk<sValues;kk++){
				I_se[cc*sValues+kk] = I_se[cc*sValues+kk]+dt*(y_se[cc*sValues+kk]-I_se[cc*sValues+kk])/tausyn_e[kk];
				y_se[cc*sValues+kk] = y_se[cc*sValues+kk]-dt*(y_se[cc*sValues+kk]/tausyn_e[kk]);
				
				
				I_sPV[cc*sValues+kk] = I_sPV[cc*sValues+kk]+dt*(y_sPV[cc*sValues+kk]-I_sPV[cc*sValues+kk])/tausyn_PV[kk];
				y_sPV[cc*sValues+kk] = y_sPV[cc*sValues+kk]-dt*(y_sPV[cc*sValues+kk]/tausyn_PV[kk]);
				
				I_sSST[cc*sValues+kk] = I_sSST[cc*sValues+kk]+dt*(y_sSST[cc*sValues+kk]-I_sSST[cc*sValues+kk])/tausyn_SST[kk];
				y_sSST[cc*sValues+kk] = y_sSST[cc*sValues+kk]-dt*(y_sSST[cc*sValues+kk]/tausyn_SST[kk]);
			
				// Total up the recurrent input (added to the EIF equations)
				I_rec += I_se[cc*sValues+kk]+I_sPV[cc*sValues+kk]+I_sSST[cc*sValues+kk];
		  	}     
			
			if(t > (lastSpike[cc]+tauref[pc])){ /* not in refractory period */
								
				v[cc] += dt*(-(v[cc]-EL[pc])+DeltaT[pc]*exp((v[cc]-VT[pc])/DeltaT[pc])
					+I_rec+mu_vec[pc])/tau_m[pc]
					+sigma_vec[pc]*randn(0,1)*sqrt(dt)/tau_m[pc]
					+shared_noise;
				
				v[cc] = fmax(v[cc],-100.0);

			} /* end if(refractory) */
		}


		/* Deal with spikes */
		for(cc=0; cc<Ntot; cc++){
			int pc = whichpop[cc];
			
			/* spike occurred */
			if((v[cc] > vth[pc]) && (ns<maxns)){

				/* update v */
				v[cc] = vre[pc];
				
				delayLoc = (delayTracker+tauDelay[pc])%maxDelay;
				spikeDelay[delayLoc+cc*maxDelay] = 1;
			}
			
			if(spikeDelay[delayTracker+cc*maxDelay]==1){
				
				spikeDelay[delayTracker+cc*maxDelay] = 0;
				
				/* record spike */
				lastSpike[cc] = t;
				ns = ns+1;
				times[ns] = t;
				tinds[ns] = cc;
				if(t > recstart){
					counts[cc] = counts[cc] + 1;
				}

				/* propagate spike */
				for(kk = wind[cc]; kk < (wind[cc+1]-1); kk++){
					int ipost = wipost[kk];
					int ppost = whichpop[ipost];
					
					for(jj=0; jj<sValues;jj++){
						
						if(ensh_realizations[kk] == jj){
							if(pc == 0 || pc == 3 || pc == 6){
								y_se[sValues*ipost+jj] += wstr[kk]/tausyn_e[jj];
							}
							else if(pc == 1 || pc == 4 || pc == 7){
								y_sPV[sValues*ipost+jj] += wstr[kk]/tausyn_PV[jj];
							}
							else if(pc == 2 || pc == 5 || pc == 8){
								y_sSST[sValues*ipost+jj] += wstr[kk]/tausyn_SST[jj];
							}
							else{
								mexErrMsgTxt("Code only works for Npop = 3, 6 or 9");
							}
						} 
					}
				}
			}
		} /* end deal with spikes */
		
		delayTracker = (delayTracker+1)%maxDelay;
	}/* end of loop over time */

	/* Issue a warning if max number of spikes reached */
	if(ns>=maxns){
		// mexPrintf("test \n");
		char warningMessage[300];
		if(tt*dt<recstart){
			sprintf(warningMessage, "Maximum number of spikes reached, simulation terminated after %.2f msec (no spikes recorded)", tt*dt);
		}else{
			sprintf(warningMessage, "Maximum number of spikes reached, simulation terminated after %.2f msec (some spikes recorded)", tt*dt);
		}
		mexWarnMsgTxt(warningMessage);
		// mexWarnMsgTxt(warningMessage2);
	}
	   
	/* Calculate the rates */
	for(i = 0; i < Ntot; i++){
		rates[i] = 1000*(double)counts[i]/ (double)(T-recstart);
		rates_official[i] = rates[i];
	}
		
	/* Free malloc'ed variables */
	mxFree(v);
	mxFree(lastSpike);
	mxFree(whichpop);
	mxFree(counts);
	mxFree(I_se);
	mxFree(I_sPV);
	mxFree(I_sSST);
	mxFree(y_se);
	mxFree(y_sPV);
	mxFree(y_sSST);
	mxFree(spikeDelay);
	
}