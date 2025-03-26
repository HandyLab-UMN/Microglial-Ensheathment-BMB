/*  
 * 
 *	rate = calc_Rate_cuEIF(E0,sigma,tau_ref,v_reset,v_t,v_th,DT,tau_m,vlb,dv);
 *
 *  
 *  Adjusted by Gregory Handy for EIF neurons from code
 *  written by James Trousdale on 08/24/2021.
 *
 *
 *  Purpose: Calculate the stationary firing rate of a white noise-driven EIF neuron.
 *
 *
 *  Parameters:
 *
 *  E0 - Effective rest potential (mV) in the absence of noise, inputs or the exponential
 *       nonlinearity (typically the leak potential + noise mean).
 *  sigma - White noise variance (mV).
 *  tau_ref - Absolute refractory period (ms).
 *  v_reset - Reset potential following the emission of a spike (mV).
 *  v_t - Soft threshold of the EIF (mV).
 *  v_th - Hard threshold of the EIF (mV). Upon reaching this threshold, the membrane
 *         potential is reset to v_reset and held there for a fixed amount of time tau_ref.
 *  DT - Spike shape parameter for the EIF (mV).
 *  tau_m - Membrane potential for the neuron (ms).
 *  vlb - Lower bound of the membrane potential (mV). Set sufficiently low, this
 *        term will not impact calculations, and is set for convenience in solving
 *        the system of ODEs for the desired statistics
 *  dv - Membrane potential step (mV) used in solving the system of ODEs for the desired
 *       statistics.
 *
 *
 *
 *  Output:
 *
 *  rate - Stationary firing rate (kHz) of the EIF.
 *
 *
 *
 *  Note:
 *
 *  Written following exactly the methods of Richardson "Spike train spectra
 *  and network response..." (2008).
 */

#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Variable declarations
    int N,i;
    
    double E0,sigma,tau_ref,v_reset,v_th,tau_m,vlb,dv,inv_sigmasq;
    double sum_p0,v_t,DT,v,psi,G;
    double *p0,*j0,*rate;
    
    mxArray *p0_mat,*j0_mat;
    
    // Input parameters
    E0=mxGetScalar(prhs[0]);
    sigma=mxGetScalar(prhs[1]);
    tau_ref=mxGetScalar(prhs[2]);
    v_reset=mxGetScalar(prhs[3]);
    v_t=mxGetScalar(prhs[4]);
    v_th=mxGetScalar(prhs[5]);
    DT=mxGetScalar(prhs[6]);
    tau_m=mxGetScalar(prhs[7]);
    vlb=mxGetScalar(prhs[8]);
    dv=mxGetScalar(prhs[9]);
    
    // Output parameters
    plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
    rate=mxGetPr(plhs[0]);
    
    N = ((int)((v_th-vlb)/dv)) + 1;
    inv_sigmasq = 1/(sigma*sigma);
  
    p0_mat = mxCreateDoubleMatrix(N, 1, mxREAL);
    p0 = mxGetPr(p0_mat);
    
    j0_mat = mxCreateDoubleMatrix(N, 1, mxREAL);
    j0 = mxGetPr(j0_mat);
    
    // Set the final conditions for the IVP. 
    p0[N-1] = 0;
    j0[N-1] = 1;
    
    v = v_th;
    
    // Solve the IVP.
    for(i = N-1; i >= 1; --i){
        psi = DT*exp((v-v_t)/(DT));
        G = (v-E0-psi)*inv_sigmasq;
        
        if(fabs(v-v_reset) < (dv/100))
            j0[i-1] = j0[i] - 1;
        else
            j0[i-1] = j0[i];
        
        p0[i-1] = p0[i]*exp(dv*G) + tau_m*j0[i]*inv_sigmasq*(exp(dv*G)-1)/G;
        
        v = v-dv;
    }
        
    
    sum_p0 = p0[0];
    for(i = 1; i < N; i++)
        sum_p0 += p0[i];
        
    *rate = 1/(dv*sum_p0 + tau_ref);
}
        
        
        
        