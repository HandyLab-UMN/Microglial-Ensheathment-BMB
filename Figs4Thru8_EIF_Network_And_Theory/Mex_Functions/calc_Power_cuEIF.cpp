/*  
 * 
 *	power = calc_Power_cuEIF(w,E0,sigma,tau_ref,v_reset,v_t,v_th,DT,tau_m,vlb,dv,rate);
 *
 *  
 *  Adjusted by Gregory Handy for EIF neurons from code
 *  written by James Trousdale on 08/24/2021.
 *
 *
 *  Purpose: Calculate the power-spectra of a white noise-driven EIF neuron.
 *
 *
 *  Parameters:
 *
 *  w - Array of temporal frequencies (kHz) at which to calculate the spectrum.
 *      Note that the theory diverges at w=0, so zero should be replaced with a
 *      small, but positive value.
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
 *  rate - Stationary firing rate (kHz) of the EIF calculated from a separate program.
 *
 *
 *
 *  Output:
 *
 *  power - Power spectrum (1/ms) of the EIF evaluated at frequencies w.
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
    int N,i,j,length_w,m;
    
    double E0,sigma,tau_ref,v_reset,v_th,tau_m,vlb,dv,inv_sigmasq,inv_tau_m;
    double v,G,psi,DT,v_t,pw_r,rate;
    double *pf_r,*pf_i,*jf_r,*jf_i,*p0_r,*p0_i,*j0_r,*j0_i,*w,*power;
    
    mxArray *pf_mat, *jf_mat, *p0_mat,*j0_mat,*power_mat;
    
    // Retrieve all input parameters.
    w=mxGetPr(prhs[0]);
    E0=mxGetScalar(prhs[1]);
    sigma=mxGetScalar(prhs[2]);
    tau_ref=mxGetScalar(prhs[3]);
    v_reset=mxGetScalar(prhs[4]);
    v_t=mxGetScalar(prhs[5]);
    v_th=mxGetScalar(prhs[6]);
    DT=mxGetScalar(prhs[7]);
    tau_m=mxGetScalar(prhs[8]);
    vlb=mxGetScalar(prhs[9]);
    dv=mxGetScalar(prhs[10]);
    rate=mxGetScalar(prhs[11]);
    
    

    
    length_w = mxGetN(prhs[0]);
    m = mxGetM(prhs[0]);
    
    if(length_w != 1 && m != 1)
    	mexErrMsgTxt("Frequency vector should be a Nx1 vector.");
    if(length_w < m)
        length_w = m;

    
    
    // Output parameters
    plhs[0] = mxCreateDoubleMatrix(length_w, 1, mxREAL);
    power = mxGetPr(plhs[0]);
    
    
    N = ((int)((v_th-vlb)/dv)) + 1;
    inv_sigmasq = 1/(sigma*sigma);
    inv_tau_m = 1/tau_m;
  
    pf_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    pf_r = mxGetPr(pf_mat);
    pf_i = mxGetPi(pf_mat);
    
    jf_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    jf_r = mxGetPr(jf_mat);
    jf_i = mxGetPi(jf_mat);
    
    p0_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    p0_r = mxGetPr(p0_mat);
    p0_i = mxGetPi(p0_mat);
    
    j0_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    j0_r = mxGetPr(j0_mat);
    j0_i = mxGetPi(j0_mat);
    
    
    // For each frequency value, solve an IVP for the spectrum at that
    // frequency.
    for(i = 0; i < length_w; i++){
        
        // Set the final conditions for the IVP.
        pf_r[N-1] = 0;
        pf_i[N-1] = 0;
        jf_r[N-1] = 1;
        jf_i[N-1] = 0;
        
        p0_r[N-1] = 0;
        p0_i[N-1] = 0;
        j0_r[N-1] = 0;
        j0_i[N-1] = 0;
        
        v = v_th;
    
        // Solve the IVP.
        for(j = N-1;j >= 1; j--){
            psi = DT*exp((v-v_t)/(DT));
            G = (v-E0-psi)*inv_sigmasq;
            
            jf_r[j-1] = jf_r[j] - dv*w[i]*2*M_PI*pf_i[j];
            jf_i[j-1] = jf_i[j] + dv*w[i]*2*M_PI*pf_r[j];
            
            pf_r[j-1] = pf_r[j]*exp(dv*G) + tau_m*jf_r[j]*inv_sigmasq*(exp(dv*G)-1)/G;
            pf_i[j-1] = pf_i[j]*exp(dv*G) + tau_m*jf_i[j]*inv_sigmasq*(exp(dv*G)-1)/G;
            
            if(fabs(v-v_reset)<dv/100){
                j0_r[j-1] = j0_r[j] - dv*w[i]*2*M_PI*p0_i[j] - cos(-w[i]*2*M_PI*tau_ref);
                j0_i[j-1] = j0_i[j] + dv*w[i]*2*M_PI*p0_r[j] - sin(-w[i]*2*M_PI*tau_ref);
            }
            else{
                j0_r[j-1] = j0_r[j] - dv*w[i]*2*M_PI*p0_i[j];
                j0_i[j-1] = j0_i[j] + dv*w[i]*2*M_PI*p0_r[j];
            }
            
            p0_r[j-1] = p0_r[j]*exp(dv*G) + tau_m*j0_r[j]*inv_sigmasq*(exp(dv*G)-1)/G;
            p0_i[j-1] = p0_i[j]*exp(dv*G) + tau_m*j0_i[j]*inv_sigmasq*(exp(dv*G)-1)/G;
            
            v = v-dv;
        }
    
        pw_r = (-(pow(j0_r[0],2) + pow(j0_i[0],2)) - j0_r[0]*jf_r[0] - j0_i[0]*jf_i[0])/(pow(j0_r[0] + jf_r[0],2) + pow(j0_i[0] + jf_i[0],2));
        
        power[i] = rate*(1 + 2*pw_r);
        
    }
}
        
        
        
        