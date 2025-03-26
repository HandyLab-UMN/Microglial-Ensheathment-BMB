/*  
 * 
 *	susc = calc_Susc_cuEIF(w,E0,sigma,tau_ref,v_reset,v_t,v_th,DT,tau_m,vlb,dv,rate);
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
 *  susc - Susceptibility (1/(mV ms)) of the EIF evaluated at frequencies w. The
 *         susceptibility is the Fourier transform of the linear response function.
 *         Units are 1/(mV ms) since this is the susceptibility to a perturbation
 *         in the membrane potential.
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
    
    double E0,sigma,tau_ref,v_reset,v_th,tau_m,vlb,dv,inv_sigmasq,G,psi,DT,v,v_t,r0,rate;
    double *p0_r,*p0_i,*j0_r,*j0_i,*pr_r,*pr_i,*jr_r,*jr_i,*pe_r,*pe_i,*je_r,*je_i;
    double *w,*susc_r,*susc_i;
    
    mxArray *p0_mat,*j0_mat,*pr_mat,*jr_mat,*pe_mat,*je_mat,*transfer_mat;
    
    
    // Input parameters
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
    plhs[0] = mxCreateDoubleMatrix(length_w, 1, mxCOMPLEX);
    susc_r = mxGetPr(plhs[0]);
    susc_i = mxGetPi(plhs[0]);
    
    
    N = ((int)((v_th-vlb)/dv)) + 1;
    inv_sigmasq = 1/(sigma*sigma);
    
    pr_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    pr_r = mxGetPr(pr_mat);
    pr_i = mxGetPi(pr_mat);
    
    jr_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    jr_r = mxGetPr(jr_mat);
    jr_i = mxGetPi(jr_mat);
    
    p0_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    p0_r = mxGetPr(p0_mat);
    p0_i = mxGetPi(p0_mat);
    
    j0_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    j0_r = mxGetPr(j0_mat);
    j0_i = mxGetPi(j0_mat);
    
    pe_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    pe_r = mxGetPr(pe_mat);
    pe_i = mxGetPi(pe_mat);
    
    je_mat = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    je_r = mxGetPr(je_mat);
    je_i = mxGetPi(je_mat);
    
    // For each frequency value, solve an IVP for the spectrum at that
    // frequency.
    for(i = 0; i < length_w; i++){
        
        // Set the final conditions for the IVP.
        pr_r[N-1] = 0;
        pr_i[N-1] = 0;
        jr_r[N-1] = 1;
        jr_i[N-1] = 0;
        
        p0_r[N-1] = 0;
        p0_i[N-1] = 0;
        j0_r[N-1] = 1;
        j0_i[N-1] = 0;
        
        pe_r[N-1] = 0;
        pe_i[N-1] = 0;
        je_r[N-1] = 0;
        je_i[N-1] = 0;
        
        v = v_th;
    
        // Solve the IVP.
        for(j = N-1;j >= 1; j--){
            psi = DT*exp((v-v_t)/(DT));
            G = (v-E0-psi)*inv_sigmasq;
            
            
            // Iterate j_r/p_r
            
            if(fabs(v-v_reset)<dv/100) {
                jr_r[j-1] = jr_r[j] - dv*w[i]*2*M_PI*pr_i[j] - cos(-w[i]*2*M_PI*tau_ref);
                jr_i[j-1] = jr_i[j] + dv*w[i]*2*M_PI*pr_r[j] - sin(-w[i]*2*M_PI*tau_ref);
            }
            else {
                jr_r[j-1] = jr_r[j] - dv*w[i]*2*M_PI*pr_i[j];
                jr_i[j-1] = jr_i[j] + dv*w[i]*2*M_PI*pr_r[j];
            }
            pr_r[j-1] = pr_r[j]*exp(dv*G) + tau_m*inv_sigmasq*jr_r[j]*(exp(dv*G)-1)/G;
            pr_i[j-1] = pr_i[j]*exp(dv*G) + tau_m*inv_sigmasq*jr_i[j]*(exp(dv*G)-1)/G;
            
            
            
            
            // Iterate j_0/p_0
             
            if(fabs(v-v_reset)<dv/100) 
                j0_r[j-1] = j0_r[j] - 1;
            else
                j0_r[j-1] = j0_r[j];
            
            p0_r[j-1] = p0_r[j]*exp(dv*G) + tau_m*inv_sigmasq*j0_r[j]*(exp(dv*G)-1)/G;
            p0_i[j-1] = p0_i[j]*exp(dv*G) + tau_m*inv_sigmasq*j0_i[j]*(exp(dv*G)-1)/G;
            
            
            // Iterate j_E/p_E
            
            je_r[j-1] = je_r[j] - dv*w[i]*2*M_PI*pe_i[j];
            je_i[j-1] = je_i[j] + dv*w[i]*2*M_PI*pe_r[j];
            
            pe_r[j-1] = pe_r[j]*exp(dv*G) + (tau_m*je_r[j]-rate*p0_r[j])*inv_sigmasq*(exp(dv*G)-1)/G;
            pe_i[j-1] = pe_i[j]*exp(dv*G) + (tau_m*je_i[j]-rate*p0_i[j])*inv_sigmasq*(exp(dv*G)-1)/G;
            
            v = v-dv;
        }
        
        susc_r[i] = -(je_r[0]*jr_r[0] + je_i[0]*jr_i[0])/(pow(jr_r[0],2) + pow(jr_i[0],2));
        susc_i[i] = (je_r[0]*jr_i[0] - je_i[0]*jr_r[0])/(pow(jr_r[0],2) + pow(jr_i[0],2));   
    }
}