/*
This runs the diffusion with recharging traps (DiRT) process
in a 2D simulation on a H-shaped domain, with the cleft, and two side domains
*/
#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"
#include <stdio.h>

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

// Receptor/Capture Region Structure
struct Receptor
{
	double start_pos;  // start of receptor
	double end_pos;    // end of receptor
	double switch_time; // time until it switches back to the absorbing state
	double status; // 0 if open, 1 if occupied
};

// Particle Structure
struct Particle
{
	//current position
	double x_pos; // x position
	double y_pos; // y position
	
	int status;   // 0 for diffusing in cleft, -1 for diffusing left of the cleft, 1 for diffusing right of the cleft, 2 for captured/escaped
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double max_time, dt;
	double D, R, L_z;
	double x_start_loc, y_start_loc, z_start_loc; 
	double tau_r;
	int N, num_trials, M, print_num, max_time_points;
	
	double domain_x_min, domain_x_max, cleft_x_min, cleft_x_max;
	double domain_y_min, domain_y_max, cleft_y_min, cleft_y_max;
	double capture_region_min, capture_region_max;
	
	double affinity;
	
	double *rec_bounds;
				
   /******
    * Import variables from matlab
    * This is messy looking and is specific to mex.
    * Ignore if you're implementing this outside of mex.
    *******/
	max_time =  mxGetScalar(prhs[0]);
	dt = mxGetScalar(prhs[1]);
	D =  mxGetScalar(prhs[2]);
	
	// Domain bounds
	domain_x_min = mxGetScalar(prhs[3]);
	domain_x_max = mxGetScalar(prhs[4]);
	cleft_x_min = mxGetScalar(prhs[5]);	
	cleft_x_max = mxGetScalar(prhs[6]);	
		
	domain_y_min = mxGetScalar(prhs[7]);
	domain_y_max = mxGetScalar(prhs[8]);
	cleft_y_min = mxGetScalar(prhs[9]);
	cleft_y_max = mxGetScalar(prhs[10]);
		
	capture_region_min = mxGetScalar(prhs[11]);
	capture_region_max = mxGetScalar(prhs[12]);
		
	x_start_loc =  mxGetScalar(prhs[13]);
	y_start_loc =  mxGetScalar(prhs[14]);
		
	N = (int)mxGetScalar(prhs[15]);
	
	// mean recharge time
	tau_r =  mxGetScalar(prhs[16]);
	M = (int)mxGetScalar(prhs[17]);
		
	int m1, m2;
	int total_prints = 0;
	
	rec_bounds = mxGetPr(prhs[18]);
	m1 = mxGetM(prhs[18]);
	m2 = mxGetN(prhs[18]);
	if(m1!=M && m2 != 2){
	    mexErrMsgTxt("rec_centers should be M x 2");
	}
		
	print_num = (int) mxGetScalar(prhs[19]);
	max_time_points = (int) mxGetScalar(prhs[20]); 
	affinity = mxGetScalar(prhs[21]);
	
	unsigned int randSeed = (int) mxGetScalar(prhs[22]);
	
	struct Receptor* receptors = mxMalloc(M*sizeof(struct Receptor));
	
	// allocate memory for particles
	struct Particle* particles = mxMalloc(N*sizeof(struct Particle));
	
	double k_constant = sqrt(2*D*dt);
	
	// trial initialization
	double current_time = 0;
	// variables of interest
	// particles_remaining: particle is removed from the domain once it is captured by a receptor
	int particles_remaining = N;
	int total_captured = 0;
	int num_recs_available = M;
	
	// helps keep track of when to store the results 
	int time_count = 1;
		
	// loop variables (could be indeclared at the start of each loop)
	int j,h,rec_loop,outer_count;
	h = 0;
		
	// used for particle reflection
	double x1, x2, y1, y2;
	double m,a,b,c;
	double x_inter1, y_inter1, x_inter2, y_inter2;
	double d1, d2, d;
	double x_inter_final, y_inter_final;
	double tangent_slope, C;
	double x_reflected, y_reflected;
	
	
	/* Allocate output vector */
	double *output_matrix;
	plhs[0] = mxCreateDoubleMatrix(max_time_points,4, mxREAL);
	output_matrix = mxGetPr(plhs[0]);
	
	
	// all receptors start in the available state (i.e., 0)
	for(j=0;j<M;j++){		
		receptors[j].start_pos = rec_bounds[j];
		receptors[j].end_pos = rec_bounds[j + M];
		
		receptors[j].switch_time = -10;
		receptors[j].status = 0;
	}
	
	// seed random number generator
	//srand(time(NULL));
	srand(randSeed);
	
	// particles starting position
	for(j=0;j<N;j++){
		particles[j].x_pos = x_start_loc;
		particles[j].y_pos = y_start_loc;
		particles[j].status = 0;
	}
				
	// main loop
	while((particles_remaining > 0 || num_recs_available < M) && current_time < max_time){
		
		// Update the status of the receptors
		for(rec_loop=0;rec_loop<M;rec_loop++)
		{
			// receptor switchs from closed (1) to open (0)
			// Particle is officially counted as "captured"
	        if(receptors[rec_loop].status == 1 && current_time>receptors[rec_loop].switch_time)
			{
				receptors[rec_loop].status = 0;
	        }
		}
		
		for(j=0; j<N; j++){
			// particle is diffusing in the cleft
			if(particles[j].status==0){
				
				// find the next point via brownian noise
				particles[j].x_pos = particles[j].x_pos+k_constant*randn(0,1);
				particles[j].y_pos = particles[j].y_pos+k_constant*randn(0,1);

				// particle left the cleft
				if(particles[j].x_pos > cleft_x_max){
					particles[j].status = 1;
				
					// particle got capture by the astrocyte
					if(particles[j].x_pos >= domain_x_max){
						particles[j].status=2;
						particles_remaining= particles_remaining - 1;	
					}// particle escaped to the ecs
					else if(particles[j].y_pos <= domain_y_min || particles[j].y_pos >= domain_y_max){
						particles[j].status=2;
						particles_remaining= particles_remaining - 1;
					}
				
				}
				else if(particles[j].x_pos < cleft_x_min){
					particles[j].status = -1;
				
					// particle got capture by the astrocyte
					if(particles[j].x_pos <= domain_x_min){
						particles[j].status=2;
						particles_remaining= particles_remaining - 1;		
					}// particle escaped to the ecs
					else if(particles[j].y_pos <= domain_y_min || particles[j].y_pos >= domain_y_max){
						particles[j].status=2;
						particles_remaining= particles_remaining - 1;
					}
				}
				/* 
				** particle hit the lower boundary in the cleft
				** three options: 1) hit receptor, 2) hit reflecting region
				*/
				else if(particles[j].y_pos <=cleft_y_min){
					// particle might have hit a receptor
					if(particles[j].x_pos >= capture_region_min && particles[j].x_pos <= capture_region_max){
						for(rec_loop=0;rec_loop<M;rec_loop++)
						{
							// particle did hit a receptor and it is avaliable
							if(particles[j].x_pos>receptors[rec_loop].start_pos
								&& particles[j].x_pos<receptors[rec_loop].end_pos
									&& current_time>=receptors[rec_loop].switch_time)
							{
								// flip a coin to see if the particle binds
								if((double)rand() / (double)RAND_MAX < affinity)
								{
									// update switching time
									receptors[rec_loop].switch_time = current_time + -log(((double) rand() / RAND_MAX))*tau_r;
									receptors[rec_loop].status = 1;

									// particle has been captured and is removed from the domain
									particles_remaining  = particles_remaining - 1;
									total_captured = total_captured + 1;
									particles[j].status = 2;

									// exits out of the receptor loop early, since particle is now captured
									rec_loop = M+1;
								}
							}
						}
						// Receptor was not hit if rec_loop==M so reflect point
						if(rec_loop==M)
						{
							particles[j].y_pos = 2*cleft_y_min-particles[j].y_pos;
						}
					}
					// particle hit the lower relecting boundary
					else{
						particles[j].y_pos = 2*cleft_y_min-particles[j].y_pos;
					}							
				}
				/*
				** particle hit upper boundary in the cleft
				** one option: 1) hit reflecting region
				*/ 
				else if(particles[j].y_pos >= cleft_y_max){
					particles[j].y_pos = 2*cleft_y_max-particles[j].y_pos;						        
				}
			
			}
			// particle is diffusing outside of the cleft
			else if(particles[j].status==1 || particles[j].status == -1){
				
				// find the next point via brownian noise
				particles[j].x_pos = particles[j].x_pos+k_constant*randn(0,1);
				particles[j].y_pos = particles[j].y_pos+k_constant*randn(0,1);
			
				// particle reentered the cleft 
				if(particles[j].x_pos <= cleft_x_max && particles[j].x_pos >= cleft_x_min && particles[j].y_pos<=cleft_y_max && particles[j].y_pos >= cleft_y_min){
					particles[j].status = 0;
				}
				// particle escapes to the astrocyte
				else if(particles[j].x_pos <= domain_x_min || particles[j].x_pos >= domain_x_max){
				
					particles[j].status=2;
					particles_remaining= particles_remaining - 1;		
				}
				// particle escapes to the ECS
				else if(particles[j].y_pos <= domain_y_min || particles[j].y_pos >= domain_y_max){
					particles[j].status=2;
					particles_remaining= particles_remaining - 1;
				}
				// particle needs to get reflected 
				else if(particles[j].x_pos >= cleft_x_min && particles[j].x_pos <= cleft_x_max)
				{
					if(particles[j].status == 1){
						particles[j].x_pos = 2*cleft_x_max-particles[j].x_pos;
					}
					else if(particles[j].status == -1){
						particles[j].x_pos = -2*cleft_x_min-particles[j].x_pos;
					}
				}
			}
			
			// if all particles are removed, end trial
			if(particles_remaining == 0){
				j = N+1;
			}
		}
		
		// Print P(t), C(t), and R(t) every 0.0001 time units
		if(time_count == print_num)
		{
			num_recs_available=0;
			for(rec_loop=0;rec_loop<M;rec_loop++)
			{
		        if(current_time>=receptors[rec_loop].switch_time)
				{
					num_recs_available = num_recs_available+1;
		        }
			}
		
			output_matrix[total_prints] = current_time;
			output_matrix[total_prints+1*max_time_points] = particles_remaining;
			output_matrix[total_prints+2*max_time_points] = total_captured;
			output_matrix[total_prints+3*max_time_points] = num_recs_available;
						
			time_count =1;
		
			total_prints = total_prints + 1;
		}else
		{
			time_count = time_count + 1;
		}

		current_time +=  dt;
	}
				
	mxFree(receptors);
	mxFree(particles);
}










