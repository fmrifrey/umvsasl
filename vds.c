#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define VDS_GAMMA 	4258.0		/* Hz/G */
#define VDS_DEBUG	0	
/* Uncomment to run main()
 * #define TESTCODE
 */

/*
%
%	VARIABLE DENSITY SPIRAL GENERATION:
%	----------------------------------
%
%	This is a general description of how the following C code
%	works.  This text is taken from a matlab script, vds.m, from
%	which the C code was derived.  However, note that the C code
%	runs considerably faster.
%
%
%	Function generates variable density spiral which traces
%	out the trajectory
%				 
%			k(t) = r(t) exp(i*q(t)), 		[1]
%
%	Where q IS THE SAME AS theta, and r IS THE SAME AS kr.
%
%		r and q are chosen to satisfy:
%
%		1) Maximum gradient amplitudes and slew rates.
%		2) Maximum gradient due to FOV, where FOV can
%		   vary with k-space radius r, as
%
%			FOV(r) = F0 + F1*r + F2*r*r 		[2]
%
%
%	INPUTS:
%	-------
%	smax = maximum slew rate G/cm/s
%	gmax = maximum gradient G/cm (limited by Gmax or FOV)
%	T = sampling period (s) for gradient AND acquisition.
%	N = number of interleaves.
%	F0,F1,F2 = FOV coefficients with respect to r - see above.
%	rmax= value of k-space radius at which to stop (cm^-1).
%		rmax = 1/(2*resolution);
%
%
%	OUTPUTS:
%	--------
%	k = k-space trajectory (kx+iky) in cm-1.
%	g = gradient waveform (Gx+iGy) in G/cm.
%	s = derivative of g (Sx+iSy) in G/cm/s.
%	time = time points corresponding to above (s).
%	r = k-space radius vs time (used to design spiral)
%	theta = atan2(ky,kx) = k-space angle vs time.
%
%
%	METHODS:
%	--------
%	Let r1 and r2 be the first derivatives of r in [1].	
%	Let q1 and q2 be the first derivatives of theta in [1].	
%	Also, r0 = r, and q0 = theta - sometimes both are used.
%	F = F(r) defined by F0,F1,F2.
%
%	Differentiating [1], we can get G = a(r0,r1,q0,q1,F)	
%	and differentiating again, we get S = b(r0,r1,r2,q0,q1,q2,F)
%
%	(functions a() and b() are reasonably easy to obtain.)
%
%	FOV limits put a constraint between r and q:
%
%		dr/dq = N/(2*pi*F)				[3]	
%
%	We can use [3] and the chain rule to give 
%
%		q1 = 2*pi*F/N * r1				[4]
%
%	and
%
%		q2 = 2*pi/N*dF/dr*r1^2 + 2*pi*F/N*r2		[5]
%
%
%
%	Now using [4] and [5], we can substitute for q1 and q2
%	in functions a() and b(), giving
%
%		G = c(r0,r1,F)
%	and 	S = d(r0,r1,r2,F,dF/dr)
%
%
%	Using the fact that the spiral should be either limited
%	by amplitude (Gradient or FOV limit) or slew rate, we can
%	solve 
%		|c(r0,r1,F)| = |Gmax|  				[6]
%
%	analytically for r1, or
%	
%	  	|d(r0,r1,r2,F,dF/dr)| = |Smax|	 		[7]
%
%	analytically for r2.
%
%	[7] is a quadratic equation in r2.  The smaller of the 
%	roots is taken, and the real part of the root is used to
%	avoid possible numeric errors - the roots should be real
%	always.
%
%	The choice of whether or not to use [6] or [7], and the
%	solving for r2 or r1 is done by calcthetadotdot().
%
%	Once the second derivative of theta(q) or r is obtained,
%	it can be integrated to give q1 and r1, and then integrated
%	again to give q and r.  The gradient waveforms follow from
%	q and r. 	
%
%	Brian Hargreaves -- Sept 2000.
%
%
*/






/* ----------------------------------------------------------------------- */
void calcthetadotdot(float slewmax, float gradmax, float kr, float krdot, float Tgsample, float Tdsample, int Ninterleaves,
				float *fov, int numfov, float *thetadotdot, float *krdotdot)
/*
 * Function calculates the 2nd derivative of kr and theta at each
 * sample point within calc_vds().  ie, this is the iterative loop
 * for calc_vds.  See the text at the top of this file for more details
 *
 *	float slewmax;			Maximum slew rate, G/cm/s		
 *	float gradmax;		 	maximum gradient amplitude, G/cm	
 *	float kr;		 	Current kr. 
 *	float krdot;			Current krdot. 
 *	float Tgsample;		Gradient Sample period (s) 	
 *	float Tdsample;		Data Sample period (s) 		
 *	int Ninterleaves;		Number of interleaves		
 *	float *fov;			FOV coefficients		
 *	int numfov;			Number of FOV coefficients	
 *	float *thetadotdot;		[output] 2nd derivative of theta.
 *	float *krdotdot;		[output] 2nd derivative of kr	
 */
{
float fovval=0;	/* FOV for this value of kr	*/
float dfovdrval=0;	/* dFOV/dkr for this value of kr	*/
float gmaxfov;		/* FOV-limited Gmax.	*/
float maxkrdot;
int count;

float tpf;	/* Used to simplify expressions. */
float tpfsq;	/* 	" 		"        */

float qdfA, qdfB, qdfC;	/* Quadratic formula coefficients */
float rootparta,rootpartb;



if (VDS_DEBUG>1)
	{
	printf("calcthetadotdot:  slewmax=%8.2f, gmax=%6.2f, \n",
			slewmax,gradmax);
	printf("        kr=%8.4f, Tg=%9.6f, N=%d, nfov=%d \n", 
			kr,Tgsample,Ninterleaves,numfov);
	}

	/* Calculate the actual FOV and dFOV/dkr for this R,
	 * based on the fact that the FOV is expressed 
	 * as a polynomial in kr.*/

for (count=0; count < numfov; count++)
	{
	fovval = fovval + fov[count]*pow(kr,count);
	if (count > 0)
		dfovdrval = dfovdrval + count*fov[count]*pow(kr,count-1);
	}

	/* Calculate FOV limit on gmax.  This is the rate of motion along
	 * a trajectory, and really should not be a limitation.  Thus,
	 * it is reasonable to comment out the following lines. */

gmaxfov = 1/VDS_GAMMA / fovval / Tdsample;	
if (gradmax > gmaxfov)
	gradmax = gmaxfov;	


	/* Maximum dkr/dt, based on gradient amplitude.  */

maxkrdot = sqrt(pow(VDS_GAMMA*gradmax,2) / (1+pow(2*M_PI*fovval*kr/Ninterleaves,2)));
if (VDS_DEBUG>1)
	printf("calcthetadotdot:  maxkrdot = %g \n",maxkrdot);

	/* These two are just to simplify expressions below */
tpf = 2*M_PI*fovval/Ninterleaves;
tpfsq = pow(tpf,2);
if (VDS_DEBUG>1)
	printf("calcthetadotdot:  tpf = %8.4f,  tpfsq = %8.4f  \n",tpf,tpfsq);




if (krdot > maxkrdot)	/* Then choose krdotdot so that krdot is in range */
	{	
	*krdotdot = (maxkrdot - krdot)/Tgsample;
	}

else			/* Choose krdotdot based on max slew rate limit. */
	{

		/* Set up for quadratic formula solution. */

	qdfA = 1+tpfsq*kr*kr;
	qdfB = 2*tpfsq*kr*krdot*krdot + 
			2*tpfsq/fovval*dfovdrval*kr*kr*krdot*krdot;
	qdfC = pow(tpfsq*kr*krdot*krdot,2) + 4*tpfsq*pow(krdot,4) +
			pow(tpf*dfovdrval/fovval*kr*krdot*krdot,2) +
			4*tpfsq*dfovdrval/fovval*kr*pow(krdot,4) -
			pow(VDS_GAMMA*slewmax,2);

	if (VDS_DEBUG>1)
		printf("calcthetadotdot:  qdfA, qdfB, qdfC = %g, %g, %g \n",
				qdfA, qdfB, qdfC);

	rootparta = -qdfB/(2*qdfA);
	rootpartb = qdfB*qdfB/(4*qdfA*qdfA) - qdfC/qdfA;
	if (VDS_DEBUG>1)
		printf("calcthetadotdot:  rootparta, rootpartb = %g, %g \n",
				rootparta, rootpartb);

	if (rootpartb < 0)	/* Safety check - if complex, take real part.*/

		*krdotdot = rootparta;

	else
		*krdotdot = rootparta + sqrt(rootpartb);


	/* Could check resulting slew rate here, as in q2r21.m. */
	}

	/* Calculate thetadotdot */

	
*thetadotdot = tpf*dfovdrval/fovval*krdot*krdot + tpf*(*krdotdot);

if (VDS_DEBUG>1)
	printf("calcthetadot:  r=%8.4f,  r'=%8.4f,  r''=%g  q''=%g \n",
		kr,krdot,*krdotdot,*thetadotdot);

}


/* ----------------------------------------------------------------------- */
void calc_vds(float slewmax, float gradmax, float Tgsample, float Tdsample, int Ninterleaves, float *fov, int numfov, float krmax,
		int ngmax, float **xgrad, float **ygrad, int *numgrad)

/*	Function designs a variable-density spiral gradient waveform
 *	that is defined by a number of interleaves, resolution (or max number
 *	of samples), and field-of-view.  
 *	The field-of-view is a polynomial function of the
 *	k-space radius, so fov is an array of coefficients so that
 *
 *	FOV = fov[0]+fov[1]*kr+fov[2]*kr^2+ ... +fov[numfov-1]*kr^(numfov-1)
 *
 * 	Gradient design is subject to a constant-slew-rate-limit model,
 * 	with maximum slew rate slewmax, and maximum gradient amplitude
 * 	of gradmax.  
 *
 * 	Tgsample is the gradient sampling rate, and Tdsample is the data
 * 	sampling rate.  It is highly recommended to OVERSAMPLE the gradient
 * 	in the design to make the integration more stable.
 *
 *
 *	float slewmax;			Maximum slew rate, G/cm/s		
 *	float gradmax;		 	maximum gradient amplitude, G/cm	
 *	float Tgsample;		Gradient Sample period (s)		
 *	float Tdsample;		Data Sample period (s)			
 *	int Ninterleaves;		Number of interleaves			
 *	float *fov;			FOV coefficients		
 *	int numfov;			Number of FOV coefficients		
 *	float krmax;			Maximum k-space extent (/cm)		
 *	int ngmax;			Maximum number of gradient samples	
 *	float **xgrad;		 	[output] X-component of gradient (G/cm) 
 *	float **ygrad;			[output] Y-component of gradient (G/cm)	
 *	int *numgrad;		 	[output] Number of gradient samples
 */
{
int gradcount=0;

float kr=0;			/* Current value of kr	*/
float krdot = 0;		/* Current value of 1st derivative of kr */
float krdotdot = 0;		/* Current value of 2nd derivative of kr */

float theta=0;			/* Current value of theta */
float thetadot=0;		/* Current value of 1st derivative of theta */
float thetadotdot=0;		/* Current value of 2nd derivative of theta */

float lastkx=0;		/* x-component of last k-location. */
float lastky=0;		/* y-component of last k-location */
float kx, ky;			/* x and y components of current k-location */

float *gxptr, *gyptr;		/* Pointers to gradient variables. */




if (VDS_DEBUG>0)
	printf("calc_vds:  First run. \n");

	/* First just find the gradient length. */

while ((kr < krmax) && (gradcount < ngmax))
	{
	calcthetadotdot(slewmax,gradmax,kr,krdot,Tgsample,Tdsample,
			Ninterleaves, fov,numfov, &thetadotdot, &krdotdot);

	/* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

	thetadot = thetadot + thetadotdot * Tgsample;
	theta = theta + thetadot * Tgsample;

	krdot = krdot + krdotdot * Tgsample;
	kr = kr + krdot * Tgsample;

	gradcount++;

	}



	/* Allocate memory for gradients. */

*numgrad = gradcount;
if (VDS_DEBUG>0)
	printf("Allocating for %d gradient points. \n",*numgrad);

*xgrad = (float *)malloc(*numgrad*sizeof(float));
*ygrad = (float *)malloc(*numgrad*sizeof(float));


	/* Reset parameters */

kr=0;
krdot=0;
theta=0;
thetadot=0;
gradcount=0;
gxptr = *xgrad;
gyptr = *ygrad;


	/* Now re-calculate gradient to find length. */

if (VDS_DEBUG>0)
	printf("calc_vds:  First run. \n");

while ((kr < krmax) && (gradcount < ngmax))
	{
	calcthetadotdot(slewmax,gradmax,kr,krdot,Tgsample,Tdsample,
			Ninterleaves, fov,numfov, &thetadotdot, &krdotdot);

	/* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

	thetadot = thetadot + thetadotdot * Tgsample;
	theta = theta + thetadot * Tgsample;

	krdot = krdot + krdotdot * Tgsample;
	kr = kr + krdot * Tgsample;

	/* Define current gradient values from kr and theta. */

	kx = kr * cos(theta);
	ky = kr * sin(theta);
	*gxptr++ = (1/VDS_GAMMA/Tgsample) * (kx-lastkx);
	*gyptr++ = (1/VDS_GAMMA/Tgsample) * (ky-lastky);
	lastkx = kx;
	lastky = ky;

	if (VDS_DEBUG>0)
		printf("Current kr is %6.3f \n",kr);

	gradcount++;
	}

}



#ifdef TESTCODE

int main(void)

{

float SLEWMAX = 17000;
float GMAX = 5;
int GRAD_UPDATE_TIME = 4;
float opfov = 240.0;
int opxres = 64;
int MAXWAVELEN = 50000;
float spalpha = 0.25;
int ntrains = 1;

float kr = opxres/(opfov/10.0)/2;

float *gx, *gy;
int n, ng;
float fov[4];
float thdd, krdd;
FILE *outfile;

fov[0] = spalpha*opfov/10.0;
fov[1] = opfov/10.0/kr*(1-spalpha);

printf("Calculating waveform.\n");
calc_vds(SLEWMAX, GMAX, GRAD_UPDATE_TIME*1e-6, GRAD_UPDATE_TIME*1e-6, ntrains, fov, 2, kr, MAXWAVELEN, &gx, &gy, &ng); 
printf("%d gradient samples \n",ng);

outfile = fopen("vdstest.txt", "wb");
for (n = 0; n < ng; n++)
	fprintf(outfile, "%f \t%f\n", gx[n], gy[n]);
fclose(outfile);
}
#endif

