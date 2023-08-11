#include <stdio.h>
#include <math.h>

/* 
 * Variable density spiral generation using the algorithm described in:
 * 	Kim DH, Adalsteinsson E, Spielman DM. Simple analytic variable density spiral design. Magn Reson Med.
 * 	2003 Jul;50(1):214-9. doi: 10.1002/mrm.10493. PMID: 12815699.
 */

/* 
 * Arguments:
 * 	kx, ky, 2D kspace trajectory arrays
 * 	D - fov (cm)
 * 	sm - maximum slew rate (G/cm/ms)
 * 	gm - maximum gradient amplitude (G/cm)
 * 	N - matrix size
 * 	N_int - number of interleaves in spiral
 * 	alpha - sampling density factor
 * 	dt - sampling rate (ms)
 */

float kimvdsp(float *kx, float *ky, int D, float sm, float gm, int N, int N_int, float alpha, float dt) {

	float t, tau;
	int i;
	int npoints;

	/* Calculate constants */
	float gam = 4.258*2*M_PI; /* gyromagnetic ratio (rad/G/ms) */
	float n = pow(1.0 - pow(1.0 - 2.0/(float)N, 1.0/alpha), - 1.0) / (float)N_int; /* # of turns in kspace */
	float omega = 2.0*M_PI*n; /* total angular displacement (rad) */
	float lambda = (float)N / (2.0 * (float)D); /* max kspace (1/cm) */
	float c1 = gam * gm / lambda / omega; /* amplitude constraint constant (Hz) */
	float c2 = sqrt(sm * gam / lambda / pow(omega, 2.0)); /* slew constraint constant (Hz) */
	float a1 = alpha + 1; /* amplitude regime sampling density factor */
	float a2 = alpha/2 + 1; /* slew regime sampling density factor */

	/* Calculate important timings */
	float T_es, T_ea, T_e, T_s2a;
	T_ea = 1.0 / (c1 * a1); /* Eq. [5] */
	T_es = 1.0 / (c2 * a2); /* Eq. [8] */
	T_s2a = pow(c1 * a1 / pow(c2 *a2, a1/a2), ((alpha + 2) / alpha)); /* Eq. [9] */
	T_e = (T_s2a > T_es) ? (T_es) : (T_ea);

	/* round T_e up to the nearest sampling interval */
	T_e = dt * ceil(T_e / dt);

	/* Loop through points in spiral */
	npoints = round(T_e / dt);
	for (i = 0; i < npoints; i++) {
		t = i*dt; /* time at point i (ms) */

		/* Calculate tau for slew & amplitude limited regimes */
		if (t < fmin(T_s2a, T_es)) /* slew limited regime */
			tau = pow(c2 * a2 * t, 1.0/a2);
		else /* amplitude limited regime */
			tau = pow(c1 * a1 * t, 1.0/a1);
		
		/* Calculate sample location at point i */
		kx[i] = lambda*pow(tau,alpha)*cos(omega*tau);
		ky[i] = lambda*pow(tau,alpha)*sin(omega*tau);
	}

	return T_e;
}
