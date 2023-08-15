#include <stdio.h>
#include <math.h>
#include "kimvdsp.h"
#include "gradtrap.h"
#include "helperfuns.h"

#define MAX_PG_WAMP 32767
#define MAXWAVELEN 10000
#define GRAD_UPDATE_TIME 4 
#define XGRAD_max 5
#define YGRAD_max 5
#define ZGRAD_max 5

/* cv's */
int opxres = 64;
float opfov = 240.0;
float spalpha = 1.0;
int ntrains = 1;
int nnav = 20;
int sptype2d = 4;
float SLEWMAX = 17000;
float GMAX = 4;

/* ipg variables */
int Gx[MAXWAVELEN], Gy[MAXWAVELEN], Gz[MAXWAVELEN];
int grad_len;

/* function prototypes */
int genspiral();
int main() {
	
	genspiral();

	return 1;
}

int genspiral() {

	/* 4 parts of the gradient waveform:
	 * 	- Z-encode (*_ze):		kZ encoding step
	 *	- Spiral (*_sp):		kXY-plane spiral trajectory (readout)
	 *	- Ramp-down (*_rd):		XY gradient ramp down
	 *	- kspace Rewinder (*_kr):	kspace rewinder (XYZ trapezoidal gradient)
	 */

	/* gradient waveforms will be in units of rad*G/ms */

	FILE* fID_ktraj = fopen("ktraj.txt", "w");
	FILE* fID_grad = fopen("grad.txt", "w");
	FILE* fID_graddac = fopen("graddac.txt", "w");

	/* declare variables */
	int n;
	float t;
	int np, np_ze, np_sp, np_rd, np_kr;
	float Tr_ze, Tp_ze, T_sp, T_rd, Tr_kr, Tp_kr;
	float t1, t2, t3;
	float gxn, gyn, gzn, kxn, kyn, kzn;
	float gx0, gy0, g0, kx0, ky0, kz0, k0;
	float h_ze, h_kr;

	/* declare waveforms */
	float kx_sp[MAXWAVELEN], ky_sp[MAXWAVELEN] = {0};
	float gx_sp[MAXWAVELEN], gy_sp[MAXWAVELEN] = {0};
	float gx[MAXWAVELEN], gy[MAXWAVELEN], gz[MAXWAVELEN] = {0};

	/* convert units */
	float dt = GRAD_UPDATE_TIME*1e-3; /* raster time (ms) */
	float D = (float)opfov / 10.0; /* fov (cm) */
	float gm = GMAX; /* gradient amplitude limit (G/cm) */
	float sm = SLEWMAX * 1e-3; /* slew limit (G/cm/ms) */
	float gam = 4.258*2*M_PI; /* gyromagnetic ratio (rad/G/ms) */
	float kmax = opxres / D / 2; /* kspace sampling radius (cm^-1) */

	/* generate the z encoding trapezoid gradient */
	gradtrap(kmax, sm, gm, dt, &h_ze, &Tr_ze, &Tp_ze);
	np_ze = round((2*Tr_ze + Tp_ze) / dt);

	/* generate the spiral trajectory */
	T_sp = kimvdsp(kx_sp, ky_sp, D, sm, gm, opxres, (sptype2d < 3) ? (ntrains) : (2*ntrains), spalpha, dt);
	np_sp = round(T_sp/dt);
	diff(kx_sp, np_sp, gam*dt, gx_sp);
	diff(ky_sp, np_sp, gam*dt, gy_sp);	

	/* calculate gradients at end of spiral */
	gx0 = gx_sp[np_sp - 1];
	gy0 = gy_sp[np_sp - 1];
	g0 = sqrt(pow(gx0,2) + pow(gy0,2));

	/* calculate gradient ramp down time and round up to nearest sampling interval */
	T_rd = g0 / sm;
	T_rd = dt * ceil(T_rd / dt);
	np_rd = round(T_rd/dt);

	/* calculate gradients at end of ramp down */
	kx0 = kx_sp[np_sp - 2] + gam * 1/2 * (T_rd + dt) * gx0;
	ky0 = ky_sp[np_sp - 2] + gam * 1/2 * (T_rd + dt) * gy0;
	kz0 = kmax;
	k0 = sqrt(pow(kx0,2) + pow(ky0,2) + pow(kz0,2));

	/* generate the kspace rewinder */
	gradtrap(k0, sm, gm, dt, &h_kr, &Tr_kr, &Tp_kr); 
	np_kr = round((2*Tr_kr + Tp_kr) / dt);

	/* calculate time markers */
	t1 = 2*Tr_ze + Tp_ze;
	t2 = t1 + T_sp;
	t3 = t2 + T_rd;

	/* calculate total number of points */
	np = np_ze + np_sp + np_rd + np_kr + 1;

	/* loop through time points */
	for (n = 0; n < np; n++) {
		t = dt * n;

		if (t <= t1) { /* Z-encode gradient */
			gxn = 0;
			gyn = 0;
			gzn = h_ze * trap(t, 0, Tr_ze, Tp_ze);
		}

		else if (t <= t2) { /* Spiral trajectory */
			gxn = gx_sp[n - np_ze - 1];
			gyn = gy_sp[n - np_ze - 1];
			gzn = 0;
		}

		else if (t <= t3) { /* Gradient ramp-down */
			gxn = gx0 * (1 - (t - t2) / T_rd);
			gyn = gy0 * (1 - (t - t2) / T_rd);
			gzn = 0;
		}

		else { /* Kspace rewinder */
			gxn = -kx0/k0 * h_kr * trap(t, t3, Tr_kr, Tp_kr);
			gyn = -ky0/k0 * h_kr * trap(t, t3, Tr_kr, Tp_kr);
			gzn = -kz0/k0 * h_kr * trap(t, t3, Tr_kr, Tp_kr);
		}


		switch (sptype2d) {
			case 1: /* spiral out */
				gx[n] = gxn;
				gy[n] = gyn;
				gz[n] = gzn;
				break;
			case 2: /* spiral in */
				gx[np - 1 - n] = gxn;
				gy[np - 1 - n] = gyn;
				gz[np - 1 - n] = gzn;
				break;
			case 3: /* spiral out-in */
				gx[n] = gxn;
				gy[n] = gyn;
				gz[n] = gzn;
				gx[2*np + nnav - 1 - n] = gxn;
				gy[2*np + nnav - 1 - n] = gyn;
				gz[2*np + nnav - 1 - n] = -gzn;
				break;
			case 4: /* spiral in-out */
				gx[np - 1 - n] = gxn;
				gy[np - 1 - n] = gyn;
				gz[np - 1 - n] = -gzn;
				gx[np + nnav + n] = gxn;
				gy[np + nnav + n] = gyn;
				gz[np + nnav + n] = gzn;
				break;
			default:
				return 0;
				break;
		}
	}

	/* calculate total number of points */
	if (sptype2d > 2)
		grad_len = 2*np + nnav;
	else
		grad_len = np;

	/* calculate kspace location, fs gradients, and write to file */
	kxn = 0.0;
	kyn = 0.0;
	kzn = 0.0;
	for (n = 0; n < grad_len; n++) {
		kxn += gam * gx[n] * dt;
		kyn += gam * gy[n] * dt;
		kzn += gam * gz[n] * dt;

		Gx[n] = 2*round(MAX_PG_WAMP/XGRAD_max * gx[n]/(2.0*M_PI) / 2.0);
		Gy[n] = 2*round(MAX_PG_WAMP/YGRAD_max * gy[n]/(2.0*M_PI) / 2.0);
		Gz[n] = 2*round(MAX_PG_WAMP/ZGRAD_max * gz[n]/(2.0*M_PI) / 2.0);
		
		fprintf(fID_ktraj, "%f \t%f \t%f\n", kxn, kyn, kzn);
		fprintf(fID_grad, "%f \t%f \t%f\n", gx[n], gy[n], gz[n]);
		fprintf(fID_graddac, "%d \t%d \t%d\n", Gx[n], Gy[n], Gz[n]);
	}

	fclose(fID_ktraj);
	fclose(fID_grad);
	fclose(fID_graddac);

	return 1;
}
