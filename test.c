#include <stdio.h>
#include <math.h>
#include "kimvdsp.h"
#include "gradtrap.h"
#include "helperfuns.h"

#define SLEWMAX 17000
#define GMAX 4

#define MAXWAVELEN 10000
#define GRAD_UPDATE_TIME 4

FILE* fID;

/* cv's */
int opxres = 64;
float opfov = 240.0;
float sp_alpha = 1.0;
int ntrains = 1;
int nnav = 20;
int sptype2d = 4;

/* function prototypes */
int genspiral();
int gradtrap(float dk, float sm, float gm, float dt, float* h, float* t_ramp, float* t_plat);

int main() {
/*
	float tr, tp, h;
	float dkx, dky, dkz, dk;
	float sm = SLEWMAX*1e-3;
	float gm = GMAX;
	float gam = 4.258*2*M_PI; */ /* gyromagnetic ratio (rad/G/ms) */
/*
	dkx = -1.0;
	dky = -0.75;
	dkz = 0.5;
	dk = sqrt(pow(dkx,2) + pow(dky,2) + pow(dkz,2));

	gradtrap(dk, sm, gm, GRAD_UPDATE_TIME*1e-3, &h, &tr, &tp);
	fprintf(stderr, "h = %f, tr = %f, tp = %f\n", h, tr, tp);

	float A = (tp + tr) / (GRAD_UPDATE_TIME*1e-3) * h;
	fprintf(stderr, "A = %f\n", A);
*/
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

	fID = fopen("ktraj.txt", "w");

	/* declare variables */
	int n;
	float t;
	int np, np_ze, np_sp, np_rd, np_kr;
	float Tr_ze, Tp_ze, T_sp, T_rd, Tr_kr, Tp_kr;
	float t1, t2, t3, t4;
	float gxn, gyn, gzn, kxn, kyn, kzn;
	float gx0, gy0, g0, kx0, ky0, kz0, k0;
	float h_ze, h_kr;

	/* declare waveforms */
	float gz_ze[MAXWAVELEN] = {0};
	float kx_sp[MAXWAVELEN], ky_sp[MAXWAVELEN] = {0};
	float gx_sp[MAXWAVELEN], gy_sp[MAXWAVELEN] = {0};
	float gx_rd[MAXWAVELEN], gy_rd[MAXWAVELEN] = {0};
	float gx_kr[MAXWAVELEN], gy_kr[MAXWAVELEN], gz_kr[MAXWAVELEN] = {0};
	float gx[MAXWAVELEN], gy[MAXWAVELEN], gz[MAXWAVELEN] = {0};

	/* convert units */
	float dt = GRAD_UPDATE_TIME*1e-3; /* raster time (ms) */
	float D = (float)opfov / 10; /* fov (cm) */
	float gm = GMAX; /* gradient amplitude limit (G/cm) */
	float sm = SLEWMAX * 1e-3; /* slew limit (G/cm/ms) */
	float gam = 4.258*2*M_PI; /* gyromagnetic ratio (rad/G/ms) */
	float kmax = opxres / D / 2; /* kspace sampling radius (cm^-1) */

	/* generate the z encoding trapezoid gradient */
	gradtrap(kmax, sm, gm, dt, &h_ze, &Tr_ze, &Tp_ze);
	np_ze = round((2*Tr_ze + Tp_ze) / dt);

	/* generate the spiral trajectory */
	T_sp = kimvdsp(kx_sp, ky_sp, D, sm, gm, opxres, (sptype2d < 3) ? (ntrains) : (2*ntrains), sp_alpha, dt);
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
	fprintf(stderr, "k0 = [%f, %f, %f]\n", kx0,ky0,kz0);

	/* generate the kspace rewinder */
	gradtrap(k0, sm, gm, dt, &h_kr, &Tr_kr, &Tp_kr); 
	np_kr = round((2*Tr_kr + Tp_kr) / dt);

	/* calculate time markers */
	t1 = 2*Tr_ze + Tp_ze;
	t2 = t1 + T_sp;
	t3 = t2 + T_rd;
	t4 = t3 + 2*Tr_kr + Tp_kr;

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
		}
	}

	/* correct the number of points */
	if (sptype2d > 2)
		np = 2*np + nnav;

	/* calculate kspace location and write to file */
	for (n = 0; n < np; n++) {
		kxn += gam * gx[n] * dt;
		kyn += gam * gy[n] * dt;
		kzn += gam * gz[n] * dt;
		fprintf(fID, "%f \t%f \t%f\n", kxn, kyn, kzn);
		if (n == np_sp + np_ze - 1)
			fprintf(stderr, "at pt %d: kx = %f, ky = %f, kz = %f\n", n,kxn,kyn,kzn);
		if (n == np_sp + np_ze + np_rd - 1)
			fprintf(stderr, "at pt %d: kx = %f, ky = %f, kz = %f\n", n,kxn,kyn,kzn);
	}

	fclose(fID);

	return 1;
}
