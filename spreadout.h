/* genspiral() function definition */
int genspiral(int N, int itr) {
	fprintf(stderr, "genspiral(): iteration %d, N = %d\n", itr, N);

	/* Declare file */
	FILE *fID;

	/* Declare/initialize variables */
	int n, m;
	float dn;
	float nturns = THETA_accel * (float)opxres / sqrt( (float)ntrains );
	int L = 4 * round( (float)N / (float)nturns / 2.0 );
	float stretch[6];
	int N_stretch;

	/* Declare waveforms */
	float kx[N], ky[N], kz[N];
	float gx[N], gy[N], gz[N];
	float sx[N], sy[N], sz[N];

	/* Declare maximum values */
	float kxy_max = (float)opxres / ((float)opfov /10.0) / 2.0;
	float kz_max = (float)nechoes / ((float)opfov /10.0) / 2.0;
	float gx_max, gy_max, gz_max;
	float sx_max, sy_max, sz_max;

	/* Declare values */
	float R[N], THETA[N], PSI[N];
	float THETA_0 = 0;
	float PSI_0 = 0;

	/* Loop through points */
	for (n = 0; n < N; n++) {
		/* Zero out ramp points */
		if (n < L || N - n < L || (sptype2d > 2 && fabs(N/2 - n) < 1 + nnav/2))
			R[n] = 0.0;
		else { 
			/* Get radius of spiral at sample */
			switch (sptype2d) {
				case 1 : /* spiral out */
					m = fabs(n);
					dn = 1.0 / (float)(N - L);
					R[n] = pow((float)fmod(m - L/2, N) * dn, R_accel);
					break;
				case 2 : /* spiral in */
					m = fabs(n);
					dn = 1.0 / (float)(N - L);
					R[n] = pow(1.0 - (float)fmod(m - L/2, N) * dn, R_accel);
					break;
				case 3 : /* spiral out-in */
					m = fabs(n - N/2);
					dn = 2.0 / (float)(N - 2*L - nnav);
					R[n] = pow(1.0 - (float)fmod(m - nnav/2, N) * dn, R_accel);
					break;
				case 4 : /* spiral in-out */
					m = fabs(n - N/2);
					dn = 2.0 / (float)(N - 2*L - nnav);
					R[n] = pow((float)fmod(m - nnav/2, N) * dn, R_accel);
					break;
				default :
					fprintf(stderr, "Invalid sptype2d: %d", sptype2d);
					return 0;
			}
		}

		/* Get phase of spiral */
		dn = 1.0 / (float)N;
		m = (sptype2d < 3) ? (n - N/2) : (fabs(n - N/2));
		if (sptype3d == 4) { /* for orbital trajectory (3d spiral) */
			THETA[n] = (float)nturns * acos(fmod(m*dn * phi1, 1.0));
			PSI[n] = 2.0 * (float)nturns * M_PI * fmod(m*dn * phi2, 1.0);
		}
		else {
			THETA[n] = M_PI * (float)nturns * (float)m*dn;
			PSI[n] = 0.0;
		}

	};

	/* Smooth waveforms using a sinusoidal convolution */
	sinsmooth(R, (sptype2d < 3) ? (N) : (N/2), L);
	reverse(R, 0, N);
	sinsmooth(R, (sptype2d < 3) ? (N) : (N/2), L);
	reverse(R, 0, N);

	/* Loop through points */
	for (n = 0; n < N; n++) {
		/* Offset phase for 2nd half of spiral for both in and out */
		if (sptype2d > 2 && n >= N/2)
			THETA_0 = M_PI;

		/* Translate polar coordinates to cartesian */
		kx[n] = kxy_max * R[n]*cos(THETA_0 + THETA[n])*cos(PSI_0 + PSI[n]);
		ky[n] = kxy_max * R[n]*sin(THETA_0 + THETA[n])*cos(PSI_0 + PSI[n]);
		if (sptype3d > 1)
			kz[n] = kxy_max * R[n]*sin(THETA_0 + THETA[n])*sin(PSI_0 + PSI[n]);
		else /* SOS case */
			kz[n] = kz_max;
	}
	if (sptype3d == 1) {
		/* Smooth kz ramp for SOS */
		sinsmooth(kz, N, L);
	}

	/* Calculate gradient waveforms by differentiating trajectory */
	diff(kx, N, GRAD_UPDATE_TIME*1e-6*GAMMA/2.0/M_PI, gx);
	diff(ky, N, GRAD_UPDATE_TIME*1e-6*GAMMA/2.0/M_PI, gy);
	diff(kz, N, GRAD_UPDATE_TIME*1e-6*GAMMA/2.0/M_PI, gz);
	gx_max = fabs(getmaxabs(gx, N));
	gy_max = fabs(getmaxabs(gy, N));
	gz_max = fabs(getmaxabs(gz, N));

	/* Calculate slew waveforms by differentiating gradient */
	diff(gx, N, GRAD_UPDATE_TIME*1e-6, sx);
	diff(gy, N, GRAD_UPDATE_TIME*1e-6, sy);
	diff(gz, N, GRAD_UPDATE_TIME*1e-6, sz);
	sx_max = fabs(getmaxabs(sx, N));
	sy_max = fabs(getmaxabs(sy, N));
	sz_max = fabs(getmaxabs(sz, N));

	/* Calculate stretch factors */
	stretch[0] = (float)(gx_max/GMAX);
	stretch[1] = (float)(gy_max/GMAX);
	stretch[2] = (float)(gz_max/GMAX);
	stretch[3] = (float)sqrt(sx_max/SLEWMAX);
	stretch[4] = (float)sqrt(sy_max/SLEWMAX);
	stretch[5] = (float)sqrt(sz_max/SLEWMAX);

	/* Determine new N if function exceeds gradient/slew limits */
	N_stretch = 4 * round(getmaxabs(stretch, 6) * N / 4.0);

	/* Determine if N can be stretched any more */
	if (N_stretch > MAXWAVELEN) {
		fprintf(stderr, "genspiral(): N_stretch = %d > MAXWAVELEN = %d, aborting...\n",
				N_stretch, MAXWAVELEN);
		return 0;
	}
	else if (itr > MAXITR) {
		fprintf(stderr, "genspiral(): itr = %d > MAXITR = %d, aborting...\n",
				itr, MAXITR);
	}
	else if (fabs(N_stretch - N) <= 16) {
		fprintf(stderr, "genspiral(): saving readout gradients with pw_g* = %dus, res_g* = %d\n", GRAD_UPDATE_TIME*N, N);
		fprintf(stderr, "genspiral(): saving readout x gradient with a_gx = %.2f G/cm\n", gx_max);
		fprintf(stderr, "genspiral(): saving readout y gradient with a_gy = %.2f G/cm\n", gy_max);
		fprintf(stderr, "genspiral(): saving readout z gradient with a_gz = %.2f G/cm\n", gz_max);
		grad_len = N;

		/* Write trajectory and gradient to files */
		fID = fopen("./ktraj.txt", "w");
		for (n = 0; n < N; n++) {
			Gx[n] = 2 * round(MAX_PG_WAMP * gx[n] / XGRAD_max / 2.0);
			Gy[n] = 2 * round(MAX_PG_WAMP * gy[n] / YGRAD_max / 2.0);
			Gz[n] = 2 * round(MAX_PG_WAMP * gz[n] / ZGRAD_max / 2.0);
			fprintf(fID, "%f \t%f \t%f\n", kx[n], ky[n], kz[n]);
		}
		fclose(fID);

		return 1;
	}
	else {
		/* Recurse */
		return genspiral(N_stretch, itr + 1);
	}

	return 0;	
};

/* genviews() function definition */
int genviews() {

	/* Declare values and matrices */
	FILE* fID = fopen("kviews.txt","w");
	int trainn, echon, n;
	float rx, ry, rz, dz;
	float Rx[9], Ry[9], Rz[9], Tz[9];
	float T_0[9], T[9];

	/* Initialize z translation to identity matrix */
	eye(Tz, 3);

	/* Get original transformation matrix in float vals */
	for (n = 0; n < 9; n++) T_0[n] = (float)rsprot[0][n] / MAX_PG_WAMP;

	/* Loop through all views */
	for (trainn = 0; trainn < ntrains; trainn++) {
		for (echon = 0; echon < nechoes; echon++) {

			/* Determine type of transformation */
			switch (sptype3d) {
				case 1 : /* Kz shifts */
					rx = 0.0;
					ry = 0.0;
					rz = (float)trainn * 2.0 * M_PI / PHI;
					dz = pow(-1, (float)trainn) * (float)trainn / (float)ntrains;
					dz += pow(-1, (float)echon) * floor((float)(echon + 1) / 2.0);
					dz *= 2.0 / (float)nechoes;
					break;
				case 2 : /* Single axis rotation */
					rx = 0.0;
					ry = (float)(trainn*nechoes + echon) * 2.0 * M_PI / PHI;
					rz = 0.0;
					dz = 1.0;
					break;
				case 3 :
				case 4 : /* Double axis rotations */
					rx = 2.0 * M_PI * (float)(trainn*nechoes + echon) / PHI;
					ry = acos(fmod(1 - 2*(trainn*nechoes + echon + 0.5) / (float)(ntrains*nechoes), 1));
					rz = 0.0;
					dz = 1.0;
					break;
				case 5: /* Debugging case */
					rx = M_PI * (float)(echon) / nechoes;
					ry = M_PI * (float)(trainn) / ntrains;
					rz = 0.0;
					dz = 1.0;
					break;
				default:
					return 0;
			}
			
			/* Calculate the transformation matrices */
			Tz[8] = dz;
			genrotmat('x', rx, Rx);
			genrotmat('y', ry, Ry);
			genrotmat('z', rz, Rz);

			/* Multiply the transformation matrices */
			multmat(3,3,3,T_0,Tz,T);
			multmat(3,3,3,Rx,T,T);
			multmat(3,3,3,Ry,T,T);
			multmat(3,3,3,Rz,T,T);

			/* Save the matrix to the table of matrices */
			fprintf(fID, "%d \t%d \t%f \t%f \t%f \t%f \t", trainn, echon, rx, ry, rz, dz);
			for (n = 0; n < 9; n++) {
				fprintf(fID, "%f \t", T[n]);
				tmtxtbl[trainn*nechoes + echon][n] = (long)round(MAX_PG_IAMP*T[n]);
			}
			fprintf(fID, "\n");
		}
	}

	/* Close the file */
	fclose(fID);

	return 1;
};

/* sinsmooth() function definition */
int sinsmooth(float *x, int N, int L) {

	/* Initialize indexing variables */
	float y[N];
	int n, l;
	float dl = 1.0 / (float)L;

	/* Make smoothing kernel */
	float k[L];
	for (l = 0; l < L; l++)
		k[l] = 0.5 * (cos(2 * M_PI * (float)l*dl - M_PI) + 1);

	/* Make a cut version of x to preserve length */
	int N_cut = N - L;
	float x_cut[N_cut];
	for (n = 0; n < N_cut; n++)
		x_cut[n] = x[n+L];

	/* Convolve waveforms */
	conv(x_cut, N_cut, k, L, y);

	/* Copy result */
	for (n = 0; n < N; n++) x[n] = y[n];

	return 1;
};
