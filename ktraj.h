/* genspiral() function definition */
int genspiral(int N, int itr) {
	fprintf(stderr, "genspiral(): iteration %d, N = %d\n", itr, N);

	/* Declare files */
	FILE *fID_ktraj;
	FILE *fID_grad;

	/* Declare/initialize variables */
	int n, m;
	float dn;
	float nturns = THETA_accel * (float)opxres / sqrt( (float)ntrains );
	int L = 4 * round( (float)N / (float)nturns / 2.0 );

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
	diff(kx, N, TSP_GRAD*1e-6*GAMMA/2.0/M_PI, gx);
	diff(ky, N, TSP_GRAD*1e-6*GAMMA/2.0/M_PI, gy);
	diff(kz, N, TSP_GRAD*1e-6*GAMMA/2.0/M_PI, gz);
	gx_max = fabs(getmaxabs(gx, N));
	gy_max = fabs(getmaxabs(gy, N));
	gz_max = fabs(getmaxabs(gz, N));

	/* Calculate slew waveforms by differentiating gradient */
	diff(gx, N, TSP_GRAD*1e-6, sx);
	diff(gy, N, TSP_GRAD*1e-6, sy);
	diff(gz, N, TSP_GRAD*1e-6, sz);
	sx_max = fabs(getmaxabs(sx, N));
	sy_max = fabs(getmaxabs(sy, N));
	sz_max = fabs(getmaxabs(sz, N));

	/* Calculate stretch factors */
	float sf_gx = (float)(gx_max/GMAX);
	float sf_gy = (float)(gy_max/GMAX);
	float sf_gz = (float)(gz_max/GMAX);
	float sf_sx = (float)sqrt(sx_max/SLEWMAX);
	float sf_sy = (float)sqrt(sy_max/SLEWMAX);
	float sf_sz = (float)sqrt(sz_max/SLEWMAX);

	/* Determine if function exceeds gradient/slew limits */
	float stretch[6] = {sf_gx, sf_gy, sf_gz, sf_sx, sf_sy, sf_sz};
	int N_stretch = 4 * round(getmaxabs(stretch, 6) * N / 4.0);

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
		/* Save values */
		res_gx = N;
		res_gy = N;
		res_gz = N;
		pw_gx = N*TSP_GRAD;
		pw_gy = N*TSP_GRAD;
		pw_gz = N*TSP_GRAD;
		a_gx = gx_max;
		a_gy = gy_max;
		a_gz = gz_max;
		ia_gx = MAX_PG_WAMP;
		ia_gy = MAX_PG_WAMP;
		ia_gz = MAX_PG_WAMP;
		fprintf(stderr, "genspiral(): saving readout x gradient with a_gx = %.2f G/cm, pw_gx = %d us, res_gx = %d\n",
			a_gx, pw_gx, res_gx);
		fprintf(stderr, "genspiral(): saving readout y gradient with a_gy = %.2f G/cm, pw_gy = %d us, res_gy = %d\n",
			a_gy, pw_gy, res_gy);
		fprintf(stderr, "genspiral(): saving readout z gradient with a_gz = %.2f G/cm, pw_gz = %d us, res_gz = %d\n",
			a_gz, pw_gz, res_gz);

		/* Write trajectory and gradient to files */
		fID_ktraj = fopen("./ktraj.txt", "w");
		fID_grad = fopen("./grad.txt","w");
		for (n = 0; n < N; n++) {
			Gx[n] = 2 * round(MAX_PG_WAMP / a_gx * gx[n] / 2.0);
			Gy[n] = 2 * round(MAX_PG_WAMP / a_gy * gy[n] / 2.0);
			Gz[n] = 2 * round(MAX_PG_WAMP / a_gz * gz[n] / 2.0);
			fprintf(fID_ktraj, "%f \t%f \t%f\n", kx[n], ky[n], kz[n]);
			fprintf(fID_grad, "%d \t%d \t%d\n", Gx[n], Gy[n], Gz[n]);
		}
		fclose(fID_ktraj);
		fclose(fID_grad);

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
	eye(T_0,3);
	eye(Tz, 3);

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
					ry = acos(fmod(1.0 - 2.0*(float)(trainn*nechoes + echon + 0.5), 1.0));
					ry *= 1.0 / (float)(ntrains * nechoes);
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
			
			/* Save the matrix to file */
			printmat(fID,3,3,T);
			fprintf(fID,"\n");

			/* Save the matrix to the table of matrices */
			for (n = 0; n < 9; n++)
				T_v[trainn*nechoes + echon][n] = (int)round(32767.0*T[n]);
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
