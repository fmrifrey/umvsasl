int readprep(int id, int *len,
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl)
{

	/* Declare variables */
	char fname[80];
	FILE *fID;
	char buff[MAXWAVELEN];
	int i, tmplen;
	double lblval, ctlval;
	
	if (id == 0) {
		/* Set all values to zero and return */
		for (i = 0; i < *len; i++) {
			rho_lbl[i] = 0;
			theta_lbl[i] = 0;
			grad_lbl[i] = 0;
			rho_ctl[i] = 0;
			theta_ctl[i] = 0;
			grad_ctl[i] = 0;
		}
		return 1;
	}

	/* Read in RF magnitude from rho file */
	sprintf(fname, "./aslprep/pulses/%05d/rho.txt", id);
	fprintf(stderr, "readprep(): opening %s...\n", fname);
	fID = fopen(fname, "r");

	/* Check if rho file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening %s\n", fname);
		return 0;
	}

	/* Loop through points in rho file */
	i = 0;
	while (fgets(buff, MAXWAVELEN, fID)) {
		sscanf(buff, "%lf %lf", &lblval, &ctlval);
		rho_lbl[i] = (int)lblval;
		rho_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);
	tmplen = i;
	
	/* Read in RF phase from theta file */
	sprintf(fname, "./aslprep/pulses/%05d/theta.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening %s\n", fname);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fgets(buff, MAXWAVELEN, fID)) {
		sscanf(buff, "%lf %lf", &lblval, &ctlval);
		theta_lbl[i] = (int)lblval;
		theta_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (tmplen != i) {
		fprintf(stderr, "readprep(): length of theta file (%d) is not consistent with rho file length (%d)\n", i, tmplen);
		return 0;
	}
	
	/* Read in RF phase from theta file */
	sprintf(fname, "./aslprep/pulses/%05d/grad.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening %s\n", fname);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fgets(buff, MAXWAVELEN, fID)) {
		sscanf(buff, "%lf %lf", &lblval, &ctlval);
		grad_lbl[i] = (int)lblval;
		grad_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (tmplen != i) {
		fprintf(stderr, "readprep(): length of grad file (%d) is not consistent with rho/theta file length (%d)\n", i, tmplen);
		return 0;
	}
	
	*len = tmplen;

	return 1;
}

int readschedule(int id, int* var, char* varname, int lines) {

	FILE* fID;
	char fname[200];
	int val;

	/* Open the schedule file */
	sprintf(fname, "./aslprep/schedules/%05d/%s.txt", id, varname);
	fprintf(stderr, "readschedule(): opening %s...\n", fname);
	fID = fopen(fname, "r");
	if (fID == 0) {
		fprintf(stderr, "File not found.\n");
		return 0;
	}

	/* Read in the array */
	int i = 0;
	while (fscanf(fID, "%d\n", &val) != EOF) {
		var[i] = val;
		i++;
	}
	fclose(fID);

	/* If only 1 line is read in (scalar --> array) */
	if (i == 1) {
		for (i = 1; i < lines; i++)
			var[i] = var[0];
	}

	/* If number of lines is less than number of lines to read */
	if (i < lines - 1)
		return -1;

	return 1;
}

int calctadjust() {
	
	int framen;

	/* Calculate adjust times */
	for (framen = 0; framen < nframes; framen++) {
		avmintr = dur_tipdowncore + nechoes * (dur_refocuscore + dur_seqcore);
		avmintr += dur_fatsatcore;
		avmintr += dur_blksatcore;
		avmintr += (prep1_id > 0) ? (dur_prep1core + prep1_pldtbl[framen]) : (0);
		avmintr += (prep2_id > 0) ? (dur_prep2core + prep2_pldtbl[framen]) : (0);
		if (optr < avmintr) {
			epic_error(use_ermes, "optr must be >= %dus", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), INT_ARG, avmintr);
			return FAILURE;
		}
		else
			tadjusttbl[framen] = optr - avmintr;
	}

	return 1;
}

int genschedule(int mod, int pld, int* lbltbl, int* pldtbl)
{
	int framen;

	/* Loop through frames */
	for (framen = 0; framen < nframes; framen++) {
		/* Set PLD */
		pldtbl[framen] = pld;

		/* Set labeling scheme */
		if (framen < nm0frames)
			lbltbl[framen] = -1;
		else {
			switch (mod) {
				case 1: /* Label, control... */
					lbltbl[framen] = (framen - nm0frames + 1) % 2; /* 1, 0, 1, 0 */
					break;
				case 2: /* Control, label... */
					lbltbl[framen] = (framen - nm0frames) % 2; /* 0, 1, 0, 1 */
					break;
				case 3: /* Label */
					lbltbl[framen] = 1;
					break;
				case 4: /* Control */
					lbltbl[framen] = 0;
					break;
			}
		}
	}	

	return 1;
}
