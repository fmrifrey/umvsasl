int readprep(int id, int *len,
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl)
{

	/* Declare variables */
	char fname[80];
	FILE *fID;
	char buff[2];
	int i;
	float lblval, ctlval;
	
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
	sprintf(fname, "./aslprep/prep%d.rho.txt", id);
	fID = fopen(fname, "r");

	/* Check if rho file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.rho.txt\n", id);
		return 0;
	}

	/* Loop through points in rho file */
	i = 0;
	while (fgets(buff, 2, fID)) {
		sscanf(buff, "%1f %1f", &lblval, &ctlval);
		fprintf(stderr, "%d: %d \t%d\n", i, (int)lblval, (int)ctlval);
		rho_lbl[i] = (int)lblval;
		rho_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);
	*len = i;
	
	/* Read in RF phase from theta file */
	sprintf(fname, "./aslprep/prep%d.theta.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.theta.txt\n", id);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fgets(buff, 2, fID)) {
		sscanf(buff, "%1f %1f", &lblval, &ctlval);
		theta_lbl[i] = (int)lblval;
		theta_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (*len != i) {
		fprintf(stderr, "readprep(): length of theta file (%d) is not consistent with rho file length (%d)\n", i, *len);
		return 0;
	}
	
	/* Read in RF phase from theta file */
	sprintf(fname, "./aslprep/prep%d.grad.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.grad.txt\n", id);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fgets(buff, 2, fID)) {
		sscanf(buff, "%1f %1f", &lblval, &ctlval);
		grad_lbl[i] = (int)lblval;
		grad_ctl[i] = (int)ctlval;
		i++;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (*len != i) {
		fprintf(stderr, "readprep(): length of grad file (%d) is not consistent with rho/theta file length (%d)\n", i, *len);
		return 0;
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
