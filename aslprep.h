int readprep(int id, int len
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl)
{

	/* Declare variables */
	char fname[80];
	FILE *fID;
	int i;
	int lblval;
	int ctlval;
	
	/* Read in RF magnitude from rho file */
	sprintf(fname, "/usr/g/bin/aslprep/prep%d.rho.txt", id);
	fID = fopen(fname, "r");

	/* Check if rho file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.rho.txt\n", id);
		return 0;
	}

	/* Loop through points in rho file */
	i = 0;
	while (fscanf(fID, "%d \t%d\n", &lblval, &ctlval)) {
		rho_lbl[i] = lblval;
		rho_ctl[i++] = ctlval;
	}
	fclose(fID);
	len = i;
	
	/* Read in RF phase from theta file */
	sprintf(fname, "/usr/g/bin/aslprep/prep%d.theta.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.theta.txt\n", id);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fscanf(fID, "%d \t%d\n", &lblval, &ctlval)) {
		theta_lbl[i] = lblval;
		theta_ctl[i++] = ctlval;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (len != i) {
		fprintf(stderr, "readprep(): length of theta file (%d) is not consistent with rho file length (%d)\n", i, len);
		return 0;
	}
	
	/* Read in RF phase from theta file */
	sprintf(fname, "/usr/g/bin/aslprep/prep%d.grad.txt", id);
	fID = fopen(fname, "r");

	/* Check if theta file was read successfully */
	if (fID == 0) {
		fprintf(stderr, "readprep(): failure opening prep%d.grad.txt\n", id);
		return 0;
	}

	/* Loop through points in theta file */
	i = 0;
	while (fscanf(fID, "%d \t%d\n", &lblval, &ctlval)) {
		grad_lbl[i] = lblval;
		grad_ctl[i++] = ctlval;
	}
	fclose(fID);

	/* Check that length is consistent */
	if (len != i) {
		fprintf(stderr, "readprep(): length of grad file (%d) is not consistent with rho/theta file length (%d)\n", i, len);
		return 0;
	}

	return 1;
}
