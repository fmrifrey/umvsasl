/*
 * GE Medical Systems
 * Copyright (C) 1996-2003 The General Electric Company
 *
 * File Name : asl3dflex.e
 * Language  : EPIC/ANSI C
 * Date      : 01-Jan-1996
 *
 * An ASL-prepped flexible spiral fse readout sequence (ASL3DFLEX),
 * built up from grass.e
 */

@inline epic.h
@inline intwave.h

@global
/*********************************************************************
 *                  ASL3DFLEX.E GLOBAL SECTION                       *
 *                                                                   *
 * Common code shared between the Host and IPG PSD processes.  This  *
 * section contains all the #define's, global variables and function *
 * declarations (prototypes).                                        *
 *********************************************************************/
#include <stdio.h>
#include <string.h>

#include "em_psd_ermes.in"
#include "grad_rf_asl3dflex.globals.h"

#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "epic_error.h"
#include "epicfuns.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#ifdef psdutil
#include "psdutil.h"
#endif
#include "psd_proto.h"
#include "epic_iopt_util.h"
#include "filter.h"

#include "asl3dflex.h"

/* Define important values */
#define MAXWAVELEN 50000 /* Maximum wave length for gradients */
#define MAXNTRAINS 50 /* Maximum number of echo trains per frame */
#define MAXNECHOES 50 /* Maximum number of echoes per echo train */
#define MAXITR 50 /* Maximum number of iterations for iterative processes */
#define GAMMA 26754 /* Gyromagnetic ratio */
#define TIMESSI 400 /* SSP instruction time */

@inline Prescan.e PSglobal
int debugstate = 1;

@ipgexport
/*********************************************************************
 *                ASL3DFLEX.E IPGEXPORT SECTION                      *
 *                                                                   *
 * Standard C variables of _any_ type common for both the Host and   *
 * IPG PSD processes. Declare here all the complex type, e.g.,       *
 * structures, arrays, files, etc.                                   *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
@inline Prescan.e PSipgexport
RF_PULSE_INFO rfpulseInfo[RF_FREE] = { {0,0} };

/* Define temporary error message string */
char tmpstr[200];

/* Declare sequencer hardware limit variables */
float XGRAD_max;
float YGRAD_max;
float ZGRAD_max;
float RHO_max;
float THETA_max;

/* Declare readout gradient waveform arrays */
int Gx[MAXWAVELEN];
int Gy[MAXWAVELEN];
int Gz[MAXWAVELEN];
int grad_len = 5000;

/* Declare table of readout gradient transformation matrices */
long tmtxtbl[MAXNTRAINS*MAXNECHOES][9];

/* Declare ASL prep waveform arrays */
int prep1_rho_lbl[MAXWAVELEN];
int prep1_theta_lbl[MAXWAVELEN];
int prep1_grad_lbl[MAXWAVELEN];

int prep1_rho_ctl[MAXWAVELEN];
int prep1_theta_ctl[MAXWAVELEN];
int prep1_grad_ctl[MAXWAVELEN];

int prep1_len = 5000;

int prep2_rho_lbl[MAXWAVELEN];
int prep2_theta_lbl[MAXWAVELEN];
int prep2_grad_lbl[MAXWAVELEN];

int prep2_rho_ctl[MAXWAVELEN];
int prep2_theta_ctl[MAXWAVELEN];
int prep2_grad_ctl[MAXWAVELEN];

int prep2_len = 5000;

/* Declare tables of asl prep pulse labeling schemes */
int prep1_lbltbl[MAXWAVELEN];
int prep1_pldtbl[MAXWAVELEN];

int prep2_lbltbl[MAXWAVELEN];
int prep2_pldtbl[MAXWAVELEN];

int tadjusttbl[MAXWAVELEN];

/* Declare core duration variables */
int dur_blksatcore;
int dur_prep1core;
int dur_prep2core;
int dur_fatsatcore;
int dur_tipdowncore;
int dur_refocuscore;
int dur_seqcore;

@cv
/*********************************************************************
 *                     ASL3DFLEX.E CV SECTION                        *
 *                                                                   *
 * Standard C variables of _limited_ types common for both the Host  *
 * and IPG PSD processes. Declare here all the simple types, e.g,    *
 * int, float, and C structures containing the min and max values,   *
 * and ID description, etc.                                          *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
@inline loadrheader.e rheadercv
@inline vmx.e SysCVs

@inline Prescan.e PScvs

int numdda = 4;			/* For Prescan: # of disdaqs ps2*/

float xmtaddScan;
int obl_debug = 0 with {0, 1, 0, INVIS, "On(=1) to print messages for obloptimize",};
int obl_method = 0 with {0, 1, 0, INVIS, "On(=1) to optimize the targets based on actual rotation matrices",};
int debug = 0 with {0,1,0,INVIS,"1 if debug is on ",};
float echo1bw = 16 with {,,,INVIS,"Echo1 filter bw.in KHz",};

/* FSE timing cvs */
int grad_buff_time = 248 with {100, , 248, INVIS, "Gradient IPG buffer time (us)",};
int trap_ramp_time = 248 with {100, , 248, INVIS, "Trapezoidal gradient ramp time (us)",};
int dolongrf = 0 with {0, 1, 0, VIS, "Option to do long (4 cycle, 6400ms rf pulses)",};
int dophasecycle = 0 with {0, 1, 0, VIS, "Option to do CPMG phase cycling (180, -180, 180...)",};
float varflipfac = 1 with {0, 1, 0, VIS, "Scaling factor for variable flip angle schedule (1 = constant fa)",};

/* Trajectory cvs */
int nechoes = 16 with {1, MAXNECHOES, 17, VIS, "Number of echoes per echo train",};
int ntrains = 1 with {1, MAXNTRAINS, 1, VIS, "Number of echo trains per frame",};
int nframes = 2 with {1, , 2, VIS, "Number of frames to acquire",};
int nnav = 20 with {0, 1000, 20, VIS, "Number of navigator points (must be even)",};
float R_accel = 0.5 with {0.05, , , VIS, "Spiral radial acceleration factor",};
float THETA_accel = 1.0 with {0, , 1, VIS, "Spiral angular acceleration factor",};
int sptype2d = 4 with {1, 4, 1, VIS, "1 = spiral out, 2 = spiral in, 3 = spiral out-in, 4 = spiral in-out",};
int sptype3d = 3 with {1, 5, 1, VIS, "1 = stack of spirals, 2 = rotating spirals (single axis), 3 = rotating spirals (2 axes), 4 = rotating orbitals (2 axes), 5 = debugging mode",};
float SLEWMAX = 17000.0 with {1000, 25000.0, 17000.0, VIS, "Maximum allowed slew rate (G/cm/s)",};
float GMAX = 4.0 with {0.5, 5.0, 4.0, VIS, "Maximum allowed gradient (G/cm)",};
int kill_grads = 0 with {0, 1, 0, VIS, "Option to turn off readout gradients",};

/* ASL prep pulse cvs */
int nm0frames = 2 with {0, , 2, VIS, "Number of M0 frames (no prep pulses are played)",};
int schedule_id = 0 with {0, , 0, VIS, "ASL labeling schedule (0 = no external schedule)",};
int doblksat = 1 with {0, 1, 1, VIS, "Option to do bulk spin saturation at end of each readout",};

int prep1_id = 0 with {0, , 0, VIS, "ASL prep pulse 1: ID number (0 = no pulse)",};
int prep1_pld = 0 with {0, , 0, VIS, "ASL prep pulse 1: post-labeling delay (us; includes background suppression)",};
int prep1_ncycles = 1 with {1, , 1, VIS, "ASL prep pulse 1: number of cycles",};
int prep1_rfmax = 234 with {0, , 0, VIS, "ASL prep pulse 1: maximum RF amplitude",};
int prep1_gmax = 3 with {0, , 3, VIS, "ASL prep pulse 1: maximum gradient amplitude",};
int prep1_mod = 1 with {1, 4, 1, VIS, "ASL prep pulse 1: labeling modulation scheme (1 = label/control, 2 = control/label, 3 = always label, 4 = always control)",};
int prep1_tbgs1 = 0 with {0, , 0, VIS, "ASL prep pulse 1: 1st background suppression delay (0 = no pulse)",};
int prep1_tbgs2 = 0 with {0, , 0, VIS, "ASL prep pulse 1: 2nd background suppression delay (0 = no pulse)",};

int prep2_id = 0 with {0, , 0, VIS, "ASL prep pulse 2: ID number (0 = no pulse)",};
int prep2_pld = 0 with {0, , 0, VIS, "ASL prep pulse 2: post-labeling delay (us; includes background suppression)",};
int prep2_ncycles = 1 with {1, , 1, VIS, "ASL prep pulse 2: number of cycles",};
int prep2_rfmax = 234 with {0, , 0, VIS, "ASL prep pulse 2: maximum RF amplitude",};
int prep2_gmax = 3 with {0, , 3, VIS, "ASL prep pulse 2: maximum gradient amplitude",};
int prep2_mod = 1 with {1, 4, 1, VIS, "ASL prep pulse 2: labeling modulation scheme (1 = label/control, 2 = control/label, 3 = always label, 4 = always control)",};
int prep2_tbgs1 = 0 with {0, , 0, VIS, "ASL prep pulse 2: 1st background suppression delay (0 = no pulse)",};
int prep2_tbgs2 = 0 with {0, , 0, VIS, "ASL prep pulse 2: 2nd background suppression delay (0 = no pulse)",};

@host
/*********************************************************************
 *                    ASL3DFLEX.E HOST SECTION                       *
 *                                                                   *
 * Write here the code unique to the Host PSD process. The following *
 * functions must be declared here: cvinit(), cveval(), cvcheck(),   *
 * and predownload().                                                *
 *                                                                   *
 *********************************************************************/
#include <math.h>
#include <stdlib.h>
#include "grad_rf_asl3dflex.h"
#include "psdopt.h"
#include "sar_pm.h"
#include "support_func.host.h"
#include "helperfuns.h"

/* fec : Field strength dependency library */
#include <sysDep.h>
#include <sysDepSupport.h>      /* FEC : fieldStrength dependency libraries */

@inline loadrheader.e rheaderhost

/** Load PSD Header **/
abstract("asl3dflex sequence");
psdname("asl3dflex");

int num_conc_grad = 3;          /* always three for grass 	*/
int entry;

/* peak B1 amplitudes */
float maxB1[MAX_ENTRY_POINTS], maxB1Seq;

/* This will point to a structure defining parameters of the filter
   used for the 1st echo */
FILTER_INFO *echo1_filt; 

/* Use real time filters, so allocate space for them instead of trying
   to point to an infinite number of structures in filter.h. */
FILTER_INFO echo1_rtfilt;

/* Golden ratio numbers */
float PHI = (1.0 + sqrt(5.0)) / 2.0; /* 1d golden ratio */
float phi1 = 0.4656; /* 2d golden ratio 1 */
float phi2 = 0.6823; /* 2d golden ratio 2 */

/* Declare function prototypes from spreadout.h */
int genspiral(int N, int itr);
int genviews();
int sinsmooth(float *x, int N, int L);

/* Declare function prototypes from aslprep.h */
int readprep(int id, int *len,
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl); 
int readschedule(int id, int* var, char* varname, int lines);
int calctadjust();
int genschedule(int mod, int pld, int* lbltbl, int* pldtbl);

/* Import functions from spreadout.h and aslprep.h (using @inline instead of #include since
 * functions reference global variables in those files)
 */
@inline spreadout.h
@inline aslprep.h

@inline Prescan.e PShostVars            /* added with new filter calcs */

static char supfailfmt[] = "Support routine %s failed";


/************************************************************************/
/*       			CVINIT    				*/
/* Invoked once (& only once) when the PSD host process	is started up.	*/
/* Code which is independent of any OPIO button operation is put here.	*/
/************************************************************************/
STATUS cvinit( void )
{

	/* turn off bandwidth option */
	cvdef(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	cvmin(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	cvmax(oprbw, 500.0 / (float)GRAD_UPDATE_TIME);
	oprbw = 500.0 / (float)GRAD_UPDATE_TIME;
	pircbnub = 0;

	/* fov */
	opfov = 240;
	pifovnub = 5;
	pifovval2 = 200;
	pifovval3 = 220;
	pifovval4 = 240;
	pifovval5 = 260;
	pifovval6 = 280;

	/* tr */
	cvdef(optr, 4500);
	optr = 4500ms;
	pitrnub = 2;
	pitrval2 = PSD_MINIMUMTR;
	pitrval3 = 4500ms;

	/* te */
	cvdef(opte, 50);
	cvmin(opte, avminte);
	opte = 50ms;
	pitrnub = 2;
	pite1val2 = PSD_MINFULLTE;
	pite1val3 = 100ms;

	/* frequency (xres) */
	cvmin(opxres, 32);
	cvmax(opxres, 128);
	cvdef(opxres, 64);
	opxres = 64;
	pixresnub = 0; /* hide option */

	/* flip angle */
	cvmin(opflip, 50);
	cvmax(opflip, 180);
	cvdef(opflip, 180);
	opflip = 180;
	pifanub = 5;
	pifaval2 = 180;
	pifaval3 = 160;
	pifaval4 = 140;
	pifaval5 = 120;
	pifaval6 = 100;

	/* hide phase (yres) option */
	piyresnub = 0;

	/* Hide inversion time */
	pitinub = 0;

	/* hide second bandwidth option */
	pircb2nub = 0;

	/* hide nex stuff */
	piechnub = 0;
	pinexnub = 0;

#ifdef ERMES_DEBUG
	use_ermes = 0;
#else /* !ERMES_DEBUG */
	use_ermes = 1;
#endif /* ERMES_DEBUG */

	configSystem();
	EpicConf();
	inittargets(&loggrd, &phygrd);

	/* Init filter slots */
	initfilter();
	
	if (_psd_rf_wait.fixedflag == 0)  { /* sets psd_grd_wait and psd_rf_wait */
		if (setsysparms() == FAILURE)  {
			epic_error(use_ermes,"Support routine setsysparams failed",
					EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms");
			return FAILURE;
		}
	}

	if( obloptimize( &loggrd, &phygrd, scan_info, exist(opslquant),
				exist(opplane), exist(opcoax), obl_method, obl_debug,
				&opnewgeo, cfsrmode ) == FAILURE )
	{
		return FAILURE;
	}
	
	/* Get sequencer hardware limits */
	gettarget(&XGRAD_max, XGRAD, &loggrd);
	gettarget(&YGRAD_max, YGRAD, &loggrd);
	gettarget(&ZGRAD_max, ZGRAD, &loggrd);
	gettarget(&RHO_max, RHO, &loggrd);
	gettarget(&THETA_max, THETA, &loggrd);


@inline Prescan.e PScvinit

#include "cvinit.in"	/* Runs the code generated by macros in preproc.*/

	return SUCCESS;
}   /* end cvinit() */

@inline InitAdvisories.e InitAdvPnlCVs

/************************************************************************/
/*       			CVEVAL    				*/
/* Called w/ every OPIO button push which has a corresponding CV. 	*/
/* CVEVAL should only contain code which impacts the advisory panel--	*/
/* put other code in cvinit or predownload				*/
/************************************************************************/
STATUS cveval( void )
{
	configSystem();
	InitAdvPnlCVs();

	pititle = 1;
	cvdesc(pititle, "Advanced pulse sequence parameters");

	piuset = use0;
	cvdesc(opuser0, "Number of frames to acquire");
	cvdef(opuser0, 2);
	opuser0 = 2;
	cvmin(opuser0, 1);
	nframes = opuser0;

	piuset += use1;
	cvdesc(opuser1, "Number of echoes in each echo train");
	cvdef(opuser1, 16);
	opuser1 = 16;
	cvmin(opuser1, 1);
	nechoes = opuser1;

	piuset += use2;
	cvdesc(opuser2, "Number of echo trains");
	cvdef(opuser2, 1);
	opuser2 = 1;
	cvmin(opuser2, 1);
	ntrains = opuser2;

	piuset += use3;
	cvdesc(opuser3, "2D spiral: 1=out 2=in 3=out-in 4=in-out");
	cvdef(opuser3, 4);
	opuser3 = 4;
	cvmin(opuser3, 1);
	cvmax(opuser3, 4);
	sptype2d = opuser3;

	piuset += use4;
	cvdesc(opuser4, "3D spiral: 1=stack 2=1-ax-rots 3=2-ax-rots 4=naut");
	cvdef(opuser4, 3);
	opuser4 = 3;
	cvmin(opuser4, 1);
	cvmax(opuser4, 4);
	sptype3d = opuser4;

	piuset += use5;
	cvdesc(opuser5, "Spiral radial acceleration factor");
	cvdef(opuser5, 0.7);
	opuser5 = 0.7;
	cvmin(opuser5, 0.1);
	cvmax(opuser5, 2);
	R_accel = opuser5;
	
	piuset += use6;
	cvdesc(opuser6, "Spiral angular acceleration factor");
	cvdef(opuser6, 0.7);
	opuser6 = 0.7;
	cvmin(opuser6, 0.1);
	THETA_accel = opuser6;

	piuset += use7;
	cvdesc(opuser7, "Variable refocuser flip angle attenuation factor");
	cvdef(opuser7, 0.6);
	cvmin(opuser7, 0.1);
	cvmax(opuser7, 1.0);
	opuser7 = 0.6;
	varflipfac = opuser7;

	piuset += use8;
	cvdesc(opuser8, "Recon script ID #");
	cvdef(opuser8, 2327);
	cvmin(opuser8, 0);
	cvmax(opuser8, 9999);
	opuser8 = 2327;
	rhrecon = opuser8;
	
	piuset += use9;
	cvdesc(opuser9, "ASL prep schedule ID #");
	cvdef(opuser9, 0);
	cvmin(opuser9, 0);
	cvmax(opuser9, 9999);
	opuser9 = 0;
	schedule_id = opuser9;
	
	/* 
	 * Calculate RF filter and update RBW:
	 *   &echo1_rtfilt: I: all the filter parameters.
	 *   exist(oprbw): I/O: desired and final allowable bw.
	 *   exist(opxres): I: output pts generated by filter.
	 *   OVERWRITE_OPRBW: oprbw will be updated.
	 */
	if( calcfilter( &echo1_rtfilt,
				exist(oprbw),
				grad_len,
				OVERWRITE_OPRBW ) == FAILURE)
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "calcfilter:echo1" );
		return FAILURE;
	}

	echo1_filt = &echo1_rtfilt;

	/* Divide by 0 protection */
	if( (echo1_filt->tdaq == 0) || 
			floatsAlmostEqualEpsilons(echo1_filt->decimation, 0.0f, 2) ) 
	{
		epic_error( use_ermes, "echo1 tdaq or decimation = 0",
				EM_PSD_BAD_FILTER, EE_ARGS(0) );
		return FAILURE;
	}

	/* For use on the RSP side */
	echo1bw = echo1_filt->bw;

@inline Prescan.e PScveval

	return SUCCESS;
}   /* end cveval() */

void getAPxParam(optval   *min,
		optval   *max,
		optdelta *delta,
		optfix   *fix,
		float    coverage,
		int      algorithm)
{
	/* Need to be filled when APx is supported in this PSD */
}

int getAPxAlgorithm(optparam *optflag, int *algorithm)
{
	return APX_CORE_NONE;
}

/************************************************************************/
/*       			CVCHECK    				*/
/* Executed on each 'next page' to ensure prescription can proceed 	*/
/* to the next page. 							*/
/************************************************************************/
STATUS cvcheck( void )
{
	/* Check if TE is valid */
	if (opte < avminte) {
		epic_error(use_ermes, "opte must be >= %dus", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), INT_ARG, avminte);
		return FAILURE;
	};

	return SUCCESS;
}   /* end cvcheck() */


/************************************************************************/
/*             		    PRE-DOWNLOAD           		        */
/* Executed prior to a download--all operations not needed for the 	*/
/* advisory panel results.  Execute the	pulsegen macro expansions for	*/
/* the predownload section here.  All internal amps, slice ordering,  	*/
/* prescan slice calc., and SAT placement calculations are performed 	*/
/* in this section.  Time anchor settings for pulsegen are done in this */
/* section too.  				 			*/
/************************************************************************/
STATUS predownload( void )
{

	/*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
	/*********************************************************************/
		
	/* Generate initial spiral trajectory */
	fprintf(stderr, "predownload(): calling genspiral()\n");
	if (genspiral(grad_len, 0) == 0) {
		epic_error(use_ermes,"failure to generate spiral waveform", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}

	/* Calculate minimum te */
	avminte = pw_rf2 + 4*psd_grd_wait + 4*trap_ramp_time + pw_gzrf2crush1 + pw_gzrf2crush2 + GRAD_UPDATE_TIME*grad_len + 2*TIMESSI + 2*psd_grd_wait;
	avminte = 1e3*ceil(avminte*1e-3 + 1);

	/* Make sure opte fits */
	if (opte < avminte || opte == PSD_MINFULLTE) {
		opte = avminte;
	}

	/* Generate view transformations */
	fprintf(stderr, "predownload(): calling genviews()\n");
	if (genviews() == 0) {
		epic_error(use_ermes,"failure to generate view transformation matrices", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	
	if (schedule_id > 0) { /* Read in schedules from file */
		
		sprintf(tmpstr, "prep1_pldtbl");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		if (readschedule(schedule_id, prep1_pldtbl, tmpstr, nframes) == -1) {
			epic_error(use_ermes, "file does not have enough lines", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
			return FAILURE;
		}		

		sprintf(tmpstr, "prep1_lbltbl");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		if (readschedule(schedule_id, prep1_lbltbl, tmpstr, nframes) == -1) {
			epic_error(use_ermes, "file does not have enough lines", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));	
			return FAILURE;
		}		
		
		sprintf(tmpstr, "prep2_pldtbl");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		if (readschedule(schedule_id, prep2_pldtbl, tmpstr, nframes) == -1) {
			epic_error(use_ermes, "file not have enough lines", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));	
			return FAILURE;
		}		
		
		sprintf(tmpstr, "prep2_lbltbl");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		if (readschedule(schedule_id, prep2_lbltbl, tmpstr, nframes) == -1) {
			epic_error(use_ermes, "file not have enough lines", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));	
			return FAILURE;
		}		
			
		sprintf(tmpstr, "doblksat");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		readschedule(schedule_id, &doblksat, tmpstr, 1);
		
		sprintf(tmpstr, "prep1_id");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		readschedule(schedule_id, &prep1_id, tmpstr, 1);
		
		sprintf(tmpstr, "prep2_id");
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		readschedule(schedule_id, &prep2_id, tmpstr, 1);

		sprintf(tmpstr, "tadjusttbl");	
		fprintf(stderr, "predownload(): reading in %s using readschedule()\n", tmpstr);	
		switch (readschedule(schedule_id, tadjusttbl, tmpstr, nframes)) {
			case -1:
				epic_error(use_ermes, "file not have enough lines", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));	
				return FAILURE;
			case 0:
				calctadjust();
				break;
		}

	}
	else { /* Generate schedules */

		fprintf(stderr, "predownload(): generating labeling schedule for prep 1 pulse\n");
		if (genschedule(prep1_mod, prep1_pld, prep1_lbltbl, prep1_pldtbl) == 0) {
			epic_error(use_ermes,"failure to generate labeling schedule for prep 1 pulse", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
			return FAILURE;
		}

		fprintf(stderr, "predownload(): generating labeling schedule for prep 2 pulse\n");
		if (genschedule(prep2_mod, prep2_pld, prep2_lbltbl, prep2_pldtbl) == 0) {
			epic_error(use_ermes,"failure to generate labeling schedule for prep 1 pulse", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
			return FAILURE;
		}
		
		/* Calculate tadjust schedule */
		calctadjust();
	}

	/* Read in asl prep pulses */
	fprintf(stderr, "predownload(): calling readprep() to read in ASL prep 1 pulse\n");
	if (readprep(prep1_id, &prep1_len,
		prep1_rho_lbl, prep1_theta_lbl, prep1_grad_lbl,
		prep1_rho_ctl, prep1_theta_ctl, prep1_grad_ctl) == 0)
	{
		epic_error(use_ermes,"failure to read in ASL prep 1 pulse", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	
	fprintf(stderr, "predownload(): calling readprep() to read in ASL prep 2 pulse\n");
	if (readprep(prep2_id, &prep2_len,
		prep2_rho_lbl, prep2_theta_lbl, prep2_grad_lbl,
		prep2_rho_ctl, prep2_theta_ctl, prep2_grad_ctl) == 0)
	{
		epic_error(use_ermes,"failure to read in ASL prep 2 pulse", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	
@inline Prescan.e PSfilter

	/* For Prescan: Inform 'Auto' Prescan about prescan parameters 	*/
	pislquant = nechoes;	/* # of 2nd pass slices */

	/* For Prescan: Declare the entry point table 	*/
	if( entrytabinit( entry_point_table, (int)ENTRY_POINT_MAX ) == FAILURE ) 
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "entrytabinit" );
		return FAILURE;
	}

	/* For Prescan: Define the entry points in the table */
	/* Scan Entry Point */
	(void)strcpy( entry_point_table[L_SCAN].epname, "scan" );
	entry_point_table[L_SCAN].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_SCAN].epprexres = grad_len;

	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	entry_point_table[L_APS2].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_APS2].epprexres = grad_len;

	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );
	entry_point_table[L_MPS2].epfilter = (unsigned char)echo1_filt->fslot;
	entry_point_table[L_MPS2].epprexres = grad_len;

	/* Turn on RF2 pulse */
	rfpulse[RF2_SLOT].activity = PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON;

	/* First, find the peak B1 for the whole sequence. */
	maxB1Seq = 0.0;
	for( entry=0; entry < MAX_ENTRY_POINTS; ++entry )
	{
		if( peakB1( &maxB1[entry], entry, RF_FREE, rfpulse ) == FAILURE )
		{
			epic_error( use_ermes, "peakB1 failed.", EM_PSD_SUPPORT_FAILURE,
					EE_ARGS(1), STRING_ARG, "peakB1" );
			return FAILURE;
		}
		if( maxB1[entry] > maxB1Seq )
		{
			maxB1Seq = maxB1[entry];
		}
	}
	maxB1Seq *= (dolongrf) ? (1) : (2);

	/* Set xmtadd according to maximum B1 and rescale for powermon,
	   adding additional (audio) scaling if xmtadd is too big.
	   Add in coilatten, too. */
	xmtaddScan = -200 * log10( maxB1[L_SCAN] / maxB1Seq ) + getCoilAtten(); 

	if( xmtaddScan > cfdbmax )
	{
		extraScale = (float)pow( 10.0, (cfdbmax - xmtaddScan) / 200.0 );
		xmtaddScan = cfdbmax;
	} 
	else
	{
		extraScale = 1.0;
	}

	if( setScale( L_SCAN, RF_FREE, rfpulse, maxB1[L_SCAN], 
				extraScale) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "setScale" );
		return FAILURE;
	}

	/* Set the parameters for the fat sat pulse */
	a_fatsatrf = 0.5 * 440 / 1250;
	ia_fatsatrf = (int)(a_fatsatrf * (float)max_pg_iamp);
	pw_fatsatrf = 4 * round(cyc_fatsatrf*1e6 / 440);
	res_fatsatrf = pw_fatsatrf / 2;	

	/* Set the parameters for the spin echo tipdown pulse */
	a_rf1 = 0.5;
	thk_rf1 = opslthick*opslquant;
	res_rf1 = 1600;
	pw_rf1 = 3200;
	cyc_rf1 = 2;
	if (dolongrf) {
		res_rf1 *= 2;
		pw_rf1 *= 2;
		cyc_rf1 *= 2;
	}
	flip_rf1 = 90;
	pw_gzrf1 = pw_rf1;
	pw_gzrf1a = trap_ramp_time;
	pw_gzrf1d = trap_ramp_time;

	/* Set the parameters for the tipdown gradient rewinder */
	pw_gzrf1r = 100;
	pw_gzrf1ra = trap_ramp_time;
	pw_gzrf1rd = trap_ramp_time;
	a_gzrf1r = -0.5*a_gzrf1 * (pw_gzrf1 + pw_gzrf1a) / (pw_gzrf1r + pw_gzrf1ra);

	/* Set the parameters for the pre-refocuser crusher */
	a_gzrf2crush1 = 1.5;
	pw_gzrf2crush1 = 100;
	pw_gzrf2crush1a = trap_ramp_time;
	pw_gzrf2crush1d = trap_ramp_time;

	/* Set the parameters for the refocuser pulse */
	a_rf2 = opflip/180.0;
	thk_rf2 = opslthick*opslquant;
	res_rf2 = 1600;
	pw_rf2 = 3200;
	cyc_rf2 = 2;
	if (dolongrf) {
		res_rf2 *= 2;
		pw_rf2 *= 2;
		cyc_rf2 *= 2;
	}
	flip_rf2 = opflip;
	pw_gzrf2 = pw_rf2;
	pw_gzrf2a = trap_ramp_time;
	pw_gzrf2d = trap_ramp_time;

	/* Set the parameters for the post-refocuser crusher */
	a_gzrf2crush2 = a_gzrf2crush1;
	pw_gzrf2crush2 = pw_gzrf2crush1;
	pw_gzrf2crush2a = pw_gzrf2crush1a;
	pw_gzrf2crush2d = pw_gzrf2crush1d;

	/* Update the readout pulse parameters */
	a_gxw = XGRAD_max;
	a_gyw = YGRAD_max;
	a_gzw = ZGRAD_max;
	ia_gxw = MAX_PG_IAMP;
	ia_gyw = MAX_PG_IAMP;
	ia_gzw = MAX_PG_IAMP;
	res_gxw = grad_len;
	res_gyw = grad_len;
	res_gzw = grad_len;
	pw_gxw = GRAD_UPDATE_TIME*grad_len;
	pw_gyw = GRAD_UPDATE_TIME*grad_len;
	pw_gzw = GRAD_UPDATE_TIME*grad_len;

	/* Update the asl prep pulse parameters */
	a_prep1rholbl = prep1_rfmax / ((float)maxB1Seq * 1e3);
	ia_prep1rholbl = (int)ceil(a_prep1rholbl * (float)max_pg_iamp);
	a_prep1rhoctl = prep1_rfmax / ((float)maxB1Seq * 1e3);
	ia_prep1rhoctl = (int)ceil(a_prep1rhoctl * (float)max_pg_iamp);
	a_prep1gradlbl = prep1_gmax / ZGRAD_max;
	ia_prep1gradlbl = (int)ceil(a_prep1gradlbl * (float)max_pg_iamp);
	a_prep1gradctl = prep1_gmax / ZGRAD_max; 
	ia_prep1gradctl = (int)ceil(a_prep1gradctl * (float)max_pg_iamp);
	
	a_prep2rholbl = prep2_rfmax / ((float)maxB1Seq * 1e3);
	ia_prep2rholbl = (int)ceil(a_prep2rholbl * (float)max_pg_iamp);
	a_prep2rhoctl = prep2_rfmax / ((float)maxB1Seq * 1e3);
	ia_prep2rhoctl = (int)ceil(a_prep2rhoctl * (float)max_pg_iamp);
	a_prep2gradlbl = prep2_gmax / ZGRAD_max;
	ia_prep2gradlbl = (int)ceil(a_prep2gradlbl * (float)max_pg_iamp);
	a_prep2gradctl = prep2_gmax / ZGRAD_max; 
	ia_prep2gradctl = (int)ceil(a_prep2gradctl * (float)max_pg_iamp);

	/* Update the bulk saturation pulse parameters */
	a_blksatrho = 1.0;
	res_blksatrho = 250;
	pw_blksatrho = GRAD_UPDATE_TIME*res_blksatrho;
	a_blksattheta = 1.0;
	res_blksattheta = res_blksatrho;
	pw_blksatrho = pw_blksatrho;
	a_blksatgrad = 0.5;
	pw_blksatgrad = 5000;
	pw_blksatgrada = trap_ramp_time;
	pw_blksatgradd = trap_ramp_time;	

	/* Set up the filter structures to be downloaded for realtime 
	   filter generation. Get the slot number of the filter in the filter rack 
	   and assign to the appropriate acquisition pulse for the right 
	   filter selection - LxMGD, RJF */
	setfilter( echo1_filt, SCAN );
	filter_echo1 = echo1_filt->fslot;

	ia_rf1 = max_pg_iamp * (*rfpulse[RF1_SLOT].amp);
	ia_rf2 = max_pg_iamp * (*rfpulse[RF2_SLOT].amp);
	entry_point_table[L_SCAN].epxmtadd = (short)rint( (double)xmtaddScan );

	/* APS2 & MPS2 */
	entry_point_table[L_APS2] = entry_point_table[L_MPS2] = entry_point_table[L_SCAN];	/* copy scan into APS2 & MPS2 */
	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );

	if( orderslice( TYPNCAT, (int)nechoes, (int)nechoes, TRIG_INTERN ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "orderslice" );
	}

	/* nex, exnex, acqs and acq_type are used in the rhheaderinit routine */
	/* -- to initialize recon header variables */
	if( floatsAlmostEqualEpsilons(opnex, 1.0, 2) )
	{
		baseline = 8;
		nex = 1;
		exnex = 1;
	}
	else
	{
		baseline = 0;
		nex = opnex;
		exnex = opnex;
	}
	acqs = nechoes;	/* Fixes the # of rhnpasses to the # of passes */
	acq_type = TYPGRAD;
@inline loadrheader.e rheaderinit   /* Recon variables */
	
	/* Set recon header variables:
	 *   rhptsize: number of bytes per data point
	 *   rhfrsize: number of data points per acquisition
	 *   rhrawsize: total number of bytes to allocate
	 *   rhrcctrl: recon image control (bitmap)
	 *   rhexecctrl: recon executive control (bitmap)
	 */ 
	cvmax(rhfrsize, 32767);
	cvmax(rhnframes, 32767);
	cvmax(rhnslices, 32767);

	rhfrsize = grad_len;
	rhnframes = ntrains*nframes + (ntrains*nframes % 2);
	rhnecho = 1;
	rhnslices = nechoes;
	rhrawsize = 2*rhptsize*rhfrsize * (rhnframes + 1) * rhnslices * rhnecho;
	
	rhrcctrl = 1; /* bit 7 (2^7 = 128) skips all recon */
	rhexecctrl = 2; /* bit 1 (2^1 = 2) sets autolock of raw files + bit 3 (2^3 = 8) transfers images to disk */

	/* Save values to rhusers */
	rhuser0 = nframes;
	rhuser1 = ntrains;
	rhuser2 = nechoes;
	rhuser3 = R_accel;
	rhuser4 = THETA_accel;

	/* Scale the transformation matrices */
	scalerotmats(tmtxtbl, &loggrd, &phygrd, ntrains*nechoes, 0);

@inline Prescan.e PSpredownload	

	return SUCCESS;
}   /* end predownload() */


@inline Prescan.e PShost


@pg
/*********************************************************************
 *                 ASL3DFLEX.E PULSEGEN SECTION                      *
 *                                                                   *
 * Write here the functional code that loads hardware sequencer      *
 * memory with data that will allow it to play out the sequence.     *
 * These functions call pulse generation macros previously defined   *
 * with @pulsedef, and must return SUCCESS or FAILURE.               *
 *********************************************************************/
#include "support_func.h"


STATUS pulsegen( void )
{
	sspinit(psd_board_type);

	/* initialize temporary time marker */
	int tmploc;

	/*********************************/
	/* Generate bulk saturation core */
	/*********************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of bulk saturation core (blksatcore)\n");
	EXTWAVE(RHO, blksatrho, psd_rf_wait, 250*GRAD_UPDATE_TIME, 1.0, 250, sech_7360.rho, , loggrd);
	EXTWAVE(THETA, blksattheta, psd_rf_wait, 250*GRAD_UPDATE_TIME, 1.0, 250, sech_7360.theta, , loggrd);
	TRAPEZOID(ZGRAD, blksatgrad, pend( &blksatrho, "blksatrho", 0) + grad_buff_time + trap_ramp_time, 0, 0, loggrd);

	/* Calculate total core length */
	dur_blksatcore = pw_blksatrho + psd_rf_wait + pw_blksatgrad + 2*trap_ramp_time + 2*grad_buff_time;
	
	fprintf(stderr, "pulsegen(): finalizing bulk saturation core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_blksatcore);
	SEQLENGTH(blksatcore, dur_blksatcore, blksatcore);
	fprintf(stderr, "\tDone.\n");


	/*****************************/
	/* Generate prep1 label core */
	/*****************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep1 label core (prep1lblcore)\n");
	INTWAVE(RHO, prep1rholbl, psd_rf_wait, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_rho_lbl, 1, loggrd); 
	INTWAVE(THETA, prep1thetalbl, psd_rf_wait, 1.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_theta_lbl, 1, loggrd); 
	INTWAVE(ZGRAD, prep1gradlbl, 0, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_grad_lbl, 1, loggrd); 
	
	/* Calculate total core length */
	dur_prep1core = GRAD_UPDATE_TIME*prep1_len + psd_rf_wait + 2*grad_buff_time;

	fprintf(stderr, "pulsegen(): finalizing prep1 label core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep1core);
	SEQLENGTH(prep1lblcore, dur_prep1core, prep1lblcore);
	fprintf(stderr, "\tDone.\n");


	/*******************************/
	/* Generate prep1 control core */
	/*******************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep1 control core (prep1ctlcore)\n");
	INTWAVE(RHO, prep1rhoctl, psd_rf_wait, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_rho_ctl, 1, loggrd); 
	INTWAVE(THETA, prep1thetactl, psd_rf_wait, 1.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_theta_ctl, 1, loggrd); 
	INTWAVE(ZGRAD, prep1gradctl, 0, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_grad_ctl, 1, loggrd); 
	
	fprintf(stderr, "pulsegen(): finalizing prep1 control core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep1core);
	SEQLENGTH(prep1ctlcore, dur_prep1core, prep1ctlcore);
	fprintf(stderr, "\tDone.\n");
	

	/*****************************/
	/* Generate prep2 label core */
	/*****************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep2 label core (prep2lblcore)\n");
	INTWAVE(RHO, prep2rholbl, psd_rf_wait, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_rho_lbl, 1, loggrd); 
	INTWAVE(THETA, prep2thetalbl, psd_rf_wait, 1.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_theta_lbl, 1, loggrd); 
	INTWAVE(ZGRAD, prep2gradlbl, 0, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_grad_lbl, 1, loggrd); 

	/* Calculate total core length */
	dur_prep2core = GRAD_UPDATE_TIME*prep2_len + psd_rf_wait + grad_buff_time;
	
	fprintf(stderr, "pulsegen(): finalizing prep2 label core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep2core);
	SEQLENGTH(prep2lblcore, dur_prep2core, prep2lblcore);
	fprintf(stderr, "\tDone.\n");


	/*******************************/
	/* Generate prep2 control core */
	/*******************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep2 control core (prep2ctlcore)\n");
	INTWAVE(RHO, prep2rhoctl, psd_rf_wait, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_rho_ctl, 1, loggrd); 
	INTWAVE(THETA, prep2thetactl, psd_rf_wait, 1.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_theta_ctl, 1, loggrd); 
	INTWAVE(ZGRAD, prep2gradctl, 0, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_grad_ctl, 1, loggrd); 
	
	fprintf(stderr, "pulsegen(): finalizing prep2 control core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep2core);
	SEQLENGTH(prep2ctlcore, dur_prep2core, prep2ctlcore);
	fprintf(stderr, "\tDone.\n");

	/*********************************/
	/* Generate fat saturation pulse */
	/*********************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of fat sat core (fatsatcore)\n");	
	SINC2(RHO, fatsatrf, psd_rf_wait, 1000, 1.0, ,0.5, , , loggrd);  
	TRAPEZOID(ZGRAD, fatsatgrad, pend(&fatsatrf, "fatsatrf", 0) + grad_buff_time + trap_ramp_time, GRAD_UPDATE_TIME*1000, 0, loggrd);

	/* Calculate total core length */
	dur_fatsatcore = pw_fatsatrf + 2*trap_ramp_time + 2*grad_buff_time + pw_fatsatgrad;

	fprintf(stderr, "pulsegen(): finalizing fat saturation core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_fatsatcore);
	SEQLENGTH(fatsatcore, dur_fatsatcore, fatsatcore);
	fprintf(stderr, "\tDone.\n");


	/***********************************/
	/* Generate spin echo tipdown core */
	/***********************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of spin echo tipdown core (tipdowncore)\n");

	fprintf(stderr, "pulsegen(): generating rf1 (90deg tipdown pulse)...\n");
	SLICESELZ(rf1, trap_ramp_time, 3200, opslthick*opslquant, 90, 1, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf1a, "gzrf1a", 0), pend( &gzrf1d, "gzrf1d", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating rf1 gradient rewinder...\n");
	TRAPEZOID(ZGRAD, gzrf1r, pend( &gzrf1d, "gzrf1d", 0 ) + grad_buff_time, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf1ra, "gzrf1ra", 0), pend( &gzrf1rd, "gzrf1rd", 0));
	fprintf(stderr, "\tDone.\n");	

	/* Calculate total core length */	
	dur_tipdowncore = opte/2 + pw_rf1/2 - pw_rf2/2 - 2*grad_buff_time - 2*trap_ramp_time - pw_gzrf2crush1 - TIMESSI;

	fprintf(stderr, "pulsegen(): finalizing spin echo tipdown core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_tipdowncore);
	SEQLENGTH(tipdowncore, dur_tipdowncore, tipdowncore);
	fprintf(stderr, "\tDone.\n");


	/*************************************/
	/* Generate spin echo refocuser core */
	/*************************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of spin echo refocuser core (refocuscore)\n");

	fprintf(stderr, "pulsegen(): generating pre-rf2 crusher...\n");
	TRAPEZOID(ZGRAD, gzrf2crush1, grad_buff_time + trap_ramp_time, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2crush1a, "gzrf2crush1a", 0), pend( &gzrf2crush1d, "gzrf2crush1d", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating rf2 (180 deg spin echo refocuser)...\n");
	SLICESELZ(rf2, pend( &gzrf2crush1d, "gzrf2crush1d", 0) + grad_buff_time + trap_ramp_time, 3200, opslthick*opslquant, opflip, 1, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2a, "gzrf2a", 0), pend( &gzrf2d, "gzrf2d", 0));
	fprintf(stderr, "\tDone.\n");	

	fprintf(stderr, "pulsegen(): generating post-rf2 crusher...\n");
	TRAPEZOID(ZGRAD, gzrf2crush2, pend( &gzrf2d, "gzrf2d", 0) + grad_buff_time + trap_ramp_time, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2crush2a, "gzrf2crush2a", 0), pend( &gzrf2crush2d, "gzrf2crush2d", 0));
	fprintf(stderr, "\tDone.\n");

	/* Calculate length of core */
	dur_refocuscore = pw_gzrf2crush1 + pw_rf2 + pw_gzrf2crush2 + 4*grad_buff_time + 6*trap_ramp_time;

	fprintf(stderr, "pulsegen(): finalizing spin echo refocuser core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_refocuscore);
	SEQLENGTH(refocuscore, dur_refocuscore, refocuscore);
	fprintf(stderr, "\tDone.\n");


	/********************************/
	/* Generate spiral readout core */
	/********************************/
	/* Calculate start of spiral */
	tmploc = (opte - GRAD_UPDATE_TIME*grad_len)/2 - pw_rf2/2 - 2*grad_buff_time - 3*trap_ramp_time - pw_gzrf2crush1 - TIMESSI;

	fprintf(stderr, "pulsegen(): generating readout x gradient using INTWAVE()...\n");
	INTWAVE(XGRAD, gxw, tmploc, XGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gx, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus, a_gxw = %f\n", pbeg( &gxw, "gxw", 0), pend( &gxw, "gxw", 0), a_gxw);
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating readout y gradient using INTWAVE()...\n");
	INTWAVE(YGRAD, gyw, tmploc, YGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gy, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus, a_gyw = %f\n", pbeg( &gyw, "gyw", 0), pend( &gyw, "gyw", 0), a_gyw);
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating readout z gradient using INTWAVE()...\n");
	INTWAVE(ZGRAD, gzw, tmploc, ZGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gz, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus, a_gzw = %f\n", pbeg( &gzw, "gzw", 0), pend( &gzw, "gzw", 0), a_gzw);
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating data acquisition instruction using ACQUIREDATA()...\n");
	ACQUIREDATA(echo1, tmploc + psd_grd_wait,,,);
	fprintf(stderr, "\tDone.\n");

	/* Calculate length of core */
	dur_seqcore = opte - pw_gzrf2 - 4*grad_buff_time - 6*trap_ramp_time - 2*pw_gzrf2crush1 - 2*TIMESSI;

	fprintf(stderr, "pulsegen(): finalizing spiral readout core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_seqcore);
	SEQLENGTH(seqcore, dur_seqcore, seqcore);
	fprintf(stderr, "\tDone.\n");


	/**********************************/
	/* Generate deadtime (empty) core */
	/**********************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of empty core (emptycore)\n");

	fprintf(stderr, "pulsegen(): finalizing empty core...\n");
	SEQLENGTH(emptycore, 1000, emptycore);
	fprintf(stderr, "\tDone.\n");


@inline Prescan.e PSpulsegen

	PASSPACK(endpass, 49ms);   /* tell Signa system we're done */
	SEQLENGTH(pass, 50ms, pass);

	buildinstr();              /* load the sequencer memory       */

	return SUCCESS;
}   /* end pulsegen() */


/* For Prescan: Pulse Generation functions */
@inline Prescan.e PSipg


@rspvar
/*********************************************************************
 *                   ASL3DFLEX.E RSPVAR SECTION                      *
 *                                                                   *
 * Declare here the real time variables that can be viewed and modi- *
 * fied while the IPG PSD process is running. Only limited standard  *
 * C types are provided: short, int, long, float, double, and 1D     *
 * arrays of those types.                                            *
 *                                                                   *
 * NOTE: Do not declare all real-time variables here because of the  *
 *       overhead required for viewing and modifying them.           *
 *********************************************************************/
extern PSD_EXIT_ARG psdexitarg;

/* Declare rsps */
int echon;
int trainn;
int framen;
int n;

/* Inherited from grass.e: */
int view;
int slice;
int dabop;
int excitation;
int rspent;
int rspdda;
int rspbas;
int rspvus;
int rspgy1;
int rspasl;
int rspesl;
int rspchp;
int rspnex;
int rspslq;
int rspsct;
short chopamp;

/* For Prescan: K */
int seqCount;

@inline Prescan.e PSrspvar 


@rsp
/*********************************************************************
 *                   ASL3DFLEX.E RSP SECTION                         *
 *                                                                   *
 * Write here the functional code for the real time processing (IPG  *
 * side). You may declare standard C variables, but of limited types *
 * short, int, long, float, double, and 1D arrays of those types.    *
 *********************************************************************/
#include <math.h>

/* For IPG Simulator: will generate the entry point list in the IPG tool */
const CHAR *entry_name_list[ENTRY_POINT_MAX] = {
	"scan", 
	"aps2",
	"mps2",
@inline Prescan.e PSeplist
};

/* Do not move the line above and do not insert any code or blank
   lines before the line above.  The code inline'd from Prescan.e
   adds more entry points and closes the list. */

/* Transmit & receive frequencies */
int *rf1_freq;
int xmitfreq;
int *receive_freq1;
int recfreq;

/* Initial transformation matrix */
long tmtx0[9];

STATUS psdinit( void )
{
	/* Initialize everything to a known state */
	setrfconfig( ENBL_RHO1 + ENBL_THETA );
	setssitime( TIMESSI/GRAD_UPDATE_TIME );
	rspqueueinit( 200 );	/* Initialize to 200 entries */
	scopeon( &seqcore );	/* Activate scope for core */
	syncon( &seqcore );		/* Activate sync for core */
	syncoff( &pass );		/* Deactivate sync during pass */
	seqCount = 0;		/* Set SPGR sequence counter */
	settriggerarray( (short)nechoes, rsptrigger );
	setrotatearray( (short)nechoes, rsprot[0] );
	setrfltrs( (int)filter_echo1, &echo1 );

	/* Store initial transformation matrix */
	getrotate(tmtx0, 0);

	return SUCCESS;
}   /* end psdinit() */


@inline Prescan.e PScore

STATUS scancore( void )
{

	float rf2fac;

	/* Determine total # of frames/trains based on entry point */
	int total_frames = (rspent == L_SCAN) ? (nframes) : (1);
	int total_trains = (rspent == L_SCAN) ? (ntrains) : (1000);

	/* Set fat sat frequency */
	setfrequency( (int)(-520 / TARDIS_FREQ_RES), &fatsatrf, 0);

	/* Set transmit frequency and phase */
	rf1_freq = (int *) AllocNode(opslquant*sizeof(int));
	setupslices(rf1_freq, rsp_info, opslquant, a_gzrf1, 1.0, opfov, TYPTRANSMIT);
	if (opslquant%2 == 1)
		xmitfreq = (int)((rf1_freq[opslquant/2] + rf1_freq[opslquant/2-1])/2);
	else
		xmitfreq = (int)rf1_freq[(opslquant+1)/2];
	setfrequency(xmitfreq, &rf1, 0);
	setphase(0.0, &rf1, 0);
	setfrequency(xmitfreq, &rf2, 0);
	setphase(M_PI/2, &rf2, 0);

	/* Set receiver frequency and phase */
	receive_freq1 = (int *) AllocNode(opslquant*sizeof(int));
	for (slice = 0; slice < opslquant; slice++) rsp_info[slice].rsprloc = 0;
	setupslices(receive_freq1, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);
	if (opslquant % 2 == 1)
		recfreq = (int)((receive_freq1[opslquant/2] + receive_freq1[opslquant/2-1])/2);
	else
		recfreq = (int)receive_freq1[(opslquant+1)/2];
	setfrequency(recfreq , &echo1, 0);
	setphase(0.0, &echo1, 0);

	if (rspent != L_SCAN || kill_grads) {
		/* Turn off the gradients */
		fprintf(stderr, "\n scaling grads to zero ");
		setiamp(0, &gxw, 0);
		setiamp(0, &gyw, 0);
		setiamp(0, &gzw, 0);
	}
	else {
		/* Restore the gradients */
		fprintf(stderr, "\n scaling grads to %d ", MAX_PG_IAMP);
		setiamp(MAX_PG_IAMP, &gxw, 0);
		setiamp(MAX_PG_IAMP, &gyw, 0);
		setiamp(MAX_PG_IAMP, &gzw, 0);
	}

	/* Loop through frames */
	for (framen = 0; framen < total_frames; framen++) {
		/* Loop through echo trains */
		for (trainn = -rspdda; trainn < total_trains; trainn++) {
		
			if (doblksat) {
				fprintf(stderr, "scancore(): playing bulk saturation core (%d us)\n", dur_blksatcore);
				/* Play bulk saturation pulse */	
				boffset(off_blksatcore);
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);
			}
			else {
				/* Set length of emptycore to duration of blksatcore */
				setperiod(dur_blksatcore, &emptycore, 0);	
			
				/* Play deadtime core (emptycore) */
				fprintf(stderr, "scancore(): playing deadtime in place of bulk saturation core (%d us)\n", dur_blksatcore);
				boffset(off_emptycore);
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);	
			}		

			if (tadjusttbl[framen] > TIMESSI) {
				/* Set length of emptycore to tadjust */
				setperiod(tadjusttbl[framen] - TIMESSI, &emptycore, 0);

				/* Play TR deadtime (emptycore) */
				fprintf(stderr, "scancore(): playing TR deadtime (tadjust) (%d us)\n", tadjusttbl[framen]);
				boffset(off_emptycore);
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);
			}
				
			
			/* Play prep1 core */
			if (prep1_id > 0) {
				if (rspent != L_SCAN || prep1_lbltbl[framen] == -1) {
					fprintf(stderr, "scancore(): playing deadtime in place of prep 1 pulse (%d us)...\n", dur_prep1core);
					setperiod(dur_prep1core, &emptycore, 0);
					boffset(off_emptycore);
				}
				else if (prep1_lbltbl[framen] == 0) { /* control */
					fprintf(stderr, "scancore(): playing prep 1 control pulse (%d us)...\n", dur_prep1core);
					boffset(off_prep1ctlcore);
				}
				else if (prep1_lbltbl[framen] == 1){ /* label */
					fprintf(stderr, "scancore(): playing prep 1 label pulse (%d us)...\n", dur_prep1core);
					boffset(off_prep1lblcore);
				}
				else {
					fprintf(stderr, "scancore(): invalid pulse 1 type: %d for frame %d...\n", prep1_lbltbl[framen], framen);
					rspexit();
				}
				
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Play pld */
				if (prep1_pldtbl[framen] > TIMESSI) {
					fprintf(stderr, "scancore(): playing prep 1 post-label delay (%d us)...\n", prep1_pldtbl[framen]);
					setperiod(prep1_pldtbl[framen], &emptycore, 0);
					boffset(off_emptycore);
					startseq(0, MAY_PAUSE);
					settrigger(TRIG_INTERN, 0);
				}

			}
			
			/* Play prep2 core */
			if (prep2_id > 0) {
				if (rspent != L_SCAN || prep2_lbltbl[framen] == -1) {
					fprintf(stderr, "scancore(): playing deadtime in place of prep 2 pulse (%d us)...\n", dur_prep2core);
					setperiod(dur_prep2core, &emptycore, 0);
					boffset(off_emptycore);
				}
				else if (prep2_lbltbl[framen] == 0) { /* control */
					fprintf(stderr, "scancore(): playing prep 2 control pulse (%d us)...\n", dur_prep2core);
					boffset(off_prep2ctlcore);
				}
				else if (prep2_lbltbl[framen] == 1){ /* label */
					fprintf(stderr, "scancore(): playing prep 2 label pulse (%d us)...\n", dur_prep2core);
					boffset(off_prep2lblcore);
				}
				else {
					fprintf(stderr, "scancore(): invalid pulse 2 type: %d for frame %d...\n", prep2_lbltbl[framen], framen);
					rspexit();
				}
				
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Play pld */
				if (prep2_pldtbl[framen] > TIMESSI) {
					fprintf(stderr, "scancore(): playing prep 2 post-label delay (%d us)...\n", prep2_pldtbl[framen]);
					setperiod(prep2_pldtbl[framen], &emptycore, 0);
					boffset(off_emptycore);
					startseq(0, MAY_PAUSE);
					settrigger(TRIG_INTERN, 0);
				}

			}
			
			/* Play fat sat core */
			fprintf(stderr, "scancore(): playing fat saturation pulse (%d us)...\n", dur_fatsatcore);
			boffset(off_fatsatcore);
			startseq(0, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);

			/* Play tipdown core */
			fprintf(stderr, "scancore(): playing tipdown (90) pulse (%d us)...\n", dur_tipdowncore);
			boffset(off_tipdowncore);
			startseq(0, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);

			/* Play readout (refocusers + spiral gradients */
			for (echon = 0; echon < nechoes; echon++) {
				
				/* Get the refocuser amplitude and set the amplitude based on variable flip angle */
				getiamp(&chopamp, &rf2, 0);
				if (echon == 0)
					rf2fac = 1.0;
				else
					rf2fac = varflipfac + (float)(echon - 1)/(float)(nechoes - 1) * (1 - varflipfac);
				setiamp((int)(rf2fac*chopamp), &rf2, 0);

				/* Play the refocuser core */
				fprintf(stderr, "scancore(): playing refocuser (180) pulse (%d us)...\n", dur_refocuscore);
				boffset(off_refocuscore);
				startseq(echon, (short)MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Reset the pulse amplitude */
				setiamp(chopamp, &rf2, 0);

				/* Determine space in memory for data storage */	
				if (trainn < 0) { /* turn DAB off for disdaqs */
					fprintf(stderr, "scancore(): playing DISDAQ echo readout %d / %d...\n",
						echon, nechoes);
					loaddab(&echo1,
						0,
						0,
						DABSTORE,
						0,
						DABOFF,
						PSD_LOAD_DAB_ALL);
				}	
				else if (rspent == L_SCAN) { /* load DAB for SCAN process */
					fprintf(stderr, "scancore(): playing scan echo readout %d / %d...\n",
						echon, nechoes);
					loaddab(&echo1,
						echon,
						0,
						DABSTORE,
						framen*ntrains + trainn + 1,
						DABON,
						PSD_LOAD_DAB_ALL);

					/* Set the transformation matrix */
					setrotate(tmtxtbl[trainn*nechoes + echon], echon);
				}
				else { /* load DAB for prescan processes */
					fprintf(stderr, "scancore(): Playing prescan echo readout %d / %d...\n",
						echon, nechoes);
					loaddab(&echo1,
						echon,
						0, /* also try trainn */
						DABSTORE,
						0,
						/* DABON for all APS, or first echo of MPS: */
						(rspent == L_APS2 || echon == 0) ? (DABON) : (DABOFF),
						PSD_LOAD_DAB_ALL);
				}

				/* Play the readout core */
				boffset(off_seqcore);
				startseq(echon, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Reset the rotation matrix */
				setrotate(tmtx0, echon);
			
				/* Negate the 180 amplitude for CPMG scheme */
				if (dophasecycle) {
					getiamp(&chopamp, &rf2, 0);
					setiamp(-chopamp, &rf2, 0);
				}

			}

			/* Reset the 180 amplitude to its absolute value */
			getiamp(&chopamp, &rf2, 0);
			setiamp((int)fabs((float)chopamp), &rf2, 0);

		}
	}
	
	fprintf(stderr, "scancore(): reached end of frame loop, sending endpass packet... \n");

	/* Send SSP packet to end scan */
	boffset( off_pass );
	setwamp(SSPD + DABPASS + DABSCAN, &endpass, 2);
	settrigger(TRIG_INTERN, 0);
	startseq(0, MAY_PAUSE);  
	
	fprintf(stderr, "Done.\n");

	return SUCCESS;
}

/* For Prescan: MPS2 Function */
STATUS mps2( void )
{
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	rspent = L_MPS2;
	rspdda = 0;
	scancore();
	rspexit();

	return SUCCESS;
}   /* end mps2() */


/* For Prescan: APS2 Function */
STATUS aps2( void )
{   
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	rspent = L_APS2;
	rspdda = 2;
	scancore();
	rspexit();

	return SUCCESS;
}   /* end aps2() */

STATUS scan( void )
{ 
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	rspent = L_SCAN;
	rspdda = 0;
	scancore();
	rspexit();

	return SUCCESS;
}


/********************************************
 * dummylinks
 *
 * This routine just pulls in routines from
 * the archive files by making a dummy call.
 ********************************************/
void dummylinks( void )
{
	epic_loadcvs( "thefile" );            /* for downloading CVs */
}

/************************ END OF ASL3DFLEX.E ******************************/

