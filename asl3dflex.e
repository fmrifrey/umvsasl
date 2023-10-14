/*
 * David Frey & Luis Hernandez-Garcia
 * University of Michigan Department of Radiology
 * Functional MRI Laboratory
 *
 * GE Medical Systems
 * Copyright (C) 1996-2003 The General Electric Company
 *
 * File Name : asl3dflex.e
 * Language  : EPIC/ANSI C
 * Date      : 10-May-2023
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
#define MAXNFRAMES 1000 /* Maximum number of temporal frames */
#define MAXITR 50 /* Maximum number of iterations for iterative processes */
#define GAMMA 26754 /* Gyromagnetic ratio (rad/s/G) */
#define TIMESSI 400 /* SSP instruction time */
#define SPOIL_SEED 21001 /* rf spoiler seed */

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
int ZGRAD_risetime;
int ZGRAD_falltime;

/* Declare readout gradient waveform arrays */
int Gx[MAXWAVELEN];
int Gy[MAXWAVELEN];
int Gz[MAXWAVELEN];
int grad_len = 5000;

/* Declare table of readout gradient transformation matrices */
long tmtxtbl[MAXNTRAINS*MAXNECHOES][9];

/* Declare table of refocuser flip angles */
float flipfactbl[MAXNECHOES];
float flipphstbl[MAXNECHOES];

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
int prep1_lbltbl[MAXNFRAMES];
int prep1_pldtbl[MAXNFRAMES];
int prep1_tbgs1tbl[MAXNFRAMES];
int prep1_tbgs2tbl[MAXNFRAMES];

int prep2_lbltbl[MAXNFRAMES];
int prep2_pldtbl[MAXNFRAMES];
int prep2_tbgs1tbl[MAXNFRAMES];
int prep2_tbgs2tbl[MAXNFRAMES];

int tadjusttbl[MAXNFRAMES];
int doblksattbl[MAXNFRAMES];

/* Declare receiver and Tx frequencies */
float recfreq;
float xmitfreq1;
float xmitfreq2;

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

float SLEWMAX = 17000.0 with {1000, 25000.0, 17000.0, VIS, "Maximum allowed slew rate (G/cm/s)",};
float GMAX = 4.0 with {0.5, 5.0, 4.0, VIS, "Maximum allowed gradient (G/cm)",};
float RFMAX = 300 with {0, 500, 300, VIS, "Maximum allowed RF amplitude (mG)",};

/* readout cvs */
int ro_mode = 0 with {0, 1, 0, VIS, "Readout mode (0 = GRE, 1 = FSE)",};
int ro_offset = 0 with { , , 0, VIS, "Readout offset time from center of echo (us)",};
int pgbuffertime = 248 with {100, , 248, INVIS, "Gradient IPG buffer time (us)",};
float phs_tip = 0.0 with { , , 0.0, VIS, "Initial transmitter phase for tipdown pulse",};
float phs_inv = M_PI/2 with { , , M_PI/2, VIS, "Transmitter phase for inversion pulse",};
float phs_rx = 0.0 with { , , 0.0, VIS, "Receiever phase",};
int phscyc_fse = 0 with {0, 1, 0, VIS, "Option to do rf phase cycling for FSE sequence",};
int rfspoil_gre = 0 with {0, 1, 0, VIS, "Option to do rf spoiling for GRE sequence",};
float varflipfac = 1 with {0, 1, 0, VIS, "Scaling factor for variable flip angle schedule (1 = constant fa)",};
int ndisdaqs = 2 with {0, , 2, VIS, "Number of disdaq echo trains at beginning of scan loop",};

/* Trajectory cvs */
int nechoes = 16 with {1, MAXNECHOES, 17, VIS, "Number of echoes per echo train",};
int ntrains = 1 with {1, MAXNTRAINS, 1, VIS, "Number of echo trains per frame",};
int nframes = 2 with {1, , 2, VIS, "Number of frames to acquire",};
int nnav = 250 with {0, 1000, 250, VIS, "Number of navigator points (must be even)",};
float spvd0 = 1.0 with {0.001, 50.0, 1.0, VIS, "Spiral center oversampling factor",};
float spvd1 = 1.0 with {0.001, 50.0, 1.0, VIS, "Spiral edge oversampling factor",};
int sptype2d = 4 with {1, 4, 1, VIS, "1 = spiral out, 2 = spiral in, 3 = spiral out-in, 4 = spiral in-out",};
int sptype3d = 3 with {1, 4, 1, VIS, "1 = stack of spirals, 2 = rotating spirals (single axis), 3 = rotating spirals (2 axes), 4 = debug mode",};
int kill_grads = 0 with {0, 1, 0, VIS, "Option to turn off readout gradients",};

/* ASL prep pulse cvs */
int nm0frames = 2 with {0, , 2, VIS, "Number of M0 frames (no prep pulses are played)",};
int schedule_id = 0 with {0, , 0, VIS, "ASL labeling schedule (0 = no external schedule)",};
int doblksat = 1 with {0, 1, 1, VIS, "Option to do bulk spin saturation at end of each readout",};

int prep1_id = 0 with {0, , 0, VIS, "ASL prep pulse 1: ID number (0 = no pulse)",};
int prep1_pld = 0 with {0, , 0, VIS, "ASL prep pulse 1: post-labeling delay (us; includes background suppression)",};
int prep1_ncycles = 1 with {1, , 1, VIS, "ASL prep pulse 1: number of cycles",};
float prep1_rfmax = 234 with {0, , 0, VIS, "ASL prep pulse 1: maximum RF amplitude",};
float prep1_gmax = 1.5 with {0, , 3, VIS, "ASL prep pulse 1: maximum gradient amplitude",};
int prep1_mod = 1 with {1, 4, 1, VIS, "ASL prep pulse 1: labeling modulation scheme (1 = label/control, 2 = control/label, 3 = always label, 4 = always control)",};
int prep1_tbgs1 = 0 with {0, , 0, VIS, "ASL prep pulse 1: 1st background suppression delay (0 = no pulse)",};
int prep1_tbgs2 = 0 with {0, , 0, VIS, "ASL prep pulse 1: 2nd background suppression delay (0 = no pulse)",};

int prep2_id = 0 with {0, , 0, VIS, "ASL prep pulse 2: ID number (0 = no pulse)",};
int prep2_pld = 0 with {0, , 0, VIS, "ASL prep pulse 2: post-labeling delay (us; includes background suppression)",};
int prep2_ncycles = 1 with {1, , 1, VIS, "ASL prep pulse 2: number of cycles",};
float prep2_rfmax = 234 with {0, , 0, VIS, "ASL prep pulse 2: maximum RF amplitude",};
float prep2_gmax = 3 with {0, , 3, VIS, "ASL prep pulse 2: maximum gradient amplitude",};
int prep2_mod = 1 with {1, 4, 1, VIS, "ASL prep pulse 2: labeling modulation scheme (1 = label/control, 2 = control/label, 3 = always label, 4 = always control)",};
int prep2_tbgs1 = 0 with {0, , 0, VIS, "ASL prep pulse 2: 1st background suppression delay (0 = no pulse)",};
int prep2_tbgs2 = 0 with {0, , 0, VIS, "ASL prep pulse 2: 2nd background suppression delay (0 = no pulse)",};

/* Declare core duration variables */
int dur_blksatcore = 0 with {0, , 0, INVIS, "Duration of the bulk saturation core (us)",};
int dur_prep1core = 0 with {0, , 0, INVIS, "Duration of the ASL prep 1 cores (us)",};
int dur_prep2core = 0 with {0, , 0, INVIS, "Duration of the ASL prep 2 cores (us)",};
int dur_bkgsupcore = 0 with {0, , 0, INVIS, "Duration of the background suppression core (us)",};
int dur_fatsatcore = 0 with {0, , 0, INVIS, "Duration of the fat saturation core (us)",};
int dur_tipdowncore = 0 with {0, , 0, INVIS, "Duration of the tipdown core (us)",};
int dur_flipcore = 0 with {0, , 0, INVIS, "Duration of the refocus core (us)",};
int dur_seqcore = 0 with {0, , 0, INVIS, "Duration of the spiral readout core (us)",};

int deadtime1_seqcore = 0 with {0, , 0, INVIS, "Pre-readout deadtime inside seqcore (us)",};
int deadtime2_seqcore = 0 with {0, , 0, INVIS, "Post-readout deadtime inside seqcore (us)",};
int deadtime_tipdowncore = 0 with {0, , 0, INVIS, "Deadtime inside tipdowncore for FSE readout (us)",};

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
#include "vds.c"

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

/* Declare trajectory generation function prototypes */
int genspiral();
int genviews();

/* Declare function prototypes from aslprep.h */
int readprep(int id, int *len,
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl); 
int readschedule(int id, int* var, char* varname, int lines);
int readschedulef(int id, float* var, char* varname, int lines);
int gentadjusttbl();
int genlbltbl(int mod, int* lbltbl);

/* Define function for calculating max B1 of a sinc pulse */
float calc_sinc_B1(float cyc_rf, int pw_rf, float flip_rf) {

	int M = 1001;
	int n;
	float w[M], x[M];
	float area = 0.0;

	/* Create an M-point symmetrical Hamming window */
	for (n = 0; n < M; n++)
		w[n] = 0.54 - 0.46*cos( 2*M_PI*n / (M-1) );
	
	/* Create a sinc pulse */
	for (n = -(M-1)/2; n < (M-1)/2; n++) {
		if (n == 0)
			x[n + (M-1)/2] = 1.0;
		else
			x[n + (M-1)/2] = sin( 4 * M_PI * cyc_rf * n / (M-1) ) / ( 4 * M_PI * cyc_rf * n / (M-1) );
	}
	
	/* Calculate the area (abswidth) */
	for (n = 0; n < M; n++) {
		area += x[n] * w[n] / M;
	}

	/* Return the B1 (derived from eq. 1 on page 2-31 in EPIC manual) */
	return (SAR_ASINC1/area * 3200/pw_rf * flip_rf/90.0 * MAX_B1_SINC1_90);
}

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
	pite1nub = 2;
	pite1val2 = 100ms;

	/* frequency (xres) */
	cvmin(opxres, 32);
	cvmax(opxres, 128);
	cvdef(opxres, 64);
	opxres = 64;
	pixresnub = 0; /* hide option */

	/* flip angle */
	cvmin(opflip, 0.0);
	cvmax(opflip, 360.0);
	cvdef(opflip, 90.0);
	opflip = 90.0;
	pifanub = 2;
	pifaval2 = 90.0;

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
	getramptime(&ZGRAD_risetime, &ZGRAD_falltime, ZGRAD, &loggrd);	

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
	cvdesc(opuser0, "Readout mode (0 = GRE, 1 = FSE)");
	cvdef(opuser0, 1);
	opuser0 = 1;
	cvmin(opuser0, 0);
	cvmax(opuser0, 1);
	ro_mode = opuser0;

	piuset += use1;
	cvdesc(opuser1, "Number of frames to acquire");
	cvdef(opuser1, 2);
	opuser1 = 2;
	cvmin(opuser1, 1);
	cvmax(opuser2, MAXNFRAMES);
	nframes = opuser1;

	piuset += use2;
	cvdesc(opuser2, "Number of echoes in each echo train");
	cvdef(opuser2, 16);
	opuser2 = 16;
	cvmin(opuser2, 1);
	cvmax(opuser2, MAXNECHOES);
	nechoes = opuser2;

	piuset += use3;
	cvdesc(opuser3, "Number of echo trains");
	cvdef(opuser3, 1);
	opuser3 = 1;
	cvmin(opuser3, 1);
	cvmax(opuser3, MAXNTRAINS);
	ntrains = opuser3;
	
	piuset += use4;
	cvdesc(opuser4, "Number of disdaq trains");
	cvdef(opuser4, 2);
	opuser4 = 2;
	cvmin(opuser4, 0);
	cvmax(opuser4, 100);
	ndisdaqs = opuser4;

	piuset += use5;
	cvdesc(opuser5, "2D spiral: 1=out 2=in 3=out-in 4=in-out");
	cvdef(opuser5, 4);
	opuser5 = 4;
	cvmin(opuser5, 1);
	cvmax(opuser5, 4);
	sptype2d = opuser5;

	piuset += use6;
	cvdesc(opuser6, "3D spiral: 1=stack 2=1-ax-rots 3=2-ax-rots");
	cvdef(opuser6, 3);
	opuser6 = 3;
	cvmin(opuser6, 1);
	cvmax(opuser6, 4);
	sptype3d = opuser6;

	piuset += use7;
	cvdesc(opuser7, "VD-spiral center oversampling factor");
	cvdef(opuser7, 1.0);
	opuser7 = 1.0;
	cvmin(opuser7, 0.001);
	cvmax(opuser7, 50.0);
	spvd0 = opuser7;

	piuset += use8;
	cvdesc(opuser8, "VD-spiral edge oversampling factor");
	cvdef(opuser8, 1.0);
	opuser8 = 1.0;
	cvmin(opuser8, 0.001);
	cvmax(opuser8, 50.0);
	spvd1 = opuser8;

	piuset += use9;
	cvdesc(opuser9, "Variable refocuser flip angle attenuation factor");
	cvdef(opuser9, 0.6);
	cvmin(opuser9, 0.1);
	cvmax(opuser9, 1.0);
	opuser9 = 0.6;
	varflipfac = opuser9;

	piuset += use10;
	cvdesc(opuser10, "Recon script ID #");
	cvdef(opuser10, 2327);
	cvmin(opuser10, 0);
	cvmax(opuser10, 9999);
	opuser10 = 2327;
	rhrecon = opuser10;
	
	piuset += use11;
	cvdesc(opuser11, "ASL prep schedule ID #");
	cvdef(opuser11, 0);
	cvmin(opuser11, 0);
	cvmax(opuser11, 9999);
	opuser11 = 0;
	schedule_id = opuser11;
	
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
	FILE* finfo;
	FILE* fsid;
	int framen, trainn, echon, ddan, slice;
	float rf1_b1, rf2_b1;
	float fatsat_b1, blksat_b1, bkgsup_b1;
	float prep1_b1, prep2_b1;
	int receive_freq[opslquant], rf1_freq[opslquant], rf2_freq[opslquant];
	int tmp_pwa, tmp_pw, tmp_pwd;
	float tmp_a, tmp_area;

	/*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
	/*********************************************************************/

	/* Create a file containing the schedule id number */
	if (schedule_id > 0) {
		fsid = fopen("asl3dflex_scheduleidnum.txt", "w");
		fprintf(fsid, "%05d", schedule_id);
		fclose(fsid);
	}

	/* Read in flip angle schedule */
	sprintf(tmpstr, "flipfactbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);
	switch (readschedulef(schedule_id, flipfactbl, tmpstr, nechoes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);
			for (echon = 0; echon < nechoes; echon++) {
				if (echon == 0)
					flipfactbl[echon] = 1.0;
				else
					flipfactbl[echon] = varflipfac + (float)(echon - 1) / (float)(nechoes - 1) * (1.0 - varflipfac);
			}
			
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in rf phase schedule */
	sprintf(tmpstr, "flipphstbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);
	switch (readschedulef(schedule_id, flipphstbl, tmpstr, nechoes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);
			for (echon = 0; echon < nechoes; echon++) {
				if (ro_mode == 0 && rfspoil_gre == 0) /* non-rf spoiled GRE */
					flipphstbl[echon] = phs_tip;
				else if (ro_mode == 0 && rfspoil_gre == 1) /* rf spoiled GRE */
					flipphstbl[echon] = phs_tip + M_PI * (fmod(pow((float)echon*SPOIL_SEED, 2.0), 2.0) - 1.0);
				else if (ro_mode == 1 && phscyc_fse == 0) /* phase cycled FSE */
					flipphstbl[echon] = phs_inv;
				else /* phase cycled FSE */
					flipphstbl[echon] = phs_inv * pow(-1.0, (float)echon);
			}
			break;
		case -1:
			return FAILURE;
	}

	/* Read in prep1_id (scalar) */
	sprintf(tmpstr, "prep1_id");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	readschedule(schedule_id, &prep1_id, tmpstr, 1);

	/* Read in prep1_pldtbl schedule */
	sprintf(tmpstr, "prep1_pldtbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep1_pldtbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep1_pldtbl[framen] = prep1_pld;
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep1_lbltbl schedule */
	sprintf(tmpstr, "prep1_lbltbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep1_lbltbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			genlbltbl(prep1_mod, prep1_lbltbl);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep1_tbgs1tbl schedule */
	sprintf(tmpstr, "prep1_tbgs1tbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep1_tbgs1tbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep1_tbgs1tbl[framen] = (framen >= nm0frames) ? (prep1_tbgs1) : (0);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep1_tbgs2tbl schedule */
	sprintf(tmpstr, "prep1_tbgs2tbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep1_tbgs2tbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep1_tbgs2tbl[framen] = (framen >= nm0frames) ? (prep1_tbgs2) : (0);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep2_id (scalar) */
	sprintf(tmpstr, "prep2_id");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	readschedule(schedule_id, &prep2_id, tmpstr, 1);
	
	/* Read in prep2_pldtbl schedule */
	sprintf(tmpstr, "prep2_pldtbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep2_pldtbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep2_pldtbl[framen] = prep2_pld;
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep2_lbltbl schedule */
	sprintf(tmpstr, "prep2_lbltbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep2_lbltbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			genlbltbl(prep2_mod, prep2_lbltbl);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep2_tbgs1tbl schedule */
	sprintf(tmpstr, "prep2_tbgs1tbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep2_tbgs1tbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep2_tbgs1tbl[framen] = (framen >= nm0frames) ? (prep2_tbgs1) : (0);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in prep2_tbgs2tbl schedule */
	sprintf(tmpstr, "prep2_tbgs2tbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, prep2_tbgs2tbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				prep2_tbgs2tbl[framen] = (framen >= nm0frames) ? (prep2_tbgs2) : (0);
			break;
		case -1:
			return FAILURE;
	}
	
	/* Read in doblksattbl schedule */
	sprintf(tmpstr, "doblksattbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, doblksattbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++)
				doblksattbl[framen] = (framen >= nm0frames) ? (1) : (0);
			break;
		case -1:
			return FAILURE;
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


	/* Update the background suppression pulse parameters */
	res_bkgsuprho = 500;
	pw_bkgsuprho = 5000;
	a_bkgsuptheta = 1.0;
	res_bkgsuptheta = res_bkgsuprho;
	pw_bkgsuptheta = pw_bkgsuprho;
	
	/* Update the bulk saturation pulse parameters */
	res_blksatrho = res_bkgsuprho/2;
	pw_blksatrho = pw_bkgsuprho/2;
	a_blksattheta = 1.0;
	res_blksattheta = res_blksatrho;
	pw_blksattheta = pw_blksatrho;

	/* Set the parameters for the fat sat pulse */
	a_fatsatrho = 0.5 * 440 / 1250;
	pw_fatsatrho = 4 * round(cyc_fatsatrho*1e6 / 440);
	res_fatsatrho = pw_fatsatrho / 2;	
	
	/* Set the parameters for the spin echo tipdown refocuser gradients */
	tmp_area = a_gzrf1 * (pw_gzrf1 + (pw_gzrf1a + pw_gzrf1d)/2.0);
	amppwgrad(tmp_area, GMAX, 0, 0, ZGRAD_risetime, 0, &tmp_a, &tmp_pwa, &tmp_pw, &tmp_pwd); 
	tmp_a *= -0.5;
	pw_gzrf1r = tmp_pw;
	pw_gzrf1ra = tmp_pwa;
	pw_gzrf1rd = tmp_pwd;
	a_gzrf1r = tmp_a;

	/* Set the parameters for the crusher gradients */
	tmp_area = opxres / (opfov/10.0); /**** figure out an appropriate area for the crushers later ******/
	amppwgrad(tmp_area, GMAX, 0, 0, ZGRAD_risetime, 0, &tmp_a, &tmp_pwa, &tmp_pw, &tmp_pwd); 	
	pw_gzblksatcrush = tmp_pw;
	pw_gzblksatcrusha = tmp_pwa;
	pw_gzblksatcrushd = tmp_pwd;
	a_gzblksatcrush = tmp_a;

	pw_gzfatsatcrush = tmp_pw;
	pw_gzfatsatcrusha = tmp_pwa;
	pw_gzfatsatcrushd = tmp_pwd;
	a_gzfatsatcrush = tmp_a;

	pw_gzrf2crush1 = tmp_pw;
	pw_gzrf2crush1a = tmp_pwa;
	pw_gzrf2crush1d = tmp_pwd;
	a_gzrf2crush1 = tmp_a;

	if (ro_mode == 0) { /* GRE - crusher acts as a refocuser for rf2 slice select gradients */
		tmp_area = a_gzrf2 * (pw_gzrf2 + (pw_gzrf2a + pw_gzrf2d)/2.0);
		amppwgrad(tmp_area, GMAX, 0, 0, ZGRAD_risetime, 0, &tmp_a, &tmp_pwa, &tmp_pw, &tmp_pwd); 	
		tmp_a *= -0.5;
	} /* otherwise, FSE - crusher is just a crusher */
	pw_gzrf2crush2 = tmp_pw;
	pw_gzrf2crush2a = tmp_pwa;
	pw_gzrf2crush2d = tmp_pwd;
	a_gzrf2crush2 = tmp_a;
		
	/* Generate initial spiral trajectory */
	fprintf(stderr, "predownload(): calling genspiral()\n");
	if (genspiral() == 0) {
		epic_error(use_ermes,"failure to generate spiral waveform", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	
	/* Generate view transformations */
	fprintf(stderr, "predownload(): calling genviews()\n");
	if (genviews() == 0) {
		epic_error(use_ermes,"failure to generate view transformation matrices", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}


	/* Calculate minimum te */
	avminte = 0;
	avminte += pw_gzrf2/2 + pw_gzrf2d/2; /* 2nd half of the rf2 pulse */
	avminte += pgbuffertime; /* buffer */
	avminte += pw_gzrf2crush2a + pw_gzrf2crush2 + pw_gzrf2crush2d; /* crush2 pulse */
	avminte += pgbuffertime; /* buffer */
	avminte += TIMESSI; /* inter-core time */
	avminte += pgbuffertime; /* buffer */
	avminte += grad_len*GRAD_UPDATE_TIME; /* readout window length */
	avminte += pgbuffertime; /* buffer */
	avminte += TIMESSI; /* inter-core time */
	avminte += pgbuffertime; /* buffer */
	avminte += pw_gzrf2crush1a + pw_gzrf2crush1 + pw_gzrf2crush1d; /* crush1 pulse */
	avminte += pgbuffertime; /* buffer */
	avminte += pw_gzrf2a + pw_gzrf2/2; /* 1st half of rf2 pulse */
	avminte = GRAD_UPDATE_TIME*ceil((float)avminte/(float)GRAD_UPDATE_TIME); /* round up to gradient sampling interval */
	cvmin(opte, avminte);

	/* Make sure opte fits */
	if (opte < avminte || opte == PSD_MINFULLTE) {
		opte = avminte;
	}

	/* Calculate tipdowncore deadtime */
	deadtime_tipdowncore = opte/2;
	deadtime_tipdowncore -= pw_gzrf1/2 + pw_gzrf1d; /* 2nd half of rf1 pulse */
	deadtime_tipdowncore -= pgbuffertime; /* buffer */
	deadtime_tipdowncore -= pw_gzrf1ra + pw_gzrf1r + pw_gzrf1rd; /* rf1 rewinder pulse */
	deadtime_tipdowncore -= pgbuffertime; /* buffer */
	deadtime_tipdowncore -= TIMESSI; /* inter-core time */
	deadtime_tipdowncore -= pgbuffertime; /* buffer */
	deadtime_tipdowncore -= pw_gzrf2crush1a + pw_gzrf2crush1 + pw_gzrf2crush1d; /* crush1 pulse */
	deadtime_tipdowncore -= pgbuffertime; /* buffer */
	deadtime_tipdowncore -= pw_gzrf2a + pw_gzrf2/2; /* 1st half of rf2 pulse */
	
	/* Calculate seqcore deadtime */
	deadtime1_seqcore = (opte - avminte)/2 + ro_offset;
	deadtime2_seqcore = (opte - avminte)/2 - ro_offset;

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
	a_prep1gradlbl = (prep1_id > 0) ? (prep1_gmax) : (0);
	ia_prep1gradlbl = (int)ceil(a_prep1gradlbl / ZGRAD_max * (float)max_pg_iamp);
	a_prep1gradctl = (prep1_id > 0) ? (prep1_gmax) : (0); 
	ia_prep1gradctl = (int)ceil(a_prep1gradctl / ZGRAD_max * (float)max_pg_iamp);
	a_prep2gradlbl = (prep2_id > 0) ? (prep2_gmax) : (0);
	ia_prep2gradlbl = (int)ceil(a_prep2gradlbl / ZGRAD_max * (float)max_pg_iamp);
	a_prep2gradctl = (prep2_id > 0) ? (prep2_gmax) : (0); 
	ia_prep2gradctl = (int)ceil(a_prep2gradctl / ZGRAD_max * (float)max_pg_iamp);

	/* First, find the peak B1 for all entry points (other than L_SCAN) */
	for( entry=0; entry < MAX_ENTRY_POINTS; ++entry )
	{
		if( peakB1( &maxB1[entry], entry, RF_FREE, rfpulse ) == FAILURE )
		{
			epic_error( use_ermes, "peakB1 failed.", EM_PSD_SUPPORT_FAILURE,
					EE_ARGS(1), STRING_ARG, "peakB1" );
			return FAILURE;
		}
	}

	/* Determine max B1 for the rest of the pulses in L_SCAN */
	rf1_b1 = calc_sinc_B1(cyc_rf1, pw_rf1, 90.0);
	fprintf(stderr, "predownload(): maximum B1 for rf1 pulse: %f\n", rf1_b1);
	if (rf1_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = rf1_b1;

	rf2_b1 = calc_sinc_B1(cyc_rf2, pw_rf2, opflip);
	fprintf(stderr, "predownload(): maximum B1 for rf2 pulse: %f\n", rf2_b1);
	if (rf2_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = rf2_b1;

	blksat_b1 = 0.03867; /* nominal max b1 of sech_7360 pulse */
	fprintf(stderr, "predownload(): maximum B1 for blksat pulse: %f\n", blksat_b1);
	if (blksat_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = blksat_b1;
	
	bkgsup_b1 = 0.03867; /* nominal max b1 of sech_7360 pulse */
	fprintf(stderr, "predownload(): maximum B1 for bkgsup pulse: %f\n", bkgsup_b1);
	if (bkgsup_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = bkgsup_b1;
	
	fatsat_b1 = calc_sinc_B1(cyc_fatsatrho, pw_fatsatrho, 90.0);
	fprintf(stderr, "predownload(): maximum B1 for fatsat pulse: %f\n", fatsat_b1);
	if (fatsat_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = fatsat_b1;
	
	prep1_b1 = (prep1_id > 0) ? (prep1_rfmax*1e-3) : (0);
	fprintf(stderr, "predownload(): maximum B1 for prep1 pulse: %f\n", prep1_b1);
	if (prep1_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = prep1_b1;

	prep2_b1 = (prep2_id > 0) ? (prep2_rfmax*1e-3) : (0);
	fprintf(stderr, "predownload(): maximum B1 for prep2 pulse: %f\n", prep2_b1);
	if (prep2_b1 > maxB1[L_SCAN]) maxB1[L_SCAN] = prep2_b1;

	/* Determine peak B1 across all entry points */
	maxB1Seq = RFMAX * 1e-3;

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

	/* Update all the rf amplitudes */
	a_rf1 = rf1_b1 / maxB1Seq;
	ia_rf1 = a_rf1 * max_pg_iamp;

	a_rf2 = rf2_b1 / maxB1Seq;
	ia_rf2 = a_rf2 * max_pg_iamp;
	
	a_blksatrho = blksat_b1 / maxB1Seq;
	ia_blksatrho = a_blksatrho * max_pg_iamp;
	
	a_bkgsuprho = bkgsup_b1 / maxB1Seq;
	ia_bkgsuprho = a_bkgsuprho * max_pg_iamp;
	
	a_fatsatrho = fatsat_b1 / maxB1Seq;
	ia_fatsatrho = a_fatsatrho * max_pg_iamp;
	
	a_prep1rholbl = prep1_b1 / maxB1Seq;
	ia_prep1rholbl = a_prep1rholbl * max_pg_iamp;
	
	a_prep1rhoctl = prep1_b1 / maxB1Seq;
	ia_prep1rhoctl = a_prep1rhoctl * max_pg_iamp;
	
	a_prep2rholbl = prep2_b1 / maxB1Seq;
	ia_prep2rholbl = a_prep2rholbl * max_pg_iamp;
	
	a_prep2rhoctl = prep2_b1 / maxB1Seq;
	ia_prep2rhoctl = a_prep2rhoctl * max_pg_iamp;

	/* Read in tadjusttbl schedule */
	sprintf(tmpstr, "tadjusttbl");
	fprintf(stderr, "predownload(): reading in %s using readschedule(), schedule_id = %d\n", tmpstr, schedule_id);	
	switch (readschedule(schedule_id, tadjusttbl, tmpstr, nframes)) {
		case 0:
			fprintf(stderr, "predownload(): generating schedule for %s...\n", tmpstr);	
			for (framen = 0; framen < nframes; framen++) {
				avmintr = 0;
				if (ro_mode == 1) /* FSE only */
					avmintr += dur_tipdowncore + TIMESSI;	
				avmintr += nechoes * (dur_flipcore + TIMESSI + dur_seqcore + TIMESSI);
				avmintr += dur_fatsatcore + TIMESSI;
				avmintr += dur_blksatcore + TIMESSI;
				if (prep1_id > 0)
					avmintr += dur_prep1core + TIMESSI + prep1_pldtbl[framen] + TIMESSI;
				if (prep2_id > 0)
					avmintr += dur_prep2core + TIMESSI + prep2_pldtbl[framen] + TIMESSI;
				if (optr < avmintr) {
					epic_error(use_ermes, "optr must be >= %dus", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), INT_ARG, avmintr);
					return FAILURE;
				}
				else
					tadjusttbl[framen] = optr - avmintr - TIMESSI;
			}
			break;
		case -1:
			return FAILURE;
	}

	/* Determine total scan time */
	pidmode = PSD_CLOCK_NORM;
	pitslice = optr;
	pitscan = 0;
	for (ddan = 0; ddan < ndisdaqs; ddan++)
		pitscan += optr;
	for (framen  = 0; framen < nframes; framen++) {
		for (trainn = 0; trainn < ntrains; trainn++) {
			pitscan = 0;
			if (ro_mode == 1) /* FSE only */
				pitscan += dur_tipdowncore + TIMESSI;
			pitscan += nechoes * (dur_flipcore + TIMESSI + dur_seqcore + TIMESSI);
			pitscan += dur_fatsatcore + TIMESSI;
			pitscan += dur_blksatcore + TIMESSI;
			if (prep1_id > 0)
				pitscan += dur_prep1core + prep1_pldtbl[framen] + TIMESSI;
			if (prep2_id > 0)
				pitscan += dur_prep2core + prep2_pldtbl[framen] + TIMESSI;
			pitscan += tadjusttbl[framen] + TIMESSI;
		}
	}
	
	/* Set up the filter structures to be downloaded for realtime 
	   filter generation. Get the slot number of the filter in the filter rack 
	   and assign to the appropriate acquisition pulse for the right 
	   filter selection - LxMGD, RJF */
	setfilter( echo1_filt, SCAN );
	filter_echo1 = echo1_filt->fslot;
	entry_point_table[L_SCAN].epxmtadd = (short)rint( (double)xmtaddScan );

	/* APS2 & MPS2 */
	entry_point_table[L_APS2] = entry_point_table[L_MPS2] = entry_point_table[L_SCAN];	/* copy scan into APS2 & MPS2 */
	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );

	/* Set up Tx/Rx frequencies */
	for (slice = 0; slice < opslquant; slice++) rsp_info[slice].rsprloc = 0;
	setupslices(rf1_freq, rsp_info, opslquant, a_gzrf1, 1.0, opfov, TYPTRANSMIT);
	setupslices(rf2_freq, rsp_info, opslquant, a_gzrf2, 1.0, opfov, TYPTRANSMIT);
	setupslices(receive_freq, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);

	/* Average together all slice frequencies */
	xmitfreq1 = 0;
	xmitfreq2 = 0;
	recfreq = 0;	
	for (slice = 0; slice < opslquant; slice++) {
		xmitfreq1 += (float)rf1_freq[slice] / (float)opslquant;
		xmitfreq2 += (float)rf2_freq[slice] / (float)opslquant;
		recfreq += (float)receive_freq[slice] / (float)opslquant;
	}

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

	/* Scale the transformation matrices */
	scalerotmats(tmtxtbl, &loggrd, &phygrd, ntrains*nechoes, 0);
	
	/* Print scan info to a file */
	finfo = fopen("scaninfo.txt", "w");

	fprintf(finfo, "Rx INFO:\n");
	fprintf(finfo, "\t%-50s%20f\n", "X-Y field of view (mm):", opfov);
	fprintf(finfo, "\t%-50s%20d\n", "Slice quantity:", opslquant);
	fprintf(finfo, "\t%-50s%20f\n", "Slice thickness (mm):", opslthick);

	fprintf(finfo, "\nTIME INFO:\n");
	fprintf(finfo, "\t%-50s%20d\n", "Number of frames:", nframes);
	fprintf(finfo, "\t%-50s%20d\n", "Number of M0 frames:", nm0frames);
	fprintf(finfo, "\t%-50s%20d\n", "Number of disdaq frames:", ndisdaqs);
	fprintf(finfo, "\t%-50s%20f\n", "Repetition time (ms):", (float)optr * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Total scan time (s):", (float)pitscan * 1e-6);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of bulk saturation core (ms):", (float)dur_blksatcore * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of prep 1 core (ms):", (float)dur_prep1core * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of prep 2 core (ms):", (float)dur_prep2core * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of background suppression core (ms):", (float)dur_bkgsupcore * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of fat saturation core (ms):", (float)dur_fatsatcore * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of FSE tipdown core (ms):", (float)dur_tipdowncore * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of flip core (ms):", (float)dur_flipcore * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Duration of readout core (ms):", (float)dur_seqcore * 1e-3);
	
	fprintf(finfo, "\nREADOUT INFO:\n");
	fprintf(finfo, "\t%-50s%20d\n", "Readout mode (0 = GRE, 1 = FSE):", ro_mode);
	fprintf(finfo, "\t%-50s%20f\n", "Echo time (ms):", (float)opte * 1e-3);
	fprintf(finfo, "\t%-50s%20d\n", "Number of echoes per echo train:", nechoes);
	fprintf(finfo, "\t%-50s%20d\n", "Number of echo trains per frame:", ntrains);
	fprintf(finfo, "\t%-50s%20d\n", "2D spiral type:", sptype2d);
	fprintf(finfo, "\t%-50s%20d\n", "3D spiral type:", sptype3d);
	fprintf(finfo, "\t%-50s%20f\n", "VD-spiral center oversampling factor:", spvd0);
	fprintf(finfo, "\t%-50s%20f\n", "VD-spiral edge oversampling factor:", spvd1);
	fprintf(finfo, "\t%-50s%20d\n", "Number of navigator points:", nnav);
	fprintf(finfo, "\t%-50s%20f\n", "Maximum slew rate (G/cm/s^2):", SLEWMAX);
	fprintf(finfo, "\t%-50s%20f\n", "Maximum gradient amplitude (G/cm/s):", GMAX);
	fprintf(finfo, "\t%-50s%20f\n", "Flip angle (deg):", opflip);
	fprintf(finfo, "\t%-50s%20f\n", "Variable flip angle attenuation factor:", varflipfac);
	fprintf(finfo, "\t%-50s%20d\n", "GRE rf spoiling (on/off):", rfspoil_gre);
	fprintf(finfo, "\t%-50s%20d\n", "FSE rf phase cycling (on/off):", phscyc_fse);
	fprintf(finfo, "\t%-50s%20f\n", "Initial transmitter phase for tipdown pulse (rad):", phs_tip);
	fprintf(finfo, "\t%-50s%20f\n", "Transmitter phase for inversion pulse (rad):", phs_inv);
	fprintf(finfo, "\t%-50s%20f\n", "Receiver phase (rad):", phs_rx);
	
	fprintf(finfo, "\nASL PREP INFO:\n");
	fprintf(finfo, "\t%-50s%20d\n", "Labeling schedule ID:", schedule_id);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 1 pulse ID:", prep1_id);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 1 modulation:", prep1_mod);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 1 number of cycles:", prep1_ncycles);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 1 gradient amplitude (G/cm/s):", prep1_gmax);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 1 RF amplitude (mG):", prep1_rfmax);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 1 post-labeling delay (ms):", (float)prep1_pld * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 1 background suppression pulse 1 time (ms):", (float)prep1_tbgs1 * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 1 background suppression pulse 2 time (ms):", (float)prep1_tbgs2 * 1e-3);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 2 pulse ID:", prep2_id);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 2 modulation:", prep2_mod);
	fprintf(finfo, "\t%-50s%20d\n", "Prep 2 number of cycles:", prep2_ncycles);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 2 gradient amplitude (G/cm/s):", prep2_gmax);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 2 RF amplitude (mG):", prep2_rfmax);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 2 post-labeling delay (ms):", (float)prep2_pld * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 2 background suppression pulse 1 time (ms):", (float)prep2_tbgs1 * 1e-3);
	fprintf(finfo, "\t%-50s%20f\n", "Prep 2 background suppression pulse 2 time (ms):", (float)prep2_tbgs2 * 1e-3);
	fprintf(finfo, "\t%-50s%20d\n", "Bulk saturation pulse (on/off):", doblksat);

	fclose(finfo);

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
	int tmploc;

	/*********************************/
	/* Generate bulk saturation core */
	/*********************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of blksatcore (bulk saturation core)\n");
	tmploc = 0;	

	fprintf(stderr, "pulsegen(): generating blksatrho & blksattheta (bulk saturation rf)...\n");
	tmploc += 0; /* start time for blksat rf pulse */
	EXTWAVE(RHO, blksatrho, tmploc, 5000, 1.0, 250, sech_7360.rho, , loggrd);
	EXTWAVE(THETA, blksattheta, tmploc, 5000, 1.0, 250, sech_7360.theta, , loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_blksatrho; /* end time for blksat rf pulse */
	fprintf(stderr, " end: %dus\n", tmploc);	

	fprintf(stderr, "pulsegen(): generating gzblksatcrush (bulk saturation crusher)...\n");
	tmploc += pgbuffertime; /* start time for blksat crusher */
	TRAPEZOID(ZGRAD, gzblksatcrush, tmploc + pw_gzblksatcrusha, 0, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzblksatcrusha + pw_gzblksatcrush + pw_gzblksatcrushd; /* end time for blksat crusher */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */
	
	fprintf(stderr, "pulsegen(): finalizing blksatcore...\n");
	dur_blksatcore = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_blksatcore);
	SEQLENGTH(blksatcore, dur_blksatcore, blksatcore);
	fprintf(stderr, "\tDone.\n");


	/*****************************/
	/* Generate prep1 label core */
	/*****************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep1lblcore (prep1 label core)\n");
	tmploc = 0;
	
	fprintf(stderr, "pulsegen(): generating prep1rholbl, prep1thetalbl & prep1gradlbl (prep1 label rf & gradients)...\n");
	tmploc += pgbuffertime; /* start time for prep1 pulse */
	INTWAVE(RHO, prep1rholbl, tmploc + psd_rf_wait, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_rho_lbl, 1, loggrd); 
	INTWAVE(THETA, prep1thetalbl, tmploc + psd_rf_wait, 1.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_theta_lbl, 1, loggrd); 
	INTWAVE(ZGRAD, prep1gradlbl, tmploc, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_grad_lbl, 1, loggrd); 
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += GRAD_UPDATE_TIME*prep1_len; /* end time for prep1 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing prep1lblcore...\n");
	dur_prep1core = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep1core);
	SEQLENGTH(prep1lblcore, dur_prep1core, prep1lblcore);
	fprintf(stderr, "\tDone.\n");


	/*******************************/
	/* Generate prep1 control core */
	/*******************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep1ctlcore (prep1 control core)\n");
	tmploc = 0;
	
	fprintf(stderr, "pulsegen(): generating prep1rhoctl, prep1thetactl & prep1gradctl (prep1 control rf & gradients)...\n");
	tmploc += pgbuffertime; /* start time for prep1 pulse */
	INTWAVE(RHO, prep1rhoctl, tmploc + psd_rf_wait, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_rho_ctl, 1, loggrd); 
	INTWAVE(THETA, prep1thetactl, tmploc + psd_rf_wait, 1.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_theta_ctl, 1, loggrd); 
	INTWAVE(ZGRAD, prep1gradctl, tmploc, 0.0, prep1_len, GRAD_UPDATE_TIME*prep1_len, prep1_grad_ctl, 1, loggrd); 
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += GRAD_UPDATE_TIME*prep1_len; /* end time for prep1 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing prep1ctlcore...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep1core);
	SEQLENGTH(prep1ctlcore, dur_prep1core, prep1ctlcore);
	fprintf(stderr, "\tDone.\n");
	

	/*****************************/
	/* Generate prep2 label core */
	/*****************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep2lblcore (prep2 label core)\n");
	tmploc = 0;
	
	fprintf(stderr, "pulsegen(): generating prep2rholbl, prep2thetalbl & prep2gradlbl (prep2 label rf & gradients)...\n");
	tmploc += pgbuffertime; /* start time for prep2 pulse */
	INTWAVE(RHO, prep2rholbl, tmploc + psd_rf_wait, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_rho_lbl, 1, loggrd); 
	INTWAVE(THETA, prep2thetalbl, tmploc + psd_rf_wait, 1.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_theta_lbl, 1, loggrd); 
	INTWAVE(ZGRAD, prep2gradlbl, tmploc, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_grad_lbl, 1, loggrd); 
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += GRAD_UPDATE_TIME*prep2_len; /* end time for prep2 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing prep2lblcore...\n");
	dur_prep2core = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep2core);
	SEQLENGTH(prep2lblcore, dur_prep2core, prep2lblcore);
	fprintf(stderr, "\tDone.\n");


	/*******************************/
	/* Generate prep2 control core */
	/*******************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of prep2ctlcore (prep2 control core)\n");
	tmploc = 0;
	
	fprintf(stderr, "pulsegen(): generating prep2rhoctl, prep2thetactl & prep2gradctl (prep2 control rf & gradients)...\n");
	tmploc += pgbuffertime; /* start time for prep2 pulse */
	INTWAVE(RHO, prep2rhoctl, tmploc + psd_rf_wait, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_rho_ctl, 1, loggrd); 
	INTWAVE(THETA, prep2thetactl, tmploc + psd_rf_wait, 1.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_theta_ctl, 1, loggrd); 
	INTWAVE(ZGRAD, prep2gradctl, tmploc, 0.0, prep2_len, GRAD_UPDATE_TIME*prep2_len, prep2_grad_ctl, 1, loggrd); 
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += GRAD_UPDATE_TIME*prep2_len; /* end time for prep2 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing prep2ctlcore...\n");
	fprintf(stderr, "\ttotal time: %dus\n", dur_prep2core);
	SEQLENGTH(prep2ctlcore, dur_prep2core, prep2ctlcore);
	fprintf(stderr, "\tDone.\n");
		
	
	/****************************************/
	/* Generate background suppression core */
	/****************************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of bkgsupcore (background suppression core)\n");
	tmploc = 0;	

	fprintf(stderr, "pulsegen(): generating bkgsuprho & bkgsuptheta (background suppression rf)...\n");
	tmploc += pgbuffertime; /* start time for bkgsup rf */
	EXTWAVE(RHO, bkgsuprho, tmploc, 5000, 1.0, 250, sech_7360.rho, , loggrd);
	EXTWAVE(THETA, bkgsuptheta, tmploc, 5000, 1.0, 250, sech_7360.theta, , loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_bkgsuprho; /* end time for bkg sup rf */
	fprintf(stderr, " end: %dus\n", tmploc);	
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing bkgsupcore...\n");
	dur_bkgsupcore = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_bkgsupcore);
	SEQLENGTH(bkgsupcore, dur_bkgsupcore, bkgsupcore);
	fprintf(stderr, "\tDone.\n");


	/*********************************/
	/* Generate fat saturation pulse */
	/*********************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of fatsatcore (fat saturation core)\n");	
	tmploc = 0;
	
	fprintf(stderr, "pulsegen(): generating fatsatrho (fat saturation rf pulse)...\n");
	tmploc += 0; /* start time for fatsatrho */
	SINC2(RHO, fatsatrho, tmploc + psd_rf_wait, 3200, 1.0, ,0.5, , , loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_fatsatrho; /* end time for fatsatrho */
	fprintf(stderr, " end: %dus\n", tmploc);	
 
	fprintf(stderr, "pulsegen(): generating gzfatsatcrush (fat saturation crusher gradients)...\n");
	tmploc += pgbuffertime; /* start time for gzfatsatcrush */
	TRAPEZOID(ZGRAD, gzfatsatcrush, tmploc + pw_gzfatsatcrusha, GRAD_UPDATE_TIME*1000, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzfatsatcrusha + pw_gzfatsatcrush + pw_gzfatsatcrushd; /* end time for gzfatsatcrush */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	fprintf(stderr, "pulsegen(): finalizing fatsatcore...\n");
	dur_fatsatcore = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_fatsatcore);
	SEQLENGTH(fatsatcore, dur_fatsatcore, fatsatcore);
	fprintf(stderr, "\tDone.\n");


	/***********************************/
	/* Generate spin echo tipdown core */
	/***********************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of tipdowncore (rf tipdown core for FSE readout)\n");
	tmploc = 0;

	fprintf(stderr, "pulsegen(): generating rf1 (90deg tipdown pulse)...\n");
	tmploc += pgbuffertime; /* start time for rf1 pulse */
	SLICESELZ(rf1, tmploc + pw_gzrf1a, 6400, (opslthick + opslspace)*opslquant, 90.0, 4, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf1a + pw_gzrf1 + pw_gzrf1d; /* end time for rf1 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);

	fprintf(stderr, "pulsegen(): generating gzrf1r (90deg tipdown gradient refocuser)...\n");
	tmploc += pgbuffertime; /* start time for gzrf1r */
	TRAPEZOID(ZGRAD, gzrf1r, tmploc + pw_gzrf1ra, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf1ra + pw_gzrf1r + pw_gzrf1rd; /* end time for gzrf1r pulse */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	tmploc += deadtime_tipdowncore; /* add deadtime to account for TE */

	fprintf(stderr, "pulsegen(): finalizing spin echo tipdown core...\n");
	dur_tipdowncore = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_tipdowncore);
	SEQLENGTH(tipdowncore, dur_tipdowncore, tipdowncore);
	fprintf(stderr, "\tDone.\n");


	/*************************************/
	/* Generate spin echo refocuser core */
	/*************************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of (flipcore) FSE inversion/GRE tipdown core\n");
	tmploc = 0;

	fprintf(stderr, "pulsegen(): generating gzrf2crush1 (pre-rf2 crusher)...\n");
	tmploc += pgbuffertime; /* start time for gzrf2crush1 */
	TRAPEZOID(ZGRAD, gzrf2crush1, tmploc + pw_gzrf2crush1a, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf2crush1a + pw_gzrf2crush1 + pw_gzrf2crush1d; /* end time for gzrf2crush1 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);

	fprintf(stderr, "pulsegen(): generating rf2 (FSE inversion/GRE tipdown)...\n");
	tmploc += pgbuffertime; /* start time for rf2 */
	SLICESELZ(rf2, tmploc + pw_gzrf2a, 3200, (opslthick + opslspace)*opslquant, 180.0, 2, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf2a + pw_gzrf2 + pw_gzrf2d; /* end time for rf2 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);

	fprintf(stderr, "pulsegen(): generating gzrf2crush2 (post-rf2 crusher for FSE/tipdown gradient rewinder for GRE)...\n");
	tmploc += pgbuffertime; /* start time for gzrf2crush2 */
	TRAPEZOID(ZGRAD, gzrf2crush2, tmploc + pw_gzrf2crush2a, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += pw_gzrf2crush2a + pw_gzrf2crush2 + pw_gzrf2crush2d; /* end time for gzrf2crush2 pulse */
	fprintf(stderr, " end: %dus\n", tmploc);	
	tmploc += pgbuffertime;

	fprintf(stderr, "pulsegen(): finalizing spin echo refocuser core...\n");
	dur_flipcore = tmploc;
	fprintf(stderr, "\ttotal time: %dus\n", dur_flipcore);
	SEQLENGTH(flipcore, dur_flipcore, flipcore);
	fprintf(stderr, "\tDone.\n");


	/*************************/
	/* Generate readout core */
	/*************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of readout core (seqcore)\n");
	tmploc = 0;

	fprintf(stderr, "pulsegen(): generating gxw, gyw, & gzw (readout gradients) and echo1 (data acquisition window)...\n");
	tmploc += deadtime1_seqcore + pgbuffertime; /* start time for readout */
	INTWAVE(XGRAD, gxw, tmploc, XGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gx, 1, loggrd);
	INTWAVE(YGRAD, gyw, tmploc, YGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gy, 1, loggrd);
	INTWAVE(ZGRAD, gzw, tmploc, ZGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gz, 1, loggrd);
	ACQUIREDATA(echo1, tmploc + psd_grd_wait,,,);
	fprintf(stderr, "\tstart: %dus, ", tmploc);
	tmploc += GRAD_UPDATE_TIME*grad_len; /* end time for readout */
	fprintf(stderr, " end: %dus\n", tmploc);
	tmploc += pgbuffertime; /* add some buffer */

	tmploc += deadtime2_seqcore; /* add post-readout deadtime */

	fprintf(stderr, "pulsegen(): finalizing spiral readout core...\n");
	dur_seqcore = tmploc;
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
int ddan;
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

	int ttmp;

	/* Determine total # of frames/trains based on entry point */
	int total_frames = (rspent == L_SCAN) ? (nframes) : (1);
	int total_trains = (rspent == L_SCAN) ? (ntrains) : (1000);
			
	/* Set rf1, rf2 tx and rx frequency */
	setfrequency((int)xmitfreq1, &rf1, 0);
	setfrequency((int)xmitfreq2, &rf2, 0);
	setfrequency((int)recfreq, &echo1, 0);

	/* Set fat sat frequency */
	setfrequency( (int)(-520 / TARDIS_FREQ_RES), &fatsatrho, 0);
		
	/* Turn off the gradients */
	fprintf(stderr, "\n scaling grads to zero ");
	setiamp(0, &gxw, 0);
	setiamp(0, &gyw, 0);
	setiamp(0, &gzw, 0);

	/* Loop through disdaqs */
	for (ddan = 0; ddan < rspdda; ddan++) {
		
		/* Play bulk saturation pulse */	
		fprintf(stderr, "scancore(): playing bulk saturation core (%d us)\n", dur_blksatcore);
		boffset(off_blksatcore);
		startseq(0, MAY_PAUSE);
		settrigger(TRIG_INTERN, 0);

		/* Play TR deadtime */
		ttmp = optr - dur_blksatcore - dur_fatsatcore - dur_tipdowncore - nechoes*(dur_flipcore + dur_seqcore);
		if (ttmp > TIMESSI) {
			fprintf(stderr, "scancore(): playing TR deadtime (%d us)\n", ttmp);
			setperiod(ttmp - TIMESSI, &emptycore, 0);
			boffset(off_emptycore);
			startseq(0, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);
		}
			
		/* Play fat sat core */
		fprintf(stderr, "scancore(): playing fat saturation pulse (%d us)...\n", dur_fatsatcore);
		boffset(off_fatsatcore);
		startseq(0, MAY_PAUSE);
		settrigger(TRIG_INTERN, 0);
	
		/* Play tipdown core */
		if (ro_mode == 1) {
			fprintf(stderr, "scancore(): playing tipdown (90) pulse (%d us) for FSE readout...\n", dur_tipdowncore);
			setphase(phs_tip, &rf1, 0);
			boffset(off_tipdowncore);
			startseq(0, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);
		}

		/* Play readout (refocusers + spiral gradients */
		for (echon = 0; echon < nechoes; echon++) {

			/* Set the refocuser flip angle & phase */
			setiamp((int)(flipfactbl[echon]*ia_rf2), &rf2, 0);
			setphase(flipphstbl[echon], &rf2, 0);
			if (ro_mode == 1)
				setphase(flipphstbl[echon], &echo1, 0);
			else
				setphase(phs_tip, &echo1, 0);

			/* Play the refocuser core */
			fprintf(stderr, "scancore(): playing flip pulse (%d us)...\n", dur_flipcore);
			boffset(off_flipcore);
			startseq(echon, (short)MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);

			/* Turn off DAQ for disdaqs */	
			fprintf(stderr, "scancore(): playing DISDAQ echo readout %d / %d...\n",
					echon, nechoes);
			loaddab(&echo1,
					0,
					0,
					DABSTORE,
					0,
					DABOFF,
					PSD_LOAD_DAB_ALL);

			/* Play the readout core */
			boffset(off_seqcore);
			startseq(echon, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);

			/* Reset the rotation matrix */
			setrotate(tmtx0, echon);

			/* Reset the rf2 amplitude */
			setiamp(ia_rf2, &rf2, 0);
		}
	}
	
	if (rspent == L_SCAN && !kill_grads) {
		/* Restore the gradients */
		fprintf(stderr, "\n scaling grads to %d ", MAX_PG_IAMP);
		setiamp(MAX_PG_IAMP, &gxw, 0);
		setiamp(MAX_PG_IAMP, &gyw, 0);
		setiamp(MAX_PG_IAMP, &gzw, 0);
	}

	/* Loop through frames */
	for (framen = 0; framen < total_frames; framen++) {
		/* Loop through echo trains */
		for (trainn = 0; trainn < total_trains; trainn++) {
		
			if (doblksattbl[framen] > 0) {
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

			if (tadjusttbl[framen] > 0) {
				/* Set length of emptycore to tadjust */
				setperiod(tadjusttbl[framen], &emptycore, 0);

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

				/* Play pld and background suppression */
				if (prep1_pldtbl[framen] > TIMESSI) {
				
					/* Initialize pld before subtracting out tbgs timing */
					ttmp = prep1_pldtbl[framen];
						
					if (prep1_tbgs1tbl[framen] > 0) {
						/* Play first background suppression delay/pulse */
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 delay (%d us)...\n", prep1_tbgs1tbl[framen]);
						setperiod(prep1_tbgs1tbl[framen], &emptycore, 0);
						ttmp -= prep1_tbgs1tbl[framen];
						boffset(off_emptycore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
						
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 (%d us)...\n", dur_bkgsupcore);
						ttmp -= dur_bkgsupcore;
						boffset(off_bkgsupcore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
					}
					
					if (prep1_tbgs2tbl[framen] > 0) {
						/* Play 2nd background suppression delay/pulse */
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 2 delay (%d us)...\n", prep1_tbgs2tbl[framen]);
						setperiod(prep1_tbgs2tbl[framen], &emptycore, 0);
						ttmp -= prep1_tbgs2tbl[framen];
						boffset(off_emptycore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
						
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 (%d us)...\n", dur_bkgsupcore);
						ttmp -= dur_bkgsupcore;
						boffset(off_bkgsupcore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
					}
					
					/* Check that ttmp is non-negative */
					if (ttmp < 0) {
						fprintf(stderr, "scancore(): invalid pld and background suppression time combination for frame %d...\n", framen);
						rspexit();
					}

					/* Play remaining PLD deadtime */
					fprintf(stderr, "scancore(): playing prep 1 post-label delay (%d us), total end delay = %d us...\n", prep1_pldtbl[framen], ttmp);
					setperiod(ttmp, &emptycore, 0);
					boffset(off_emptycore);
					startseq(0, MAY_PAUSE);
					settrigger(TRIG_INTERN, 0);
				}

			}
			
			/* Play prep2 core */
			if (prep2_id > 0) {
				if (rspent != L_SCAN || prep2_lbltbl[framen] == -1) {
					fprintf(stderr, "scancore(): playing deadtime in place of prep 1 pulse (%d us)...\n", dur_prep2core);
					setperiod(dur_prep2core, &emptycore, 0);
					boffset(off_emptycore);
				}
				else if (prep2_lbltbl[framen] == 0) { /* control */
					fprintf(stderr, "scancore(): playing prep 1 control pulse (%d us)...\n", dur_prep2core);
					boffset(off_prep2ctlcore);
				}
				else if (prep2_lbltbl[framen] == 1){ /* label */
					fprintf(stderr, "scancore(): playing prep 1 label pulse (%d us)...\n", dur_prep2core);
					boffset(off_prep2lblcore);
				}
				else {
					fprintf(stderr, "scancore(): invalid pulse 1 type: %d for frame %d...\n", prep2_lbltbl[framen], framen);
					rspexit();
				}
				
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Play pld and background suppression */
				if (prep2_pldtbl[framen] > TIMESSI) {
				
					/* Initialize pld before subtracting out tbgs timing */
					ttmp = prep2_pldtbl[framen];
						
					if (prep2_tbgs1tbl[framen] > 0) {
						/* Play first background suppression delay/pulse */
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 delay (%d us)...\n", prep2_tbgs1tbl[framen]);
						setperiod(prep2_tbgs1tbl[framen], &emptycore, 0);
						ttmp -= prep2_tbgs1tbl[framen];
						boffset(off_emptycore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
						
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 (%d us)...\n", dur_bkgsupcore);
						ttmp -= dur_bkgsupcore;
						boffset(off_bkgsupcore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
					}
					
					if (prep2_tbgs2tbl[framen] > 0) {
						/* Play 2nd background suppression delay/pulse */
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 2 delay (%d us)...\n", prep2_tbgs2tbl[framen]);
						setperiod(prep2_tbgs2tbl[framen], &emptycore, 0);
						ttmp -= prep2_tbgs2tbl[framen];
						boffset(off_emptycore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
						
						fprintf(stderr, "scancore(): playing prep 1 bkg suppression pulse 1 (%d us)...\n", dur_bkgsupcore);
						ttmp -= dur_bkgsupcore;
						boffset(off_bkgsupcore);
						startseq(0, MAY_PAUSE);
						settrigger(TRIG_INTERN, 0);
					}
					
					/* Check that ttmp is non-negative */
					if (ttmp < 0) {
						fprintf(stderr, "scancore(): invalid pld and background suppression time combination for frame %d...\n", framen);
						rspexit();
					}

					/* Play remaining PLD deadtime */
					fprintf(stderr, "scancore(): playing prep 1 post-label delay (%d us), total end delay = %d us...\n", prep2_pldtbl[framen], ttmp);
					setperiod(ttmp, &emptycore, 0);
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
			if (ro_mode == 1) {
				fprintf(stderr, "scancore(): playing tipdown (90) pulse (%d us) for FSE readout...\n", dur_tipdowncore);
				setphase(phs_tip, &rf1, 0);
				boffset(off_tipdowncore);
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);
			}

			/* Play readout (refocusers + spiral gradients */
			for (echon = 0; echon < nechoes; echon++) {
				
				/* Set the refocuser flip angle & phase */
				setiamp((int)(flipfactbl[echon]*ia_rf2), &rf2, 0);
				setphase(flipphstbl[echon], &rf2, 0);
				if (ro_mode == 1)
					setphase(flipphstbl[echon], &echo1, 0);
				else
					setphase(phs_tip, &echo1, 0);
				
				/* Play the refocuser core */
				fprintf(stderr, "scancore(): playing flip pulse (%d us)...\n", dur_flipcore);
				boffset(off_flipcore);
				startseq(echon, (short)MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Determine space in memory for data storage */	
				if (rspent == L_SCAN) { /* load DAB for SCAN process */
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
				
				/* Reset the rf2 amplitude */
				setiamp(ia_rf2, &rf2, 0);
			}

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
	rspdda = ndisdaqs;
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


@host

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
	int tmp_pwa, tmp_pw, tmp_pwd;
	float Tr_ze, Tp_ze, T_sp, T_rd, Tr_kr, Tp_kr;
	float t1, t2, t3;
	float gxn, gyn, gzn, kxn, kyn, kzn;
	float gx0, gy0, g0, kx0, ky0, kz0, k0;
	float h_ze, h_kr;
	float F[2];

	/* declare waveforms */
	float *gx_sp, *gy_sp;
	float gx[MAXWAVELEN], gy[MAXWAVELEN], gz[MAXWAVELEN];

	/* convert units */
	float dt = GRAD_UPDATE_TIME*1e-6; /* raster time (s) */
	float D = (float)opfov / 10.0; /* fov (cm) */
	float gm = GMAX; /* gradient amplitude limit (G/cm) */
	float sm = SLEWMAX; /* slew limit (G/cm/s) */
	float gam = 4258; /* gyromagnetic ratio (Hz/G) */
	float kxymax = opxres / D / 2.0; /* kspace xy sampling radius (cm^-1) */
	float kzmax = nechoes / D / 2.0; /* kspace z sampling radius for SOS (cm^-1) */

	/* generate the z encoding trapezoid gradient */
	amppwgrad(kzmax/gam*1e6, gm, 0, 0, ZGRAD_risetime, 0, &h_ze, &tmp_pwa, &tmp_pw, &tmp_pwd);
	Tr_ze = (tmp_pwa + tmp_pwd) / 2.0 * 1e-6;
	Tp_ze = tmp_pw * 1e-6;
	np_ze = round((2*Tr_ze + Tp_ze) / dt);

	/* generate the spiral trajectory */
	F[0] = spvd0 * D;
	F[1] = (D*spvd1 - F[0]) / kxymax;
	calc_vds(sm, gm, dt, dt, ((sptype2d<3)?(1):(2))*ntrains, F, 2, kxymax, MAXWAVELEN, &gx_sp, &gy_sp, &np_sp);
	T_sp = dt * np_sp;

	/* calculate gradients at end of spiral */
	gx0 = gx_sp[np_sp - 1];
	gy0 = gy_sp[np_sp - 1];
	g0 = sqrt(pow(gx0,2) + pow(gy0,2));

	/* calculate gradient ramp down time and round up to nearest sampling interval */
	T_rd = g0 / sm;
	T_rd = dt * ceil(T_rd / dt);
	np_rd = round(T_rd/dt);

	/* calculate gradients at end of ramp down */
	kx0 = gam * dt * fsumarr(gx_sp, np_sp - 1) + gam * 1/2 * (T_rd + dt) * gx0;
	ky0 = gam * dt * fsumarr(gy_sp, np_sp - 1) + gam * 1/2 * (T_rd + dt) * gy0;
	kz0 = kzmax;
	k0 = sqrt(pow(kx0,2) + pow(ky0,2) + pow(kz0,2));

	/* generate the kspace rewinder */
	amppwgrad(k0/gam*1e6, gm, 0, 0, ZGRAD_risetime, 0, &h_kr, &tmp_pwa, &tmp_pw, &tmp_pwd);
	Tr_kr = (tmp_pwa + tmp_pwd) / 2.0 * 1e-6;
	Tp_kr = tmp_pw * 1e-6;
	np_kr = round((2*Tr_kr + Tp_kr) / dt);

	/* calculate time markers */
	t1 = 2*Tr_ze + Tp_ze;
	t2 = t1 + T_sp;
	t3 = t2 + T_rd;

	/* calculate total number of points */
	np = np_ze + np_sp + np_rd + np_kr + 1;

	/* loop through time points */
	memset(gx, 0, sizeof gx);
	memset(gy, 0, sizeof gx);
	memset(gz, 0, sizeof gx);
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

		Gx[n] = 2*round(MAX_PG_WAMP/XGRAD_max * gx[n] / 2.0);
		Gy[n] = 2*round(MAX_PG_WAMP/YGRAD_max * gy[n] / 2.0);
		Gz[n] = 2*round(MAX_PG_WAMP/ZGRAD_max * gz[n] / 2.0);
		
		fprintf(fID_ktraj, "%f \t%f \t%f\n", kxn, kyn, kzn);
		fprintf(fID_grad, "%f \t%f \t%f\n", gx[n], gy[n], gz[n]);
		fprintf(fID_graddac, "%d \t%d \t%d\n", Gx[n], Gy[n], Gz[n]);
	}

	fclose(fID_ktraj);
	fclose(fID_grad);
	fclose(fID_graddac);

	return 1;
}

int genviews() {

	/* Declare values and matrices */
	FILE* fID = fopen("kviews.txt","w");
	int trainn, echon, n;
	float rx, ry, rz, dz;
	float Rx[9], Ry[9], Rz[9], Tz[9];
	float T_0[9], T[9];

	/* Initialize z translation to identity matrix */
	eye(Tz, 3);

	/* Get original transformation matrix */
	for (n = 0; n < 9; n++) T_0[n] = (float)rsprot[0][n] / MAX_PG_WAMP;
	orthonormalize(T_0, 3, 3);

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
					dz = 0.0;
					break;
				case 3 : /* Double axis rotations */
					rx = 2.0 * M_PI * (float)(trainn*nechoes + echon) / PHI;
					ry = acos(fmod(1 - 2*(trainn*nechoes + echon + 0.5) / (float)(ntrains*nechoes), 1));
					rz = 0.0;
					dz = 0.0;
					break;
				case 4: /* Debugging case */
					rx = M_PI * (float)(echon) / nechoes;
					ry = M_PI * (float)(trainn) / ntrains;
					rz = 0.0;
					dz = 0.0;
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
			multmat(3,3,3,Rz,T,T);
			multmat(3,3,3,Rx,T,T);
			multmat(3,3,3,Ry,T,T);

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

int readprep(int id, int *len,
		int *rho_lbl, int *theta_lbl, int *grad_lbl,
		int *rho_ctl, int *theta_ctl, int *grad_ctl)
{

	/* Declare variables */
	char fname[80];
	FILE *fID;
	char buff[200];
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
	while (fgets(buff, 200, fID)) {
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
	while (fgets(buff, 200, fID)) {
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
	while (fgets(buff, 200, fID)) {
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
	
	if (id == 0)
		return 0;

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

int readschedulef(int id, float* var, char* varname, int lines) {

	FILE* fID;
	char fname[200];
	float val;
	
	if (id == 0)
		return 0;

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
	while (fscanf(fID, "%f\n", &val) != EOF) {
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

int genlbltbl(int mod, int* lbltbl)
{
	int framen;

	/* Loop through frames */
	for (framen = 0; framen < nframes; framen++) {

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

/************************ END OF ASL3DFLEX.E ******************************/

