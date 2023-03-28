/*
 * GE Medical Systems
 * Copyright (C) 1996-2003 The General Electric Company
 *
 * File Name : grass.e
 * Language  : EPIC/ANSI C
 * Date      : 01-Jan-1996
 *
 * This is the solution to the PSD Generation Tutorial: grass.e
 * A Gradient Recalled Acquisition Steady State (GRASS)
 */
/* do not modify anything above this line */

/* MRIhc15728 Added psdiopt.h file*/
/* 14.5   10/05/2006 TS    MRIhc15304 - Coil info related changes */
/* PX25  06/Jun/2014 YT    HCSDM00289004 - add APx functions */

@inline epic.h
@inline intwave.h

@global
/*********************************************************************
 *                    GRASS.E GLOBAL SECTION                         *
 *                                                                   *
 * Common code shared between the Host and IPG PSD processes.  This  *
 * section contains all the #define's, global variables and function *
 * declarations (prototypes).                                        *
 *********************************************************************/
#include <stdio.h>
#include <string.h>

#include "em_psd_ermes.in"
#include "grad_rf_grass.globals.h"

#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "epic_error.h"
#include "epicfuns.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#include "psdutil.h"
#include "psd_proto.h"
#include "epic_iopt_util.h"
#include "filter.h"

#include "grass.h"

/* Define important values */
#define MAXWAVELEN 50000 /* Maximum wave length for gradients */
#define MAXNTRAINS 50 /* Maximum number of echo trains per frame */
#define MAXNECHOES 50 /* Maximum number of echoes per echo train */
#define MAXITR 50 /* Maximum number of iterations for iterative processes */
#define GAMMA 26754 /* Gyromagnetic ratio */
#define TSP_GRAD 4 /* Scanner gradient transmitter sampling rate */
#define TSP_RF 2 /* Scanner rf transmitter sampling rate */

@inline Prescan.e PSglobal
int debugstate = 1;

@ipgexport
/*********************************************************************
 *                  GRASS.E IPGEXPORT SECTION                        *
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

/* Declare gradient waveform arrays */
int Gx[MAXWAVELEN];
int Gy[MAXWAVELEN];
int Gz[MAXWAVELEN];

@cv
/*********************************************************************
 *                       GRASS.E CV SECTION                          *
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

int numdda = 4;			/* For Prescan: # of disdaqs ps2*/

float xmtaddScan;

int tlead = 25us with {
    0, , 25us, INVIS, "Init deadtime",
};
float echo1bw = 16 with {
    , , , INVIS, "Echo1 filter bw.in KHz",
};
int obl_debug = 0 with {
    0, 1, 0, INVIS, "On(=1) to print messages for obloptimize",
};
int obl_method = 0 with {
    0, 1, 0, INVIS, "On(=1) to optimize the targets based on actual rotation matrices",
};

int debug = 0 with {0,1,0,INVIS,"1 if debug is on ",};

int pw_gy1_tot;       /* temp time accumulation */
float yfov_aspect = 1.0 with {0,,,INVIS, "acquired Y FOV aspect ratio to X",};

/* Trajectory cvs */
int nechoes = 17 with {1, MAXNECHOES, 17, VIS, "Number of echoes per echo train",};
int ntrains = 1 with {1, MAXNTRAINS, 1, VIS, "Number of echo trains per frame",};
int nnav = 20 with {0, 1000, 20, VIS, "Number of navigator points (must be even)",};
float R_accel = 0.5 with {0.05, , , VIS, "Spiral radial acceleration factor",};
float THETA_accel = 1.0 with {0, , 1, VIS, "Spiral angular acceleration factor",};
int sptype2d = 1 with {1, 4, 1, VIS, "1 = spiral out, 2 = spiral in, 3 = spiral out-in, 4 = spiral in-out",};
int sptype3d = 1 with {1, 4, 1, VIS, "1 = stack of spirals, 2 = rotating spirals (single axis), 3 = rotating spirals (2 axises), 4 = rotating orbitals (2 axises)",};
float SLEWMAX = 17000.0 with {, 25000.0, 17000.0, VIS, "Maximum allowed slew rate (G/cm/s)",};
float GMAX = 4.0 with {0.5, 5.0, 4.0, VIS, "Maximum allowed gradient (G/cm)",};

@inline Prescan.e PScvs


@host
/*********************************************************************
 *                      GRASS.E HOST SECTION                         *
 *                                                                   *
 * Write here the code unique to the Host PSD process. The following *
 * functions must be declared here: cvinit(), cveval(), cvcheck(),   *
 * and predownload().                                                *
 *                                                                   *
 *********************************************************************/
#include <math.h>
#include <stdlib.h>
#include "grad_rf_grass.h"

#include "psdopt.h"
#include "sar_pm.h"
#include "support_func.host.h"

/* fec : Field strength dependency library */
#include <sysDep.h>
#include <sysDepSupport.h>      /* FEC : fieldStrength dependency libraries */

#include "helperfuns.h"

@inline loadrheader.e rheaderhost

/** Load PSD Header **/
abstract("grass sequence");
psdname("grass");

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

/* Declare function prototypes from ktraj.h */
int genspiral(int N, int itr);
int genviews();
int sinsmooth(float *x, int N, int L);

/* Import functions from ktraj.h (using @inline instead of #include since
 * functions reference global variables in this file)
 */
@inline ktraj.h

@inline Prescan.e PShostVars            /* added with new filter calcs */

static char supfailfmt[] = "Support routine %s failed";


/************************************************************************/
/*       			CVINIT    				*/
/* Invoked once (& only once) when the PSD host process	is started up.	*/
/* Code which is independent of any OPIO button operation is put here.	*/
/************************************************************************/
STATUS cvinit( void )
{
	/* turn off bandwidth option - fixed! */
	oprbw = 500.0 / (float)TSP_GRAD;
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
	optr = 4500ms;
	pitrnub = 2;
	pitrval2 = PSD_MINIMUMTR;
	pitrval3 = 4500ms;

	/* te */
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

	/* hide phase (yres) option */
	piyresnub = 0;

	/* show flip angle menu */
	pifanub = 2;

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

	if( obloptimize( &loggrd, &phygrd, scan_info, exist(opslquant),
				exist(opplane), exist(opcoax), obl_method, obl_debug,
				&opnewgeo, cfsrmode ) == FAILURE )
	{
		return FAILURE;
	}


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

	/* 
	 * Calculate RF filter and update RBW:
	 *   &echo1_rtfilt: I: all the filter parameters.
	 *   exist(oprbw): I/O: desired and final allowable bw.
	 *   exist(opxres): I: output pts generated by filter.
	 *   OVERWRITE_OPRBW: oprbw will be updated.
	 */
	if( calcfilter( &echo1_rtfilt,
				exist(oprbw),
				exist(opxres),
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

	/* Calculate minimum tr */	
	avmintr = 2s;

@inline Prescan.e PScveval

	return SUCCESS;
}   /* end cveval() */

	void
getAPxParam(optval   *min,
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

	/* Set the defaults for the excitation pulse */
	a_rf1 = opflip/180.0;
	thk_rf1 = opslthick;
	res_rf1 = 1600;
	pw_rf1 = 3200;
	flip_rf1 = opflip;
	pw_gzrf1 = 3200;
	pw_gzrf1a = 300;
	pw_gzrf1d = 300;
	
	/* Generate 2d spiral */
	if (genspiral(5000, 0) == 0) {
		epic_error( use_ermes, "Error: failure to generate spiral waveform", EE_ARGS(1), EE_ARGS(0));
		return FAILURE;
	}

	/* Generate view transformations */
	if (genviews() == 0) {
		epic_error( use_ermes, "Error: failure to generate view transformation matrices", EE_ARGS(1), EE_ARGS(0));
		return FAILURE;
	}

	
	/* Set up the filter structures to be downloaded for realtime 
	   filter generation. Get the slot number of the filter in the filter rack 
	   and assign to the appropriate acquisition pulse for the right 
	   filter selection - LxMGD, RJF */
	setfilter( echo1_filt, SCAN );
/*
	filter_echo1 = echo1_filt->fslot;
*/

@inline Prescan.e PSfilter
	
	/* For Prescan: Inform 'Auto' Prescan about prescan parameters 	*/
	pislquant = opslquant;	/* # of 2nd pass slices */
	/* slquant1 = max # of locations in 1 pass */

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

	entry_point_table[L_SCAN].epxmtadd = (short)rint( (double)xmtaddScan );

	/* APS2 & MPS2 */
	entry_point_table[L_APS2] = entry_point_table[L_MPS2] = entry_point_table[L_SCAN];	/* copy scan into APS2 & MPS2 */
	(void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
	(void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );

	if( orderslice( TYPNCAT, (int)opslquant, (int)1, TRIG_INTERN ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "orderslice" );
		return FAILURE;
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
	acqs = opslquant;	/* Fixes the # of rhnpasses to the # of passes */
	acq_type = TYPGRAD;
@inline loadrheader.e rheaderinit   /* Recon variables */

	scalerotmats( rsprot, &loggrd, &phygrd, (int)(opslquant), obl_debug );

@inline Prescan.e PSpredownload

	return SUCCESS;
}   /* end predownload() */


@inline Prescan.e PShost


@pg
/*********************************************************************
 *                   GRASS.E PULSEGEN SECTION                        *
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
	
	fprintf(stderr, "Generating RF1 pulse (90deg tipdown) with a_rf1 = %.2f, pw_rf1 = %d... ",
		a_rf1, pw_rf1);
	SLICESELZ(rf1, 30ms, 3200, opslthick, opflip, 1, 1, loggrd);
	fprintf(stderr, " done.\n");
	
	fprintf(stderr, "Generating readout x gradient with a_gx = %.2f, pw_gx = %d, res_gx = %d... ",
		a_gx, pw_gx, res_gx);
	INTWAVE(XGRAD, gx, pend( &gzrf1d, "gzrf1d", 0), 1.0, 1000, 4000, Gx, 0, loggrd);
	fprintf(stderr, " done.\n");

	fprintf(stderr, "Generating readout y gradient with a_gy = %.2f, pw_gy = %d, res_gy = %d... ",
		a_gy, pw_gy, res_gy);
	INTWAVE(YGRAD, gy, pend( &gzrf1d, "gzrf1d", 0), 1.0, 1000, 4000, Gy, 0, loggrd);
	fprintf(stderr, " done.\n");

	fprintf(stderr, "Generating readout z gradient with a_gz = %.2f, pw_gz = %d, res_gz = %d... ",
		a_gz, pw_gz, res_gz);
	INTWAVE(ZGRAD, gz, pend( &gzrf1d, "gzrf1d", 0), 1.0, 1000, 4000, Gz, 0, loggrd);
	fprintf(stderr, " done.\n");
	
	SEQLENGTH(seqcore, optr, seqcore); /* set the sequence length to optr */

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
 *                     GRASS.E RSPVAR SECTION                        *
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

int echon;
int trainn;
int view;
int excitation;
int dabop;
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
 *                     GRASS.E RSP SECTION                           *
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

int *rf1_freq;
int *receive_freq1;


STATUS psdinit( void )
{
	/* Initialize everything to a known state */
	echon = trainn = 0;
	setrfconfig( (short)5 );	/* Only activate rho1 */
	setssitime( 100 );		/* Set ssi counter to 400 us. */
	rspqueueinit( 200 );	/* Initialize to 200 entries */
	scopeon( &seqcore );	/* Activate scope for core */
	syncon( &seqcore );		/* Activate sync for core */
	syncoff( &pass );		/* Deactivate sync during pass */
	seqCount = 0;		/* Set SPGR sequence counter */
	settriggerarray( (short)1, rsptrigger );
	setrotatearray( (short)ntrains*nechoes, rsprot[0] );
/*
	setrfltrs( (int)filter_echo1, &echo1 );
*/
	return SUCCESS;
}   /* end psdinit() */


@inline Prescan.e PScore


/* For Prescan: MPS2 Function */
STATUS mps2( void )
{
	boffset(off_seqcore);
	/* Initialize RSP parameters */
	rspent = L_MPS2;	
	rspdda = 2;
	rspbas = 0;
	rspvus = 30000;
	rspgy1 = 0;
	rspnex = 2;
	rspesl = -1;
	rspasl = pre_slice;
	rspslq = ntrains*nechoes;
	rspsct = 0;

	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	prescanCore();
	rspexit();

	return SUCCESS;
}   /* end mps2() */


/* For Prescan: APS2 Function */
STATUS aps2( void )
{
	boffset(off_seqcore);
	/* Initialize RSP parameters */
	rspent = L_APS2;	
	rspdda = 2;
	rspbas = 0;
	rspvus = 1024;
	rspgy1 = 0;
	rspnex = 2;
	rspesl = -1;
	rspasl = -1;
	rspslq = ntrains*nechoes;
	rspsct = 0;

	if( psdinit() == FAILURE )
	{
		return rspexit();
	}

	prescanCore();
	rspexit();

	return SUCCESS;
}   /* end aps2() */


/************************************************************************/
/* 				SCAN FUNCTION				*/
/************************************************************************/
STATUS scan( void )
{

	int echon, trainn;

	if( psdinit() == FAILURE )
	{
		return rspexit();
	}
	boffset( off_seqcore );	/* start the hardware in the 'core' sequence */
	setrotatearray( ntrains*nechoes, rsprot[0] );
	settriggerarray( (short)1, rsptrigger );
	setssitime( 250 );	/* allow time to update sequencer memory */
    
	/* Calculate the RF & slice frequencies */
	rf1_freq = (int *)AllocNode( opslquant * sizeof(int) );
	receive_freq1 = (int *)AllocNode( opslquant * sizeof(int) );

	/* Set the Slice Frequency */
	setupslices( rf1_freq, rsp_info, opslquant, a_gzrf1, (float)1, opfov,
			TYPTRANSMIT );
	setupslices( receive_freq1, rsp_info, opslquant,(float)0, echo1bw, opfov,
			TYPREC);

	setiamp( ia_rf1,  &rf1, 0 );

	for (trainn = 0; trainn < ntrains; trainn++) {
		for (echon = 0; echon < nechoes; echon++) {
			fprintf(stderr, "scan(): playing train %d/%d, echo %d/%d... ",
				trainn + 1, ntrains, echon + 1, nechoes);
			
			settrigger(TRIG_INTERN,0);
			setrotate((s32 *)rsprot, trainn*nechoes + echon);
			
			boffset( off_seqcore ); /* reset the hardware in the 'core' sequence */
			startseq( 0, (short)MAY_PAUSE );
			boffset( off_seqcore ); /* reset the hardware in the 'core' sequence */
			fprintf(stderr, "done.\n");

		}
	}

	rspexit();

	return SUCCESS;
}   /* end scan() */


/*************************************************/
/* Runtime core section executed during prescan. */
/*************************************************/
STATUS prescanCore( void )
{
	setiamp( 0, &gx, 0 );
	setiamp( 0, &gy, 0 );
	setiamp( 0, &gz, 0 );
	if( scan() == FAILURE )
	{
		return rspexit();
	}
	setiamp( ia_gx, &gx, 0 );
	setiamp( ia_gy, &gy, 0 );
	setiamp( ia_gz, &gx, 0 );

	return SUCCESS;
}   /* end prescanCore() */


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

/************************ END OF GRASS.E ******************************/

