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

int pos_start = 0 with {
    0, , , INVIS, "Start time for sequence. ",
};
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
int endview_iamp; /* last instruction phase amp */
float endview_scale; /* ratio of last instruction amp to maximum value */

/* Trajectory cvs */
int nechoes = 17 with {1, MAXNECHOES, 17, VIS, "Number of echoes per echo train",};
int ntrains = 1 with {1, MAXNTRAINS, 1, VIS, "Number of echo trains per frame",};
int nnav = 20 with {0, 1000, 20, VIS, "Number of navigator points (must be even)",};
float R_accel = 0.5 with {0.05, , , VIS, "Spiral radial acceleration factor",};
float THETA_accel = 1.0 with {0, , 1, VIS, "Spiral angular acceleration factor",};
int sptype2d = 1 with {1, 4, 1, VIS, "1 = spiral out, 2 = spiral in, 3 = spiral out-in, 4 = spiral in-out",};
int sptype3d = 1 with {1, 4, 1, VIS, "1 = stack of spirals, 2 = rotating spirals (single axis), 3 = rotating spirals (2 axises), 4 = rotating orbitals (2 axises)",};
float SLEWMAX = 17000.0 with {5000.0, 25000.0, 17000.0, VIS, "Maximum allowed slew rate (G/cm/s)",};
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

/* Declare gradient waveform arrays */
int Gx[MAXWAVELEN];
int Gy[MAXWAVELEN];
int Gz[MAXWAVELEN];
float a_Gx;
float a_Gy;
float a_Gz;
int pw_G = 5000;

/* Declare view transformation table */
int T_v[MAXNTRAINS*MAXNECHOES][9];

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
	opte = PSD_MINFULLTE;
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

	/* Generate 2d spiral */
	if (genspiral(pw_G, 0) == 0) {
		epic_error( use_ermes, "failure to generate spiral waveform", EE_ARGS(1), EE_ARGS(0));
		return FAILURE;
	}

	/* Generate view transformations */
	genviews();

	/*
	 * The minimum TR is based on the time before the RF pulse +
	 * half the RF pulse + the TE time + the last half of the
	 * readout + the time for the end of sequence killers
	 */
	avmintr = 1ms + pw_rf1 / 2 + exist(opte) + echo1_rtfilt.tdaq / 2 + 2ms;

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
	STATUS
cvcheck( void )
{
	if( exist(optr) < avmintr )
	{
		int min_tr = (int)ceil( (double)avmintr / 1ms );

		epic_error( use_ermes,
				"The TR needs to be increased to %d ms for the current prescription.",
				EM_PSD_TR_OUT_OF_RANGE, EE_ARGS(1), INT_ARG, min_tr );
		return ADVISORY_FAILURE;
	}

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
	STATUS
predownload( void )
{
	/* Set the defaults for the excitation pulse */
	a_rf1 = opflip/180.0;
	thk_rf1 = opslthick;
	res_rf1 = 320;
	flip_rf1 = opflip;

	/* Set the phase encode amplitude*/

	yfov_aspect = nop*exist(opphasefov);
	rhnframes = opyres*fn*yfov_aspect;

	if( endview( (int)(rhnframes / fn), &endview_iamp ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "endview" );
		return FAILURE;
	}

	endview_scale = (float)max_pg_iamp / (float)endview_iamp;

	if( amppwencode( &grady[GY1_SLOT], &pw_gy1_tot,
				FMin( 2, loggrd.ty_xyz,loggrd.ty / endview_scale ),
				(int)loggrd.yrt,
				(float)(nop * opfov * opphasefov), (int)(rhnframes / fn ),
				(float)0.0 /* offset area */ ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "amppwencode:gy1" );
		return FAILURE;
	}

	grady[GY1_SLOT].num = 1;

	/* Set the Read Out gradient amplitude */
	if( ampfov( &a_gxw, echo1_filt->bw, opfov ) == FAILURE )
	{
		epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
				EE_ARGS(1), STRING_ARG, "ampfov" );
		return FAILURE;
	}

	/* Duration of read lobe to match acquisition interval */
	pw_gxw = echo1_filt->tdaq;

	a_gy1a = -a_gy1a;
	a_gy1b = -a_gy1b;
	a_gy1 = -0.9;     /* to make it visible in plotpulse */

	a_gyr1a = -a_gy1a;
	a_gyr1b = -a_gy1b;
	a_gyr1 = a_gy1;   /* to make it visible in plotpulse */

	pw_gyr1 = pw_gy1;
	pw_gyr1a = pw_gy1a;
	pw_gyr1d = pw_gy1d;

/*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
/*********************************************************************/

	/* Set up the filter structures to be downloaded for realtime 
	   filter generation. Get the slot number of the filter in the filter rack 
	   and assign to the appropriate acquisition pulse for the right 
	   filter selection - LxMGD, RJF */
	setfilter( echo1_filt, SCAN );
	filter_echo1 = echo1_filt->fslot;

@inline Prescan.e PSfilter

		/* Set the Slope of the Read Out window's leading edge */
		if( optramp( &pw_gxwa, a_gxw, loggrd.tx, loggrd.xrt, TYPDEF ) == FAILURE )
		{
			epic_error( use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
					EE_ARGS(1), STRING_ARG, "optramp" );
			return FAILURE;
		}
	pw_gxwd = pw_gxwa;		/* Set trailing edge ramp to same duration. */

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

	ia_rf1 = max_pg_iamp * (*rfpulse[RF1_SLOT].amp);
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

	/***************************** Gradient Structure Initializations ******************************/

	gradx[GX1_SLOT].num = 1;
	gradx[GXW2_SLOT].num = 1;
	grady[GY1_SLOT].num = 1;
	gradz[GZRF1_SLOT].num = 1;
	gradz[GZ1_SLOT].num = 1;

	avepepowscale(&(grady[GY1_SLOT].scale), rhnframes, rhnframes/2);

	gradx[GX1_SLOT].powscale = 1.0;
	gradx[GXW2_SLOT].powscale = 1.0;
	grady[GY1_SLOT].powscale = 1.0;
	gradz[GZRF1_SLOT].powscale = 1.0;
	gradz[GZ1_SLOT].powscale = 1.0;





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


	STATUS
pulsegen( void )
{
	sspinit(psd_board_type);

	/* RF wave */
	SLICESELZ(rf1, 1ms, 3200us, opslthick, opflip, 1, , loggrd); 

	/* Z Dephaser */
	TRAPEZOID(ZGRAD, gz1, pend( &gzrf1d, "gzrf1d", 0 ) + pw_gz1a, (int)(-0.5 * a_gzrf1 * (pw_rf1 + pw_gzrf1d)), , loggrd);

	/* X Readout */
	TRAPEZOID(XGRAD, gxw, RUP_GRD(pmid( &gzrf1, "gzrf1", 0 ) + opte - pw_gxw / 2), 0, TYPNDEF, loggrd);

	/* Frequency Dephaser */
	TRAPEZOID(XGRAD, gx1, pbeg( &gxwa, "gxwa", 0 ) - pw_gx1 - pw_gx1d, (int)(-0.5 * a_gxw * (pw_gxw + pw_gxwa)), , loggrd);

	/* Phase Encoding */
	TRAPEZOID2(YGRAD, gy1, RDN_GRD(pend(&rf1,"rf1",0) + rfupd), TRAP_ALL_SLOPED, , , endview_scale, loggrd);

	/* Data Acquisition */
	ACQUIREDATA(echo1, pbeg( &gxw, "gxw", 0 ), , , );

	TRAPEZOID2(YGRAD, gyr1, pend( &gxw, "gxw", 0), TRAP_ALL_SLOPED, , , endview_scale, loggrd);

	/* Z & X Killers */
	TRAPEZOID(ZGRAD, gzk, pend( &gxwd, "gxwd", 0 ) + pw_gzka, 980, , loggrd);
	TRAPEZOID(XGRAD, gxk, pend( &gxwd, "gxwd", 0 ) + pw_gxka, 980, , loggrd);

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

short viewtable[513];
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
	view = slice = excitation = 0;
	setrfconfig( (short)5 );	/* Only activate rho1 */
	setssitime( 100 );		/* Set ssi counter to 400 us. */
	rspqueueinit( 200 );	/* Initialize to 200 entries */
	scopeon( &seqcore );	/* Activate scope for core */
	syncon( &seqcore );		/* Activate sync for core */
	syncoff( &pass );		/* Deactivate sync during pass */
	seqCount = 0;		/* Set SPGR sequence counter */
	settriggerarray( (short)slquant1, rsptrigger );
	setrotatearray( (short)slquant1, rsprot[0] );
	setrfltrs( (int)filter_echo1, &echo1 );

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
	rspslq = slquant1;
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
	rspslq = slquant1;
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
	if( psdinit() == FAILURE )
	{
		return rspexit();
	}
	boffset( off_seqcore );	/* start the hardware in the 'core' sequence */
	setrotatearray( opslquant, rsprot[0] );
	settriggerarray( opslquant, rsptrigger );
	setssitime( 250 );	/* allow time to update sequencer memory */

	/* Calculate the RF & slice frequencies */
	rf1_freq = (int *)AllocNode( opslquant * sizeof(int) );
	receive_freq1 = (int *)AllocNode( opslquant * sizeof(int) );

	/* Set the Slice Frequency */
	setupslices( rf1_freq, rsp_info, opslquant, a_gzrf1, (float)1, opfov,
			TYPTRANSMIT );
	setupslices( receive_freq1, rsp_info, opslquant,(float)0, echo1bw, opfov,
			TYPREC);

	setiamp( ia_rf1, &rf1, 0 );
	setupphasetable( viewtable, TYPNORM, opyres );

	/* The SLICE loop */
	for( slice = 0; slice < opslquant; ++slice )
	{
		setfrequency( rf1_freq[slice], &rf1, 0 );
		setfrequency( receive_freq1[slice], &echo1, 0 );

		/* equilibrium views */
		dabop = 0;
		loaddab( &echo1, 0, 0, dabop, (int)0, DABOFF, PSD_LOAD_DAB_ALL );  /* each slice is a pass, slice index in each pass should be 0 */
		setiampt( viewtable[1], &gy1, 0 );
		setiampt( viewtable[1], &gyr1, 0 );  

		for( view = 0; view < 4; ++view )
		{
			startseq( 0, (short)MAY_PAUSE );
			getiamp( &chopamp, &rf1, 0 );
			setiamp( -chopamp, &rf1, 0 );
		}

		for( view = 0; view < rhbline; ++view )
		{
			loaddab( &echo1, 0, 0, dabop, (int)0, DABON, PSD_LOAD_DAB_ALL );  /* each slice is a pass, slice index in each pass should be 0 */
			startseq( 0, (short)MAY_PAUSE );
			getiamp( &chopamp, &rf1, 0 );
			setiamp( -chopamp, &rf1, 0 );
			dabop = 1; /* accumulate base views */
		}

		/* The VIEW loop */
		for( view = 1; view < (opyres + 1); ++view )
		{
			for( excitation = 1; excitation <= opnex; ++excitation )
			{
				if( excitation == 1 )
				{
					dabop = 0;
				}
				else
				{
					dabop = 3 - 2 * (excitation % 2);
				}

				setiampt( viewtable[view], &gy1, 0 );
				setiampt( viewtable[view], &gyr1, 0 );
				/* loaddab loads SSP packet used by data acquisition */
				loaddab( &echo1, 0, 0, dabop, view, DABON, PSD_LOAD_DAB_ALL );  /* each slice is a pass, slice index in each pass should be 0 */

				startseq( 0, (short)MAY_PAUSE );
				getiamp( &chopamp, &rf1, 0 );
				setiamp( -chopamp, &rf1, 0 );

			} /* end-of-excitation loop */
		}  /* end-of-view loop */

		boffset( off_pass );
		if( slice == (opslquant - 1) ) /* Last pass */
		{
			/* Set DAB pass packet to end of scan */
			setwamp( SSPD + DABPASS + DABSCAN, &endpass, 2 );
		}
		else
		{
			/* Set DAB pass packet to end of pass */
			setwamp( SSPD + DABPASS, &endpass, 2 );
		}
		startseq( 0, (short)MAY_PAUSE );

		boffset( off_seqcore ); /* reset the hardware in the 'core' sequence */
	} /* end-of-slice loop */

	rspexit();

	return SUCCESS;
}   /* end scan() */


/*************************************************/
/* Runtime core section executed during prescan. */
/*************************************************/
STATUS prescanCore( void )
{

	/*
	 * Core rsp routine for prescan entry points. Same as scan, only
	 *  no PE gradients or chopping. 
	 */
	boffset( off_seqcore );

	if( psdinit() == FAILURE )
	{
		rspexit();
	}

	setrotatearray( 1, rsprot[0] );
	settriggerarray( 1, rsptrigger );

	rf1_freq = (int *)AllocNode( sizeof(int) );
	receive_freq1 = (int *)AllocNode( sizeof(int) );

	(*rf1_freq) = (int)(GAM * a_gzrf1 * rsp_info[rspesl].rsptloc / (10 * TARDIS_FREQ_RES));
	(*receive_freq1) = (int)((float)cfreceiveroffsetfreq / TARDIS_FREQ_RES);

	setiamp( ia_rf1, &rf1, 0 );

	setfrequency( (*rf1_freq), &rf1, 0 );
	setfrequency( (*receive_freq1), &echo1, 0 );

	dabop = 0;

	loaddab( &echo1, 0, 0, dabop, (int) 0, DABOFF, PSD_LOAD_DAB_ALL );

	setiampt( 0, &gy1, 0 );

	for( view = 0; view < numdda; ++view ) 
	{
		startseq( 0, (short)MAY_PAUSE );
		getiamp( &chopamp, &rf1, 0 );
		setiamp( -chopamp, &rf1, 0 );
	}

	for( view = 1; view < rspvus + 1; ++view ) 
	{
		loaddab( &echo1, 0, 0, dabop, view, DABON, PSD_LOAD_DAB_ALL );
		startseq( 0, (short)MAY_PAUSE );
		getiamp( &chopamp, &rf1, 0 );
		setiamp( -chopamp, &rf1, 0 );
	} /* views */

	rspexit();

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

