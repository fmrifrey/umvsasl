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
#ifdef psdutil
#include "psdutil.h"
#endif
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
#define TIMESSI 400 /* SSP instruction time */

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

/* Define temporary error message string */
char *tmpstr;

/* Declare sequencer hardware limit variables */
float XGRAD_max;
float YGRAD_max;
float ZGRAD_max;
float THETA_max;

/* Declare gradient waveform arrays */
int Gx[MAXWAVELEN];
int Gy[MAXWAVELEN];
int Gz[MAXWAVELEN];
int grad_len = 5000;

/* Declare table of rotation matrices */
long tmtxtbl[MAXNTRAINS*MAXNECHOES][9];

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

@inline Prescan.e PScvs

int numdda = 4;			/* For Prescan: # of disdaqs ps2*/

float xmtaddScan;
int obl_debug = 0 with {0, 1, 0, INVIS, "On(=1) to print messages for obloptimize",};
int obl_method = 0 with {0, 1, 0, INVIS, "On(=1) to optimize the targets based on actual rotation matrices",};
int debug = 0 with {0,1,0,INVIS,"1 if debug is on ",};
float echo1bw = 16 with {,,,INVIS,"Echo1 filter bw.in KHz",};

/* FSE timing cvs */
int trapramptime = 100 with {100, , 100, INVIS, "Trapezoidal gradient ramp time (us)",};
int gradbufftime = TIMESSI with {TIMESSI, , TIMESSI, INVIS, "Gradient IPG buffer space (us)",};
int tadjust = 0 with {0, , 0, INVIS, "Deadtime after readout to fill the TR (us)",};

/* Trajectory cvs */
int nechoes = 16 with {1, MAXNECHOES, 17, VIS, "Number of echoes per echo train",};
int ntrains = 1 with {1, MAXNTRAINS, 1, VIS, "Number of echo trains per frame",};
int nframes = 1 with {1, , 1, VIS, "Number of frames to acquire",};
int nnav = 20 with {0, 1000, 20, VIS, "Number of navigator points (must be even)",};
float R_accel = 0.5 with {0.05, , , VIS, "Spiral radial acceleration factor",};
float THETA_accel = 1.0 with {0, , 1, VIS, "Spiral angular acceleration factor",};
int sptype2d = 4 with {1, 4, 1, VIS, "1 = spiral out, 2 = spiral in, 3 = spiral out-in, 4 = spiral in-out",};
int sptype3d = 3 with {1, 4, 1, VIS, "1 = stack of spirals, 2 = rotating spirals (single axis), 3 = rotating spirals (2 axises), 4 = rotating orbitals (2 axises)",};
float SLEWMAX = 17000.0 with {1000, 25000.0, 17000.0, VIS, "Maximum allowed slew rate (G/cm/s)",};
float GMAX = 4.0 with {0.5, 5.0, 4.0, VIS, "Maximum allowed gradient (G/cm)",};
int kill_grads = 0 with {0, 1, 0, VIS, "Option to turn off readout gradients",};

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
#include "helperfuns.h"

/* fec : Field strength dependency library */
#include <sysDep.h>
#include <sysDepSupport.h>      /* FEC : fieldStrength dependency libraries */

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
STATUS
cvinit( void )
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
	cvmin(optr, avmintr);
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
	cvmax(opflip, 130);
	cvdef(opflip, 90);
	opflip = 90;
	pitinub = 0;

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
STATUS
cveval( void )
{

    configSystem();
    InitAdvPnlCVs();

	pititle = 1;
	cvdesc(pititle, "Advanced pulse sequence parameters");
	
	piuset = use0;
	cvdesc(opuser0, "Number of frames to acquire");
	cvdef(opuser0, 1);
	opuser0 = 1;
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
	cvdesc(opuser3, "Initial spiral trajectory type");
	cvdef(opuser3, 1);
	opuser3 = 4;
	cvmin(opuser3, 1);
	cvmax(opuser3, 4);
	sptype2d = opuser3;

	piuset += use4;
	cvdesc(opuser4, "3d trajectory transformations type");
	cvdef(opuser4, 3);
	opuser4 = 1;
	cvmin(opuser4, 1);
	cvmax(opuser4, 4);
	sptype3d = opuser4;

	/* Get sequencer hardware limits */
	gettarget(&XGRAD_max, XGRAD, &loggrd);
	gettarget(&YGRAD_max, YGRAD, &loggrd);
	gettarget(&ZGRAD_max, ZGRAD, &loggrd);
	gettarget(&THETA_max, THETA, &loggrd);

	cvmax(rhfrsize, 32767);
	rhfrsize = grad_len;
	cvmax(rhnframes, 32767);
	rhnframes = 2*ceil((nframes + 1)/2);
	rhrawsize = 2*rhptsize*rhfrsize * nframes * nechoes * ntrains;

	rhrcctrl = 1;
	rhexecctrl = 11;
	autolock = 1;
	
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

	a_rf1 = (float)opflip / (float)180;
	a_rf2 = 2.0 * a_rf1;

	ia_rf1 = a_rf1 * MAX_PG_IAMP;
	ia_rf2 = a_rf2 * MAX_PG_IAMP;

	/* Calculate minimum te */
	avminte = pw_rf2 + 4*gradbufftime + 4*trapramptime + pw_gzrf2crush1 + pw_gzrf2crush2 + GRAD_UPDATE_TIME*grad_len;

	/* Calculate minimum tr */
	avmintr = ((float)nechoes + 0.5) * opte;
	
	/* Calculate adjust time */
	tadjust = optr - avmintr;

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
	/* Check if TE is valid */
	if (opte < avminte) {
		epic_error(use_ermes, "opte must be >= %dus", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), INT_ARG, avminte);
		return FAILURE;
	};
	
	/* Check if TR is valid */
	if (optr < avmintr) {
		epic_error(use_ermes, "optr must be >= %dus", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), INT_ARG, avmintr);
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
STATUS
predownload( void )
{
    /*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
    /*********************************************************************/

	/* Set the parameters for the spin echo tipdown pulse */
	a_rf1 = opflip/180.0;
	thk_rf1 = opslthick*opslquant;
	res_rf1 = 1600;
	pw_rf1 = 3200;
	flip_rf1 = opflip;
	pw_gzrf1 = 3200;
	pw_gzrf1a = trapramptime;
	pw_gzrf1d = trapramptime;
	
	/* Set the parameters for the tipdown crusher */
	a_gzrf1crush = 1.5;
	pw_gzrf1crush = 2*trapramptime + 4;
	pw_gzrf1crusha = trapramptime;
	pw_gzrf1crushd = trapramptime;
	
	/* Set the parameters for the pre-refocuser crusher */
	a_gzrf2crush1 = 1.5;
	pw_gzrf2crush1 = 2*trapramptime + 4;
	pw_gzrf2crush1a = trapramptime;
	pw_gzrf2crush1d = trapramptime;
	
	/* Set the parameters for the refocuser pulse */
	a_rf2 = 2.0*opflip/180.0;
	thk_rf2 = opslthick*opslquant;
	res_rf2 = 1600;
	pw_rf2 = 3200;
	flip_rf2 = opflip;
	pw_gzrf2 = 3200;
	pw_gzrf2a = trapramptime;
	pw_gzrf2d = trapramptime;
	
	/* Set the parameters for the post-refocuser crusher */
	a_gzrf2crush2 = 1.5;
	pw_gzrf2crush2 = 2*trapramptime + 4;
	pw_gzrf2crush2a = trapramptime;
	pw_gzrf2crush2d = trapramptime;
	
	/* Update the pulse parameters */
	a_gxw = XGRAD_max;
	a_gyw = YGRAD_max;
	a_gzw = ZGRAD_max;
	res_gxw = grad_len;
	res_gyw = grad_len;
	res_gzw = grad_len;
	pw_gxw = GRAD_UPDATE_TIME*grad_len;
	pw_gyw = GRAD_UPDATE_TIME*grad_len;
	pw_gzw = GRAD_UPDATE_TIME*grad_len;

    /* Set up the filter structures to be downloaded for realtime 
       filter generation. Get the slot number of the filter in the filter rack 
       and assign to the appropriate acquisition pulse for the right 
       filter selection - LxMGD, RJF */
    setfilter( echo1_filt, SCAN );
    filter_echo1 = echo1_filt->fslot;
    
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
    entry_point_table[L_SCAN].epprexres = grad_len;

    (void)strcpy( entry_point_table[L_APS2].epname, "aps2" );
    entry_point_table[L_APS2].epfilter = (unsigned char)echo1_filt->fslot;
    entry_point_table[L_APS2].epprexres = grad_len;

    (void)strcpy( entry_point_table[L_MPS2].epname, "mps2" );
    entry_point_table[L_MPS2].epfilter = (unsigned char)echo1_filt->fslot;
    entry_point_table[L_MPS2].epprexres = grad_len;

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
	
	/* Generate initial spiral trajectory */
	fprintf(stderr, "cveval(): calling genspiral()\n");
	if (genspiral(grad_len, 0) == 0) {
		epic_error(use_ermes,"failure to generate spiral waveform", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}
	
	/* Generate view transformations */
	fprintf(stderr, "cveval(): calling genviews()\n");
	if (genviews() == 0) {
		epic_error(use_ermes,"failure to generate view transformation matrices", EM_PSD_SUPPORT_FAILURE, EE_ARGS(0));
		return FAILURE;
	}

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
	
	/* initialize temporary time marker */
	int tmploc;


	/********************************/
	/* Generate spiral readout core */
	/********************************/
	/* Calculate start of spiral */
	tmploc = (opte - (pw_gzrf2crush1 + pw_rf2 + pw_gzrf2crush2 + 4*gradbufftime + 4*trapramptime) - GRAD_UPDATE_TIME*grad_len)/2;

	fprintf(stderr, "pulsegen(): generating readout x gradient using INTWAVE()...\n");
	INTWAVE(XGRAD, gxw, tmploc, XGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gx, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gxw, "gxw", 0), pend( &gxw, "gxw", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating readout y gradient using INTWAVE()...\n");
	INTWAVE(YGRAD, gyw, tmploc, YGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gy, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gyw, "gyw", 0), pend( &gyw, "gyw", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating readout z gradient using INTWAVE()...\n");
	INTWAVE(ZGRAD, gzw, tmploc, ZGRAD_max, grad_len, GRAD_UPDATE_TIME*grad_len, Gz, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzw, "gzw", 0), pend( &gzw, "gzw", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating data acquisition instruction using ACQUIREDATA()...\n");
	ACQUIREDATA(echo1, tmploc + psd_grd_wait,,,);
	fprintf(stderr, "\tDone.\n");
	
	/* Calculate length of core */
	tmploc = opte - (pw_gzrf2crush1 + pw_rf2 + pw_gzrf2crush2 + 4*gradbufftime);

	fprintf(stderr, "pulsegen(): finalizing spiral readout core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", tmploc);
	SEQLENGTH(seqcore, tmploc, seqcore);
	fprintf(stderr, "\tDone.\n");


	/***********************************/
	/* Generate spin echo tipdown core */
	/***********************************/	
	fprintf(stderr, "pulsegen(): beginning pulse generation of spin echo tipdown core (tipdowncore)\n");

	fprintf(stderr, "pulsegen(): generating rf1 (90deg tipdown pulse)...\n");
	SLICESELZ(rf1, trapramptime, 3200, opslthick, opflip, 1, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf1a, "gzrf1a", 0), pend( &gzrf1d, "gzrf1d", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating rf1 crusher...\n");
	TRAPEZOID(ZGRAD, gzrf1crush, pend( &gzrf1d, "gzrf1d", 0 ) + gradbufftime, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf1crusha, "gzrf1crusha", 0), pend( &gzrf1crushd, "gzrf1crushd", 0));
	fprintf(stderr, "\tDone.\n");	

	/* Calculate length of core */	
	tmploc = opte/2 - pw_rf2/2 - 2*gradbufftime - pw_gzrf2crush1;

	fprintf(stderr, "pulsegen(): finalizing spin echo tipdown core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", tmploc);
	SEQLENGTH(tipdowncore, tmploc, tipdowncore);
	fprintf(stderr, "\tDone.\n");
	

	/*************************************/
	/* Generate spin echo refocuser core */
	/*************************************/
	fprintf(stderr, "pulsegen(): beginning pulse generation of spin echo refocuser and readout core (seqcore)\n");

	fprintf(stderr, "pulsegen(): generating pre-rf2 crusher...\n");
	TRAPEZOID(ZGRAD, gzrf2crush1, gradbufftime + trapramptime, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2crush1a, "gzrf2crush1a", 0), pend( &gzrf2crush1d, "gzrf2crush1d", 0));
	fprintf(stderr, "\tDone.\n");

	fprintf(stderr, "pulsegen(): generating rf2 (180 deg spin echo refocuser)...\n");
	SLICESELZ(rf2, pend( &gzrf2crush1d, "gzrf2crush1d", 0) + gradbufftime + trapramptime, 3200, opslthick, opflip, 1, 1, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2a, "gzrf2a", 0), pend( &gzrf2d, "gzrf2d", 0));
	fprintf(stderr, "\tDone.\n");	

	fprintf(stderr, "pulsegen(): generating post-rf2 crusher...\n");
	TRAPEZOID(ZGRAD, gzrf2crush2, pend( &gzrf2d, "gzrf2d", 0) + gradbufftime + trapramptime, 3200, 0, loggrd);
	fprintf(stderr, "\tstart: %dus, end: %dus\n", pbeg( &gzrf2crush2a, "gzrf2crush2a", 0), pend( &gzrf2crush2d, "gzrf2crush2d", 0));
	fprintf(stderr, "\tDone.\n");

	/* Calculate length of core */
	tmploc = pw_gzrf2crush1 + pw_rf2 + pw_gzrf2crush2 + 4*gradbufftime + 4*trapramptime;

	fprintf(stderr, "pulsegen(): finalizing spin echo refocuser core...\n");
	fprintf(stderr, "\ttotal time: %dus\n", tmploc);
	SEQLENGTH(refocuscore, tmploc, refocuscore);
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

/* Declare rsps */
int ispre;
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

/* Transmit & receive frequencies */
int *rf1_freq;
int xmitfreq;
int *receive_freq1;
int recfreq;

/* Initial transformation matrix */
long tmtx0[9];

STATUS
psdinit( void )
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
    scalerotmats(tmtxtbl, &loggrd, &phygrd, ntrains*nechoes, 0);

    /* Store initial transformation matrix */
    getrotate(tmtx0, 0);

    return SUCCESS;
}   /* end psdinit() */


@inline Prescan.e PScore

STATUS
scancore( void )
{
	
	/* Set transmit frequency and phase */
	rf1_freq = (int *) AllocNode(opslquant*sizeof(int));
	setupslices(rf1_freq, rsp_info, opslquant, a_gzrf1, 1.0, opfov, TYPTRANSMIT);
	xmitfreq = (int)((rf1_freq[opslquant/2] + rf1_freq[opslquant/2-1])/2);
	setfrequency(xmitfreq, &rf2, 0);
	setiphase(0, &rf1, 0);

	/* Set receiver frequency and phase */
	receive_freq1 = (int *) AllocNode(opslquant*sizeof(int));
	setupslices(receive_freq1, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);
	recfreq = (int)((receive_freq1[opslquant/2] + receive_freq1[opslquant/2-1])/2);
	setfrequency(recfreq , &echo1, 0);
	setphase(0.0, &echo1, 0);
	
	if (ispre || kill_grads) {
		/* Turn off the gradients */
		setiamp(0, &gxw, 0);
		setiamp(0, &gyw, 0);
		setiamp(0, &gzw, 0);
	}
	else {
		/* Restore the gradients */
		setiamp(MAX_PG_IAMP, &gxw, 0);
		setiamp(MAX_PG_IAMP, &gyw, 0);
		setiamp(MAX_PG_IAMP, &gzw, 0);
	}

	/* Loop through frames */
	for (framen = 0; framen < ((ispre) ? (1) : (nframes)); framen++) {
		/* Loop through echo trains */
		for (trainn = -rspdda; trainn < ((ispre) ? (1000) : (ntrains)); trainn++) {
			/* Play tipdown core */
			boffset(off_tipdowncore);
			startseq(0, MAY_PAUSE);
			settrigger(TRIG_INTERN, 0);

			/* Play readout (refocusers + spiral gradients */
			for (echon = 0; echon < nechoes; echon++) {
				/* Play the refocuser core */
				boffset(off_refocuscore);
				startseq(echon, (short)MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Allocate space in memory for data using loaddab */			
				if (trainn < 0) /* DISDAQs */
					loaddab(&echo1, 0, 0, DABSTORE, 0, DABOFF, PSD_LOAD_DAB_ALL);
				else if (ispre) /* Prescans (after DISDAQ) */
					loaddab(&echo1, 0, 0, DABSTORE, 0, DABON, PSD_LOAD_DAB_ALL);
				else /* Scan */
					loaddab(&echo1, echon, 0, DABSTORE, framen*ntrains + trainn + 1, DABON, PSD_LOAD_DAB_ALL);

				/* Set the transformation matrix */
				setrotate(tmtxtbl[(trainn < 0) ? (0) : (trainn*nechoes + echon)], echon);

				/* Play the readout core */
				boffset(off_seqcore);
				startseq(echon, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);

				/* Reset the rotation matrix */
				setrotate(tmtx0, echon);
			}

			if (tadjust > 0) {
				/* Set length of emptycore to tadjust */
				setperiod(tadjust, &emptycore, 0);

				/* Play TR deadtime (emptycore) */
				boffset(off_emptycore);
				startseq(0, MAY_PAUSE);
				settrigger(TRIG_INTERN, 0);
			}

		}
	}

	/* Send SSP packet to end scan */
	boffset( off_pass );
	setwamp(SSPD + DABPASS + DABSCAN, &endpass, 2);
	startseq(0, MAY_PAUSE);
	settrigger(TRIG_INTERN, 0);

	rspexit();

    return SUCCESS;
}

/* For Prescan: MPS2 Function */
STATUS
mps2( void )
{
    if( psdinit() == FAILURE )
    {
        return rspexit();
    }

    ispre = 1;
    rspdda = 0;
    scancore();
    rspexit();

    return SUCCESS;
}   /* end mps2() */


/* For Prescan: APS2 Function */
STATUS
aps2( void )
{   
    if( psdinit() == FAILURE )
    {
        return rspexit();
    }

    ispre = 1;
    rspdda = 2;
    scancore();
    rspexit();

    return SUCCESS;
}   /* end aps2() */

STATUS
scan( void )
{ 
    if( psdinit() == FAILURE )
    {
        return rspexit();
    }

    ispre = 0;
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
void
dummylinks( void )
{
    epic_loadcvs( "thefile" );            /* for downloading CVs */
}

/************************ END OF GRASS.E ******************************/

