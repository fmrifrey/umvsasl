/*@Start***********************************************************/
/* GEMSBG C source File
* Copyright (C) 1990 The General Electric Company
*
*    File Name:  grassTutorial.e
*    Developer:  General Electric Company
*/

/*@Synopsis 
* PSD Generation Tutorial: A Gradient Recalled Acquisition Steady State
*  This sequence is a basic shell of a PSD. Refer to the 
*  EPIC Software Reference Manual for information on how to use
*  this PSD shell to generate a working PSD.  The solution is
*  called grass.e
*/     
/*@End*********************************************************/

@inline epic.h
@global
#include "em_psd_ermes.in"
#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "filter.h"

@cv
int num_filter_slots;

/************************************************************************/
/*				HOST SECTION 				*/
/************************************************************************/
@host

abstract("Tutorial GRASS sequence");
psdname("grass");

/************************************************************************/
/*       			CVINIT    				*/
/* Invoked once (& only once) when the PSD host process	is started up.	*/
/* Code which is independent of any OPIO button operation is put here.	*/
/************************************************************************/
int cvinit()
{
#include "cvinit.in"
pircbnub=pircb2nub=0;
pite2nub=6;

EpicConf();
inittargets(&loggrd, &phygrd);
  
if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant),
		exist(opplane), exist(opcoax), 0, 0,
		&opnewgeo, cfsrmode)==FAILURE)
  return FAILURE;

TARDIS_FREQ_OFFSET=RCV_FREQ_CERD;

return SUCCESS;
}

/************************************************************************/
/*       			CVEVAL    				*/
/* Called w/ every OPIO button push which has a corresponding CV. 	*/
/* CVEVAL should only contain code which impacts the advisory panel--	*/
/* put other code in cvinit or predownload				*/
/************************************************************************/
int cveval()
{
return SUCCESS;
}

/************************************************************************/
/*       			CVCHECK    				*/
/* Executed on each 'next page' to ensure prescription can proceed 	*/
/* to the next page. 							*/
/************************************************************************/
int cvcheck()
{
return SUCCESS;
}

/************************************************************************/
/*             		    PRE-DOWN LOAD           		        */
/* Executed prior to a download--all operations not needed for the 	*/
/* advisory panel results.  Execute the	pulsegen macro expansions for	*/
/* the predownload section here.  All internal amps, slice ordering,  	*/
/* prescan slice calc., and SAT placement calculations are performed 	*/
/* in this section.  Time anchor settings for pulsegen are done in this */
/* section too.  				 			*/
/************************************************************************/
int predownload()
{
/*********************************************************************/
#include "predownload.in"	/* include 'canned' predownload code */
/*********************************************************************/

num_filter_slots = 0;
psd_board_type = PSDCERD;
/* @inline Prescan.e PSfilter  */

return SUCCESS;
}

/************************************************************************/
/*             		    PULSE GENERATION        		        */
/************************************************************************/
@pg
pulsegen()
{

psd_board_type = PSDCERD;
sspinit(psd_board_type);

SEQLENGTH(seqcore,optr,seqcore); /* set the sequence length to optr */
buildinstr();              /* load the sequencer memory       */
}

/************************************************************************/
/*             		 REAL-TIME SEQ. PROCESS        		        */
/*		(The signals for the actual IPG hardware)		*/
/************************************************************************/
@rsp

CHAR *entry_name_list[] = {"scan", 0};



