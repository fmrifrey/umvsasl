/* *************************************
 * grad_rf_grass.h
 * This structure is used to track the 
 * rf heating, SAR heating, Grad coil heating,
 * grad amplifier heating.
 *
 * COMMENTS
 * sccs    date	   initials	comments
 *       03/06/99    RAK        Added one more element to the gradient structure.
 *                              This element is defined in epic.h for SGD heating. 
 * ********************************** */

#include "epic_pulse_structs.h"
#include "sar_pm.h"

/* only do this once in any given compilation.*/
#ifndef  grad_rf_grass_INCL
#define  grad_rf_grass_INCL

RF_PULSE rfpulse[MAX_RFPULSE] = {
  /* RFPULSE 0 - RF1 Pulse */
  {  (int *)&pw_rf1,
     (float *)&a_rf1, 
     SAR_ABS_SINC1,
     SAR_PSINC1,
     SAR_ASINC1,
     SAR_DTYCYC_SINC1,
     SAR_MAXPW_SINC1,
     1,
     MAX_B1_SINC1_90,
     MAX_INT_B1_SQ_SINC1_90,
     MAX_RMS_B1_SINC1_90,
     90.0,
     &flip_rf1,
     3200.0,
     1250.0,
     PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,
     0,
     0,
     1.0,
     (int *)&res_rf1,
     0,
     (int *)&wg_rf1,
     1
  },

#include "rf_Prescan.h"
};

GRAD_PULSE gradx[MAX_GRADX] = {
  {
  }
};

GRAD_PULSE grady[MAX_GRADY] = {
  {
  }
};

GRAD_PULSE gradz[MAX_GRADZ] = {
  {
  }
};

#endif  /* grad_rf_grass_INCL */
