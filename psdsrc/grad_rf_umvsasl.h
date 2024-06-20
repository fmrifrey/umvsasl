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
#ifndef  grad_rf_umvsasl_INCL
#define  grad_rf_umvsasl_INCL

RF_PULSE rfpulse[MAX_RFPULSE] = {
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
