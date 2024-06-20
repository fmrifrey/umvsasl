/*
 *  GE Medical Systems
 *  Copyright (C) 1997 The General Electric Company
 *  
 *  grass.h
 *  
 *  This file contains the prototypes declarations for all functions
 *  in grass.e
 *  
 *  Language : ANSI C
 *  Author   : Gabriel Fernandez
 *  Date     : 21/Mar/1999
 */
/* do not edit anything above this line */

/*
   Version    Author     Date       Comment
----------------------------------------------------------------------
     1.0      GFN     21-Mar-1999   Initial version.
 */


#ifndef umvsasl_h
#define umvsasl_h

/*
 * @host section
 */

/*
 * @pg section
 */

/*
 * @rsp section
 */
STATUS psdinit( void );
STATUS mps2( void );
STATUS aps2( void );
STATUS scan( void );
STATUS scancore( void );
void dummylinks( void );

#endif /* grass_h */

