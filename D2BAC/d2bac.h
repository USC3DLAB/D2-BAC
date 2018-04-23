/****************************************************************************
 *  Program d2alg.h is the header file for d2alg.c a program that 			*
 *  implements the d2algorithm for 2-stage SMIP problems. 					*
 * 																			*
 *        Author:    Lewis Ntaimo 											*
 *         Date :    Feb 6, 2003 											*
 *	    Revised :    Feb 17, 2003 											*
 * 																			*
 * 																			*
 ****************************************************************************/


/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */
#include <ilcplex/cplex.h>
/* Bring in the declarations for the string and character functions and malloc */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


#ifndef  CPX_PROTOTYPE_MIN

static void
   free_and_null (char **ptr),
   usage         (char *progname);

#else

static void
   free_and_null (),
   usage         ();

#endif


 
 
 



