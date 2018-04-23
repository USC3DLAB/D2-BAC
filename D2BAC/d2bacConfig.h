/***************************************************************************************
 *  Program dbacConfig.h is the configurable header file for dbacMain.c a program that *
 *  implements the D^2 algorithm for 2-stage SMIP problems.			       *  
 *              								       *
 *         Author:  Lewis Ntaimo 						       *
 *         Date :  March 19, 2003 						       *
 *	   Revised :  								       *
 *										       *
 * 										       *
 ***************************************************************************************/


/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions and malloc */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/**
 * Define the smallest number greater than zero. Any solution value
 * greater than this number is considered nonzero
 */
#define NONZERO_LB 1e-4

/**
 * Define the allowed deviation from integrality for integral solution 
 */
#define INT_PRECISION 1e-6

/**
 * Define the allowed deviation between obj ub and master prob solution 
 */
//#define EPISILON 1e-6
#define EPISILON 3.0

/**
 * Define the allowed % optimality gap [100*(UB-LB)/UB]
 */
#define PERCENT_GAP 0.0100

/**
 * Define the allowed difference % gap between current and previous one 
 * before envoking subproblem MIP solves
 */
#define PERCENT_GAP_THRES 0.0010

/**
 * Define the allowed minimum percent gap 
 * before envoking subproblem MIP solves
 */
#define PERCENT_GAP_LIM 20.0

/**
 * Define the allowed gap between current master obj value and previous obj value
 */
#define MASTER_OBJVAL_GAP 0.01

/**
 * Define the defaulty penalty on the feas obj column for the subproblem 
 */
#define PENALTY 1e6 //CPX_INFBOUND 




