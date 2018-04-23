/******************************************************************************
 *  Program d2baclp.h is the header file for d2alg.c a program that           *
 *  implements the D^2 algorithm for 2-stage SMIP problems.                   *
 *              							      *
 *         Author:  Lewis Ntaimo 			  	              *
 *          Date :  March 4, 2003 					      *
 *	     Revised : Jan 7, 2003   						      *
 *									      *
 * 									      *
 ******************************************************************************/


/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions and malloc */
//#include <ctype.h>
//#include <stdlib.h>
//#include <string.h>


#define STAGES 2
#define LENGTH 400
#define FACTOR 100		// Just for convenience in malloc
#define NAMELEN 50
#define FIELDLEN 500
#define NUMSCENS 2010		// Number of scenarios: this is the size initially
                                // allocated for scenario arrays
#define ROWSIZE_A 2000		// Initial size of row (constraint) arrays for A
#define ROWSIZE_W 24000		// Initial size of row (constraint) arrays for W


/**
 * Structure to hold master problem info
 */
typedef struct {
   int nrows;      	// an integer indicating the number of rows in matrix A
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg_A, cmatind_A
   double *obj;    	// objective vector for the subproblem
   char *ctype;    	// Stores ctypes for for subproblem lp for determining integer variables
   char *sense;    	// Temp array for storing the contraint senses read from the core file
   double *rhs;   	// RHS vector for subproblem
   char **colnames;	// Array of pointers to column names in colnamestore array
   char *colnamestore;	// Column name array
   //char **rownames;	// Row or constraint name array

   /* Constraint matrix A storage variables */
   // Row sparse matrix format
   int nzcnt_A;		// number of nonzeros in matrix W
   int *rmatbeg_A;	// an array containing indices of where each col begins in the array cmatval
                        // and cmatind
   int *rmatcnt_A;     // an array containing the number of entries in each column of W
   int *rmatind_A;	// an array containing the row indices associated with the elements of cmatval
   double *rmatval_A;	// an array containing nonzero coefficients of the specified columns
   int rmatspace_A;	// an integer indicating the size of matrix A
   int spacesize_A;	// Current allocated space for A
   
   double rhsCoef;	// Rhs coefs in optimality cut
   double *cutCoefs;	// Array to store optimality cut coefs
   
} masterproblem_t;

/**
 * Structure to hold subproblem info
 */
typedef struct  {
   int nrowsW;      	// an integer indicating the number of rows in matrix W (constant)
   int nrowsWmip;      // an integer indicating the number of rows in matrix W for mip (constant)
   int nrows_mip;      	// an integer indicating the number of rows in matrix W^k for mip
   int nrows_total;     // an integer indicating the number of rows in matrix W^k for C^3 lp
   int nrows;      	// an integer indicating the number of rows in matrix W^k
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg_W, cmatind_W
   double *obj;    	// objective vector for the subproblem
   char *ctype;    	// Stores ctypes for for subproblem lp for determining integer variables
   char *ctype_ctns;	// Array for storing ctypes 'C'=`- continuous for subproblem lp
   char *sense;    	// Temp array for storing the contraint senses read from the core file
   char *sense_geq;    	// Temp array for storing the contraint senses >= as req'd by the D2 algorithm
   double *rhs;   	// RHS vector for subproblem r as read from the core file
   double **rhsRho;   	// RHS vector for subproblem scenario omega rh0(w) = r(w) - T(w)x
   char **colnames;	// Column or variable name array
   char **rownames;	// Row or constraint name array
   char *colnamestore;	// column name array
   char *rownamestore;	// row name array
   double *lb;		// Array of lower bounds on the cols
   double *ub;		// Array of upper bounds on the cols
   int *indices_ctype;	// Array of length ncols containing the indices for the ctypes
   int *indices_row;	// Array of length nrows containing the numerical indices
     			// of the rows corresponding to the constraints for which 
     			// the rhs coefs are to be changed for the subproblem lp
   int *indices_rowmip;	// Array of length nrows containing the numerical indices
     			// of the rows corresponding to the constraints for which 
     			// the rhs coefs are to be changed for the subproblem mip			
    double *duals;	// Array of dual multipliers

   /* Constraint matrix W storage variables : Initialy it is a column-by-column sparse
      matrix amd then it is changed to row-by-row to allow for addition of D^2 cut pi
      coefs (this is more convenient)
   */
   int nzcnt_W;		// number of nonzeros in matrix W
   int *cmatbeg_W;	// an array containing indices of where each col begins in the array cmatval
                        // and cmatind
   int *cmatcnt_W;     // an array containing the number of entries in each column of W
   int *cmatind_W;	// an array containing the row indices associated with the elements of cmatval
   double *cmatval_W;	// an array containing nonzero coefficients of the specified columns
   int cmatspace_W;	// an integer indicating the size of matrix W
   int spacesize_W;	// Current allocated space for W
   
   int *disjVarWindex; // For each row of pi coefs this array contains indices indicating 
                       // indicating the disjunction variable used to generate that row
                       // Need this when forming the C3-LP. For a given disjunction var
    int nd2cuts;        // Number of D^2 cuts added thus far
   
   // Constant technology matrix T(w) = T storage variables 
   int nrows_T;	   	// an integer indicating the number of rows in the matrix T 
   int nzcnt_T;	// Number of nonzeros in matrix T for scenario w
   int *cmatbeg_T;	// An array containing indices of where each col begins in the array cmatval
                       	// and cmatind 
   int *cmatcnt_T;     	// Array containing the number of entries in each column of T
   int *cmatind_T;	// Array containing the row indices associated with the elements of cmatval
   double*cmatval_T;  	// Array containing nonzero coefficients of the specified columns 
   int cmatspace_T;    // An integers indicating the size of matrix T for each scenario w
   int spacesize_T;	// Current allocated space for T
} subproblem_t;


/**
 * Structure to hold info from the STOCH file: stochastic data
 * r(w) vector, T(w) matrix and prob (support) vector for all the scenarios
 */
typedef struct  {
    char *probname;       	// The name of the problem
    char *content;        	// File content: PERIODS
    char *distn;		// Support distrn type: DSCRETE
    char **scenName;		// Scenario name array
    double *coef;		// coef for variable 'col in contraint row
    double **rhs;		// rhs elements for each scenario
                                // corresponding to the exact location in the subprobPtr->rhs
    double *scenProb;		// Array of probabilities for each scenario
    double *scenCondProb;	// Conditional probability for each scenario
    int nscens;			// total number of scenarios read from stoch file
    int nrows;           	// Number of rows in the T(w) matrix
    int ncols;	   	        // an integer indicating the length of the arrays
                                // cmatbeg_T(w), cmatind_T(w)
                                
    // Random Objective: Added Dec 26, 2003
    int obj_cnt;	// Number of scenario subproblem random obj coefs	
    double **obj;	// Random obj elements for each scenario subproblem
    int *obj_index;	// Array of indices of random obj elements for scenario subproblem                            
                                
                                
                                
    // Technology matrix T(w) storage variables read from the stoch file 
    // This is a col by col sparse matrix 
    int *cnzcnt_T;	 	// vector of number of nonzeros in matrix T for scenario w
    int **cmatbeg_T;		// an array containing indices of where each col begins in the array cmatval
                       		// and cmatind for each scenario
    int **cmatcnt_T;     	// an array containing the number of entries in each column of T
    int **cmatind_T;		// an array containing the row indices associated with the elements of cmatval
   				// for each scenario w
    double **cmatval_T;  	// an array containing nonzero coefficients of the specified columns for each
	                        // scenario w
    int *cmatspace_T;        	// an array of integers indicating the size of matrix T
                                // for each scenario w
    int nscenspace_T;		// Current allocated space for this T in terms of scenarios                         
                                
    // Technology matrix T(w) for newly added rows 
    // This is a row by row sparse matrix 
    int rnrows;           	// Number of rows in the rowT(w) matrix
                                // cmatbeg_T(w), cmatind_T(w)
    int *rmatspace_T;
    int *rnzcnt_T;	 	// vector of number of nonzeros in matrix T for scenario w
    int **rmatbeg_T;		// an array containing indices of where each col begins in the array cmatval
                       		// and cmatind for each scenario
    int **rmatcnt_T;     	// an array containing the number of entries in each row of T
    int **rmatind_T;		// an array containing the row indices associated with the elements of cmatval
   				// for each scenario w
    double **rmatval_T;  	// an array containing nonzero coefficients of the specified columns for each
	                        // scenario w  	                                             
} stochfile_info;


/**
 * Structure to hold solution for all scenario subproblems 
 * in row sparse matrix format
 */
typedef struct  {
   int nrows;      	// an integer indicating the number of rows in matrix S
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg_S, cmatind_S
   double obj;		// incumbent optimal objective value
   double *solnX;       // Incumbent solution
   /* Solution matrix S storage variables */
   int nzcnt_S;		// number of nonzeros in matrix 
   int *rmatbeg_S;	// an array containing indices of where each row begins in the array cmatval
                        // and cmatind
   int *rmatcnt_S;     // an array containing the number of entries in each column of W
   int *rmatind_S;	// an array containing the row indices associated with the elements of cmatval
   double *rmatval_S;	// an array containing nonzero coefficients of the specified columns
   int rmatspace_S;	// an integer indicating the size of matrix 
   
   int ncondScens;      // Number of scenarios used in computing the conditional
    	                // average scenario solutions
   int scenToDrop;     // Next scenario to drop from the C3 lp obj
   int disj_scen_index; // Disjunctive scenario in the condition prob array
   int *condScens;	// Array containing indices of scenarios used in computing
   			// the conditional expectation of the scenario solutions
   			// Scenarion must have a fraction component of the disjunction
   			// variable
   double *condScenProbs;	// Array containing the corresponding cond. scenario probabilities

   /* Incumbent subproblem solution sparse matrix for all scenarios 
      This is a row sparse matrix with each row corresponding to a scenario */
   int nzcnt;		// number of nonzeros in matrix 
   int *rmatbeg;	// an array containing indices of where each row begins in the array cmatval
                        // and cmatind
   int *rmatcnt;        // an array containing the number of entries in each row
   int *rmatind;	// an array containing the row indices associated with the elements of cmatval
   double *rmatval;	// an array containing nonzero coefficients of the specified rows
   int rmatspace;	// an integer indicating the size of matrix 
   
   double best_ub;
   double best_lb;
} solnMatrix_t;

/**
 * Structure to hold the data for the  dual lp to the  C3-LP
 * Convenient to work with the dual
 */
typedef struct  {
   int nrows;      	// an integer indicating the number of rows in constraint matrix
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg,
                        // cmatind
   int nrowsWk;      	// an integer indicating the number of rows in constraint matrix 
                        // W^K for this disjunction variable
                        
   //int nrows_sub;	// Number of rows in second stage subproblem -- for convenience
   //int ncols_sub;	// Number of rows in second stage subproblem -- -- for convenience
   int nscens;	        // Number of scenario subproblems -- for convenience
   int *indices;	// Array to contain column indices for the obj function
   
   double *obj;    	// objective vector for the lp
   char *ctype;	      	// Array for storing  variable  types 
   char *sense;    	// Temp array for storing the contraint senses
   double *rhs;   	// RHS vector for problem
   double *lb;   	// An array containing the lower bound on each of the variables
                        // in the subproblem
   double *ub;   	// An array containing the ub bound on each of the variables
                        // in the subproblem  
   // constraint matrix
   int nzcnt;		// number of nonzeros in matrix
   int *cmatbeg;	// an array containing indices of where each col begins in the array cmatval
                        // and cmatind
   int *cmatcnt;	// an array containing the number of entries in each column of W
   int *cmatind;	// an array containing the row indices associated with the elements of cmatval
   double *cmatval;	// an array containing nonzero coefficients of the specified columns
   int cmatspace;	// an integer indicating the max size of the constraint matrix
} c3lpProb_t;

/**
 * Structure to hold the data for the  RHS LP
 */
typedef struct  {
   int nrows;      	// an integer indicating the number of rows in constraint matrix
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg,
                        // cmatind
   double *obj;    	// objective vector for the lp
   char *ctype;	      	// Array for storing  variable  types 
   char *sense;    	// Temp array for storing the contraint senses
   double *rhs;   	// RHS vector for problem
   double *lb;   	// An array containing the lower bound on each of the variables
                        // in the subproblem
   double *ub;   	// An array containing the ub bound on each of the variables
                        // in the subproblem
   // constraint matrix
   int nzcnt;		// number of nonzeros in matrix
   int *cmatbeg;	// an array containing indices of where each col begins in the array cmatval
                        // and cmatind
   int *cmatcnt;	// an array containing the number of entries in each column of W
   int *cmatind;	// an array containing the row indices associated with the elements of cmatval
   double *cmatval;	// an array containing nonzero coefficients of the specified columns
   int cmatspace;	// an integer indicating the max size of the constraint matrix
} rhslpProb_t;

/**
 * Structure to hold the data for the  B&B reverse polar LP
 */
typedef struct  {
   int nrows;      	// an integer indicating the number of rows in constraint matrix
   int ncols;	   	// an integer indicating the length of the arrays cmatbeg,
                        // cmatind
   double *obj;    	// objective vector for the lp
   char *ctype;	      	// Array for storing  variable  types 
   char *sense;    	// Temp array for storing the contraint senses
   double *rhs;   	// RHS vector for problem
   double *lb;   	// An array containing the lower bound on each of the variables
                        // in the subproblem
   double *ub;   	// An array containing the ub bound on each of the variables
                        // in the subproblem
   // constraint matrix
   int nzcnt;		// number of nonzeros in matrix
   int *cmatbeg;	// an array containing indices of where each col begins in the array cmatval
                        // and cmatind
   int *cmatcnt;	// an array containing the number of entries in each column of W
   int *cmatind;	// an array containing the row indices associated with the elements of cmatval
   double *cmatval;	// an array containing nonzero coefficients of the specified columns
   int cmatspace;	// an integer indicating the max size of the constraint matrix
} bblpProb_t;

