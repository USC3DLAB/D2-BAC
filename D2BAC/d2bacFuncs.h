/******************************************************************************
 *  Program d2bacFuncs.h is the header file for d2bacMain.c a program that    *
 *  reads SMPS format file for 2-stage SMIP problems.						  *
 *              						              *
 *         Author:  Yang Yuan 											  *
 *          Date :  Oct 18, 2006 											  *
 *	     Revised :  														  *
 *			Based on the previous work by Lewis Ntaimo						  *
 * 																			  *
 ******************************************************************************/

/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions and malloc */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "d2baclp.h" 	        /* prototypes for lp data structures    */

/***** Function Declarations *******/

static void
freeAndNull (double *ptr);
/**
 * Frees the memory allocated to a double array pointed to by ptr
 *
 */

double myabs(double num);
/**
 * Returns absolute value of num
 */

void printSparseMatrix(int ncols, int nzcnt, int *cmatbeg, int *cmatcnt, int *cmatind, double *cmatval, FILE *fpout);
/**
 * Reads the TIME file data and puts the data in the structure t_info
 * @param ncols an integer indicating the length of the arrays cmatbeg
 * @param nzcnt number of nonzeros in the matrix
 * @param cmatbeg an array containing indices of where each col begins in the array cmatval
 *                and cmatind
 * @param cmatcnt an array containing the number of entries in each column of the matrix
 * @param cmatind an array containing the row indices associated with the elements of cmatval
 * @param cmatval an array containing nonzero coefficients of the specified columns
 */

void printMatrix(int ncols, int nzcnt, int *cmatbeg, int *cmatind, double *cmatval, FILE *fpout);
/**
 * Reads the TIME file data and puts the data in the structure t_info
 * @param ncols an integer indicating the length of the arrays cmatbeg
 * @param nzcnt number of nonzeros in the matrix
 * @param cmatbeg an array containing indices of where each col begins in the array cmatval
 *                and cmatind
 * @param cmatind an array containing the row indices associated with the elements of cmatval
 * @param cmatval an array containing nonzero coefficients of the specified columns
 */

void printWMatrix(int ncols_rows, int nzcnt, int *cmatbeg, int *cmatind,
                       double *cmatval, FILE *fpout);
/**
 * Prints a given column/row-wise sparse matrix
 * @param ncols_rows an integer indicating the length of the arrays cmatbeg
 * @param nzcnt number of nonzeros in the matrix
 * @param cmatbeg an array containing indices of where each col begins in the array cmatval
 *                and cmatind
 * @param cmatcnt an array containing the number of entries in each column of the matrix
 * @param cmatind an array containing the row indices associated with the elements of cmatval
 * @param cmatval an array containing nonzero coefficients of the specified columns
 * @param fpout output file pointer
 */

void
printSoln(int ncols_master, char **colnames_master, solnMatrix_t *solnPtr,
          int nrows, int ncols, char **colnames_sub, FILE *fpSolnOut);
/**
 * This function prints out the optimal solution to a specified file
 * @param solnX array of first-stage solution
 * @param ncols_master number of columns in master problem
 * @param colnames_master pointer to master problem column names
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param nrows   	integer indicating number of rows of the matrix
 *                      This is equal to the number of scenario subproblems
 * @param ncols   	integer indicating number of columns of the matrix
 * @param colnames_sub pointer to sub problem column names
 * @param fpSolnOut output file pointer
 */

double
minimum(double val1, double val2);
/**
 * Returns the minimum of the two values val1 and val2
 * @param val1 first value
 * @param val2 second value
 * @return the minimum of the two input values
 */

double
maximum(double val1, double val2);
/**
 * Returns the maximum of the two values val1 and val2
 * @param val1 first value
 * @param val2 second value
 * @return the maximum of the two input values
 */

double
getSubProbLowerBound(subproblem_t *subprobPtr);
/**
 * Computes and returns the value of the lower bound on the second stage problem.
 * This is required by the D^2 algorithm to keep the subproblem obj positive
 * for all scenario subproblems
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @return the lower bound SUBPROB_LB_L on second stage problem
 */

void
setIncumbent(double *solnX, double objMaster, double expObjlb, int ncolsX, solnMatrix_t *solnPtr);
/**
 * Sets the incumbent solution
 * @param solnX stage one solution array
 * @param objMaster stage one obj
 * @param expObjlb lower bound on expected obj value
 * @param ncolsX number of stage 1 columns (vars)
 * @param solnPtr pointer to solution data structure
 */

int
isFractionalSoln(solnMatrix_t *solnPtr, char* ctype);
/**
 * This function finds the first fractional component in the sub problem
 * solution for all the scenarios and returns a nonzero value.
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */

int
isFractionalVector(double *vec, int len);
/**
 * This function checks if any of the components of the vector is fractional
 * @param vec the vector
 * @param len the length of the vector
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */

int
isEqual(double *vecOne, double *vecTwo, int len);
/**
 * This function checks if the two vectors are the same: Assumes the two vectors
 * are of the same length
 * @param vecOne first the vector
 * @param vecTwo second vector
 * @param len the length of the vectors
 *@return status = 0 if the vectors are the same. Returns a nonzero value otherwise
 */


int
memAllocStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub);
/**
 * Allocates memory to SMPS Format STOCH file data structure variables
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @param nrows_sub  integer indicating number of rows in subproblem lp
 * @param ncols_master  integer indicating number of columns in master lp
 * @param ncols_sub  integer indicating number of columns in subproblem lp
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */

/*	Memory deallocation script. Written on May 18, 2015	*/
int freeStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub);


int
memReAllocStochMatrixT(stochfile_info *stochdataPtr, int random_T);
/**
 * Allocates memory to SMPS Format STOCH file data structure variables
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */

int
memAllocMasterProblemStructs(masterproblem_t *masterprobPtr, int nrows_master,
                             int ncols_master);
/**
 * Allocates memory to subproblem data structure variables
 * @param masterprobPtr  pointer to the master problem data structure
 * @param nrows_master  integer indicating number of stage 1 rows (constraints)
 * @param ncols_master  integer indicating number of stage 1 columns (vars)
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */

/*	Memory deallocation script. Written on May 18, 2015	*/
int freeMasterProblemStructs(masterproblem_t *masterprobPtr, int ncols_master);


int
memAllocSubProblemStruct(subproblem_t *subprobPtr, int nrows_sub, int ncols_sub, int ncols_master,
                         int numscens);
/**
 * Allocates memory to subproblem data structure variables
 * @param subprobPtr  pointer to the subproblem structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param ncols_master   integer indicating number of stage 1 columns (vars)
 * @param numscens   number of scenario subproblems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */

/*	Memory deallocation script. Written on May 18, 2015	*/
int freeSubProblemStruct(subproblem_t *subprobPtr, int numscens);

int
memAllocSolnMatrix(solnMatrix_t *solnPtr, int nrows, int ncols, int ncolsX);
/**
 * Allocates memory to master and all scenario subproblem solution data structure
 * This is a row by row sparse matrix format
 * @param solnMatrix_t  pointer to the the soln data struct with sparse matrix structure
 * @param nrows   	integer indicating number of rows of the matrix
 *                      This is equal to the number of scenario subproblems
 * @param ncolsX   	integer indicating number of columns in the master (x's)
 * @param ncols   	integer indicating number of columns in subproblem
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */


int
memAllocC3parameter(c3lpProb_t *c3lpPtr,stochfile_info *stochdataPtr, subproblem_t *subprobPtr, double *lambda01, double *lambda11,
					double *lambda_02_12, double *pi, double *pi_0, double *grad01, double *grad11, double *grad02_12,
                                        int *ind, double *val, double *secstage, double *secobj);

int
memAllocC3LPproblemStruct(c3lpProb_t *c3LPdualPtr, int nrows_sub, int ncols_sub, int nscenarios);
/**
 * Allocates memory to c3 dual lp data structure variables
 * @param c3LPdualPtr  pointer to the C3 Dual lp structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param nscenarios   integer indicating total number of stage 2 scenario problems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */

int
memAllocRHSLPproblemStruct(rhslpProb_t *rhslpPtr, int nrows_master, int ncols_master);
/**
 * Allocates memory to RHS lp data structure variables
 * @param rhslpPtr  pointer to the C3 Dual lp structure
 * @param nrows_master 	integer indicating number of stage 1 rows (constraints)
 * @param ncols_master 	integer indicating number of stage 1 columns (vars)
 * @param nscenarios   integer indicating total number of stage 2 scenario problems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */

int
loadTimeFile(char *rowname_start, char *colname_start, char *filename, FILE *fpout);
/**
 * Reads the TIME file data and puts the data in the structure timedataPtr
 * @param rowname_start subproblem first row name
 * @param colname_start subproblem first column name
 * @param filename  TIME file name
 * @param fpout output file pointer
 * @return 0 if TIME file is successfully read, otherwise return a nonzero integer
 */

int
loadStochFile(stochfile_info *stochdataPtr, subproblem_t *subprobPtr, char *filename,
              int *random_T, int *random_obj, int ncols_master, FILE *fpout);
/**
 * Reads the STOCH file data and puts the data in the structure timedataPtr
 * @param timedataPtr pointer to data structure to hold STOCH file data
 * @param subprobPtr pointer to master problem data structure
 * @param filename  STOCH file name
 * @param probname  problem name as read from the TIME file
 * @param random_T a pointer to an integer value: o indicates that the technology
 *                 matrix T(w) = T is constant, 1 indicates that it is random T(W)
 * @param random_obj a pointer to an integer value: o indicates that the scenario subprob
 *                 objective is constant, 1 indicates that it is random
 * @param ncols_master number of columns in the master problem
 * @param fpout output file pointer
 * @return 0 if STOCH file is successfully read, otherwise return a nonzero integer
 */


int
loadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,
                      int ncols_master, int nrows_core, int ncols_core, FILE *fpout);
/**
 * Extracts  master problem data from the lp read from the core file and puts it in a data structure
 * It is assumed that the master problem data structure has already been allocated memory
 * @param lp_core pointer to the lp model read from the core file
 * @param master prob_t  pointer to the master problem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_master   integer indicating number of master problem columns (vars)
 * @param nrows_core   integer indicating number of core file lp problem rows (constraints)
 * @param ncols_core   integer indicating number of core file lp problem columns (vars)
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */

int
addOptColToMasterProbLP(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr);
/**
 * This function adds an optimality (theta) column to master problem CPLEX lp object
 * for adding optimality cuts as required by the D^2 algorithm
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr a pointer to the master problem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */

int
setupSubProbMip(CPXENVptr env, CPXLPptr lp_submip, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub,  FILE *fpout);
/**
 * Sets up the subproblem mip CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP. This LP object will be used for solving subproblem mips for
 * upper bounding. This function also updates the subproblem data structure number
 * of rows in the mip obj. D^2 cuts will sequentially be added to this model in main.
 * @param lp_submip pointer to the subproblem mip model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise returns a nonzero integer
 */

int
loadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr, FILE *fpout);
/**
 * Sets up the subproblem lp CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP. This LP object will be used for solving subproblem lps.
 * The function also extracts subproblem data from the lp and puts it in subproblem data
 * structure to be used in creating D^2 cuts. It is assumed that the subproblem data structure
 * has already been allocated memory
 * @param lp pointer to the lp model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */

int
loadSubProblemDataLL(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr, FILE *fpout);
/**
 * Sets up the subproblem lp CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP. This LP object will be used for solving subproblem lps.
 * The function also extracts subproblem data from the lp and puts it in subproblem data
 * structure to be used in creating D^2 cuts. It is assumed that the subproblem data structure
 * has already been allocated memory
 * @param lp pointer to the lp model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */
void
resetSubprobRhsAndT(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
                    double *solnX, int random_T, FILE *fpout);
/**
 * Resets the subproblem RHS vector r(w) and T(w) or T matrix for the constraints that
 * have been reset to 'G' sense from 'L' by multiplying throught by -1
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param subprobPtr  pointer to the sub problem lp data structure
 * @param solnX array of first-stage solution
 * @param random_T a nonzero value indicates random T(w) matrix, a zero value indicates
 *        that is it constant T(w) = T
 * @param fpout output file pointer
 */

void
createsubprobRHS(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                 double *solnX, int random_T, FILE *fpout);
/**
 * Creates the subproblem RHS vector rho(w) = r(w)-T(w)x for a given scenario w
 * random technology matrix T(w). The resulting RHS rho(w) is stored in
 * subprobPtr->rhsrho[scenario]. If T is random then all the T(w) coefs are stored
 * in stochdataPtr structure. Otherwise, the constant values are stored in the
 * T structure in subprobPtr pointer and the newly appended rows (from the D^2-Cut)
 * are stored in the stochdataPtr structure.
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param subprobPtr  pointer to the sub problem lp data structure
 * @param scenario current scenario
 * @param solnX array of first-stage solution
 * @param random_T a zero value indicates constant T matrix
 *        while a nonzero value indicates random T(w) matrix
 * @param fpout output file pointer
 */

void
createsubprobMipRHS(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                 double *solnX, int random_T, FILE *fpout);
/**
 * Creates the subproblem mip RHS vector rho(w) = r(w)-T(w)x for a given scenario w
 * random technology matrix T(w) for Laporte and Louveaux alg. The resulting RHS rho(w) is stored in
 * subprobPtr->rhsrho[scenario]. If T is random then all the T(w) coefs are stored
 * in stochdataPtr structure. Otherwise, the constant values are stored in the
 * T structure in subprobPtr pointer and the newly appended rows (from the D^2-Cut)
 * are stored in the stochdataPtr structure.
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param subprobPtr  pointer to the sub problem lp data structure
 * @param scenario current scenario
 * @param solnX array of first-stage solution
 * @param random_T a zero value indicates constant T matrix
 *        while a nonzero value indicates random T(w) matrix
 * @param fpout output file pointer
 */

void
storeSubProbSoln(solnMatrix_t *solnPtr, double *solnY, int scenario);
/**
 * Allocates memory to all scenario subproblem solution data structure
 * This is a row by row sparse matrix format
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param solnY  array of scenario subproblem solution
 * @param scenario   scenario index (starting from 0)
 */

void
storeSubProbIncumbSoln(solnMatrix_t *solnPtr, double *solnY, int scenario);
/**
 * Allocates memory to all scenario subproblem solution data structure
 * This is a row by row sparse matrix format
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param solnY  array of scenario subproblem solution
 * @param scenario   scenario index (starting from 0)
 */

int
getDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                  int *disj_ind, char *ctype);
/**
 * This function finds the first fractional component in the sub problem
 * solution for all the scenarios and sets its index to disj_var.
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param disj_var pointer to an integer var for storing the col (var) index
 *                with the first fractional component -- this is set as the
 *                disjunction variable
 *@param disj_scen pointer to an integer var for storing the scenario with this
 *					fractional solution component
 *@param disj_ind pointer to an integer for storing the disj_var index in the
 *                soln matrix
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */


int
getNextDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                      int *disj_ind, char *ctype);
/**
 * This function finds the next first fractional component in the sub problem
 * solution for all the scenarios and sets its index to disj_var.
 * The previous disjunction var and scenario are disj_var and disj_scen
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param disj_var pointer to an integer var for storing the col (var) index
 *                with the next first fractional component -- this is set as the
 *                disjunction variable
 *@param disj_scen pointer to an integer var for storing the scenario with this
 *					fractional solution component
 *@param disj_ind pointer to an integer for storing the disj_var index in the
 *                soln sparse matrix : i.e cmatval[disj_ind] = this var soln
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */


int
getMaxDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                     int *disj_ind, char *ctype);
/**
 * This function finds the fractional component in the sub problem
 * solution for all the scenarios with the LARGEST fractional component
 * and sets its index to disj_var. Scans all nonzero solutions!!!
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param disj_var pointer to an integer var for storing the col (var) index
 *                with the first fractional component -- this is set as the
 *                disjunction variable
 *@param disj_scen pointer to an integer var for storing the scenario with this
 *					fractional solution component
 *@param disj_ind pointer to an integer for storing the disj_var index in the
 *                soln matrix
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */


void
getAverageSubProblemsSoln(solnMatrix_t *solnPtr, double *scenProb, double *averSoln);
/**
 * This function averages all the scenario solutions based on the probability
 * of outcome for each scenario
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param scenProb pointer to array containing the outcome prob for each scenario
 *@param averSoln pointer to array for storing the averaged solution
 */

void
getCondAverageSubProblemsSoln(solnMatrix_t *solnPtr, double *scenProb, int disj_var_index,
                              int disj_scen, double *averSoln);
/**
 * This function computes an averaged solutions for all scenario subproblems
 * with a fractional "disj_var_index" solution component. This is used as the object for the
 * C^3 LP for
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param scenProb pointer to array containing the outcome prob for each scenario
 *@param disj_var_index disjunction variable index
 *@param disj_scen disjunction scenario index
 *@param averSoln pointer to array for storing the averaged solution
 */


void
getConditionC3objCoefs(solnMatrix_t *solnPtr, double *scenProb, double *condC3Obj,
                       int scen);
/**
 * This function computes an averaged solutions for all scenario subproblems
 * with a fractional "disj_var_index" solution component. This is used as the object for the
 * C^3 LP for
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param scenProb pointer to array containing the outcome prob for each scenario
 *@param disj_var_index disjunction variable index
 *@param averSoln pointer to array for storing the averaged solution
 */

void
loadRHSLPdata(rhslpProb_t *rhslpPtr, masterproblem_t *masterprobPtr, FILE *fpout);
/**
 * Loads data into RHS lp data structure variables without the t_00 and t_01 columns
 * These columns will be added and deleted later for each scenario as the coefs for
 * these cols become available.
 * @param rhslpPt  pointer to the RHS lp structure
 * @param masterprobPtr  pointer to the master lp data structure
 * @param fpout output file pointer
 */

void
loadC3LPdata(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout);
/**
 * Loads data into C^33 lp data structure variables for the current iteration
 * The columns are arranged in the order: pi_i's, pi_0's, lambda_01's, lambda_11's,
 * lambda_02, and lambda_12
 * @param c3lpPtr  pointer to the C^3 lp data structure
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param stochdataPtr  pointer to the subproblem stochastic data structure
 * @param aveSubprobSoln double pointer to averaged solution of all the scenario subproblems
 * @param solnX array of first-stage solution
 * @param disj_var_index disjunction variable index
 * @param disj_scen disjunction scenario
 * @param disj_var_floor floor value of the disjunction var current soln
 * @param disj_var_ceil ceiling value of the disjunction var current soln
 * @param random_T a pointer to an integer value: o indicates that the technology
 *                 matrix T(w) = T is constant, 1 indicates that it is random T(W)
 * @param fpout output file pointer
 */

void
loadC3LPdataNew(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout);
/**
 * Loads data into C^33 lp data structure variables for the current iteration
 * The columns are arranged in the order: pi_i's, pi_0's, lambda_01's, lambda_11's,
 * lambda_02, and lambda_12
 * @param c3lpPtr  pointer to the C^3 lp data structure
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param stochdataPtr  pointer to the subproblem stochastic data structure
 * @param aveSubprobSoln double pointer to averaged solution of all the scenario subproblems
 * @param solnX array of first-stage solution
 * @param disj_var_index disjunction variable index
 * @param disj_scen disjunction scenario
 * @param disj_var_floor floor value of the disjunction var current soln
 * @param disj_var_ceil ceiling value of the disjunction var current soln
 * @param random_T a pointer to an integer value: o indicates that the technology
 *                 matrix T(w) = T is constant, 1 indicates that it is random T(W)
 * @param fpout output file pointer
 */

void
loadC3LP(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout);

void
updateConstrMatrixWk(subproblem_t *subprobPtr, double *solnPi, int disjVar, int kIter);
/**
 * This function appends the pi vector to the W^k matrix in the k-th iteration
 * and updates the disjVarWIndex array by storing the disjvar index for this pi.
 * Wk is the subproblem constraint matrix in scenario k
 * and is in row sparse matrix format
 * @param subprobPtr  pointer to the subproblem structure
 * @param solnPi array of pi values
 * @param disjVar disjunction variable index
 * @param kIter integer indicating k-th iteration
 */

 int
addNewRowToSubProbWmat(CPXENVptr env, CPXLPptr lp, double *solnPi, int len);
/**
 * This function appends the pi vector to the W^k matrix in the k-th iteration
 * Wk is the subproblem constraint matrix in scenario k
 * and is in row sparse matrix format
 * @param env  a pointer to the CPLEX environment
 * @param lp  a pointer to the CPLEX LP problem object
 * @param solnPi array of pi values
 * @param len length of the solnPi array
 * @return returns 0 on success, otherwise returns a nonzero value
 */

double
scalProd(double *array1, double *array2, int len);
/**
 * This function computes the scalar product of two vectors given by the two arrays
 * @param array1 first array
 * @param array2 second array
 * @param len the length of the two arrays
 * @return returns the scala product of the two arrays
 */

double
scalProd2(double *array1, double *array2, int start, int len);

/**
 * This function computes the scalar product of two vectors of different sizes
 * given by the two arrays.
 * @param array1 first array
 * @param array2 second array
 * @param start the begining index for array1
 * @param len the length of array2, which starts at index 0
 * @return returns the scala product of the two arrays
 */

double
getNu01(subproblem_t *subprobPtr, double *lambda_1, double *rhs,
        double lambda_2, double integral_bd, int disj_var);
/**
 * This function computes the nu_0 or nu_1 scalars
 * Must pass the -ve of lambda_02 for computing nu_0
 * @param array1 first array
 * @param array2 second array
 * @param len the length of the two arrays
 * @return returns the scala product of the two arrays
 */

void
getGammaH(double *lambda_h, stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
          double *gamma_h, int scenario, int random_T, int disj_var);
/**
 * This function computes the vector product between a vector and a matrix
 * The matrix is in column sparse matrix format
 * @param lambda_h array representing lambda_01 or lambda_11 from C^3 LP soln
 * @param stochdataPtr pointer to the data structure for the sparse matrix
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param gamma_h array for storing the gamma_0 or gamma_1 values
 * @param scenario  subproblem scenario
 * @param random_T a zero value indicates constant T matrix
 *        while a nonzero value indicates random T(w) matrix
 */

int
addNewColsToRHSlp(CPXENVptr env, CPXLPptr lp_rhs, double *gamma_0, double gamma_1[],
                  double nu_0, double nu_1, int len);
/**
 * This function adds tau_00 and tau_01 cols to the RHS lp for this scenario.
 * The other cols have already been created
 * @param env  a pointer to the CPLEX environment
 * @param lp_rhs  a pointer to the CPLEX LP problem object for the RHS lp
 * @param gamma_0 array of gamma_o values as computed from the C^3 lp
 * @param gamma_0 array of gamma_1 values as computed from the C^3 lp
 * @param nu_0 nu_0 value as computed from the C^3 lp
 * @param nu_1 nu_1 value as computed from the C^3 lp
 * @param len length of the gamma arrays
 * @return returns 0 on success, otherwise returns a nonzero value
 */

void
updateRHSrAndMatT(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
                  double nu, double *gamma, int scenario, double *solnX);
/**
 * Updates the rhs vector r(w) and technology matrix T(w) for a given scenario
 * as well as the rhs rh0(w) = r(w) - T(w)x for the updated scenario lp new row
 * The T(w) matrix is a row-byrow sparse matrix
 * @param stochdataPtr pointer to data structure holding the stochastic data
 * @param subprobPtr pointer to subproblem lp data structure
 * @param delta the "delta" value from the RHS LP solution
 * @param sigma_0 the "sigma_0" value from the RHS LP solution
 * @param sigma_i the "sigma" array from the RHS LP solution
 * @param solnX pointer to master problem solution
 */

int
addNewRowToMaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr,
                  double *solnX, double expObjval, double SUBPROB_LB);
/**
 * This function appends an optimality cut to the master problem
 * The optimality cut is generated using Laporte and Louveaux method
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr  pointer to the master problem data structure
 * @param solnX array of master x solution values
 * @param expObjval_t expected objective value of the scenario subproblems
 * @return returns 0 on success, otherwise returns a nonzero value
 */

int
addBendersCutToMaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr,
                      double SUBPROB_LB);
/**
 * This function appends an optimality cut to the master problem
 * The optimality cut is generated using Laporte and Louveaux method
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr  pointer to the master problem data structure
 * @param cutCoefs array of master x solution values
 * @param expObjval_t expected objective value of the scenario subproblems
 * @return returns 0 on success, otherwise returns a nonzero value
 */

int
dropScenSolnFromC3obj(solnMatrix_t *solnPtr, double *c3objcoefs,
                      double* scenProb, int disj_scen, double * c3secondobj);
/**
 * This function recomputes the condition expectation of scenario
 * subproblem solutions with fractional disjunctive var component.
 * by dropping off one scenario at a time excluding the disjunction
 * scenario
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param c3objcoefs array of C^3 lp obj coefs
 * @param scenProb scenario original probability
 * @param disj_scen integer indicating the index of the disjunction scenario
 * @param scenToDrop scenario whose soln to drop from C^3 LP obj
 */

int
freeC3lpModelAndData(CPXENVptr env, CPXLPptr lp_c3, c3lpProb_t *c3lpPtr);
/**
 * This frees up the C^3 lp CPLEX object and data arrays
 * @param env  a pointer to the CPLEX environment object
 * @param lp_c3  a pointer to the CPLEX LP problem object
 * @param c3lpPtr  a pointer to the C^3 data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */

void
computeBendersCutCoefs(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                       int curr_nrows, double *duals, double *cutCoefs, double *rhsCoef,
                       int random_T, double total_redcosts, FILE *fpout);
/**
 * Creates the subproblem RHS vector rho(w) = r(w)-T(w)x for a given scenario w
 * random technology matrix T(w). The resulting RHS rho(w) is stored in
 * subprobPtr->rhsrho[scenario]. If T is random then all the T(w) coefs are stored
 * in stochdataPtr structure. Otherwise, the constant values are stored in the
 * T structure in subprobPtr pointer and the newly appended rows (from the D^2-Cut)
 * are stored in the stochdataPtr structure.
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param subprobPtr  pointer to the sub problem lp data structure
 * @param scenario current scenario
 * @param solnX array of first-stage solution
 * @param random_T a zero value indicates constant T matrix
 *        while a nonzero value indicates random T(w) matrix
 * @param total_redcosts sum of reduced costs due to binary variables
 * @param fpout output file pointer
 */

 int
addBinaryConstrsSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr);
/**
 * This function appends expliticty binary constraints to the
 * W matrix in the subprobl CPLEX lp object
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param subprobPtr a pointer to the subproblem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */

int
addFeasColToSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr);
/**
 * This function adds a column to subprobl CPLEX lp object to enable
 * complete resourse as required by the D^2 algorithm
 * This column is penalized in the object function
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param subprobPtr a pointer to the subproblem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */

double
getTotalDualsDueToReducedCosts(subproblem_t *subprobPtr, double *redcosts, double *solnY);
/**
 * This function sums the total duals for binary constraints extracted as
 * the negative of the reduced costs
 * @param subprobPtr a pointer to the subproblem data structure
 * @param redcost a pointer to the subproblem reduced costs
 * @param solnY a pointer to the subproblem solution
 * @return returns sum  of the duals (-ve of the reduced costs for the binary vars
 * whose solution value is one (tight)
 */

 int
DEPloadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,
                      int ncols_master, int nrows_core, int ncols_core, FILE *fpout);
/**
 * Extracts  master problem data from the lp read from the core file and puts it in a data structure
 * It is assumed that the master problem data structure has already been allocated memory
 * @param lp_core pointer to the lp model read from the core file
 * @param master prob_t  pointer to the master problem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_master   integer indicating number of master problem columns (vars)
 * @param nrows_core   integer indicating number of core file lp problem rows (constraints)
 * @param ncols_core   integer indicating number of core file lp problem columns (vars)
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */

 int
DEPloadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
                   FILE *fpout);
/**
 * Sets up the subproblem lp CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP and then extracts subproblem data ( q, T, W, r) from the lp and
 * puts it in subproblem data structure to be used in creating the DEP lp. It is assumed
 * that the subproblem data structure has already been allocated memory.
 * @param lp pointer to the lp model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */

int
DEPsetupDEP(CPXENVptr env, CPXLPptr lp_dep, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
            double random_T, int nrows_master, int ncols_master,
            char** depcolnames, FILE *fpout);
/**
 * Sets up the D.E.P. Assumes all necessary lp data subproblem has been extracted
 * @param env the CPLEX environment
 * @param lp_dep pointer to the lp model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param random_T a nonzero value indicates that the technology matrix (T) is random
 * @param nrows_master number of first stage rows
 * @param ncols_master number of first stage columns
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */


int
getnumscenarios(char *filename);
/**
 * Reads the STOCH file data and determines the total number of scenarios
 * @param filename  STOCH file name
 * @return 0 number of scenarios
 */

double
getrandobjsubproblowerbound(stochfile_info *stochdataPtr, subproblem_t *subprobPtr);
/**
 * Computes and returns the value of the lower bound on the second stage problem.
 * This is required by the D^2 algorithm to keep the subproblem obj positive
 * for all scenario subproblems
 * @param stochdataPtr pointer to the subproblem lp stoch data structure
 * @param subprobPtr   pointer to the subproblem lp data structure
 * @return the lower bound SUBPROB_LB_L on second stage problem
 */

int
resetconstraints(CPXENVptr env, CPXLPptr lp);
/**
 * Resets <= to <= constraints in the CORE FILE lp. This is required by the D2 Algorithm.
 * This is done row by row by multiplying both sides of the constraint by -1.
 * The scenario RHS must also be multiplied by -1 using the function resetScenarioRHS
 * @param env CPLEX environment
 * @param lp pointer to the lp model
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */


/************************************************************************
 *			BRANCH-AND-BOUND FUNCS PROTOTYPES  		*
 *									*									*
 * This part contains functions for an implementation of a best-node   	*
 * strategy branch-and-bound algorithm. The size of the tree is limited	*
 * to the user input maximum number of nodes to explore.		*
 *  									*
 *   uses: orderedlist functions above					*
 * 									*
 *  Author: Lewis Ntaimo						*
 *  Date:   December 21, 2003						*
 *          January 7, 2003						*
 ************************************************************************/

 /**
 * This function finds the "largest" fractional solution index in the sub problem
 * solution for a give scenario solution.
 *@param solnY subproblem solution
 *@param numcols length of solnY array
 *@param index pointer to an integer to return the max fractional index
 *@param value pointer to a double to return the max fractional soln component
 *@param ctype the subproblem array containing the variable types
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */
int
BBgetMaxFractionIndex(double *solnY, int numcols, int *index, double *value, char* ctype);


/**
 * This function appends  branching constraints to the subproblem CPLEX lp object
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param numcontrs number of variable constraints to add
 * @param varind array of branching variable indices
 * @param varval array of branching variable coefs
 * @param varbnds array of contraint rhs values (variable bounds)
 * @return returns 0 on success, otherwise returns a nonzero value
 */
int
BBaddbranchconstrs(CPXENVptr env, CPXLPptr lp_sub, int numconstrs, int *varind, double *varval,
                   double *varbnds);

/**
 * This function performs a branch-and-bound (B&B) alg that uses the best  node
 * strategy. The size of the B&B tree is based on the maximum # of nodes to
 * explore passed to this function.
 * @param env CPLEX environment pointer
 * @param lp_sub subproblem lp pointer
 * @param subprobPtr   pointer to the subproblem lp data structure
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @return status zero for success and non-zero otherwise
 */
int
BBdobranchbound(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
                int max_numnodes, int scenario, double *bb_nu, double **bb_gamma, int *num_nodesexplored,
                int *num_termnodes, double *cputime_sublp, double *bestbound, double *solnY_sub,
                int *fract_sub, int random_T);



/**
 * Allocates memory to B&B reverse polar lp data structure variables
 * @param bblpPtr  pointer to data structure
 * @param nrows_A integer indicating number of stage 1 rows (constraints)
 * @param ncols_A integer indicating number of stage 1 columns (vars)
 * @param num_nodes integer indicating number of explored nodes in the TB&B tree
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
int
BBmalloclpstruct(bblpProb_t *bblpPtr, int nrows_A, int ncols_A, int num_nodes);

/**
 * Loads data into TB&B lp data structure variables
 * @param bblpPt  pointer to the TB&B lp structure
 * @param masterprobPtr  pointer to the master lp data structure
 * @param num_nodesexplored number of terminal nodes
 * @param nu pointer to nu coefs for each terminal node
 * @param gamma pointer to gamma arrays for each terminal node
 * @param solnX master problem X solution
 * @param eta_k current (iteration k) value of master problem eta variable
 */
void
BBloadreversepolarLP(bblpProb_t *bblpPtr, masterproblem_t *masterprobPtr, int num_nodesexplored,
                     double *nu, double **gamma, double *solnX, double eta_k);


/**
 * This frees up the BB lp CPLEX object and data arrays
 * @param env  a pointer to the CPLEX environment object
 * @param lp_bb  a pointer to the CPLEX LP problem object
 * @param bblpPtr  a pointer to the BB data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
int
BBfreelpmodel(CPXENVptr env, CPXLPptr lp_bb, bblpProb_t *bblpPtr);


/**
 * Computes the TB&B node LP nu coef and gamma values for creating the TB&B
 * reverse polar lp.
 * @param stochdataPtr pointer to sub problem stochastic data structure
 * @param subprobPtr  pointer to the sub problem lp data structure
 * @param scenario current scenario
 * @param solnX array of first-stage solution
 * @param random_T a zero value indicates constant T matrix
 *        while a nonzero value indicates random T(w) matrix
 * @param total_redcosts sum of reduced costs due to binary variables
 */
void
BBcomputegammanu (stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                  int curr_nrows, double *duals, double *cutCoefs, double *rhsCoef,
                  int random_T, double total_redcosts);

/**
 * This function appends a "D2 optimality cut" to the master problem
 * The optimality cut is generated using Laporte and Louveaux method
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr  pointer to the master problem data structure
 * @param rhs RHS cut coefficient
 * @param coefs array of lhs cut coefs for x variables
 * @return returns 0 on success, otherwise returns a nonzero value
 */
int
BBaddD2optcuttomaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr, double rhs,
                      double *coefs);









