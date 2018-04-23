 /***********************************************************************************************
 *  Program d2bacFuncs.c contains the functions for d2bacMain.c to solve					*
 *  two stage SMIP problems														.													*
 *																																								*
 *         Author: 			Yang Yuan																										*
 *         Date : 			Sep 10, 2006																									*
 *	   		Finished :  	Oct 30, 2006																									*
 *         Revised :  	Feb 20, 2007																									*
 *																																								*
 *			Based on the preivous work by Lewis Ntaimo															*
 ************************************************************************************************/


/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions and malloc */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "time.h"
#include <math.h>
#include "d2bacFuncs.h"  /* Prototypes for SMPS reader functions */
#include "d2bacConfig.h"  	/* Configurable header file */
#include "orderedlist.h"        /* header file for ordered linked list c file  */

//#define DEBUG_FUNC  // Debugging mode for functions if uncommented


double myabs(double num)
/**
 * Returns absolute value of num
 */
{
  if (num >= 0)
     return num;
  else
     return -num;
}

void printSparseMatrix(int ncols, int nzcnt, int *cmatbeg, int *cmatcnt, int *cmatind,
                       double *cmatval, FILE *fpout)
/**
 * Prints a given column-wise sparse matrix
 * @param ncols an integer indicating the length of the arrays cmatbeg
 * @param nzcnt number of nonzeros in the matrix
 * @param cmatbeg an array containing indices of where each col begins in the array cmatval
 *                and cmatind
 * @param cmatcnt an array containing the number of entries in each column of the matrix
 * @param cmatind an array containing the row indices associated with the elements of cmatval
 * @param cmatval an array containing nonzero coefficients of the specified columns
 * @param fpout output file pointer
 */
{
  int j; // counter

  fprintf(fpout, "\nMatrix nonzeros: %d \n",nzcnt);

  fprintf(fpout, "\nMatrix cmatbeg: \n");
  for (j = 0; j < ncols; j++){
      fprintf(fpout, "%d ", cmatbeg[j]);
  }

  fprintf(fpout, "\nMatrix cmatcnt: \n");
  for (j = 0; j < ncols; j++){
      fprintf(fpout, "%d ", cmatcnt[j]);
  }

  fprintf(fpout, "\nMatrix cmatind: \n");
  for (j = 0; j < nzcnt; j++){
      fprintf(fpout, "%d ", cmatind[j]);
  }

  fprintf(fpout, "\nMatrix cmatval: \n");
  for (j = 0; j < nzcnt; j++){
      fprintf(fpout, "%6.6g ", cmatval[j]);
  }
  fprintf(fpout, "\n\n");

}
/************************** End printSparseMatrix() function *****************************/

void printMatrix(int ncols_rows, int nzcnt, int *cmatbeg, int *cmatind,
                       double *cmatval, FILE *fpout)
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
{
  int i, j; // counter

  fprintf(fpout, "\nMatrix nonzeros: %d\n", nzcnt);

  fprintf(fpout, "ncols_nrows: %d \n\n", ncols_rows);


  for (i = 0; i < ncols_rows-1; i++ ){
     fprintf(fpout, "ROW/COL: %d\n", i);
     for (j = cmatbeg[i]; j < cmatbeg[i+1]; j++){
     	 fprintf(fpout, "Index = %d \t Val = %f \n", cmatind[j], cmatval[j]);
     }
     fprintf(fpout, "\n");
  }
  fprintf(fpout, "ROW/COL: %d\n", ncols_rows-1);;
  for (j = cmatbeg[ncols_rows-1]; j < nzcnt; j++ ){
      fprintf(fpout, "Index = %d \t Val = %f \n", cmatind[j], cmatval[j]);
  }
  fprintf(fpout, "\n");

}
/************************** End printMatrix() function *****************************/


void printWMatrix(int ncols_rows, int nzcnt, int *cmatbeg, int *cmatind,
                       double *cmatval, FILE *fpout)
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
{
  int i, j; // counter
  int count;

  fprintf(fpout, "\nMatrix nonzeros: %d\n", nzcnt);

  fprintf(fpout, "ncols_nrows: %d \n\n", ncols_rows);

  fprintf(fpout, "\nCOLUMNS:\n");
  for (i = 0; i < ncols_rows-1; i++ ){
     for (j = cmatbeg[i]; j < cmatbeg[i+1]; j++){
     	 fprintf(fpout, "%4d ", cmatind[j]);
     }
     fprintf(fpout, "\n");
  }
  for (j = cmatbeg[ncols_rows-1]; j < nzcnt; j++ ){
      fprintf(fpout, "%4d ", cmatind[j]);
  }
  fprintf(fpout, "\n");


  fprintf(fpout, "\nVALUES \t NONZEROS:\n");
  for (i = 0; i < ncols_rows-1; i++ ){
     count = 0;
     for (j = cmatbeg[i]; j < cmatbeg[i+1]; j++){
     	 fprintf(fpout, " %f ", cmatval[j]);
     	 count++;
     }
     fprintf(fpout, "\t%4d \t%4g\n", count, 1.0*nzcnt/count);
  }
  count = 0;
  for (j = cmatbeg[ncols_rows-1]; j < nzcnt; j++ ){
      fprintf(fpout, " %f ", cmatval[j]);
      count++;
  }
  fprintf(fpout, "\t%4d \t%4g\n", count, 1.0*nzcnt/count);

}

void
printSoln(int ncols_master, char **colnames_master, solnMatrix_t *solnPtr,
          int nrows, int ncols, char **colnames_sub, FILE *fpSolnOut)
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
 {
   int i, j;
   int start, stop;
   int col;

   fprintf(fpSolnOut, "\n********* OPTIMAL SOLUTION ********** \n");
   fprintf(fpSolnOut, "\n[lb, ub] = [%.6g, %.6g]\n", solnPtr->best_lb, solnPtr->best_ub);

   fprintf(fpSolnOut, "\nObjective: %6.6g\n", solnPtr->obj );

   fprintf(fpSolnOut, "\nMaster Problem: \n");
   for (i = 0; i < ncols_master; i++) {
      if (solnPtr->solnX[i] > NONZERO_LB)
          fprintf(fpSolnOut," %s = %6.6g\n", colnames_master[i], solnPtr->solnX[i]);
   }
   fprintf(fpSolnOut, "\nScenario Subproblem: \n");
   for (i = 0; i < nrows; i++) {
      fprintf(fpSolnOut, "scenario %d: \n", i);
      if (i == nrows-1) {
          start = solnPtr->rmatbeg[i];
          stop  = solnPtr->nzcnt;
       } else {
          start = solnPtr->rmatbeg[i];
          stop  = solnPtr->rmatbeg[i+1];
       }
      for (j = start; j < stop; j++) {
          col = solnPtr->rmatind[j];
          fprintf(fpSolnOut,"  %s = %6.6g\n", colnames_sub[col], solnPtr->rmatval[j]);
      }
  } // End outer for loop
  fprintf(fpSolnOut, "\n");

  //fprintf(fpSolnOut, "\nAll other variables are equal to zero. \n");

 } // ***************End printSoln ***********************//



double
minimum(double val1, double val2)
/**
 * Returns the minimum of the two values val1 and val2
 * @param val1 first value
 * @param val2 second value
 * @return the minimum of the two input values
 */
{
  if (val1 < val2)
      return val1;
  else
      return val2;
} //*************************** End minimum() ****************************************


double
maximum(double val1, double val2)
/**
 * Returns the maximum of the two values val1 and val2
 * @param val1 first value
 * @param val2 second value
 * @return the maximum of the two input values
 */
{
  if (val1 > val2)
      return val1;
  else
      return val2;
} //*************************** End minimum() ****************************************


double
getSubProbLowerBound(subproblem_t *subprobPtr)
/**
 * Computes and returns the value of the lower bound on the second stage problem.
 * This is required by the D^2 algorithm to keep the subproblem obj positive
 * for all scenario subproblems
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @return the lower bound SUBPROB_LB_L on second stage problem
 */
{
 int i;
 double lb_L = 0;

 for (i = 0; i < subprobPtr->ncols; i++) {
    if (subprobPtr->ctype[i] == 'B' || subprobPtr->ctype[i] == 'I') {
       if (subprobPtr->obj[i] < 0)
   	  lb_L += subprobPtr->obj[i];
    } //else { // if ctype[i] is 'G' need to opt the subproblems!!
      // Will take care of this case later!!
      //}

 } // End for loop

  return -1*lb_L;

} // ******************* getSubProbLowerBound() ************************** //


void
setIncumbent(double *solnX, double objMaster, double expObjlb, int ncolsX, solnMatrix_t *solnPtr)
/**
 * Sets the incumbent solution
 * @param solnX stage one solution array
 * @param objMaster stage one obj
 * @param expObjlb current lower bound on expected obj value
 * @param ncolsX number of stage 1 columns (vars)
 * @param solnPtr pointer to solution data structure
 */
{
   int i, j;

   //printf("setIncumbent():\n");
   //printf("objMaster = %f\n expObjlb = %f\n solnPtr->solnX[ncolsX-1] = %f \n", objMaster, expObjlb, solnPtr->solnX[ncolsX-1]);

   solnPtr->obj = objMaster + expObjlb - solnX[ncolsX-1];
    //printf("solnPtr->obj: %f\n", solnPtr->obj);

   // Copy stage 1 solution
   for (i = 0; i < ncolsX; i++) {
      solnPtr->solnX[i] = solnX[i];
      //printf("\n solnX[%d] = %f\n", i, solnX[i]);
   }

   // Copy stage 2 solution
   solnPtr->nzcnt =  solnPtr->nzcnt_S;

   for (j = 0; j < solnPtr->nrows; j++){
       solnPtr->rmatbeg[j] =  solnPtr->rmatbeg_S[j];
       solnPtr->rmatcnt[j] =  solnPtr->rmatcnt_S[j];
   }

   for (j = 0; j < solnPtr->nzcnt_S; j++){
       solnPtr->rmatind[j] =  solnPtr->rmatind_S[j];
       solnPtr->rmatval[j] =  solnPtr->rmatval_S[j];
   }


} // ******************* setIncumbent() ************************** //


int
isFractionalSoln(solnMatrix_t *solnPtr, char* ctype)
/**
 * This function finds the first fractional component in the sub problem
 * solution for all the scenarios and returns a nonzero value.
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */
{
  int i, j;
  double lval;
  double uval;

  for (i = 0; i < solnPtr->nrows; i++){ // do row by row (scenario) scanning
        //printf("scenario %d\n", i);
        //printf("beg[i] = %d\n", solnPtr->rmatbeg_S[i]);
        //printf("beg[i+1] = %d\n", solnPtr->rmatbeg_S[i+1]);
	for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
	{
	  if (ctype[j] != 'C') {
	     //printf("val[%d] = %f", j, solnPtr->rmatval_S[j]);
   	     //lval = solnPtr->rmatval_S[j] - floor(solnPtr->rmatval_S[j]);
	     //uval = ceil(solnPtr->rmatval_S[j]) - solnPtr->rmatval_S[j];
	     //printf("lval = %f	uval = %f\n", lval, uval);
	     //if (lval > INT_PRECISION || uval > INT_PRECISION){ // Fractional component
	     if (solnPtr->rmatval_S[j] > INT_PRECISION &&
	          solnPtr->rmatval_S[j] < ceil(solnPtr->rmatval_S[j])-INT_PRECISION) { // Fractional component
	         //printf("disj_scen  = %d\n", i);
	         //printf("disj_var  = %d\n", solnPtr->rmatind_S[j]);
	        //*disj_var = solnPtr->rmatind_S[j];     // column index
	        //*disj_scen = i; // row index
	        //*disj_ind = j;

	        return 1;
	      }
	   } // End if
	} // End inner for loop
   } // end outer for loop

   return 0;

} //****************** isFractionalSoln ***************************


int
isFractionalVector(double *vec, int len)
/**
 * This function checks if any of the components of the vector is fractional
 * @param vec the vector
 * @param len the length of the vector
 *@return status = 0 if no fractional component is present: integer soln
 *                   found. Otherwise, a nonzero value is returned
 */
{
  int i;


  for (i = 0; i < len; i++){
      if (vec[i] > INT_PRECISION && vec[i] < ceil(vec[i])-INT_PRECISION)
	    return 1;
   }

   return 0;

} //****************** isFractionalSoln ***************************

int
isEqual(double *vecOne, double *vecTwo, int len)
/**
 * This function checks if the two vectors are the same: Assumes the two vectors
 * are of the same length and have non-negative values
 * @param vecOne first the vector
 * @param vecTwo second vector
 * @param len the length of the vectors
 *@return status = 0 if the vectors are the same. Returns a nonzero value otherwise
 */
{
  int i;
  double temp;


  for (i = 0; i < len; i++){
      temp = vecOne[i] - vecTwo[i];
      if (temp < 0)
         temp *= -1;
      //if (vecOne[i] != vecTwo[i])
      if (temp >= NONZERO_LB)
	    return 1;
  }

  return 0;

} //****************** isEqual ***************************



/*	Memory deallocation scripts. Written on May 18, 2015	*/

int freeStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub)
{
	// already null
	if (stochdataPtr == NULL) 	return 0;

	free(stochdataPtr->probname	);
	free(stochdataPtr->content	);
	free(stochdataPtr->distn	);
	free(stochdataPtr->scenProb	);
	free(stochdataPtr->scenCondProb	);
	free(stochdataPtr->	obj_index	);
	
	int numscens = stochdataPtr->nscens;
	
	for (int j = 0; j < numscens; j++)
	{
		free(stochdataPtr->scenName[j]);
		free(stochdataPtr->rhs[j]);
		free(stochdataPtr->obj[j]);
	}
	free(stochdataPtr->scenName	);
	free(stochdataPtr->rhs		);
	free(stochdataPtr->obj	);
	
	for (int j = 0; j < numscens; j++){
		free(stochdataPtr->cmatbeg_T[j]);
		free(stochdataPtr->cmatcnt_T[j]);
		free(stochdataPtr->cmatind_T[j]);
		free(stochdataPtr->cmatval_T[j]);
	}
	free(stochdataPtr->cnzcnt_T);
	free(stochdataPtr->cmatbeg_T);
	free(stochdataPtr->cmatcnt_T);
	free(stochdataPtr->cmatind_T);
	free(stochdataPtr->cmatval_T);

	for (int j = 0; j < numscens; j++){
		free(stochdataPtr->rmatbeg_T[j]);
		free(stochdataPtr->rmatind_T[j]);
		free(stochdataPtr->rmatval_T[j]);
	}
	free(stochdataPtr->rnzcnt_T);
	free(stochdataPtr->rmatbeg_T);
	free(stochdataPtr->rmatind_T);
	free(stochdataPtr->rmatval_T);
	
	free(stochdataPtr->coef);

	//free(stochdataPtr);		// SMH: I believe no freeing necessary, as this pointer points to a global variable, which will be killed after the program ends.
	stochdataPtr = NULL;

	return 0;
}

int
memAllocStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub)
/**
 * Allocates memory to SMPS Format STOCH file data structure variables
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @param nrows_sub  integer indicating number of rows in subproblem lp
 * @param ncols_master  integer indicating number of columns in master lp
 * @param ncols_sub  integer indicating number of columns in subproblem lp
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */
 {
   int j;

   // STOCH file data structure variables
   //stochdataPtr->nrows  = nrows_sub;
   stochdataPtr->ncols  = ncols_master;
   int rmatspace_T;
   int space;
    int numscens = stochdataPtr->nscens;
   //printf("nrows_sub = %d\n", nrows_sub);
   //printf("numscens = %d\n", numscens);
   //printf("ncols_sub = %d\n", ncols_sub);

   stochdataPtr->probname = (char*)malloc(NAMELEN*sizeof(char));
   stochdataPtr->content  = (char*)malloc(FIELDLEN*sizeof(char));
   stochdataPtr->distn = (char*)malloc(FIELDLEN*sizeof(char));
   stochdataPtr->scenName  = (char**)malloc(numscens*sizeof(char*));
   stochdataPtr->rhs      = (double**)malloc(numscens*sizeof(double*));
   stochdataPtr->scenProb = (double*)malloc(numscens*sizeof(double));
   stochdataPtr->scenCondProb = (double*)malloc(numscens*sizeof(double));

   //Added DEC 26, 2003
   stochdataPtr->obj = (double**)malloc(numscens*sizeof(double*));
   stochdataPtr->obj_index = (int*)malloc(ncols_sub*sizeof(int));

   if (stochdataPtr->probname == NULL || stochdataPtr->content   == NULL ||
       stochdataPtr->distn    == NULL || stochdataPtr->scenName  == NULL ||
       stochdataPtr->rhs       == NULL || stochdataPtr->scenProb == NULL ||
       stochdataPtr->scenCondProb == NULL || stochdataPtr->obj == NULL   ||
       stochdataPtr->obj_index == NULL)
   {
   	fprintf (stderr, "memAllocStochFileStruct(...):\n");
   	fprintf(stderr, "Failure to allocate memory to stochfile pointers\n");
     	fprintf(stderr, "Exiting...\n");
     	return(1);
  }

  for (j = 0; j < numscens; j++) {
     	stochdataPtr->scenName[j] = (char*)malloc(FIELDLEN*sizeof(char));
     	stochdataPtr->rhs[j]      = (double*)malloc((ROWSIZE_W)*sizeof(double));
     	stochdataPtr->obj[j]      = (double*)malloc(ncols_sub*sizeof(double));
     	if(stochdataPtr->scenName[j] == NULL || stochdataPtr->rhs[j] == NULL ||
     	   stochdataPtr->obj[j] == NULL) {
     		fprintf (stderr, "memAllocStochFileStructs(...):\n");
   		fprintf(stderr, "\nFailure to allocate memory to stochfile scenName[]  and rhs[] \n");
     		fprintf(stderr, "Exiting...\n");
     		return(1);
     	}
  }// End for loop

     // printf("Number 9b\n");

  // Random Technology Matrix as read from the stoch file if available
  // This is a col by col sparse matrix

  stochdataPtr->rnrows  = 0;	// Initialize to zero

  stochdataPtr->cnzcnt_T = (int*)malloc(numscens*sizeof(int));
  stochdataPtr->cmatbeg_T = (int**)malloc(numscens*sizeof(int*));
  stochdataPtr->cmatcnt_T = (int**)malloc(numscens*sizeof(int*));
  stochdataPtr->cmatind_T = (int**)malloc(numscens*sizeof(int*));
  stochdataPtr->cmatval_T = (double**)malloc(numscens*sizeof(double*));

  if (stochdataPtr->cnzcnt_T == NULL || stochdataPtr->rhs == NULL ||
      stochdataPtr->cmatbeg_T == NULL || stochdataPtr->cmatbeg_T == NULL ||
      stochdataPtr->cmatind_T == NULL || stochdataPtr->cmatval_T == NULL) //||       stochdataPtr->cmatspace_T == NULL)
    {
        fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
        fprintf(stderr, "Failure to allocate memory to subproblem T(w) data pointers\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

    space = ncols_master*ROWSIZE_W;
    for (j = 0; j < numscens; j++){
    	stochdataPtr->cmatbeg_T[j]  = (int*)malloc(stochdataPtr->ncols*sizeof(int));
    	stochdataPtr->cmatcnt_T[j]  = (int*)malloc(stochdataPtr->ncols*sizeof(int));
    	stochdataPtr->cmatind_T[j]  = (int*)malloc(space*sizeof(int));
    	stochdataPtr->cmatval_T[j]  = (double*)malloc(space*sizeof(double));

    	if (stochdataPtr->cmatbeg_T[j] == NULL || stochdataPtr->cmatcnt_T[j] == NULL ||
    	    stochdataPtr->cmatind_T[j] == NULL || stochdataPtr->cmatval_T[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to subproblem stochastic T(w)[] data arrays\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop




  // Random Technology Matrix for update by the D2 algorithm
  // This is a row by row sparse matrix
  stochdataPtr->nscenspace_T = numscens;

  stochdataPtr->rnzcnt_T = (int*)calloc(numscens, sizeof(int));
  stochdataPtr->rmatbeg_T = (int**)malloc(numscens*sizeof(int*));
  stochdataPtr->rmatind_T = (int**)calloc(numscens, sizeof(int*));
  stochdataPtr->rmatval_T = (double**)malloc(numscens*sizeof(double*));

  if (stochdataPtr->rnzcnt_T  == NULL || stochdataPtr->rmatbeg_T == NULL ||
      stochdataPtr->rmatbeg_T == NULL || stochdataPtr->rmatind_T == NULL ||
      stochdataPtr->rmatval_T == NULL ) // stochdataPtr->rmatspace_T == NULL) //|| stochdataPtr->rmatcnt_T == NULL)
    {
        fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
        fprintf(stderr, "Failure to allocate memory to subproblem T(w) row-by-row data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

   // printf("Number 10\n");
    rmatspace_T  = ncols_master*ROWSIZE_W;      // Initial assumed size
    for (j = 0; j < numscens; j++){
    	stochdataPtr->rmatbeg_T[j]  = (int*)malloc(ROWSIZE_W*sizeof(int));
    	stochdataPtr->rmatind_T[j]  = (int*)malloc(rmatspace_T*sizeof(int));
    	stochdataPtr->rmatval_T[j]  = (double*)malloc(rmatspace_T*sizeof(double));

    	if (stochdataPtr->rmatbeg_T[j] == NULL || stochdataPtr->rmatind_T[j] == NULL
    	    || stochdataPtr->rmatval_T[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to subproblem stochastic T(w)[] row-by-row data arrays\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop


  return(0); // successfull return

}//************************** End memAllocStochFileStruct function *****************************


int
memReAllocStochMatrixT(stochfile_info *stochdataPtr, int random_T)
/**
 * Allocates memory to SMPS Format STOCH file data structure variables
 * @param stochdataPtr  pointer to the STOCH file data structure
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */
 {
   int j;
   int oldsize = stochdataPtr->nscenspace_T;
   int numscens = stochdataPtr->nscens;
   int newsize = oldsize+numscens;

   // Reset allocated size
   stochdataPtr->nscenspace_T = newsize;

   // Reallocate more memory for the rhs and prob vectors
   stochdataPtr->rhs      = (double**)realloc(stochdataPtr->rhs, newsize*sizeof(double*));
   stochdataPtr->scenProb = (double*)realloc(stochdataPtr->scenProb, newsize*sizeof(double));

   if (stochdataPtr->rhs  == NULL || stochdataPtr->scenProb  == NULL) {
          fprintf(stderr, "\n memReallocStochMatrixT(): \n");
          fprintf(stderr, "Failure to allocate memory to stoch rhs and prob data pointers\n");
          fprintf(stderr, "Exiting...\n");
          return(1);
   }

   for (j = oldsize; j < newsize; j++) {
   	stochdataPtr->rhs[j] = (double*)realloc(stochdataPtr->rhs[j], newsize*sizeof(double));
   	if (stochdataPtr->rhs[j]  == NULL ) {
          fprintf(stderr, "\n memReallocStochMatrixT(): \n");
          fprintf(stderr, "Failure to allocate memory to stochdataPtr->rhs[%d] pointer\n", j);
          fprintf(stderr, "Exiting...\n");
          return(1);
   }
   }

   // Reallocate memory for the random Technology matrix T(w)
   if (random_T) {
       stochdataPtr->cnzcnt_T = (int*)realloc(stochdataPtr->cnzcnt_T, newsize*sizeof(int));
       stochdataPtr->cmatbeg_T = (int**)realloc(stochdataPtr->cmatbeg_T, newsize*sizeof(int*));
       stochdataPtr->cmatcnt_T = (int**)realloc(stochdataPtr->cmatcnt_T, newsize*sizeof(int*));
       stochdataPtr->cmatind_T = (int**)realloc(stochdataPtr->cmatind_T, newsize*sizeof(int*));
       stochdataPtr->cmatval_T = (double**)realloc(stochdataPtr->cmatval_T, newsize*sizeof(double*));
       stochdataPtr->cmatspace_T = (int*)realloc(stochdataPtr->cmatspace_T, newsize*sizeof(int));

      if (stochdataPtr->cnzcnt_T == NULL  || stochdataPtr->rhs == NULL ||
          stochdataPtr->cmatbeg_T == NULL || stochdataPtr->cmatbeg_T == NULL ||
          stochdataPtr->cmatind_T == NULL || stochdataPtr->cmatval_T == NULL ||
          stochdataPtr->cmatspace_T == NULL)
      {
          fprintf(stderr, "\n memReallocStochMatrixT(): \n");
          fprintf(stderr, "Failure to allocate memory to subproblem T(w) data pointers\n");
          fprintf(stderr, "Exiting...\n");
          return(1);
      }
      // Initialize the counter
      for (j = oldsize; j < newsize; j++)
          stochdataPtr->cnzcnt_T[j] = 0;
   }// End if statement

    return(0); // successfull return

}//************************** End memReallocStochMatrixT()function *****************************



/*	Memory deallocation script. Written on May 18, 2015	*/

int freeMasterProblemStructs(masterproblem_t *masterprobPtr, int ncols_master)
{
	// already null
	if (masterprobPtr == NULL) 	return 0;

	free(masterprobPtr->obj);
	free(masterprobPtr->ctype);
	free(masterprobPtr->sense);
	free(masterprobPtr->rhs);
	free(masterprobPtr->cutCoefs);
	free(masterprobPtr->colnamestore);
	free(masterprobPtr->colnames);
	
	free(masterprobPtr->rmatbeg_A);
	free(masterprobPtr->rmatcnt_A);
	free(masterprobPtr->rmatind_A);
	free(masterprobPtr->rmatval_A);

	return 0;
}

int
memAllocMasterProblemStructs(masterproblem_t *masterprobPtr, int nrows_master, int ncols_master)
/**
 * Allocates memory to subproblem data structure variables
 * @param masterprobPtr  pointer to the master problem data structure
 * @param nrows_master  integer indicating number of stage 1 rows (constraints)
 * @param ncols_master  integer indicating number of stage 1 columns (vars)
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */
 {
    // initializations
    int status;
    int j;

    // Initialize subproblem data structure vars
    // Will add a column for the theta variable, i.e, optimality variable
    masterprobPtr->nrows  = nrows_master;
    masterprobPtr->ncols  = ncols_master;
    masterprobPtr->rhsCoef = 0;
    masterprobPtr->rmatspace_A  = (ncols_master+1)*ROWSIZE_A; // assumed max number of nonzeros in A

    //printf("masterprobPtr->nrows = %d\n masterprobPtr->ncols = %d\n", masterprobPtr->nrows, masterprobPtr->ncols);

    masterprobPtr->obj = (double*)malloc((ncols_master+1)*sizeof(double));
    masterprobPtr->ctype = (char*)malloc((ncols_master+1)*sizeof(char));
    masterprobPtr->sense = (char*)malloc(ROWSIZE_A*sizeof(char));
    masterprobPtr->rhs = (double*)malloc(ROWSIZE_A*sizeof(double));
    masterprobPtr->colnames = (char**)malloc((ncols_master+1)*sizeof(char*));
    masterprobPtr->colnamestore = (char*)malloc((ncols_master+1)*FIELDLEN*sizeof(char));
    masterprobPtr->cutCoefs = (double*)malloc((ncols_master+1)*sizeof(double));

    masterprobPtr->rmatbeg_A = (int*)malloc((ncols_master+1)*sizeof(int));
    masterprobPtr->rmatcnt_A = (int*)malloc((ncols_master+1)*sizeof(int));
    masterprobPtr->rmatind_A = (int*)malloc(masterprobPtr->rmatspace_A*sizeof(int));
    masterprobPtr->rmatval_A = (double*)malloc(masterprobPtr->rmatspace_A*sizeof(double));

    if (masterprobPtr->obj   == NULL || masterprobPtr->ctype == NULL || masterprobPtr->sense == NULL ||
        masterprobPtr->rhs   == NULL || masterprobPtr->rmatbeg_A  == NULL ||
        masterprobPtr->rmatind_A == NULL || masterprobPtr->rmatind_A == NULL ||
        masterprobPtr->rmatval_A == NULL || masterprobPtr->colnames == NULL ||
        masterprobPtr->cutCoefs == NULL || masterprobPtr->colnamestore == NULL)
    {
        fprintf(stderr, "memAllocMasterProblemStructs(...): \n");
        fprintf(stderr, "Failure to allocate memory to masterproblem data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }
    for (j = 0; j < (ncols_master+1); j++){
    	masterprobPtr->colnames[j] = (char*)malloc(NAMELEN*sizeof(char));
    	if (masterprobPtr->colnames[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocMasterProblemStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to master problem colnames[].\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop

    return(0); // successfull return

}//************************** End memAllocMasterProblemStruct function *****************************


int freeSubProblemStruct(subproblem_t *subprobPtr, int numscens)
{
	free(subprobPtr->obj);
	free(subprobPtr->ctype);
	free(subprobPtr->ctype_ctns);
	free(subprobPtr->sense);
	free(subprobPtr->sense_geq);
	free(subprobPtr->rhs);
	
//	for (int j = 0; j < subprobPtr->nrows; j++){
//		free(subprobPtr->rownames[j]);
//	} // end for loop

	free(subprobPtr->rownamestore);
	free(subprobPtr->colnamestore);
	free(subprobPtr->rownames);
	free(subprobPtr->colnames);
	free(subprobPtr->lb);
	free(subprobPtr->ub);
	free(subprobPtr->indices_row);
	free(subprobPtr->indices_rowmip);
	free(subprobPtr->indices_ctype);
	free(subprobPtr->disjVarWindex);
	free(subprobPtr->duals);
	
	free(subprobPtr->cmatbeg_W);
	free(subprobPtr->cmatcnt_W);
	free(subprobPtr->cmatind_W);
	free(subprobPtr->cmatval_W);
	
	
	free(subprobPtr->cmatbeg_T);
	free(subprobPtr->cmatcnt_T);
	free(subprobPtr->cmatind_T);
	free(subprobPtr->cmatval_T);

	
	for (int j=0; j < numscens; j++){
		free(subprobPtr->rhsRho[j]);
	}
	free(subprobPtr->rhsRho);

	subprobPtr = NULL;
	
	return 0;
}

int
memAllocSubProblemStruct(subproblem_t *subprobPtr, int nrows_sub, int ncols_sub, int ncols_master,
                         int numscens)
/**
 * Allocates memory to subproblem data structure variables
 * @param subprobPtr  pointer to the subproblem structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param ncols_master   integer indicating number of stage 1 columns (vars)
 * @param numscens   number of scenario subproblems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
 {
    // initializations
    int j;

    // Initialize the subproblem data structure
    nrows_sub += ncols_sub; // Will add binary constraints for each var

    //subprobPtr->nrowsW   = nrows_sub;
    //subprobPtr->nrows   = nrows_sub;
    //subprobPtr->nrows_T = nrows_sub;
    subprobPtr->ncols   = ncols_sub;
    subprobPtr->cmatspace_W  = (ncols_sub+1)*ROWSIZE_W; // assumed max number of nonzeros in W
    subprobPtr->nd2cuts   = 0;

    subprobPtr->cmatspace_T = (nrows_sub+1)*ncols_master;
    subprobPtr->nzcnt_T = nrows_sub*ncols_master;

    //printf("subprobPtr->ncols = %d\n", subprobPtr->ncols);



    subprobPtr->obj = (double*)malloc((ncols_sub+1)*sizeof(double));
    subprobPtr->ctype = (char*)malloc((ncols_sub+1)*sizeof(char));
    subprobPtr->ctype_ctns = (char*)malloc((ncols_sub+1)*sizeof(char));
    subprobPtr->sense = (char*)malloc(ROWSIZE_W*sizeof(char));
    subprobPtr->sense_geq = (char*)malloc(ROWSIZE_W*sizeof(char));
    subprobPtr->rhs = (double*)malloc(nrows_sub*sizeof(double));
    subprobPtr->rhsRho = (double**)malloc(numscens*sizeof(double*));
    subprobPtr->rownames      = (char**)malloc(nrows_sub*sizeof(char*));
    subprobPtr->colnames      = (char**)malloc((subprobPtr->ncols+1)*sizeof(char*));
    subprobPtr->rownamestore  = (char*)malloc(nrows_sub*FIELDLEN*sizeof(char));
    subprobPtr->colnamestore  = (char*)malloc((subprobPtr->ncols+1)*FIELDLEN*sizeof(char));
    subprobPtr->lb = (double*)malloc((ncols_sub+1)*sizeof(double));
    subprobPtr->ub = (double*)malloc((ncols_sub+1)*sizeof(double));
    subprobPtr->indices_row = (int*)malloc(ROWSIZE_W*sizeof(int));
    subprobPtr->indices_rowmip = (int*)malloc(ROWSIZE_W*sizeof(int));
    subprobPtr->indices_ctype = (int*)malloc((ncols_sub+1)*sizeof(int));
    subprobPtr->disjVarWindex = (int*)malloc(ROWSIZE_W*sizeof(int));
    subprobPtr->duals = (double*)malloc(ROWSIZE_W*sizeof(double));

    subprobPtr->cmatbeg_W = (int*)malloc((ncols_sub+1)*nrows_sub*sizeof(int));
    subprobPtr->cmatcnt_W = (int*)malloc((ncols_sub+1)*nrows_sub*sizeof(int));
    subprobPtr->cmatind_W = (int*)malloc(subprobPtr->cmatspace_W*sizeof(int));
    subprobPtr->cmatval_W = (double*)malloc(subprobPtr->cmatspace_W*sizeof(double));


    subprobPtr->cmatbeg_T = (int*)malloc((ncols_master)*sizeof(int));
    subprobPtr->cmatcnt_T = (int*)malloc((ncols_master)*sizeof(int));
    subprobPtr->cmatind_T = (int*)malloc(subprobPtr->cmatspace_T*sizeof(int));
    subprobPtr->cmatval_T = (double*)malloc(subprobPtr->cmatspace_T*sizeof(double));

    if (subprobPtr->obj   == NULL || subprobPtr->ctype == NULL || subprobPtr->ctype_ctns == NULL ||
        subprobPtr->sense == NULL || subprobPtr->rhsRho   == NULL || subprobPtr->rownames   == NULL ||
        subprobPtr->colnames   == NULL || subprobPtr->cmatbeg_W  == NULL ||
        subprobPtr->cmatind_W == NULL || subprobPtr->cmatind_W == NULL ||
        subprobPtr->cmatval_W == NULL || subprobPtr->sense_geq == NULL ||
        subprobPtr->lb == NULL || subprobPtr->ub == NULL || subprobPtr->indices_row == NULL ||
        subprobPtr->indices_ctype == NULL || subprobPtr->rhs   == NULL ||
        subprobPtr->disjVarWindex == NULL || subprobPtr->indices_rowmip == NULL)
    {
        fprintf(stderr, "Function memAllocSubProblemStruct(...): \n");
        fprintf(stderr, "Failure to allocate memory to subproblem data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }
 //printf("DONE\n");

    for (j = 0; j < nrows_sub; j++){
    	subprobPtr->rownames[j] = (char*)malloc(NAMELEN*sizeof(char));
    	if (subprobPtr->rownames[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to sub problem rownames[].\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop

    for (j = 0; j < subprobPtr->ncols+1; j++){
    	subprobPtr->colnames[j] = (char*)malloc(NAMELEN*sizeof(char));
    	if (subprobPtr->colnames[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to sub problem colnames[].\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop

    for (j = 0; j < numscens; j++){
    	subprobPtr->rhsRho[j] = (double*)malloc(ROWSIZE_W*sizeof(double));
    	if (subprobPtr->rhsRho[j] == NULL)
    	{
        	fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
        	fprintf(stderr, "Failure to allocate memory to sub problem subprobPtr->rhsRho[j].\n");
        	fprintf(stderr, "Exiting...\n");
        	return(1);
    	}
    } // end for loop

    return(0); // successfull return

}//************************** End memAllocSubProblemStruct function *****************************



int
memAllocSolnMatrix(solnMatrix_t *solnPtr, int nrows, int ncols, int ncolsX)
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
 {
    // Initialize the subproblem data structure
    solnPtr->nrows  = nrows;
    solnPtr->ncols  = ncols;
    solnPtr->rmatspace_S  = nrows*ncols;
    solnPtr->nzcnt_S = 0;

    //printf("nrows = %d\t ncols = %d\t solnPtr->rmatspace_S = %d", nrows, ncols, solnPtr->rmatspace_S);


    solnPtr->solnX = (double*)malloc((ncolsX)*sizeof(double));
    solnPtr->rmatbeg_S = (int*)malloc((nrows+1)*sizeof(int));
    solnPtr->rmatcnt_S = (int*)malloc((nrows+1)*sizeof(int));
    solnPtr->rmatind_S = (int*)malloc(solnPtr->rmatspace_S*sizeof(int));
    solnPtr->rmatval_S = (double*)malloc(solnPtr->rmatspace_S*sizeof(double));
    solnPtr->condScens = (int*)malloc(nrows*sizeof(int));
    solnPtr->condScenProbs = (double*)malloc(nrows*sizeof(double));

    if (solnPtr->rmatbeg_S == NULL || solnPtr->rmatind_S == NULL ||
        solnPtr->rmatind_S == NULL || solnPtr->rmatval_S == NULL ||
        solnPtr->solnX     == NULL || solnPtr->condScens == NULL ||
        solnPtr->condScenProbs == NULL)
    {
        fprintf(stderr, "Function memAllocSubProbSolnMatrix: \n");
        fprintf(stderr, "Failure to allocate memory to scenario subproblems soln array\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

    // Storage vars for the incumbent solution
    solnPtr->rmatspace  = nrows*ncols; // assumed max number of nonzeros in S
    solnPtr->nzcnt = 0;
    solnPtr->rmatbeg = (int*)malloc((nrows+1)*sizeof(int));
    solnPtr->rmatcnt = (int*)malloc((nrows+1)*sizeof(int));
    solnPtr->rmatind = (int*)malloc(solnPtr->rmatspace*sizeof(int));
    solnPtr->rmatval = (double*)malloc(solnPtr->rmatspace*sizeof(double));

    if (solnPtr->rmatbeg == NULL || solnPtr->rmatind == NULL ||
        solnPtr->rmatind == NULL || solnPtr->rmatval == NULL )
    {
        fprintf(stderr, "Function memAllocSubProbSolnMatrix: \n");
        fprintf(stderr, "Failure to allocate memory to scenario subproblems incumbent soln array\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

    return(0); // successfull return

} //************************** End memAllocSubProbSolnMatrix function *****************************

int
memAllocC3parameter(c3lpProb_t *c3lpPtr,stochfile_info *stochdataPtr, subproblem_t *subprobPtr, double *lambda01, double *lambda11,
					double *lambda_02_12, double *pi, double *pi_0, double *grad01, double *grad11, double *grad02_12, int *ind,
					double *val, double *secstage, double *secobj)

{
	int i, j, k;

	i = subprobPtr->ncols + 2;
	j = stochdataPtr->nscens;
	k = c3lpPtr->nrowsWk + 2;

    lambda01 = (double *) malloc (k*sizeof(double));
    lambda11 = (double *) malloc (k*sizeof(double));
    lambda_02_12 = (double *) malloc (2*sizeof(double));
    pi = (double *) malloc (i*sizeof(double));
    pi_0 = (double *) malloc (j*sizeof(double));
    grad01   = (double *) malloc (k*sizeof(double));
    grad11   = (double *) malloc (k*sizeof(double));
    grad02_12= (double *) malloc (2*sizeof(double));
    ind = (int *) malloc ((2*k)*sizeof(int));
    val = (double *) malloc ((2*k)*sizeof(double));
    secstage = (double *) malloc (j * sizeof(double));
    secobj        = (double *) malloc (j * sizeof(double));

    if (lambda01      == NULL || lambda11   == NULL || lambda_02_12 == NULL ||
        pi            == NULL || pi_0       == NULL || grad01       == NULL ||
		grad11        == NULL || grad02_12  == NULL || ind          == NULL ||
        val           == NULL || secstage   == NULL || secobj       == NULL)
    {
        fprintf(stderr, "Function memAllocC3parameter(): \n");
        fprintf(stderr, "Failure to allocate memory to c3 parameter data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

	return (0);

}

int
memAllocC3LPproblemStruct(c3lpProb_t *c3lpPtr, int nrows_sub, int ncols_sub, int nscenarios)
/**
 * Allocates memory to c3 dual lp data structure variables
 * @param c3lpPtr  pointer to the C3 Dual lp structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param nscenarios   integer indicating total number of stage 2 scenario problems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
 {
    // initializations
    int status;
    int j;
    int ncols, nrows;

    // Initializations
    c3lpPtr->nzcnt = 0;
    //c3lpPtr->ncols  = ncols_sub + nscenarios + 2*nrows_sub + 2;
    c3lpPtr->ncols  =  ncols_sub + 2*ROWSIZE_W* + 3;
    c3lpPtr->nrows  = 2*ncols_sub + nscenarios;  //need deeper considerations
    //c3lpPtr->nrows_sub  = nrows_sub;
    //c3lpPtr->ncols_sub  = ncols_sub;
    c3lpPtr->nscens  = nscenarios;

    nrows = c3lpPtr->nrows;
    ncols = c3lpPtr->ncols;
    //c3lpPtr->cmatspace  = nrows*ncols;  // assumed max number of nonzeros in c3 matrix
	 //yarra
    c3lpPtr->cmatspace  = nrows*(ncols / 100);

     // printf("c3lpPtr->nrows = %d\n c3lpPtr->ncols = %d\n", nrows, ncols);

    c3lpPtr->obj = (double*)malloc(2*ncols*sizeof(double));
    c3lpPtr->ctype = (char*)malloc(2*ncols*sizeof(char));
    c3lpPtr->sense = (char*)malloc(2*nrows*sizeof(char));
    c3lpPtr->rhs = (double*)malloc(2*nrows*sizeof(double));
    c3lpPtr->lb = (double*)malloc(2*ncols*sizeof(double));
    c3lpPtr->ub = (double*)malloc(2*ncols*sizeof(double));

    c3lpPtr->cmatbeg = (int*)malloc(2*ncols*sizeof(int));
    c3lpPtr->cmatcnt = (int*)malloc(2*ncols*sizeof(int));
    c3lpPtr->cmatind = (int*)malloc(c3lpPtr->cmatspace*sizeof(int));
    c3lpPtr->cmatval = (double*)malloc(c3lpPtr->cmatspace*sizeof(double));
    c3lpPtr->indices = (int*)malloc(ncols*sizeof(int));

    if (c3lpPtr->obj      == NULL)
        printf("c3lpPtr->obj == NULL\n");
    else if(c3lpPtr->lb   == NULL)
        printf("c3lpPtr->ctype == NULL\n");
    else if(c3lpPtr->cmatind   == NULL)
        printf("c3lpPtr->cmatind == NULL\n");
    else if(c3lpPtr->indices   == NULL)
        printf("c3lpPtr->indices == NULL\n");
    else if(c3lpPtr->sense   == NULL)
        printf("c3lpPtr->sense == NULL\n");
    else if(c3lpPtr->cmatval   == NULL)
        printf("c3lpPtr->cmatval == NULL\n");



    if (c3lpPtr->obj      == NULL || c3lpPtr->ctype   == NULL ||
        c3lpPtr->sense    == NULL || c3lpPtr->rhs     == NULL ||
        c3lpPtr->cmatbeg  == NULL || c3lpPtr->cmatind == NULL ||
        c3lpPtr->cmatind  == NULL || c3lpPtr->cmatval == NULL ||
        c3lpPtr->lb       == NULL || c3lpPtr->ub      == NULL)
    {
        fprintf(stderr, "Function memAllocC3LPproblemStruct(): \n");
        fprintf(stderr, "Failure to allocate memory to c3 lp data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }


    // Initialize column indices
    for (j = 0; j < ncols; j++){
        // Skip freasibility column
        //if (j != ncols-1)
            c3lpPtr->indices[j] = j;
    }

    return(0); // successfull return

} //************************** End memAllocC3LPproblemStruct function ********************

void
loadC3LP(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout)

{
   int i, j, k;
   int m2;        // Number of rows for this W^k matrix
   int n2;        // Number of cols for this W^k matrix including the feasibility variable
   int nscens;    // Number of subproblem scenarios

   int nrows; 	  // Number of rows for this C^3 LP = c3lpPtr->nrows;
   int ncols;	  // Number of cols for this C^3 LP = c3lpPtr->ncols;

   int col_start, row_start1, row_start2;
   int col_count;
   int row_start, row_stop;
   int mynd2cuts = 0;		// counter for number of d^2 cuts added for this disj var
   int curr_row;
   int addthis_row;             // A zero value indicates to extract row from W^k matrix
    				// Otherwise, skip the row, it's not for this disj var
   int row_index;
   int cnt;
   int my_index;

   // Compute the number of rows added for the current disjunction variable thus far
//   printf("\nsubprobPtr->nd2cuts = %d\n",subprobPtr->nd2cuts);
   for(i = 0; i < subprobPtr->nd2cuts; i++){
      if (subprobPtr->disjVarWindex[i] == disj_var_index)
         mynd2cuts++;
   }

   m2 = subprobPtr->nrowsW + mynd2cuts;
   //m2 = subprobPtr->nrowsW;
   c3lpPtr->nrowsWk = m2; // update number of rows in W^k matrix for this disj_var

//   printf ("\n c3lpPtr->nrowsWk = %d \n", c3lpPtr->nrowsWk);
   n2 = subprobPtr->ncols;
   nscens =  stochdataPtr->nscens;
   nrows = 2*n2; 	//INITIALIZE THE NO# OF COLUMNS FOR THE ORIGINAL C3SLP
   ncols = n2 + 2*m2 + 3; //LAMBDA 01, 11, 02, 12, PI, THETA,
   c3lpPtr->nrows = nrows;
   c3lpPtr->ncols = ncols;

   #ifdef DEBUG_FUNC
      fprintf(stdout, " **nrowsW = %d**\n", subprobPtr->nrowsW);
      fprintf(stdout, " **mynd2cuts = %d**\n", mynd2cuts);
      fprintf(stdout, " **c3 rows = %d**\n", m2);
      fprintf(stdout, " **c3lpPtr->nrows = %d**\n", c3lpPtr->nrows);
      fprintf(stdout, " **c3lpPtr->ncols = %d**\n", c3lpPtr->ncols);
      fprintf(stdout, " **subprobPtr->nrows = %d**\n", subprobPtr->nrows_total);
      fprintf(stdout, " **subprobPtr->ncols = %d**\n", subprobPtr->ncols);
   #endif

   #ifdef DEBUG_MAIN
       fprintf (stdout, "\nd2algMain():\n");
       fprintf (stdout, "Average subproblem solution is > \n");
       for (i = 0; i < subprobPtr->ncols; i++)
             fprintf(stdout, "aveSubprobSoln[%d] = %f\n", i, aveSubprobSoln[i] );
   #endif


   //******* Set object coefs ************//
   // pi_i's column
   for (i = 0; i < n2-1; i++) {
       c3lpPtr->obj[i] = aveSubprobSoln[i];
       //fprintf(stdout, "c3lpPtr->obj[%d] = %f\n", i, c3lpPtr->obj[i]);
   }

   // pi feasibility column
   c3lpPtr->obj[n2-1] = 0;

   for (i = n2; i < ncols-1; i++)
   {
	   c3lpPtr->obj[i] = 0;
   }

   c3lpPtr->obj[ncols-1] = 1;

   //*********** Set the ctype array continuous vars ********//
   for (i = 0; i < ncols; i++)      // Without t_00 and t_01 cols
       c3lpPtr->ctype[i] = 'C';

   //************* Set the sense array to >= ***************//
   for (i = 0; i < nrows; i++)      // Without t_00 and t_01 cols
       c3lpPtr->sense[i] = 'G';

   //************** Set rhs coefs **************************//
   for (i = 0; i < nrows; i++)
        c3lpPtr->rhs[i] = 0;

   //*********** Set the bounds on each variable ************//
   // pi_i's and pi_0's variables
   for (i = 0; i < n2; i++){
       c3lpPtr->lb[i] = -1;
       c3lpPtr->ub[i] =  1;
   }

   // lambda's columns
   for (i = n2; i < ncols-1; i++) {
       if (i< subprobPtr->ncols) {
          if ( subprobPtr->sense[i] == 'E') // Make multiplier free
              c3lpPtr->lb[i] = -CPX_INFBOUND;
           else
               c3lpPtr->lb[i] = 0;
       } else {
          c3lpPtr->lb[i] = 0;
       }
       c3lpPtr->ub[i] =  CPX_INFBOUND;
   }

   c3lpPtr->lb[ncols-1] = -1;
   c3lpPtr->ub[ncols-1] = 1;

   //************** Set the constraint matrix coefs in sparse format ***********//

   // Set the column count for cols pi_i's
   for (i = 0; i < n2; i++) {
       c3lpPtr->cmatcnt[i] = 2;
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d\n", i, c3lpPtr->cmatcnt[i]);
   }
   // Initialize to zero for the lambda columns.
   for (i = n2; i < ncols; i++) {
       c3lpPtr->cmatcnt[i] = 0;   // Initialize to 0
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d \n", i, c3lpPtr->cmatcnt[i]);
   }
   //fprintf(stdout, "Insert pi's coefs in C^3 LP");

   // pi_i's columns
   c3lpPtr->nzcnt = 0; // Initialize nonzero counter
   row_start1 = 0;
   row_start2 = n2;
   for (i = 0; i < n2; i++) {
   	c3lpPtr->cmatbeg[i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }
   //fprintf(stdout, "Insert pi_0's coefs in C^3 LP");
/*****
   // pi_0's columns
   row_start1 = 2*n2;
   row_start2 = 2*n2+nscens;
   for (i = 0; i < nscens; i++) {
   	c3lpPtr->cmatbeg[n2+i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }
****/
  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);

   //lambda_01's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens) to c < (n2+nscens+m2)
   //rols start at r = 0 to r < n2 and r = 2*n2 to r < (2*n2 + nscens)

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   // Just extract the original W rows and then the rows for 	//
   // this disjunction variable					//

   //fprintf(stdout, "Insert lambda_01 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable
   for (i = 0; i < subprobPtr->nrows_total; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
                 // printf ("Skip i            = %d \n", i);
             } else {// extract this previously added pi row for this disj var
                  addthis_row = 0;
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                 // printf ("Add i            = %d \n", i);
                 // printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }
       //printf ("c3 lp W^k row i            = %d \n", i);
      // printf ("c3 lp W^k pi row_index     = %d \n", row_index);

       if (addthis_row == 0) {
          // fprintf(stdout, "          W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2] = c3lpPtr->nzcnt;
          //printf("\n\n*******Lambda11: c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row in W^k
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
       	  } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
       	  }
          //fprintf(stdout, "W^k row_start = %d \n", row_start);
          //fprintf(stdout, "W^k row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
              //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                 c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];
                 //fprintf(stdout, "Row = %d  val[%d] = %f\n", row_index,
                 //                 c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                 c3lpPtr->cmatind[c3lpPtr->nzcnt] = subprobPtr->cmatind_W[j];
                 c3lpPtr->nzcnt++;      			// count nonzeros
                 c3lpPtr->cmatcnt[row_index+n2]++;  // count nonzeros in this column
             }
           }
/****
           // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
           for (k = 0; k < nscens; k++) {               // k = scenario
               // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
               ///////////// 6-9-03 ////////////////
             if (i < subprobPtr->nrows_T ) {
                  if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                      c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
                      //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][i]);
        	      c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	      c3lpPtr->nzcnt++;			// count nonzeros
        	      c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
                  }
              } else if (i < subprobPtr->nrowsW) {  // Add the rhs binary bounds
                   c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        	   c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	   c3lpPtr->nzcnt++;			         // count nonzeros
        	   c3lpPtr->cmatcnt[row_index+n2+nscens]++;      // count nonzeros in this column
              } else { // Added rows to W^k
                    my_index = i - subprobPtr->nrowsW + subprobPtr->nrows_T;
                   if (subprobPtr->rhsRho[k][my_index] > NONZERO_LB || subprobPtr->rhsRho[k][my_index] < -NONZERO_LB) {
                      c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][my_index];
                      //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][my_index]);
        	      c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	      c3lpPtr->nzcnt++;			// count nonzeros
        	      c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
                  }
              }

           }
 ****/
       } // End if (addthis_row)

   } // End outer for loop

   //lambda_11's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens+m2) to c < ncols
   //rols start at r = n2 to r < 2*n2 and r = 2*n2+S to r < nrows

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   //fprintf(stdout, "Insert lambda_11 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable

   for (i = 0; i < subprobPtr->nrows_total; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
             } else {// extract this previously added pi row
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                  //printf ("Add i            = %d \n", i);
                  //printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }

      // printf ("c3 lp W^k row i   = %d\n", i);
      // printf ("c3 lp W^k pi row_index = %d\n", row_index);

       if (addthis_row == 0) {
       	  //fprintf(stdout, "          W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2+m2] = c3lpPtr->nzcnt;
          //printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
          } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
          }

         // fprintf(stdout, "row_start = %d \n", row_start);
         // fprintf(stdout, "row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
          //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                 c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];
                // fprintf(stdout, "Row = %d  val[%d] = %f\n",
                //         row_index, c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2 + subprobPtr->cmatind_W[j];
                c3lpPtr->nzcnt++;      			// count nonzeros
                c3lpPtr->cmatcnt[row_index+n2+m2]++;          // count nonzeros in this column
             }
          }
   /***
          // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
          for (k = 0; k < nscens; k++) {               // k = scenario
             // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
             /////////////// 6-9-03 ////////////////
             if (i < subprobPtr->nrows_T) {
                if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                     c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
                     //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][i]);
        	     c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	     c3lpPtr->nzcnt++;			// count nonzeros
        	     c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;  // count nonzeros in this column
                }
              } else if (i < subprobPtr->nrowsW){ // Add the rhs binary bounds
                    c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        	    c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	    c3lpPtr->nzcnt++;			                // count nonzeros
        	    c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;         // count nonzeros in this column
              } else {
                   my_index = i - subprobPtr->nrowsW + subprobPtr->nrows_T;
                   if (subprobPtr->rhsRho[k][my_index] > NONZERO_LB || subprobPtr->rhsRho[k][my_index] < -NONZERO_LB) {
                     c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][my_index];
                     //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][my_index]);
        	     c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	     c3lpPtr->nzcnt++;			// count nonzeros
        	     c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;  // count nonzeros in this column
                   }
              }
          }
		  ***/
       } // End if(addthis_row)
   } // End outer for loop

  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
 //  printf("*******ncols-2 = %d\n", ncols-2);

   // lambda_02's column: add the disjunction var floor coefs
   c3lpPtr->cmatbeg[ncols-3] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-3]++; 				// count nonzeros in this column

   //fprintf(stdout, "disj_var_floor = %f \n", i, disj_var_floor);
   // Add the disj_var_floor coef to this column
   /***
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_floor > NONZERO_LB) {
            fprintf(stdout, "******** disj_var_floor = %f \n", i, disj_var_floor);
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*disj_var_floor;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-2]++; 			// count nonzeros in this column
       }
    }
   ***/

   // lambda_12's column: add the disjunction var ceil coefs
   c3lpPtr->cmatbeg[ncols-2] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2+disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-2]++; 				// count nonzeros in this column

  // fprintf(stdout, "disj_var_ceil = %f \n", i, disj_var_ceil);
   // Add the disj_var_ceil coef to this column
   /****
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_ceil > NONZERO_LB) {
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = disj_var_ceil;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + nscens + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-1]++; 			// count nonzeros in this column
       }
    }
    ****/

	c3lpPtr->cmatbeg[ncols-1] = c3lpPtr->nzcnt;

    c3lpPtr->cmatbeg[ncols] = c3lpPtr->nzcnt;

 } //******************************* loadC3LPdataNew ************************//

/***********************************
int
memAllocC3LPproblemStruct(c3lpProb_t *c3lpPtr, int nrows_sub, int ncols_sub, int nscenarios)
 /*
 * Allocates memory to c3 dual lp data structure variables
 * @param c3lpPtr  pointer to the C3 Dual lp structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param nscenarios   integer indicating total number of stage 2 scenario problems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer

 {
    // initializations
    int status;
    int j;
    int ncols, nrows;

    // Initializations
    c3lpPtr->nzcnt = 0;
    //c3lpPtr->ncols  = ncols_sub + nscenarios + 2*nrows_sub + 2;
    c3lpPtr->ncols  =  ncols_sub + nscenarios + 2*ROWSIZE_W* + 2;
    c3lpPtr->nrows  = 20*(ncols_sub + nscenarios);
    //c3lpPtr->nrows_sub  = nrows_sub;
    //c3lpPtr->ncols_sub  = ncols_sub;
    c3lpPtr->nscens  = nscenarios;

    nrows = c3lpPtr->nrows;
    ncols = c3lpPtr->ncols;
    //c3lpPtr->cmatspace  = nrows*ncols;  // assumed max number of nonzeros in c3 matrix
    c3lpPtr->cmatspace  = nrows*(ncols_sub + nscenarios + 2*nrows_sub + 2);

     // printf("c3lpPtr->nrows = %d\n c3lpPtr->ncols = %d\n", nrows, ncols);

    c3lpPtr->obj = (double*)malloc(2*ncols*sizeof(double));
    c3lpPtr->ctype = (char*)malloc(2*ncols*sizeof(char));
    c3lpPtr->sense = (char*)malloc(2*nrows*sizeof(char));
    c3lpPtr->rhs = (double*)malloc(2*nrows*sizeof(double));
    c3lpPtr->lb = (double*)malloc(2*ncols*sizeof(double));
    c3lpPtr->ub = (double*)malloc(2*ncols*sizeof(double));

    c3lpPtr->cmatbeg = (int*)malloc(2*ncols*sizeof(int));
    c3lpPtr->cmatcnt = (int*)malloc(2*ncols*sizeof(int));
    c3lpPtr->cmatind = (int*)malloc(c3lpPtr->cmatspace*sizeof(int));
    c3lpPtr->cmatval = (double*)malloc(c3lpPtr->cmatspace*sizeof(double));
    c3lpPtr->indices = (int*)malloc(2*(ncols+nscenarios)*sizeof(int));

    if (c3lpPtr->obj      == NULL)
        printf("c3lpPtr->obj == NULL\n");
    else if(c3lpPtr->lb   == NULL)
        printf("c3lpPtr->ctype == NULL\n");
    else if(c3lpPtr->cmatind   == NULL)
        printf("c3lpPtr->cmatind == NULL\n");
    else if(c3lpPtr->indices   == NULL)
        printf("c3lpPtr->indices == NULL\n");
    else if(c3lpPtr->sense   == NULL)
        printf("c3lpPtr->sense == NULL\n");
    else if(c3lpPtr->cmatval   == NULL)
        printf("c3lpPtr->cmatval == NULL\n");



    if (c3lpPtr->obj      == NULL || c3lpPtr->ctype   == NULL ||
        c3lpPtr->sense    == NULL || c3lpPtr->rhs     == NULL ||
        c3lpPtr->cmatbeg  == NULL || c3lpPtr->cmatind == NULL ||
        c3lpPtr->cmatind  == NULL || c3lpPtr->cmatval == NULL ||
        c3lpPtr->lb       == NULL || c3lpPtr->ub      == NULL)
    {
        fprintf(stderr, "Function memAllocC3LPproblemStruct(): \n");
        fprintf(stderr, "Failure to allocate memory to c3 lp data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }


    // Initialize column indices
    for (j = 0; j < ncols_sub+nscenarios; j++){
        // Skip freasibility column
        //if (j != ncols-1)
            c3lpPtr->indices[j] = j;
    }

    return(0); // successfull return

} //************************** End memAllocC3LPproblemStruct function ********************
****/

int
memAllocRHSLPproblemStruct(rhslpProb_t *rhslpPtr, int nrows_master, int ncols_master)
/**
 * Allocates memory to RHS lp data structure variables
 * @param rhslpPtr  pointer to the C3 Dual lp structure
 * @param nrows_master 	integer indicating number of stage 1 rows (constraints)
 * @param ncols_master 	integer indicating number of stage 1 columns (vars)
 * @param nscenarios   integer indicating total number of stage 2 scenario problems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
 {
    // initializations
    int status;
    int j;
    int ncols, nrows;

    // Assignments
    rhslpPtr->nrows  = 2*ncols_master + 5;
    rhslpPtr->ncols  = 2 + ncols_master + 2*nrows_master; // Does not include cols t_00 & t_01

    nrows = rhslpPtr->nrows;
    ncols = rhslpPtr->ncols;
    rhslpPtr->cmatspace  = nrows*ncols;  // assumed max number of nonzeros in W


    //printf("rhslpPtr->nrows = %d\n rhslpPtr->ncols = %d\n", nrows, ncols);

    rhslpPtr->obj = (double*)malloc(ncols*sizeof(double));
    rhslpPtr->ctype = (char*)malloc(ncols*sizeof(char));
    rhslpPtr->sense = (char*)malloc(nrows*sizeof(char));
    rhslpPtr->rhs = (double*)malloc(nrows*sizeof(double));
    rhslpPtr->lb = (double*)malloc(ncols*sizeof(double));
    rhslpPtr->ub = (double*)malloc(ncols*sizeof(double));

    rhslpPtr->cmatbeg = (int*)malloc(ncols*sizeof(int));
    rhslpPtr->cmatcnt = (int*)malloc(ncols*sizeof(int));
    rhslpPtr->cmatind = (int*)malloc(rhslpPtr->cmatspace*sizeof(int));
    rhslpPtr->cmatval = (double*)malloc(rhslpPtr->cmatspace*sizeof(double));


    if (rhslpPtr->obj      == NULL || rhslpPtr->ctype   == NULL ||
        rhslpPtr->sense    == NULL || rhslpPtr->rhs     == NULL ||
        rhslpPtr->cmatbeg  == NULL || rhslpPtr->cmatind == NULL ||
        rhslpPtr->cmatind  == NULL || rhslpPtr->cmatval == NULL ||
        rhslpPtr->lb       == NULL || rhslpPtr->ub      == NULL)
    {
        fprintf(stderr, "Function memAllocC3DualLPproblemStruct(...): \n");
        fprintf(stderr, "Failure to allocate memory to subproblem data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

    return(0); // successfull return

} //************************** End memAllocRHSLPproblemStruct function ********************



int loadTimeFile(char *rowname_start, char *colname_start, char *filename, FILE *fpout)
/**
 * Reads the TIME file data and puts the data in the structure timedataPtr
 * @param rowname_start subproblem first row name
 * @param colname_start subproblem first column name
 * @param filename  TIME file name
 * @param fpout output file pointer
 * @return 0 if TIME file is successfully read, otherwise return a nonzero integer
 */
{
   FILE *time;
   char field1[NAMELEN], field2[NAMELEN];
   char buffer[LENGTH];

   int j;

   // Open the TIME file
   time = fopen(filename,"r" );
   if(time == NULL) {
   	fprintf (stderr, "\n loadTimeFile():\n");
     	fprintf(stderr, "Could not open the TIME file %s for reading!\n", filename);
        return(1);
    }

   // Process the TIME file
   if (fgets (buffer, LENGTH, time) != NULL) {
        //fprintf(stdout, "%s \n", buffer);
        sscanf(buffer, "%s %s", field1, field2);
        //fprintf(stdout, "field1: %s \n", field1);
        //fprintf(stdout, "field2: %s \n", field2);
        if ( strcmp(field1, "TIME") != 0) {
           fprintf (stderr, "\n loadTimeFile():\n");
           fprintf(stderr, "The first line of file %s must start with ""TIME"" and not %s \n!", filename, field1);
           fprintf(stderr, "Exiting ...\n");
           return(1);
        }
   }

   if (fgets (buffer, LENGTH, time) != NULL){
        sscanf(buffer, "%s%s", field1, field2);
   }

    if (fgets (buffer, LENGTH, time) == NULL){
         fprintf (stderr, "\n loadTimeFile():\n");
         fprintf(stderr, "The time file %s must have a description line for the first colname and row names", filename);
         fprintf(stderr, "Exiting...\n");
         return(1);
    }
    if (fgets (buffer, LENGTH, time) == NULL){
         fprintf (stderr, "\n loadTimeFile():\n");
         fprintf(stderr, "The time file %s must have a description line for the first subproblem colname and row names", filename);
         fprintf(stderr, "Exiting...\n");
         return(1);
    }
    sscanf(buffer, "%s%s", colname_start, rowname_start);
    //fprintf(fpout, "%s  %s \n", colname_start, rowname_start);;


   if (fgets (buffer, LENGTH, time) != NULL){
       sscanf(buffer, "%s", field1);
       if ( strcmp(field1, "ENDATA") != 0) {
           fprintf (stderr, "\n loadTimeFile():\n");
           fprintf(stderr, "The last line of file %s must start with ""ENDATA"" \n!", filename);
           fprintf(stderr, "Exiting ...\n");
           return(0);
        }
    } else {
        fprintf (stderr, "\n loadTimeFile():\n");
    	fprintf(stderr, "The last line of file %s must start with ""ENDATA"" \n!", filename);
        fprintf(stderr, "Exiting ...\n");
        return(1);
    }

   // close the TIME file
   fclose(time);
   return (0);

} /************************** End loadTimeFile() function *****************************/



int
loadStochFile(stochfile_info *stochdataPtr, subproblem_t *subprobPtr, char *filename,
              int *random_T, int *random_obj, int ncols_master, FILE *fpout)
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

 {

   FILE *stoch;
   char field1[NAMELEN], field2[NAMELEN], field3[NAMELEN], field4[NAMELEN];
   char buffer[LENGTH];
   char myfilename[LENGTH];
   int i, j, col, scenario;  // counters
   double tempDouble;
   int rowcount = 0;
   int status;
   double sum;
   int numscens = stochdataPtr->nscens;
   int scen_index;


   // Open the TIME file
   stoch = fopen(filename,"r" );
   if(stoch == NULL) {
   	fprintf(stderr, "\nloadStochFile(...): \n");
     	fprintf(stderr, " Could not open the TIME file %s for reading!\n", filename);
     	return(1);
    }


   //************************* Read the TIME file ************************************
   // Read first line: e.g. STOCH	example
   if (fgets (buffer, LENGTH, stoch) != NULL) {
        sscanf(buffer, "%s%s", myfilename, stochdataPtr->probname);
        #ifdef DEBUG_FUNC
            fprintf(fpout, "\nloadStochFile(...): \n");
            fprintf(fpout, "\n %s  %s\n", myfilename, stochdataPtr->probname);
       #endif
        if ( strcmp(myfilename, "STOCH") != 0) {
           fprintf(stderr, "\nloadStochFile(...): \n");
           fprintf(stderr, " The first line of the STOCH file %s must start with ""STOCH"" \n!", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
   }


    // Read the second line: e.g SCENARIOS	DISCRETE
   if (fgets (buffer, LENGTH, stoch) != NULL){
        sscanf(buffer, "%s%s", stochdataPtr->content, stochdataPtr->distn);
        #ifdef DEBUG_FUNC
            fprintf(fpout, "\nloadStochFile(...): \n");
            fprintf(fpout, "%s  %s\n", stochdataPtr->content, stochdataPtr->distn);
        #endif
       if ( strcmp(stochdataPtr->content, "SCENARIOS") != 0){
       	   fprintf(stderr, "\nloadStochFile(...): \n");
           fprintf(stderr, " The second line of file %s must start with ""SCENARIOS""!\n", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
        if ( strcmp(stochdataPtr->distn, "DISCRETE") != 0){
       	   fprintf(stderr, "\nloadStochFile(...): \n");
           fprintf(stderr, " The second word in the second line of file %s must be ""DISCRETE""!\n", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
   }


   // Initialize the counters first
   *random_T = 0; 			// First assume constant technology matrix;
   *random_obj = 0; 			// First assume constant scenario subproblem objective;
   scen_index = -1;   	                // Scenario index


   for (j = 0; j < numscens; j++) {
       stochdataPtr->cnzcnt_T[j] = 0;
       for (i = 0; i < ncols_master-1; i++)
            stochdataPtr->cmatcnt_T[j][i] = 0;
   }

   // Read the rest of the file
   while (fgets (buffer, LENGTH, stoch) != NULL) {
      //fprintf(stderr, "BUFFER: %s\n", buffer);
      sscanf(buffer, "%s", field1);

     // printf("Number 10b\n");

      if (( strcmp(field1, "ENDATA") == 0)) { // End of file
          scen_index++; // Total number of scenarios

      	 // Probabilites must add to one
      	 sum = 0;
      	 for (i = 0; i < scen_index; i++)
      	     sum +=  stochdataPtr->scenProb[i];

      	 if (sum < 1-NONZERO_LB || sum > 1+NONZERO_LB){
      	      fprintf (stderr, "loadStochFile: \n");
      	      fprintf(stderr, "Scenario probabities do NOT add to one!\n");
      	      fprintf(stderr, "Please check your STOCH file %s\n", filename);
      	      fprintf(stderr, "Exiting...\n");
      	      //exit(0);
      	      return (1);
      	 }

      	 // Make copy of probabilities for conditional prob
      	for (i = 0; i < scen_index; i++)
      	     stochdataPtr->scenCondProb[i] =  stochdataPtr->scenProb[i];
      	#ifdef DEBUG_FUNC
      	    for (i = 0; i < scen_index; i++)
      	         printf("scenCondProb[%d] = %f \n", i, stochdataPtr->scenCondProb[i]);
      	#endif

         fclose(stoch);   // close the STOCH file
         return(0);      // successfull return
      } else if (strcmp(field1, "SC") == 0){  // Read in new scenario data T(w) amd r(w)
      	 rowcount = 0; 			// reset the row counter
      	 scen_index++; 			// count scenario
      	 stochdataPtr->obj_cnt = 0;	// Re-initialize scenario subproblem random objective
      	 				// coefficients count

      	/* if (scen_index >= stochdataPtr->nscenspace_T){
      	    // Allocate more space for scenarios
      	    fprintf(stderr, "NEED TO ALLOCATE MORE SPACE FOR SCENARIOS!!\n");
      	    fprintf(stderr, "scen_index = %d: scen_index = %d\n",
      	            scen_index, stochdataPtr->nscenspace_T);

      	    status = memReAllocStochMatrixT(stochdataPtr, *random_T);
      	    if (status) {
      	        fprintf(stderr, "loadStochFile(...):\n");
      	        fprintf(stderr, "Failed to reallocate more memory to stoch data structure\n");
      	        fprintf(stderr, "scen_index = %d: scen_index = %d\n",
      	            scen_index, stochdataPtr->nscenspace_T);
      	    }

      	 }
     */


      	 sscanf(buffer, "%s%s%s%lf%s", field1, field2, field3, &tempDouble, field4);
      	 // Store the scenario name
      	 strcpy(stochdataPtr->scenName[scen_index], field2);

      	 // Read in the prob for this scenario
      	 if (strcmp(field3, "'ROOT'") == 0 || strcmp(field3, "ROOT") == 0){
      	 	stochdataPtr->scenProb[scen_index] = tempDouble;
      	 	#ifdef DEBUG_FUNC
      	 		fprintf(fpout, "loadStochFile(...):\n");
      	 		fprintf(fpout, "ECHO: = %s  %s  %s  %f\n", field1, field2, field3,tempDouble);
      	        #endif
      	 } else {
      	    fprintf(stderr, "\nloadStochFile(...): \n");
      	    fprintf(stderr, "The third word of line: %s\n", buffer);
      	    fprintf(stderr, "Must be ""'ROOT'"" and not ""%s""\n", field3);
      	    fprintf(stderr, "Exiting...\n");
      	    return (1);
      	 }
      } else if (strcmp(field1, "RHS") == 0){ // Read the RHS

      	 sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);
      	 // Get the row index for this rhs
		  
		  printf("%s %s\n", field1, field2);
      	 for (i = 0; i < subprobPtr->nrows; i++){
      	    if (strcmp(subprobPtr->rownames[i], field2) == 0){
      	      break;
      	    }
      	 } // End for loop

      	 // Multiply rhs by -1 for constraints originally with a <= sense
      	 if (subprobPtr->sense[i] == 'L')
      	      stochdataPtr->rhs[scen_index][i] = -1*tempDouble;
      	 else
      	      stochdataPtr->rhs[scen_index][i] = tempDouble;

	 #ifdef DEBUG_FUNC
	      printf("buffer: %s\n", buffer);
 	      fprintf(stderr, "loadStochFile(...):\n");
      	      fprintf(stderr, "stochdataPtr->rhs[%d][%d] = %f\n", scen_index,
      	 	             rowcount, stochdataPtr->rhs[scen_index][rowcount]);
      	 #endif
      	 rowcount++; // Increment row counter

      } else if (strcmp(field1, "OBJ") == 0) { // Stochastic Technology matrix T(w)

         //printf("Random objective function for scenario: %d\n", scen_index);
         *random_obj = 1;

         sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);

      	 //fprintf(stdout, "field1 = %s \t field2 = %s\t  obj = %f \n", field1, field2, tempDouble);

      	 // Get the row index for this rhs
      	 //fprintf(stdout, "subprobPtr->ncols = %d \n", subprobPtr->ncols);
      	 for (i = 0; i < subprobPtr->ncols; i++) {
      	    if (strcmp(subprobPtr->colnames[i], field2) == 0){
      	      break;
      	    }
      	 } // End for loop

      	 //fprintf(stdout, "obj_index[%d] = %d \n", i, stochdataPtr->obj_index[i]);

      	 stochdataPtr->obj[scen_index][stochdataPtr->obj_cnt] = tempDouble;

      	 if (scen_index == 0) { // Store random objective indices only once
             stochdataPtr->obj_index[stochdataPtr->obj_cnt] = i;
         }
         stochdataPtr->obj_cnt++;

      } else { // Stochastic Technology matrix T(w)

        printf("Stochastic Technology matrix T(w)\n");
        *random_T = 1; 		      // Technology matrix is random. More work!!;
        scenario = scen_index;     // Current scenario
      	j = stochdataPtr->cnzcnt_T[scenario]; // Current nonzero count
      	//printf("scenario = %d\n", scenario);
      	//printf("nzcnt_T[%d] = %d\n", scenario, j);

      	// Read in col name, row name and T(w) nonzero
      	sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);
      	//printf("colnames[j] = %s\n", field1);
      	//printf("rownames[j] = %s\n", field2);

      	// Store this nonzero T(w) value
      	stochdataPtr->cmatval_T[scenario][j] = tempDouble;

      	// Set start index for this colname read in
      	if (j == 0) { // start of a column
      	    col = 0;
      	    stochdataPtr->cmatbeg_T[scenario][col] = 0;
      	    stochdataPtr->cmatcnt_T[scenario][col]++; // count nonzeros in this column
      	    //printf("cmatbeg_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatbeg_T[scenario][col]);
      	    //printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
      	} else {
      	    // count nonzeros in this column
      	    if (strcmp(subprobPtr->colnames[j-1], subprobPtr->colnames[j]) == 0) {
      	    	stochdataPtr->cmatcnt_T[scenario][col]++;
      	    	//printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
      	     } else { // new column
      	        col++;
      	        stochdataPtr->cmatbeg_T[scenario][col] = j; // set beginning of this new col read
      	        stochdataPtr->cmatcnt_T[scenario][col]++;   // count nonzeros in this col
      	        //printf("cmatbeg_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatbeg_T[scenario][col]);
      	        //printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
      	     }

      	}


      	// Find row index for current rowname read and store it
      	for (i = 0; i < subprobPtr->nrows; i++) {
      	    if (strcmp(subprobPtr->rownames[i], field2) == 0){
      	         stochdataPtr->cmatind_T[scenario][j] = i;
      	         //printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
      	         break;
      	     }
      	}


      	stochdataPtr->cnzcnt_T[scenario]++; // Increment nonzeros counter

      	#ifdef DEBUG_FUNC
      	    printf("cmatval_T[%d][%d] = %f\n", scenario, j, stochdataPtr->cmatval_T[scenario][j]);
      	    printf("cnzcnt_T[%d] = %d\n", scenario, stochdataPtr->cnzcnt_T[scenario]);
      	    fprintf(fpout, "loadStochFile(...):\n");
      	    fprintf(fpout, "ECHO: = %s  %s    %f\n", field1,
      	            field2, stochdataPtr->cmatval_T[scenario][j]);
      	#endif

      }


   } // end while loop

   return 1;

   fprintf(stderr, "STOCH file %s must end with ""ENDATA""!\n", filename);
   fprintf(stderr, "Exiting...\n");
   // close the TIME file
   fclose(stoch);


   return (1);

} //************************** End loadStochFile() function *****************************



int
loadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,
                      int ncols_master, int nrows_core, int ncols_core, FILE *fpout)
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
 {

    int status;
    int surplus; // To check for sufficiency of the aloocated array size
    int i, j;

    #ifdef DEBUG_FUNC
        fprintf (stderr, "loadMasterProblemData(...):\n");
        fprintf(stderr, "masterprobPtr->nrows = %d\n masterprobPtr->ncols = %d\n", masterprobPtr->nrows, masterprobPtr->ncols);
    #endif


    //***** Create master problem lp by deleting stage_2 cols and rows from core file lp  *******
    status = CPXdelcols(env, lp_core, ncols_master, ncols_core-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to delete stage-2 cols from lp_core, error %d.\n", status);
       return status;
    }

    /* Delete stage-2 rows (constraints) */
    status = CPXdelrows(env, lp_core, nrows_master, nrows_core-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to delete stage-2 rows from lp_core rows, error %d.\n", status);
       return status;
    }
    //////////////////// ADDED Jan 7, 2002 //////////////////
    //            Reset the <= contraints to >=            //
    /////////////////////////////////////////////////////////
    status = resetconstraints(env, lp_core);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
       return status;
    }
    ////////////////////////////////////////////////////////


    // Write master to file in lp format
    #ifdef DEBUG_FUNC
        fprintf (stdout, "\nloadMasterProblemData():\n");
   	status = CPXwriteprob(env, lp_core, "master.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write lp_core problem lp to file, error %d.\n", status);
      		return status;
   	}
    #endif


    //************* Get and store the ctype (var type) array for the master problem lp ********
    status = CPXgetctype(env, lp_core, masterprobPtr->ctype, 0, ncols_master-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to get ctype from lp_core for the master lp, error %d.\n", status);
       return (status);
    }
    #ifdef DEBUG_FUNC
        fprintf(fpout,"loadMasterProblemData(...): \n Master problem ctype: \n");
        for (j = 0; j < ncols_master; j++)
            fprintf(fpout, "col %d: %c\n", j, masterprobPtr->ctype[j]);
    #endif

   //************* Extract the obj from the lp_core from the core file ****************
   status = CPXgetobj(env, lp_core, masterprobPtr->obj, 0, ncols_master-1);
   if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to write get obj from lp_core for the master lp, error %d.\n", status);
       return(status);
   }

   //******************* Extract the A matrix from lp_core from the core file ****************
   // in row sparse format: more convenient for the formation of the RHS LP later
   if (nrows_master > 0) {
         status = CPXgetrows(env, lp_core, &masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatind_A,
                      masterprobPtr->rmatval_A, masterprobPtr->rmatspace_A, &surplus, 0, nrows_master-1);
   }
   if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to get rows from lp_core for A matrix, error %d.\n", status);
       fprintf (stderr, "Surplus value is: %d.\n", surplus);
       return(status);
   }

  //******************** Populate the cmatcnt_A matrix: **************************************
  // an array containing the number of entries in each column of A
  //masterprobPtr->rmatcnt_A[0] = masterprobPtr->rmatbeg_A[1];
 // for (i = 1; i < ncols_master-1; i++){
 //    masterprobPtr->rmatcnt_A[i] = masterprobPtr->rmatbeg_A[i+1]-masterprobPtr->rmatbeg_A[i];
 // }
 // masterprobPtr->rmatcnt_A[ncols_master-1] = masterprobPtr->nzcnt_A - masterprobPtr->rmatbeg_A[ncols_master-1];

  //****** Get the sense vector for the range of contraints for the subproblem ******/
  status = CPXgetsense(env, lp_core, masterprobPtr->sense, 0, nrows_master-1);
  if ( status ) {
     fprintf (stderr, "loadMasterProblemData(...):\n");
     fprintf (stderr, "Failure to get the senses from lp_core problem, error %d.\n", status);
     return(status);
  }

  //****** Get the rhs (b) vector for the master problem *******************************
  status = CPXgetrhs(env, lp_core, masterprobPtr->rhs, 0, nrows_master-1);
  if ( status ) {
     fprintf (stderr, "loadMasterProblemData(...):\n");
     fprintf (stderr, "Failure to get the rhs vector (b) from lp_core problem, error %d.\n", status);
     return(status);
  }


  //**************************** Print the A matrix, senses, and rhs ********************
  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadMasterProblemData(...):\n");
     fprintf(stdout, "masterprobPtr->rmatspace_A = %d\n", masterprobPtr->rmatspace_A);
     fprintf(stdout, "masterprobPtr->nzcnt_A = %d\n", masterprobPtr->nzcnt_A);
     fprintf(stdout, "The A matrix is: \n");
     printSparseMatrix(nrows_master, masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatcnt_A, masterprobPtr->rmatind_A,
                       masterprobPtr->rmatval_A, stdout);

     fprintf(stdout, "master problem senses > \n");
     for (j = 0; j < masterprobPtr->nrows; j++)
           fprintf(stdout, "Row %d: %c\n", j, masterprobPtr->sense[j]);

     fprintf(stdout, "master problem rhs (b) > \n");
     for (j = 0; j < masterprobPtr->nrows; j++)
           fprintf(stdout, "Row %d: %f\n", j, masterprobPtr->rhs[j]);
  #endif

  return (0); // successful return

 } //************************** End loadMasterProblemData function *****************************

int
setupSubProbMip(CPXENVptr env, CPXLPptr lp_submip, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub,  FILE *fpout)
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
 {

    int status;
    int surplus; // To check for sufficiency of the aloocated array size
    int i, j;
    char mipname[20];
    int curr_row;

    strcpy(mipname, "subprobmip1.lp");

   //Create subproblem lp by deleting stage_1 cols and rows from core file lp

   //***************** Delete stage-1 rows (constraints) *********************************
   if (nrows_master > 0) {
        status = CPXdelrows(env, lp_submip, 0, nrows_master-1);
       if ( status ) {
           fprintf(stderr,"\nFunction setupSubProbMip(...): \n");
           fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
           return(status);
        }
   }

  //********************** Delete stage-1 columns (vars) **********************************

  status = CPXdelcols(env, lp_submip, 0, ncols_master-1);
  if ( status ) {
       fprintf(fpout, "\nFunction setupSubProbMip(...): \n");
       fprintf (stderr, "Failure to delete lp_core cols, error %d.\n", status);
       return(status);
  }

  subprobPtr->nrows_mip = CPXgetnumrows (env, lp_submip);
  subprobPtr->nrowsWmip = subprobPtr->nrows_mip;
  subprobPtr->nrows = subprobPtr->nrows_mip;

   //////////////////// ADDED Jan 7, 2002 //////////////////
   //            Reset the <= contraints to >=            //
   /////////////////////////////////////////////////////////
   status = resetconstraints(env, lp_submip);
   if ( status ) {
       fprintf (stderr, "setupSubProbMip(...):\n");
       fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
       return status;
   }
   ////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // Add explicit feasibility column to the model to enable complete recourse //
  // as require by the D^2 algorithm                                          //
  //////////////////////////////////////////////////////////////////////////////

   status = addFeasColToSubProbLP(env, lp_submip, subprobPtr);
   if ( status ) {
  	fprintf (stderr, "setupSubProbMip:\n");
      	fprintf (stderr, "Failure to add binary contraints to lp_submip, error %d.\n", status);
      	return(status);
   }

 //********************** Write prob to file in lp format **************************
 #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_submip, mipname, NULL);
  	if ( status ) {
   		fprintf (stderr, "Function setupSubProbMip():\n");
      		fprintf (stderr, "Failure to write lp_submip to file, error %d.\n", status);
      		return(status);
 	}
  #endif

  //****** Initialize the indices array corresponding to the constraints *********//
 //************* for which the rhs coefs are to be changed *********************//
   for (j = 0; j < subprobPtr->nrows_mip; j++) {
     subprobPtr->indices_rowmip[j] = j;
   }
   // Remember to update this array every time a cut is added to the subproblem!!!
  #ifdef DEBUG_FUNC
      fprintf (stdout, "Function setupSubProbMip:\n");
      for(i = 0; i < subprobPtr->nrows_mip; i++){
          fprintf(stdout, "subprobPtr->indices_rowmip[%d] = %d \n", i, subprobPtr->indices_rowmip[i]);
      }
  #endif

   subprobPtr->ncols = CPXgetnumcols (env, lp_submip);

   return (0);
 } // End function setupSubProbMip

int
loadSubProblemDataLL(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
                   FILE *fpout)
/**
 * Sets up the subproblem lp CPLEX LP for L & L alg by deleting the first stage rows and cols
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
 {

    int status;
    int surplus; // To check for sufficiency of the alocated array size
    int i, j;
    char lpname[20];
    int curr_row;
    int numscens = stochdataPtr->nscens;


    strcpy(lpname, "subproblp.lp");;


   //Create subproblem lp by deleting stage_1 cols and rows from core file lp

   //***************** Delete stage-1 rows (constraints) *********************************
   if (nrows_master > 0) {
         status = CPXdelrows(env, lp_sub, 0, nrows_master-1);
         if ( status ) {
              fprintf(stderr,"\nFunction loadSubProblemData(...): \n");
              fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
              return(status);
         }
   }

   // Get the sense vector for the range of contraints for the subproblem       //
  status = CPXgetsense(env, lp_sub, subprobPtr->sense, 0, subprobPtr->nrows-1);
  if ( status ) {
     fprintf (stderr, "loadSubProblemData:\n");
     fprintf (stderr, "Failure to get the senses from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem senses\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->sense[%d:] %c\n", j, subprobPtr->sense[j]);
       }
  #endif

   //////////////////// ADDED Jan 7, 2002 //////////////////
   //            Reset the <= contraints to >=            //
   /////////////////////////////////////////////////////////
   status = resetconstraints(env, lp_sub);
   if ( status ) {
       fprintf (stderr, "loadSubProblemDataLL(...):\n");
       fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
       return status;
   }
   ////////////////////////////////////////////////////////

   //**********Extract the T matrix from the core file and store it in T *****************
   status = CPXgetcols(env, lp_sub, &subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatind_T, subprobPtr->cmatval_T,
                       subprobPtr->cmatspace_T, &surplus, 0, ncols_master-1);
   if ( status ) {
      fprintf (stderr, "Failure to extraxt T from lp_sub from T(w) matrix, error %d.\n", status);
      fprintf (stderr, "Surplus value is: %d.\n", surplus);
      return(status);
   }

  //******************** Print the T matrix ******************************************
  #ifdef DEBUG_FUNC
     fprintf(stderr, "\nFunction loadSubProblemData(...): \n");
     fprintf(stderr, "subprobPtr->nzcnt_T = %d\n", subprobPtr->nzcnt_T);
     fprintf(stderr, "*******The T matrix is********: \n");
     printSparseMatrix(ncols_master, subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatbeg_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, stderr);
  #endif


  //********************** Delete stage-1 columns (vars) **********************************
  status = CPXdelcols(env, lp_sub, 0, ncols_master-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "Failure to delete lp_sub cols, error %d.\n", status);
     return(status);
  }

  //************* Get and store the ctype (var type) array from the subproblem lp ********
  status = CPXgetctype(env, lp_sub, subprobPtr->ctype, 0, ncols_sub-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "**Failure to get ctype from lp_sub, error %d.\n", status);
     return (status);
  }
  #ifdef DEBUG_FUNC
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf(stdout,"Subproblem ctype: \n");
     for (j = 0; j < ncols_sub; j++)
         fprintf(stdout, "col %d: %c\n", j, subprobPtr->ctype[j]);
  #endif

  // Store a vector of continuous ctypes for the linear relaxation
  // and ctype indices
  // Include feasibility column to be added soon
  subprobPtr->ctype[ncols_sub] = 'C';
  for (j = 0; j < ncols_sub+1; j++){
     subprobPtr->ctype_ctns[j] = 'C';
     subprobPtr->indices_ctype[j] = j;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Add explicit feasibility column to the model to enable complete recourse //
  // as require by the D^2 algorithm                                          //
  //////////////////////////////////////////////////////////////////////////////

   status = addFeasColToSubProbLP(env, lp_sub, subprobPtr);
   if ( status ) {
  	fprintf (stderr, "loadSubProblemData:\n");
      	fprintf (stderr, "Failure to add binary contraints to lp_sub, error %d.\n", status);
      	return(status);
   }



   // Make updates
   nrows_sub          = CPXgetnumrows (env, lp_sub);
   subprobPtr->nrows  = nrows_sub;
   subprobPtr->nrowsW = nrows_sub;;
   subprobPtr->nrows_T = nrows_sub;
   ncols_sub          = CPXgetnumcols (env, lp_sub);
   subprobPtr->ncols  = ncols_sub;

   stochdataPtr->nrows  = nrows_sub;	// Initialize: VERY IMPORTANT!!!!!!



  //********************** Write prob to file in lp format **************************
 #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_sub, lpname, NULL);
  	if ( status ) {
   		fprintf (stderr, "loadSubProblemData:\n");
      		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
      		return(status);
 	}
  #endif


  //************* Extract the obj from the lp_sub from the core file ****************
  status = CPXgetobj(env, lp_sub, subprobPtr->obj, 0, ncols_sub-1);
  if ( status ) {
      fprintf (stderr, "loadSubProblemData:\n");
      fprintf (stderr, "Failure to get obj from lp_sub, error %d.\n", status);
      return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem obj\n");
       for (j = 0; j < subprobPtr->ncols; j++) {
           fprintf(stdout, "subprobPtr->obj[%d]: %f\n", j, subprobPtr->obj[j]);
       }
  #endif


  // Get the rhs vector for the range of contraints for the subproblem  //
  // and initialize each scenario rhs with these values                 //
  status = CPXgetrhs(env, lp_sub, subprobPtr->rhs, 0, nrows_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get the rhs from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem rhs\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->rhs[%d]: %f\n", j, subprobPtr->rhs[j]);
       }
  #endif



  // Initialize each scenario rhs with these values: Need to rewrite this section for efficiency !
  // Don't know the number of scenarios yet!!

  // Initialize the scenario subproblem rhs(w)
  for (i = 0; i < numscens; i++){
      for (j = 0; j < subprobPtr->nrows; j++)
          stochdataPtr->rhs[i][j] = subprobPtr->rhs[j];
  }



  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadSubProblemData:\n");
     for (i = 0; i < 3; i++){
        fprintf(stdout, "Subproblem rhs for scenario %d\n", i);
        for (j = 0; j < subprobPtr->nrows; j++) {
            fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", i, j, stochdataPtr->rhs[i][j]);
         }
     }
  #endif

  //Get the lower bounds on the subproblem cols (vars)
  status = CPXgetlb(env, lp_sub, subprobPtr->lb, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get lower bounds lb from lp_sub, error %d.\n", status);
     return (status);
  }
  //Get the upper bounds on the subproblem cols (vars)
  status = CPXgetub(env, lp_sub, subprobPtr->ub, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Function loadSubProblemData \n");
     fprintf (stderr, "Failure to get upper bounds ub from lp_sub, error %d.\n", status);
     return (status);
  }


  #ifdef DEBUG_FUNC
      fprintf (stdout, "loadSubProblemData:\n");
      fprintf(stdout, "\n\n\n  Subproblem lp column ub and lb:\n");
      for(i = 0; i < ncols_sub; i++){
         fprintf(stdout, "lb[%d] = %f	ub[%d] = %f \n", i,subprobPtr->lb[i],i,subprobPtr->ub[i]);
      }
  #endif

 //****** Initialize the indices array corresponding to the constraints *********//
 //************* for which the rhs coefs are to be changed *********************//
   for (j = 0; j < subprobPtr->nrows; j++) {
     subprobPtr->indices_row[j] = j;
   }
   // Remember to update this array every time a cut is added to the subproblem!!!
  #ifdef DEBUG_FUNC
      fprintf (stdout, "loadSubProblemData:\n");
      for(i = 0; i < subprobPtr->nrows; i++){
          fprintf(stdout, "subprobPtr->indices_row[%d] = %d \n", i, subprobPtr->indices_row[i]);
      }
  #endif

  return (0); // Successful return

 } //************************** End loadSubProblemDataLL function *****************************



int
loadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
                   FILE *fpout)
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
 {

    int status;
    int surplus; // To check for sufficiency of the aloocated array size
    int i, j;
    char lpname[20];
    int curr_row;
    char *ub = NULL;
    int *indices = NULL;
    int numscens = stochdataPtr->nscens;

    strcpy(lpname, "subproblp.lp");;


   //Create subproblem lp by deleting stage_1 cols and rows from core file lp

   //***************** Delete stage-1 rows (constraints) *********************************
   if (nrows_master > 0) {
         status = CPXdelrows(env, lp_sub, 0, nrows_master-1);
         if ( status ) {
               fprintf(stderr,"\nFunction loadSubProblemData(...): \n");
               fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
               return(status);
         }
   }
   // Get the sense vector for the range of contraints for the subproblem       //
  status = CPXgetsense(env, lp_sub, subprobPtr->sense, 0, subprobPtr->nrows-1);
  if ( status ) {
     fprintf (stderr, "loadSubProblemData:\n");
     fprintf (stderr, "Failure to get the senses from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem senses\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->sense[%d:] %c\n", j, subprobPtr->sense[j]);
       }
  #endif

   //////////////////// ADDED Jan 7, 2002 //////////////////
   //            Reset the <= contraints to >=            //
   /////////////////////////////////////////////////////////
   status = resetconstraints(env, lp_sub);
   if ( status ) {
       fprintf (stderr, "loadSubProblemData(...):\n");
       fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
       return status;
   }

   #ifdef DEBUG_FUNC
      fprintf (stdout, "writing subprob lp:\n");
        // Write sub prob to file in lp format
   	status = CPXwriteprob(env, lp_sub, "prob3.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "resetconstraints:\n");
      		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		return status;
   	}
   #endif
   ////////////////////////////////////////////////////////





   //**********Extract the T matrix from the core file and store it in T *****************
   status = CPXgetcols(env, lp_sub, &subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatind_T, subprobPtr->cmatval_T,
                       subprobPtr->cmatspace_T, &surplus, 0, ncols_master-1);
   if ( status ) {
      fprintf (stderr, "Failure to extraxt T from lp_sub from T(w) matrix, error %d.\n", status);
      fprintf (stderr, "Surplus value is: %d.\n", surplus);
      return(status);
   }

  //******************** Print the T matrix ******************************************
  #ifdef DEBUG_FUNC
     fprintf(stderr, "\nFunction loadSubProblemData(...): \n");
     fprintf(stderr, "subprobPtr->nzcnt_T = %d\n", subprobPtr->nzcnt_T);
     fprintf(stderr, "*******The T matrix is********: \n");
     printSparseMatrix(ncols_master, subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatbeg_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, stderr);
  #endif


  //********************** Delete stage-1 columns (vars) **********************************
  status = CPXdelcols(env, lp_sub, 0, ncols_master-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "Failure to delete lp_sub cols, error %d.\n", status);
     return(status);
  }

  //************* Get and store the ctype (var type) array from the subproblem lp ********
  status = CPXgetctype(env, lp_sub, subprobPtr->ctype, 0, ncols_sub-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "*Failure to get ctype from lp_sub, error %d.\n", status);
     return (status);
  }
  #ifdef DEBUG_FUNC
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf(stdout,"Subproblem ctype: \n");
     for (j = 0; j < ncols_sub; j++)
         fprintf(stdout, "col %d: %c\n", j, subprobPtr->ctype[j]);
  #endif



  // Store a vector of continuous ctypes for the linear relaxation
  // and ctype indices
  // Include feasibility column to be added soon
  subprobPtr->ctype[ncols_sub] = 'C';
  for (j = 0; j < ncols_sub+1; j++){
     subprobPtr->ctype_ctns[j] = 'C';
     subprobPtr->indices_ctype[j] = j;
  }


  /////////////////// Updates

  subprobPtr->nrows   = CPXgetnumrows (env, lp_sub);
  subprobPtr->nrows_T = CPXgetnumrows (env, lp_sub);
  stochdataPtr->nrows = CPXgetnumrows (env, lp_sub);;	// Initialize: VERY IMPORTANT!!!!!!
  #ifdef DEBUG_FUNC
       printf("subprobPtr->nrows = %d \n", subprobPtr->nrows);
       printf("subprobPtr->nrows_T = %d \n", subprobPtr->nrows_T);
       printf("stochdataPtr->nrows = %d \n", stochdataPtr->nrows);

        //Write sub prob to file in lp format
   	status = CPXwriteprob(env, lp_sub, "prob2.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "resetconstraints:\n");
      		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		return status;
   	}
    #endif




  //////////////////////////////////////////////////////////////////////////////
  // Add explicit feasibility column to the model to enable complete recourse //
  // as require by the D^2 algorithm                                          //
  //////////////////////////////////////////////////////////////////////////////

   status = addFeasColToSubProbLP(env, lp_sub, subprobPtr);

   if ( status ) {
  	fprintf (stderr, "loadSubProblemData:\n");
      	fprintf (stderr, "Failure to add binary contraints to lp_sub, error %d.\n", status);
      	return(status);
   }

  ///////////////////////////////////////////////////////////////////////
  // Add explicit binary contraints to the model for the D^2 algorithm //
  // cut generation.                                                   //
  ///////////////////////////////////////////////////////////////////////

   status = addBinaryConstrsSubProbLP(env, lp_sub, subprobPtr);
   if ( status ) {
  	fprintf (stderr, "loadSubProblemData:\n");
      	fprintf (stderr, "Failure to add binary contraints to lp_sub, error %d.\n", status);
      	return(status);
   }


   // Make updates
   subprobPtr->nrows_total  = CPXgetnumrows (env, lp_sub);;
   subprobPtr->nrowsW       = CPXgetnumrows (env, lp_sub);;
   ncols_sub          = CPXgetnumcols (env, lp_sub);
   subprobPtr->ncols  = ncols_sub;


  //********************** Write prob to file in lp format **************************
  #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_sub, lpname, NULL);
  	if ( status ) {
   		fprintf (stderr, "loadSubProblemData:\n");
      		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
      		return(status);
 	}
  #endif

  //************* Extract the obj from the lp_subfrom the core file ****************
  status = CPXgetobj(env, lp_sub, subprobPtr->obj, 0, ncols_sub-1);
  if ( status ) {
      fprintf (stderr, "loadSubProblemData:\n");
      fprintf (stderr, "Failure to get obj from lp_sub, error %d.\n", status);
      return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem obj\n");
       for (j = 0; j < subprobPtr->ncols; j++) {
           fprintf(stdout, "subprobPtr->obj[%d]: %f\n", j, subprobPtr->obj[j]);
       }
  #endif

  // Get the rhs vector for the range of contraints for the subproblem  //
  // and initialize each scenario rhs with these values                 //
  status = CPXgetrhs(env, lp_sub, subprobPtr->rhs, 0, nrows_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get the rhs from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem rhs\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->rhs[%d]: %f\n", j, subprobPtr->rhs[j]);
       }
  #endif

  // Initialize each scenario rhs with these values:
  // Initialize the scenario subproblem rhs(w)
  for (i = 0; i < numscens; i++){
      for (j = 0; j < subprobPtr->nrows; j++)
          stochdataPtr->rhs[i][j] = subprobPtr->rhs[j];
  }



  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadSubProblemData:\n");
     for (i = 0; i < numscens; i++){
        fprintf(stdout, "Subproblem rhs for scenario %d\n", i);
        for (j = 0; j < subprobPtr->nrows; j++) {
            fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", i, j, stochdataPtr->rhs[i][j]);
         }
     }
  #endif

  //Get the lower bounds on the subproblem cols (vars)
  status = CPXgetlb(env, lp_sub, subprobPtr->lb, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get lower bounds lb from lp_sub, error %d.\n", status);
     return (status);
  }
  //Get the upper bounds on the subproblem cols (vars)
  status = CPXgetub(env, lp_sub, subprobPtr->ub, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get upper bounds ub from lp_sub, error %d.\n", status);
     return (status);
  }


  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadSubProblemData:\n");
     fprintf(stdout, "\n\n\n  Subproblem lp column ub and lb:\n");
     for(i = 0; i < ncols_sub; i++){
        fprintf(stdout, "lb[%d] = %f	ub[%d] = %f \n", i,subprobPtr->lb[i],i,subprobPtr->ub[i]);
     }
  #endif


 //****** Initialize the indices array corresponding to the constraints *********//
 //************* for which the rhs coefs are to be changed *********************//
   for (j = 0; j < subprobPtr->nrows; j++) {
     subprobPtr->indices_row[j] = j;
   }
   // Remember to update this array every time a cut is added to the subproblem!!!
  #ifdef DEBUG_FUNC
      fprintf (stdout, "loadSubProblemData:\n");
      for(i = 0; i < subprobPtr->nrows; i++){
          fprintf(stdout, "subprobPtr->indices_row[%d] = %d \n", i, subprobPtr->indices_row[i]);
      }
  #endif

  return (0); // Successful return

 } //************************** End loadSubProblemData function *****************************


void
resetSubprobRhsAndT(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
                    double *solnX, int random_T, FILE *fpout)
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
 {
   int i, j, curr_row;
   int ncols = stochdataPtr->ncols;
   int scenario;
   int start, stop;

   //Reset the r(w) vector for all the scenarios
   for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
       for (i = 0; i < subprobPtr->nrows_T; i++){
            if (subprobPtr->sense[i] == 'L') // this row  was reset to 'G' by multiplying by -1
       	         stochdataPtr->rhs[scenario][i] = -1*stochdataPtr->rhs[scenario][i];
       }
   }

   if (random_T == 0 ) { // Constant T matrix
      // Reset the constant T matrix
      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->cmatbeg_T[i+1];
          } else {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->nzcnt_T;
          }

          for (j = start; j < stop; j++){
              curr_row = subprobPtr->cmatind_T[j];
              if (subprobPtr->sense[curr_row] == 'L')	// this row  was reset to 'G' by multiplying by -1
                   subprobPtr->cmatval_T[j] = -1*subprobPtr->cmatval_T[j];
          }
      } // End outer for loop

   } else { // Random T(w) matrix

      // Reset the random matrix T(w) for all scenarios
      for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
          for (i = 0; i < ncols; i++){
              if (i == ncols-1) {
                  start = stochdataPtr->cmatbeg_T[scenario][i];
                  stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
              } else {
                  start = stochdataPtr->cmatbeg_T[scenario][i];
                  stop  = stochdataPtr->cnzcnt_T[scenario];
              } // End if/else statement

              for (j = start; j < stop; j++){
                  curr_row = stochdataPtr->cmatind_T[scenario][j];
                  if (subprobPtr->sense[curr_row] == 'L')	// this row  was reset to 'G' by multiplying by -1
                       stochdataPtr->cmatval_T[scenario][j] = -1*stochdataPtr->cmatval_T[scenario][j];
              }
          }
      } // End outer scenario for loop

  } // end if/else statement


   //****** Print the Rhs for this scenario ******/
   #ifdef DEBUG_FUNC
        fprintf(stderr, "resetSubprobRhsAndT: \n");
        for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
             fprintf(stderr, "\n\nRhs for scenario %d \n", scenario);
  	     for (i = 0; i < subprobPtr->nrows_T; i++)
  	          fprintf(stderr, "r[%d] = %f\n", i, stochdataPtr->rhs[scenario][i]);
  	}

   #endif


   //****** Print the T(w) or T matrix ******/
   if (random_T == 0) {
      #ifdef DEBUG_FUNC
          fprintf(stderr, "\n\nresetSubprobRhsAndT: \n");
          fprintf(stderr, "Constant matrix T: \n");
  	     for (i = 0; i < subprobPtr->nzcnt_T; i++)
  	        fprintf(stderr, "cmatval_T[%d] = %f\n", i, subprobPtr->cmatval_T[i]);
       #endif

   } else {

       #ifdef DEBUG_FUNC
          fprintf(stderr, "\n\nresetSubprobRhsAndT(): \n");
          fprintf(stderr, "Random T(w): \n");
          for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
             fprintf(stderr, "\n\nT(w) for scenario %d \n", scenario);
  	     for (i = 0; i < stochdataPtr->cnzcnt_T[scenario]; i++)
  	          fprintf(stderr, "cmatval_T[%d][%d] = %f\n", scenario, i, stochdataPtr->cmatval_T[scenario][i]);
  	  }
       #endif
   }

 } //******************************* End resetSubprobRhsAndT *************************************



void
createsubprobRHS(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                 double *solnX, int random_T, FILE *fpout)
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
 {
   int i, j;
   int curr_row;
   int curr_col;
   int nrows;
   int ncols = stochdataPtr->ncols;
   int start, stop;


     //fprintf(stdout, "\n\n createsubprobRHS(): SCENARIO %d\n\n", scenario);

     //fprintf(stdout, " Initialization: \n");
   //Initialize subproblem RHS h(w) = r(w)
   for (i = 0; i < subprobPtr->nrows; i++){
        subprobPtr->rhsRho[scenario][i] = stochdataPtr->rhs[scenario][i];
        //fprintf(stdout, "subprobPtr->rhsRho[%d][%d] = %f \n", scenario, i, stochdataPtr->rhs[scenario][i]);
   }

   //fprintf(stdout, " End initialization: \n");

   if (random_T == 0) {  // Constant technology matrix T
      // Compute rho(w) col by col : T matrix is column by column sparse matrix

      //for (i = 0; i < ncols; i++)
      //      printf("subprobPtr->cmatbeg_T[%d] = %d  \n", i, subprobPtr->cmatbeg_T[i]);


      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->nzcnt_T;
          } else {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->cmatbeg_T[i+1];
          }
          //fprintf(stdout, "\n\n");
          for (j = start; j < stop; j++){
             //fprintf(stdout, "start = %d stop = %d\n", start, stop);
             curr_row = subprobPtr->cmatind_T[j];
             subprobPtr->rhsRho[scenario][curr_row] -= subprobPtr->cmatval_T[j]*solnX[i];
             //fprintf(stdout, "col = %d	cmatval_T[%d] = %f: ", i, j,
             //        subprobPtr->cmatval_T[j]);
             //fprintf(stdout, "solnX[%d]	= %6.6f\n", i, solnX[i]);
             //fprintf(stdout, "@subprobPtr->rhsRho[%d][%d] = %6.6f\n", scenario, curr_row,
             //                 subprobPtr->rhsRho[scenario][curr_row]);
          }
      }

   //exit(0);

   } else { //Random Technology matrix T(w)

      // Compute rho(w) col by col : T(w) matrix is column by column sparse matrix
      // This is the T(w) read from the stoch file
      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
          } else {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cnzcnt_T[scenario];
          }
          for (j = start; j < stop; j++){
             curr_row = stochdataPtr->cmatind_T[scenario][j];
             subprobPtr->rhsRho[scenario][curr_row] -= stochdataPtr->cmatval_T[scenario][j]*solnX[i];
             fprintf(stdout, "col = %d	val[%d] = %f: \n", i, j, stochdataPtr->cmatval_T[scenario][j]);
          }
      }

   } // End if/else statement


   // Compute rho(w) row by row for the newly added rows
   // The T(w) matrix for the newly added rows is in row-by-row sparse matrix format
   //fprintf(stdout, "\n\n\n subprobPtr->nrows > subprobPtr->nrows_T: %d > %d \n\n\n", subprobPtr->nrows, subprobPtr->nrows_T);
   if (subprobPtr->nrows > subprobPtr->nrows_T) { // New rows have been added
          #ifdef DEBUG_FUNC
              fprintf(stdout, "\n Number of rows added = %d \n", stochdataPtr->rnrows );

             printf("The T matrix is: \n");
             printMatrix(stochdataPtr->rnrows, stochdataPtr->rnzcnt_T[scenario], stochdataPtr->rmatbeg_T[scenario],
                      stochdataPtr->rmatind_T[scenario], stochdataPtr->rmatval_T[scenario], stdout) ;
          #endif
         // Continue computing rho(w) row by row due to the T(w) matrix, which is
         // in sparse format
         nrows = stochdataPtr->rnrows;
         for (i = 0; i < nrows; i++) {
             if (i == nrows-1) {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rnzcnt_T[scenario];
              } else {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rmatbeg_T[scenario][i+1];
              }

              curr_row = subprobPtr->nrows_T + i;
              //fprintf(stdout, "subprobPtr->nrows_T = %d\n", subprobPtr->nrows_T);
              //fprintf(stdout, "curr_row            = %d\n", curr_row);
              for (j = start; j < stop; j++){
                  curr_col = stochdataPtr->rmatind_T[scenario][j];
                  //fprintf(stdout, "stochdataPtr->rmatind_T[%d][%d] = %d\n",
                  //                    scenario, j, stochdataPtr->rmatind_T[scenario][j]);
                  if (curr_col >= stochdataPtr->ncols){
                     fprintf(stderr, "createsubprobRHS():\n");
                     fprintf(stderr, "*LARGE curr_col = %d, exiting...\n", curr_col);
                     fprintf(stdout, "Before: subprobPtr->rhsRho[%d][%d] = %f \n", scenario, curr_row, subprobPtr->rhsRho[scenario][curr_row]);
                     fprintf(stdout, "stochdataPtr->rmatval_T[%d][%d] = %f \n", scenario, j, stochdataPtr->rmatval_T[scenario][j]);
                     fprintf(stdout, "solnX[%d] = %f\n",curr_col, solnX[curr_col] );
                     fprintf(stdout, "stochdataPtr->rmatind_T[%d][%d] = %d\n",
                                      scenario, j, stochdataPtr->rmatind_T[scenario][j]);
                    fprintf(stdout, "stochdataPtr->rmatind_T[%d][%d] = %d\n",
                                      scenario, j+1, stochdataPtr->rmatind_T[scenario][j+1]);
                     exit(1);
                  }

                  subprobPtr->rhsRho[scenario][curr_row] -=
                              stochdataPtr->rmatval_T[scenario][j]*solnX[curr_col];
                  //fprintf(stdout, "subprobPtr->rhsRho[%d][%d] = %f \n", scenario, curr_row, subprobPtr->rhsRho[scenario][curr_row]);
              }

         } // end outer for loop

    } //End if statement


     //****** Print the RHS for this scenario ******/
     #ifdef DEBUG_FUNC
        fprintf(stdout, "\ncreatesubprobRHS(): \n");
        fprintf(stdout, "RHS for scenario %d \n", scenario);
  	for (i = 0; i < subprobPtr->nrows; i++)
  	    fprintf(stdout, "subprobPtr->rhsRho[%d][%d] = %f ", scenario, i, subprobPtr->rhsRho[scenario][i]);
     #endif


 } //******************************* End createsubprobRHS *************************************//


void
createsubprobMipRHS(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                 double *solnX, int random_T, FILE *fpout)
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
 {
   int i, j;
   int curr_row;
   int curr_col;
   int nrows;
   int ncols;
   int start, stop;

   ncols = stochdataPtr->ncols;

   //Initialize subproblem RHS h(w) = r(w)
   for (i = 0; i < subprobPtr->nrows; i++){
        subprobPtr->rhsRho[scenario][i] = stochdataPtr->rhs[scenario][i];
        //fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f \n", scenario, i, stochdataPtr->rhs[scenario][i]);
   }


   if (random_T == 0) {  // Constant technology matrix T
      // Compute rho(w) col by col : T matrix is column by column sparse matrix

      //for (i = 0; i < ncols; i++)
      //      printf("subprobPtr->cmatbeg_T[%d] = %d  \n", i, subprobPtr->cmatbeg_T[i]);


      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->nzcnt_T;
          } else {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->cmatbeg_T[i+1];
          }
          //fprintf(stdout, "\n\n");
          for (j = start; j < stop; j++){
             //fprintf(stdout, "start = %d stop = %d\n", start, stop);
             curr_row = subprobPtr->cmatind_T[j];
             subprobPtr->rhsRho[scenario][curr_row] -= subprobPtr->cmatval_T[j]*solnX[i];
             //fprintf(stdout, "col = %d	cmatval_T[%d] = %f: ", i, j,
             //        subprobPtr->cmatval_T[j]);
            // fprintf(stdout, "solnX[%d]	= %6.6f\n", i, solnX[i]);
            // fprintf(stdout, "subprobPtr->rhsRho[%d][%d] = %6.6f\n", scenario, curr_row,
            //                  subprobPtr->rhsRho[scenario][curr_row]);
          }
      }

   //exit(0);

   } else { //Random Technology matrix T(w)

      // Compute rho(w) col by col : T(w) matrix is column by column sparse matrix
      // This is the T(w) read from the stoch file
      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
          } else {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cnzcnt_T[scenario];
          }
          for (j = start; j < stop; j++){
             curr_row = stochdataPtr->cmatind_T[scenario][j];
             subprobPtr->rhsRho[scenario][curr_row] -= stochdataPtr->cmatval_T[scenario][j]*solnX[i];
             //fprintf(fpout, "col = %d	val[%d] = %f: \n", i, j, stochdataPtr->cmatval_T[scenario][j]);
          }
      }

   } // End if/else statement


     //****** Print the RHS for this scenario ******/
     #ifdef DEBUG_FUNC
        fprintf(fpout, "createsubprobRHS(): \n");
        fprintf(fpout, "RHS for scenario %d \n", scenario);
  	for (i = 0; i < subprobPtr->nrows; i++)
  	    fprintf(stdout, "rho(w)[%d] = %f ", i, subprobPtr->rhsRho[scenario][i]);
   #endif


 } //******************************* End createsubprobMipRHS *************************************//


void
storeSubProbSoln(solnMatrix_t *solnPtr, double *solnY, int scenario)
/**
 * Allocates memory to all scenario subproblem solution data structure
 * This is a row by row sparse matrix format
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param solnY  array of scenario subproblem solution
 * @param scenario   scenario index (starting from 0)
 */
 {
   int i, j;
   int nzcnt = 0; // nonzeros for this scenario soln


   if (scenario == 0) {             // Reset the solution matrix
   	solnPtr->nzcnt_S = 0;
   }

   solnPtr->rmatbeg_S[scenario] = solnPtr->nzcnt_S;

   //printf("ncols = %d \n", solnPtr->ncols);
  // printf("nrows = %d \n", solnPtr->nrows);
   // Loop thru the current scenario subproblem soln and store nonzeros
   for (j = 0; j < solnPtr->ncols; j++){
   	if (solnY[j] > NONZERO_LB) {
   	    i = solnPtr->nzcnt_S;
   	    solnPtr->rmatval_S[i] = solnY[j];
   	    solnPtr->rmatind_S[i] = j;
   	    solnPtr->nzcnt_S++;
   	    nzcnt++;
   	    //printf("nz[%d] = %f\n", j, solnPtr->rmatval_S[i] );
   	}
   }
   // Set the number of nonzeros for this scenario soln
   solnPtr->rmatcnt_S[scenario] = nzcnt;

   // Set the rmatbeg for scenario+1 to nzcnt: for convenience only
   solnPtr->rmatbeg_S[scenario+1] = solnPtr->nzcnt_S;

 }


void
storeSubProbIncumbSoln(solnMatrix_t *solnPtr, double *solnY, int scenario)
/**
 * Allocates memory to all scenario subproblem solution data structure
 * This is a row by row sparse matrix format
 * @param solnMatrix_t  pointer to the the soln sparse matrix structure
 * @param solnY  array of scenario subproblem solution
 * @param scenario   scenario index (starting from 0)
 */
 {
   int i, j;
   int nzcnt = 0; // nonzeros for this scenario soln


   if (scenario == 0) {             // Reset the solution matrix
   	solnPtr->nzcnt = 0;
   }

   solnPtr->rmatbeg[scenario] = solnPtr->nzcnt;

   //printf("ncols = %d \n", solnPtr->ncols);
  // printf("nrows = %d \n", solnPtr->nrows);
   // Loop thru the current scenario subproblem soln and store nonzeros
   for (j = 0; j < solnPtr->ncols; j++){
   	if (solnY[j] > NONZERO_LB) {
   	    i = solnPtr->nzcnt;
   	    solnPtr->rmatval[i] = solnY[j];
   	    solnPtr->rmatind[i] = j;
   	    solnPtr->nzcnt++;
   	    nzcnt++;
   	    //printf("nz[%d] = %f\n", j, solnPtr->rmatval[i] );
   	}
   }
   // Set the number of nonzeros for this scenario soln
   solnPtr->rmatcnt[scenario] = nzcnt;

   // Set the rmatbeg for scenario+1 to nzcnt: for convenience only
   solnPtr->rmatbeg[scenario+1] = solnPtr->nzcnt;

 } //*********** End storeSubProbIncumbSoln ************


int
getDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                  int *disj_ind, char* ctype)
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
{
  int i, j;
  double lval;
  double uval;

  for (i = 0; i < solnPtr->nrows; i++){ // do row by row (scenario) scanning
        //printf("scenario %d\n", i);
        //printf("beg[i] = %d\n", solnPtr->rmatbeg_S[i]);
        //printf("beg[i+1] = %d\n", solnPtr->rmatbeg_S[i+1]);
	for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
	{
	   if (ctype[j] != 'C') {
	      //printf("val[%d] = %f", j, solnPtr->rmatval_S[j]);
   	      lval = solnPtr->rmatval_S[j] - floor(solnPtr->rmatval_S[j]);
	      uval = ceil(solnPtr->rmatval_S[j]) - solnPtr->rmatval_S[j];
	      //printf("lval = %f	uval = %f\n", lval, uval);
	      if (lval > INT_PRECISION || uval > INT_PRECISION){ // Fractional component
	         //printf("disj_scen  = %d\n", i);
	         //printf("disj_var  = %d\n", solnPtr->rmatind_S[j]);
	        *disj_var = solnPtr->rmatind_S[j];     // column index
	        *disj_scen = i; // row index
	        *disj_ind = j;

	        return 1;
	      }
	    } // End outer if
	}
   } // end outer for loop

   return 0;

} //****************** End getDisjunctionVar() ***************************




int
getNextDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                      int *disj_ind, char *ctype)
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
{
  int i, j;
  double lval;
  double uval;
  int start;

  for (i = *disj_scen; i < solnPtr->nrows; i++){ // do row by row (scenario) scanning
	if (i == *disj_scen)
	    start =  *disj_ind+1;
	else
            start = solnPtr->rmatbeg_S[i];

	//printf("scenario %d\n", i);
       // printf("beg[i] = %d\n", start);
       // printf("beg[i+1] = %d\n", solnPtr->rmatbeg_S[i+1]);

	for (j = start; j < solnPtr->rmatbeg_S[i+1]; j++)
	{
	   if (ctype[j] != 'C') {
	      lval = solnPtr->rmatval_S[j] - floor(solnPtr->rmatval_S[j]);
	      uval = ceil(solnPtr->rmatval_S[j]) - solnPtr->rmatval_S[j];
	      if (lval> INT_PRECISION || uval > INT_PRECISION){ // Fractional component
	         *disj_var = solnPtr->rmatind_S[j];     // column index
	         *disj_scen = i; // row index
	         *disj_ind = j;

	          return 1;
	       }
	    } // End outer if
         }
   } // end outer for loop

	return 0;

} //****************** End getNextDisjunctionVar() ***************************


int
getMaxDisjunctionVar(solnMatrix_t *solnPtr, int *disj_var, int *disj_scen,
                     int *disj_ind, char* ctype)
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
{
  int i, j;
  int status = 0;
  double lval, uval;
  double curr_max;
  double max = 0;

  for (i = 0; i < solnPtr->nrows; i++) // do row by row (scenario) scanning
  {
	for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
	{
	   if (ctype[solnPtr->rmatind_S[j]] != 'C') {
	      lval = solnPtr->rmatval_S[j] - floor(solnPtr->rmatval_S[j]);
	      uval = ceil(solnPtr->rmatval_S[j]) - solnPtr->rmatval_S[j];
              //if (lval > INT_PRECISION || uval > INT_PRECISION){ // Fractional component
              if (solnPtr->rmatval_S[j] > INT_PRECISION && solnPtr->rmatval_S[j] <1- INT_PRECISION){ // Fractional component
		  curr_max = minimum(lval, uval);
		  if (curr_max > max){
		      max = curr_max;
		      *disj_var = solnPtr->rmatind_S[j];     // column index
		      *disj_scen = i; // row index
		      *disj_ind = j;
		      status = 1;
	           }
             } // end outer if
           } // End outer if
        }
   } // end outer for loop

	return status;

} //****************** End getMaxDisjunctionVar() ***************************


void
getAverageSubProblemsSoln(solnMatrix_t *solnPtr, double *scenProb, double *averSoln)
/**
 * This function computes an averaged solutions for all scenario subproblems
 * based on the probability of outcome for each scenario.
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param scenProb pointer to array containing the outcome prob for each scenario
 *@param averSoln pointer to array for storing the averaged solution
 */
{
  int i, j;
  int curr_col;


  // Initialize the averSoln array with zeros
  for (i = 0; i < solnPtr->ncols; i++)
      averSoln[i] = 0;

  // Average the solutions scenario by scenario
  for (i = 0; i < solnPtr->nrows; i++) // do row by row (scenario) scanning
  {

     for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
     {
	curr_col = solnPtr->rmatind_S[j];
	averSoln[curr_col] += scenProb[i]*solnPtr->rmatval_S[j];
	//printf("scenProb[%d] = %f\n", i, scenProb[i]);
	//printf("val[%d] = %f\n", j, solnPtr->rmatval_S[j]);

     }
  } // end outer for loop

} //****************** End getAverageSoln **********************************//



void
getCondAverageSubProblemsSoln(solnMatrix_t *solnPtr, double *scenProb, int disj_var_index,
                              int disj_scen, double *averSoln)
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
{

  int i, j;
  int curr_col;
  double sumprobFrac;
  int frac;
  double lval;
  double uval;


  //printf("Computing Conditional Average Soln \n");
  // Initialize the averSoln array with zeros
  for (i = 0; i < solnPtr->ncols; i++)
      averSoln[i] = 0;

  // Initialize counter for scenarions with fractiona disj_var_index component
  solnPtr->ncondScens = 0;

  sumprobFrac = 0;

  // average the solutions scenario by scenario for which the disj_var_index
  // component is fractional
  for (i = 0; i < solnPtr->nrows; i++) // do row by row (scenario) scanning
  {
     //printf("scenario %d \n", i);
     frac = 0;
     // Check if disj_var_index is fractional for this scenario
     for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
     {
	curr_col = solnPtr->rmatind_S[j];
	//printf("disj_var = %d \n", disj_var_index);
	//printf("curr_col = %d \n", curr_col);
	if (disj_var_index == curr_col) {
	    //printf("curr_col == disj_var \n");
	    lval = solnPtr->rmatval_S[j] - floor(solnPtr->rmatval_S[j]);
	    uval = ceil(solnPtr->rmatval_S[j]) - solnPtr->rmatval_S[j];
	    if (lval > INT_PRECISION || uval > INT_PRECISION){ // Fractional component
	         frac = 1;
	         solnPtr->condScens[solnPtr->ncondScens] = i;
	         sumprobFrac += scenProb[i];
	         //printf("sumprobFrac = %f\n", sumprobFrac);
	        // printf("scenProb[%d] = %f\n", i, scenProb[i]);
	        // printf("Scenario [%d] fractional component %f\n",
	        //         solnPtr->condScens[solnPtr->ncondScens], solnPtr->rmatval_S[j]);
	         if (i == disj_scen)
	             solnPtr->disj_scen_index = solnPtr->ncondScens;
	         solnPtr->ncondScens++;
	    }
	    break;
	} // End outer if
     } // End for loop

     if (sumprobFrac == 1){ // All scenarios have frac disj_var component!
         sumprobFrac = 0;
     }
     if (frac == 1) { // Average this scenario solution
         for (j = solnPtr->rmatbeg_S[i]; j < solnPtr->rmatbeg_S[i+1]; j++)
         {
	   curr_col = solnPtr->rmatind_S[j];
	   averSoln[curr_col] += scenProb[i]*solnPtr->rmatval_S[j];
	   //printf("averCondSoln[%d] = %f\n", curr_col, averSoln[curr_col]);
         }
     } // end if statement

  } // End outer for loop

  //printf("sumprobFrac = %f\n", sumprobFrac);
 // printf("solnPtr->ncondScens = %d \n", solnPtr->ncondScens);
 // printf("averCondSoln before conditioning \n");
 // for (j = 0; j < solnPtr->ncols ; j++){
 //     printf("averCondSoln[%d] = %f\n", j, averSoln[j]);
 // }

  // Compute conditional average values
  for (j = 0; j < solnPtr->ncols ; j++){
      averSoln[j] = averSoln[j]/sumprobFrac;
      //printf("averCondSoln[%d] = %f\n", j, averSoln[j]);
  }

 // printf("solnPtr->ncondScens = %d\n", solnPtr->ncondScens);
  // Compute conditional probs for the frac (disj_var_index) scenarios
  for (i = 0; i < solnPtr->ncondScens ; i++){
      j = solnPtr->condScens[i];
      solnPtr->condScenProbs[i] = scenProb[j]/sumprobFrac;
      //printf("condScenProb[%d] = %f\n", j, solnPtr->condScenProbs[i]);
  }

  // Initialize scenario to drop in case the current fractional soln
  // is not cut off
  solnPtr->scenToDrop = 0;

  //printf("Done Computing Conditional Average Soln \n");

} //****************** End getCondAverageSoln **********************************//



void
getConditionC3objCoefs(solnMatrix_t *solnPtr, double *scenProb, double *condC3Obj,
                       int scen)
/**
 * This function computes an averaged solutions for all scenario subproblems
 * based on the probability of outcome for each scenario.
 * @param solnPtr pointer to solution sparse matrix for all the
 *                       scenarios
 *@param scenProb pointer to array containing the outcome prob for each scenario
 *@param condC3Obj pointer to array for storing the conditional solution
 *@param scen scenario whose soln will be used as c^3 lp obj
 */
{
  int j;
  int curr_col;
  int start, stop;

  // Initialize the averSoln array with zeros
  for (j = 0; j < solnPtr->ncols; j++)
      condC3Obj[j] = 0;

  if (scen == solnPtr->nrows-1) {
       start = solnPtr->rmatbeg_S[scen];
       stop  = solnPtr->nzcnt_S;
   } else {
       start = solnPtr->rmatbeg_S[scen];
       stop  = solnPtr->rmatbeg_S[scen+1];
   }

  // Get conditional solutions of scenario subproblem
  for (j = start; j < stop; j++) {
      curr_col = solnPtr->rmatind_S[j];
      condC3Obj[curr_col] += scenProb[scen]*solnPtr->rmatval_S[j];
	//printf("scenProb[%d] = %f\n", scen, scenProb[scen]);
	//printf("val[%d] = %f\n", j, solnPtr->rmatval_S[j]);
  }

} //****************** End getAverageSoln **********************************//


void
loadC3LPdata(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout)
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
 {
   int i, j, k;
   int m2;        // Number of rows for this W^k matrix
   int n2;        // Number of cols for this W^k matrix including the feasibility variable
   int nscens;    // Number of subproblem scenarios

   int nrows; 	  // Number of rows for this C^3 LP = c3lpPtr->nrows;
   int ncols;	  // Number of cols for this C^3 LP = c3lpPtr->ncols;

   int col_start, row_start1, row_start2;
   int col_count;
   int row_start, row_stop;
   int mynd2cuts = 0;		// counter for number of d^2 cuts added for this disj var
   int curr_row;
   int addthis_row;             // A zero value indicates to extract row from W^k matrix
    				// Otherwise, skip the row, it's not for this disj var
   int row_index;
   int cnt;

   // Compute the number of rows added for the current disjunction variable thus far
   //printf("\nsubprobPtr->nd2cuts = %d\n",subprobPtr->nd2cuts);
   for(i = 0; i < subprobPtr->nd2cuts; i++){
      if (subprobPtr->disjVarWindex[i] == disj_var_index)
         mynd2cuts++;
   }

   m2 = subprobPtr->nrowsW + mynd2cuts;
   //m2 = subprobPtr->nrowsW;
   c3lpPtr->nrowsWk = m2; // update number of rows in W^k matrix for this disj_var
   n2 = subprobPtr->ncols;
   nscens =  stochdataPtr->nscens;
   nrows = 2*n2 + 2*nscens;
   ncols = n2 + nscens + 2*m2 + 2;
   c3lpPtr->nrows = nrows;
   c3lpPtr->ncols = ncols;

   #ifdef DEBUG_FUNC
      fprintf(stdout, " **nrowsW = %d**\n", subprobPtr->nrowsW);
      fprintf(stdout, " **mynd2cuts = %d**\n", mynd2cuts);
      fprintf(stdout, " **c3 rows = %d**\n", m2);
      fprintf(stdout, " **c3lpPtr->nrows = %d**\n", c3lpPtr->nrows);
      fprintf(stdout, " **c3lpPtr->ncols = %d**\n", c3lpPtr->ncols);
      fprintf(stdout, " **subprobPtr->nrows = %d**\n", subprobPtr->nrows);
      fprintf(stdout, " **subprobPtr->ncols = %d**\n", subprobPtr->ncols);
   #endif

   #ifdef DEBUG_MAIN
        fprintf (stdout, "\nd2algMain():\n");
        fprintf (stdout, "Average subproblem solution is > \n");
        for (i = 0; i < subprobPtr->ncols; i++)
            fprintf(stdout, "aveSubprobSoln[%d] = %f\n", i, aveSubprobSoln[i] );
    #endif
        //fprintf(stdout, "ECHO: \n" );

   //******* Set object coefs ************//
   // pi_i's column
   for (i = 0; i < n2-1; i++) {
       c3lpPtr->obj[i] = aveSubprobSoln[i];
       //fprintf(stdout, "c3lpPtr->obj[%d] = %f\n", i, c3lpPtr->obj[i]);
   }

   // pi feasibility column
   c3lpPtr->obj[n2-1] = 0;

   // pi_0's column
   for (i = n2; i < n2+nscens; i++) {
       //if (i == disj_scen+n2)
           c3lpPtr->obj[i] = -1*stochdataPtr->scenProb[i-n2];
       //else
       //    c3lpPtr->obj[i] = 0;
       //fprintf(stderr, "\n\n**c3lpPtr->obj[i] = %d**\n",i, c3lpPtr->obj[i]);
   }

   // lambda's column
   for (i = n2+nscens; i < ncols; i++) {
       c3lpPtr->obj[i] = 0;
   }

   //*********** Set the ctype array continuous vars ********//
   for (i = 0; i < ncols; i++)      // Without t_00 and t_01 cols
       c3lpPtr->ctype[i] = 'C';

   //************* Set the sense array to >= ***************//
   for (i = 0; i < nrows; i++)       /// Without t_00 and t_01 cols
       c3lpPtr->sense[i] = 'G';

   //************** Set rhs coefs **************************//
   for (i = 0; i < nrows; i++)
        c3lpPtr->rhs[i] = 0;

   //*********** Set the bounds on each variable ************//
   // pi_i's and pi_0's variables
   for (i = 0; i < n2+nscens; i++){
       c3lpPtr->lb[i] = -1;
       c3lpPtr->ub[i] =  1;
   }

   // lambda's columns
   for (i = n2+nscens; i < ncols; i++) {
       if (i-nscens < subprobPtr->ncols) {
          if ( subprobPtr->sense[i-nscens] == 'E') // Make multiplier free
              c3lpPtr->lb[i] = -CPX_INFBOUND;
           else
               c3lpPtr->lb[i] = 0;
       } else {
          c3lpPtr->lb[i] = 0;
       }
       c3lpPtr->ub[i] =  CPX_INFBOUND;
   }

   //************** Set the constraint matrix coefs in sparse format ***********//

   // Set the column count for cols pi_i's and pi_0's
   for (i = 0; i < n2+nscens; i++) {
       c3lpPtr->cmatcnt[i] = 2;
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d\n", i, c3lpPtr->cmatcnt[i]);
   }
   // Initialize to zero for the lambda columns.
   for (i = n2+nscens; i < ncols; i++) {
       c3lpPtr->cmatcnt[i] = 0;   // Initialize to 0
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d \n", i, c3lpPtr->cmatcnt[i]);
   }
   //fprintf(stdout, "Insert pi's coefs in C^3 LP");

   // pi_i's columns
   c3lpPtr->nzcnt = 0; // Initialize nonzero counter
   row_start1 = 0;
   row_start2 = n2;
   for (i = 0; i < n2; i++) {
   	c3lpPtr->cmatbeg[i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }
   //fprintf(stdout, "Insert pi_0's coefs in C^3 LP");

   // pi_0's columns
   row_start1 = 2*n2;
   row_start2 = 2*n2+nscens;
   for (i = 0; i < nscens; i++) {
   	c3lpPtr->cmatbeg[n2+i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }

  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);

   //lambda_01's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens) to c < (n2+nscens+m2)
   //rols start at r = 0 to r < n2 and r = 2*n2 to r < (2*n2 + nscens)

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   // Just extract the original W rows and then the rows for 	//
   // this disjunction variable					//

   //fprintf(stdout, "Insert lambda_01 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable
   for (i = 0; i < subprobPtr->nrows; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
                 // printf ("Skip i            = %d \n", i);
             } else {// extract this previously added pi row for this disj var
                  addthis_row = 0;
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                 // printf ("Add i            = %d \n", i);
                 // printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }
       //printf ("c3 lp W^k row i            = %d \n", i);
      // printf ("c3 lp W^k pi row_index     = %d \n", row_index);

       if (addthis_row == 0) {
           //fprintf(stdout, "W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2+nscens] = c3lpPtr->nzcnt;
          //printf("\n\n*******Lambda11: c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row in W^k
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
       	  } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
       	  }
          //fprintf(stdout, "W^k row_start = %d \n", row_start);
          //fprintf(stdout, "W^k row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
             //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];

                //fprintf(stdout, "Row = %d  val[%d] = %f\n", row_index,
                //                 c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                c3lpPtr->cmatind[c3lpPtr->nzcnt] = subprobPtr->cmatind_W[j];
                c3lpPtr->nzcnt++;      			// count nonzeros
                c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
             }
           }

           // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
           for (k = 0; k < nscens; k++) {               // k = scenario
               // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
               if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                   c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
                   //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][i]);
        	   c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	   c3lpPtr->nzcnt++;			// count nonzeros
        	   c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
               }
           }
       } // End if (addthis_row)

   } // End outer for loop

   //lambda_11's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens+m2) to c < ncols
   //rols start at r = n2 to r < 2*n2 and r = 2*n2+S to r < nrows

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   //fprintf(stdout, "Insert lambda_11 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable

   for (i = 0; i < subprobPtr->nrows_total; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
             } else {// extract this previously added pi row
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                  //printf ("Add i            = %d \n", i);
                  //printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }

      // printf ("c3 lp W^k row i   = %d\n", i);
      // printf ("c3 lp W^k pi row_index = %d\n", row_index);

       if (addthis_row == 0) {
       	  // fprintf(stdout, "W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2+nscens+m2] = c3lpPtr->nzcnt;
          //printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
          } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
          }

         // fprintf(stdout, "row_start = %d \n", row_start);
         // fprintf(stdout, "row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
          //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                 c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];
                // fprintf(stdout, "Row = %d  val[%d] = %f\n",
                //         row_index, c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2 + subprobPtr->cmatind_W[j];
                c3lpPtr->nzcnt++;      			// count nonzeros
                c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;          // count nonzeros in this column
             }
          }

          // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
          for (k = 0; k < nscens; k++) {               // k = scenario
             // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
             if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                   c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
        	   c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	   c3lpPtr->nzcnt++;			// count nonzeros
        	   c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;  // count nonzeros in this column
             }
          }
       } // End if(addthis_row)
   } // End outer for loop

  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
 //  printf("*******ncols-2 = %d\n", ncols-2);

   // lambda_02's column: add the disjunction var floor coefs
   c3lpPtr->cmatbeg[ncols-2] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-2]++; 				// count nonzeros in this column

   //fprintf(stdout, "disj_var_floor = %f \n", i, disj_var_floor);
   // Add the disj_var_floor coef to this column
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_floor > NONZERO_LB) {
            fprintf(stdout, "******** disj_var_floor = %f \n", i, disj_var_floor);
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*disj_var_floor;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-2]++; 			// count nonzeros in this column
       }
    }

   // lambda_12's column: add the disjunction var ceil coefs
   c3lpPtr->cmatbeg[ncols-1] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2+disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-1]++; 				// count nonzeros in this column

  // fprintf(stdout, "disj_var_ceil = %f \n", i, disj_var_ceil);
   // Add the disj_var_ceil coef to this column
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_ceil > NONZERO_LB) {
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = disj_var_ceil;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + nscens + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-1]++; 			// count nonzeros in this column
       }
    }

    c3lpPtr->cmatbeg[ncols] = c3lpPtr->nzcnt;

 } //******************************* loadC3LPdata ************************//





void
loadC3LPdataNew(c3lpProb_t *c3lpPtr, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
             double *aveSubprobSoln, double *solnX, int disj_var_index, int disj_scen,
             int disj_var_floor, int disj_var_ceil, int random_T, FILE *fpout)
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
 {
   int i, j, k;
   int m2;        // Number of rows for this W^k matrix
   int n2;        // Number of cols for this W^k matrix including the feasibility variable
   int nscens;    // Number of subproblem scenarios

   int nrows; 	  // Number of rows for this C^3 LP = c3lpPtr->nrows;
   int ncols;	  // Number of cols for this C^3 LP = c3lpPtr->ncols;

   int col_start, row_start1, row_start2;
   int col_count;
   int row_start, row_stop;
   int mynd2cuts = 0;		// counter for number of d^2 cuts added for this disj var
   int curr_row;
   int addthis_row;             // A zero value indicates to extract row from W^k matrix
    				// Otherwise, skip the row, it's not for this disj var
   int row_index;
   int cnt;
   int my_index;

   // Compute the number of rows added for the current disjunction variable thus far
   //printf("\nsubprobPtr->nd2cuts = %d\n",subprobPtr->nd2cuts);
   for(i = 0; i < subprobPtr->nd2cuts; i++){
      if (subprobPtr->disjVarWindex[i] == disj_var_index)
         mynd2cuts++;
   }

   m2 = subprobPtr->nrowsW + mynd2cuts;
   //m2 = subprobPtr->nrowsW;
   c3lpPtr->nrowsWk = m2; // update number of rows in W^k matrix for this disj_var
   n2 = subprobPtr->ncols;
   nscens =  stochdataPtr->nscens;
   nrows = 2*n2 + 2*nscens;
   ncols = n2 + nscens + 2*m2 + 2;
   c3lpPtr->nrows = nrows;
   c3lpPtr->ncols = ncols;

   #ifdef DEBUG_FUNC
      fprintf(stdout, " **nrowsW = %d**\n", subprobPtr->nrowsW);
      fprintf(stdout, " **mynd2cuts = %d**\n", mynd2cuts);
      fprintf(stdout, " **c3 rows = %d**\n", m2);
      fprintf(stdout, " **c3lpPtr->nrows = %d**\n", c3lpPtr->nrows);
      fprintf(stdout, " **c3lpPtr->ncols = %d**\n", c3lpPtr->ncols);
      fprintf(stdout, " **subprobPtr->nrows = %d**\n", subprobPtr->nrows_total);
      fprintf(stdout, " **subprobPtr->ncols = %d**\n", subprobPtr->ncols);
   #endif

   #ifdef DEBUG_MAIN
       fprintf (stdout, "\nd2algMain():\n");
       fprintf (stdout, "Average subproblem solution is > \n");
       for (i = 0; i < subprobPtr->ncols; i++)
             fprintf(stdout, "aveSubprobSoln[%d] = %f\n", i, aveSubprobSoln[i] );
   #endif


   //******* Set object coefs ************//
   // pi_i's column
   for (i = 0; i < n2-1; i++) {
       c3lpPtr->obj[i] = aveSubprobSoln[i];
       //fprintf(stdout, "c3lpPtr->obj[%d] = %f\n", i, c3lpPtr->obj[i]);
   }

   // pi feasibility column
   c3lpPtr->obj[n2-1] = 0;

   // pi_0's column
   for (i = n2; i < n2+nscens; i++) {
       //if (i == disj_scen+n2)
           c3lpPtr->obj[i] = -1*stochdataPtr->scenProb[i-n2];
       //else
       //    c3lpPtr->obj[i] = 0;
       //fprintf(stderr, "\n\n**c3lpPtr->obj[i] = %d**\n",i, c3lpPtr->obj[i]);
   }

   // lambda's column
   for (i = n2+nscens; i < ncols; i++) {
       c3lpPtr->obj[i] = 0;
   }

   //*********** Set the ctype array continuous vars ********//
   for (i = 0; i < ncols; i++)      // Without t_00 and t_01 cols
       c3lpPtr->ctype[i] = 'C';

   //************* Set the sense array to >= ***************//
   for (i = 0; i < nrows; i++)      // Without t_00 and t_01 cols
       c3lpPtr->sense[i] = 'G';

   //************** Set rhs coefs **************************//
   for (i = 0; i < nrows; i++)
        c3lpPtr->rhs[i] = 0;

   //*********** Set the bounds on each variable ************//
   // pi_i's and pi_0's variables
   for (i = 0; i < n2+nscens; i++){
       c3lpPtr->lb[i] = -1;
       c3lpPtr->ub[i] =  1;
   }

   // lambda's columns
   for (i = n2+nscens; i < ncols; i++) {
       if (i-nscens < subprobPtr->ncols) {
          if ( subprobPtr->sense[i-nscens] == 'E') // Make multiplier free
              c3lpPtr->lb[i] = -CPX_INFBOUND;
           else
               c3lpPtr->lb[i] = 0;
       } else {
          c3lpPtr->lb[i] = 0;
       }
       c3lpPtr->ub[i] =  CPX_INFBOUND;
   }

   //************** Set the constraint matrix coefs in sparse format ***********//

   // Set the column count for cols pi_i's and pi_0's
   for (i = 0; i < n2+nscens; i++) {
       c3lpPtr->cmatcnt[i] = 2;
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d\n", i, c3lpPtr->cmatcnt[i]);
   }
   // Initialize to zero for the lambda columns.
   for (i = n2+nscens; i < ncols; i++) {
       c3lpPtr->cmatcnt[i] = 0;   // Initialize to 0
       //fprintf(stdout, "c3lpPtr->cmatcnt[%d] = %d \n", i, c3lpPtr->cmatcnt[i]);
   }
   //fprintf(stdout, "Insert pi's coefs in C^3 LP");

   // pi_i's columns
   c3lpPtr->nzcnt = 0; // Initialize nonzero counter
   row_start1 = 0;
   row_start2 = n2;
   for (i = 0; i < n2; i++) {
   	c3lpPtr->cmatbeg[i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }
   //fprintf(stdout, "Insert pi_0's coefs in C^3 LP");

   // pi_0's columns
   row_start1 = 2*n2;
   row_start2 = 2*n2+nscens;
   for (i = 0; i < nscens; i++) {
   	c3lpPtr->cmatbeg[n2+i] = c3lpPtr->nzcnt;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start1;
        c3lpPtr->nzcnt++;
        c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        c3lpPtr->cmatind[c3lpPtr->nzcnt] = row_start2;
        c3lpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }

  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);

   //lambda_01's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens) to c < (n2+nscens+m2)
   //rols start at r = 0 to r < n2 and r = 2*n2 to r < (2*n2 + nscens)

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   // Just extract the original W rows and then the rows for 	//
   // this disjunction variable					//

   //fprintf(stdout, "Insert lambda_01 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable
   for (i = 0; i < subprobPtr->nrows_total; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
                 // printf ("Skip i            = %d \n", i);
             } else {// extract this previously added pi row for this disj var
                  addthis_row = 0;
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                 // printf ("Add i            = %d \n", i);
                 // printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }
       //printf ("c3 lp W^k row i            = %d \n", i);
      // printf ("c3 lp W^k pi row_index     = %d \n", row_index);

       if (addthis_row == 0) {
          // fprintf(stdout, "          W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2+nscens] = c3lpPtr->nzcnt;
          //printf("\n\n*******Lambda11: c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row in W^k
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
       	  } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
       	  }
          //fprintf(stdout, "W^k row_start = %d \n", row_start);
          //fprintf(stdout, "W^k row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
              //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                 c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];
                 //fprintf(stdout, "Row = %d  val[%d] = %f\n", row_index,
                 //                 c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                 c3lpPtr->cmatind[c3lpPtr->nzcnt] = subprobPtr->cmatind_W[j];
                 c3lpPtr->nzcnt++;      			// count nonzeros
                 c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
             }
           }

           // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
           for (k = 0; k < nscens; k++) {               // k = scenario
               // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
               ///////////// 6-9-03 ////////////////
             if (i < subprobPtr->nrows_T ) {
                  if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                      c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
                      //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][i]);
        	      c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	      c3lpPtr->nzcnt++;			// count nonzeros
        	      c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
                  }
              } else if (i < subprobPtr->nrowsW) {  // Add the rhs binary bounds
                   c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        	   c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	   c3lpPtr->nzcnt++;			         // count nonzeros
        	   c3lpPtr->cmatcnt[row_index+n2+nscens]++;      // count nonzeros in this column
              } else { // Added rows to W^k
                    my_index = i - subprobPtr->nrowsW + subprobPtr->nrows_T;
                   if (subprobPtr->rhsRho[k][my_index] > NONZERO_LB || subprobPtr->rhsRho[k][my_index] < -NONZERO_LB) {
                      c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][my_index];
                      //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][my_index]);
        	      c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
        	      c3lpPtr->nzcnt++;			// count nonzeros
        	      c3lpPtr->cmatcnt[row_index+n2+nscens]++;  // count nonzeros in this column
                  }
              }
               /////////////////////////////////////
           }
       } // End if (addthis_row)

   } // End outer for loop

   //lambda_11's columns: add the matrix W coefs and rhs rho(w) coefs
   //cols start at c = (n2+nscens+m2) to c < ncols
   //rols start at r = n2 to r < 2*n2 and r = 2*n2+S to r < nrows

   // *** Loop over W^k nonzeros and extract each row data ***//
   // Remember W^k is in row sparse matrix format //
   //fprintf(stdout, "Insert lambda_11 coefs in C^3 LP");

   cnt = 0; // start row index counter for added pi rows for this
            // disjunction variable

   for (i = 0; i < subprobPtr->nrows_total; i++) {   // i = row of W^k
       addthis_row = 0;
       // Extract the appropriate row for this disj_var
       if (i >= subprobPtr->nrowsW){ // done extracting the original W rows
             curr_row = i - subprobPtr->nrowsW;
             if (subprobPtr->disjVarWindex[curr_row] != disj_var_index) {
                  addthis_row = 1; // skip this row
             } else {// extract this previously added pi row
                  row_index = subprobPtr->nrowsW + cnt;
                  cnt++;
                  //printf ("Add i            = %d \n", i);
                  //printf ("Add i            = %d \n", i);
             }
       } else {
         row_index = i;
       }

      // printf ("c3 lp W^k row i   = %d\n", i);
      // printf ("c3 lp W^k pi row_index = %d\n", row_index);

       if (addthis_row == 0) {
       	  //fprintf(stdout, "          W^k row %d added to c^3 lp matrix\n", i);
          // Set col  beg for this C^3 constraint matrix col
          c3lpPtr->cmatbeg[row_index+n2+nscens+m2] = c3lpPtr->nzcnt;
          //printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
          if (i == subprobPtr->nrows_total-1) { // last row
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->nzcnt_W;
          } else {
             row_start = subprobPtr->cmatbeg_W[i];
             row_stop  = subprobPtr->cmatbeg_W[i+1];
          }

         // fprintf(stdout, "row_start = %d \n", row_start);
         // fprintf(stdout, "row_stop = %d \n", row_stop);

          for (j = row_start; j < row_stop; j++){
          //////// ADDED JUNE 11 /////////////////
             if (subprobPtr->cmatval_W[j] > NONZERO_LB || subprobPtr->cmatval_W[j] < -NONZERO_LB){
                 c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*subprobPtr->cmatval_W[j];
                // fprintf(stdout, "Row = %d  val[%d] = %f\n",
                //         row_index, c3lpPtr->nzcnt, c3lpPtr->cmatval[c3lpPtr->nzcnt]);
                c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2 + subprobPtr->cmatind_W[j];
                c3lpPtr->nzcnt++;      			// count nonzeros
                c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;          // count nonzeros in this column
             }
          }

          // Add the subproblem rhsRho[scenario][row] coef for each scenario omega
          for (k = 0; k < nscens; k++) {               // k = scenario
             // Get the rhs h(w) = r(w) - T(w)x = subprobPtr->rhsRho[k][i] for this scenario
             /////////////// 6-9-03 ////////////////
             if (i < subprobPtr->nrows_T) {
                if (subprobPtr->rhsRho[k][i] > NONZERO_LB || subprobPtr->rhsRho[k][i] < -NONZERO_LB) {
                     c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][i];
                     //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][i]);
        	     c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	     c3lpPtr->nzcnt++;			// count nonzeros
        	     c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;  // count nonzeros in this column
                }
              } else if (i < subprobPtr->nrowsW){ // Add the rhs binary bounds
                    c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
        	    c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	    c3lpPtr->nzcnt++;			                // count nonzeros
        	    c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;         // count nonzeros in this column
              } else {
                   my_index = i - subprobPtr->nrowsW + subprobPtr->nrows_T;
                   if (subprobPtr->rhsRho[k][my_index] > NONZERO_LB || subprobPtr->rhsRho[k][my_index] < -NONZERO_LB) {
                     c3lpPtr->cmatval[c3lpPtr->nzcnt] = subprobPtr->rhsRho[k][my_index];
                     //fprintf(stdout, ">>>>>> rhsRho[%d][%d] = %f\n", k, i, subprobPtr->rhsRho[k][my_index]);
        	     c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 +nscens + k;
        	     c3lpPtr->nzcnt++;			// count nonzeros
        	     c3lpPtr->cmatcnt[row_index+n2+nscens+m2]++;  // count nonzeros in this column
                   }
              }
          }
       } // End if(addthis_row)
   } // End outer for loop

  // printf("\n\n*******c3lpPtr->nzcnt = %d\n", c3lpPtr->nzcnt);
 //  printf("*******ncols-2 = %d\n", ncols-2);

   // lambda_02's column: add the disjunction var floor coefs
   c3lpPtr->cmatbeg[ncols-2] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = 1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-2]++; 				// count nonzeros in this column

   //fprintf(stdout, "disj_var_floor = %f \n", i, disj_var_floor);
   // Add the disj_var_floor coef to this column
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_floor > NONZERO_LB) {
            fprintf(stdout, "******** disj_var_floor = %f \n", i, disj_var_floor);
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1*disj_var_floor;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-2]++; 			// count nonzeros in this column
       }
    }

   // lambda_12's column: add the disjunction var ceil coefs
   c3lpPtr->cmatbeg[ncols-1] = c3lpPtr->nzcnt;
   c3lpPtr->cmatval[c3lpPtr->nzcnt] = -1;
   c3lpPtr->cmatind[c3lpPtr->nzcnt] = n2+disj_var_index;
   c3lpPtr->nzcnt++;
   c3lpPtr->cmatcnt[ncols-1]++; 				// count nonzeros in this column

  // fprintf(stdout, "disj_var_ceil = %f \n", i, disj_var_ceil);
   // Add the disj_var_ceil coef to this column
   for (k = 0; k < nscens; k++) {                   // k = scenario
       if (disj_var_ceil > NONZERO_LB) {
            c3lpPtr->cmatval[c3lpPtr->nzcnt] = disj_var_ceil;
            c3lpPtr->cmatind[c3lpPtr->nzcnt] = 2*n2 + nscens + k;
            c3lpPtr->nzcnt++;
            c3lpPtr->cmatcnt[ncols-1]++; 			// count nonzeros in this column
       }
    }

    c3lpPtr->cmatbeg[ncols] = c3lpPtr->nzcnt;

 } //******************************* loadC3LPdataNew ************************//



void
loadRHSLPdata(rhslpProb_t *rhslpPtr, masterproblem_t *masterprobPtr, FILE *fpout)
/**
 * Loads data into RHS lp data structure variables without the t_00 and t_01 columns
 * These columns will be added and deleted later for each scenario as the coefs for
 * these cols become available.
 * @param rhslpPt  pointer to the RHS lp structure
 * @param masterprobPtr  pointer to the master lp data structure
 * @param fpout output file pointer
 */
 {
   int i, j;
   int nrows = rhslpPtr->nrows;
   int ncols = rhslpPtr->ncols;
   int m1    = masterprobPtr->nrows;
   int n1    = masterprobPtr->ncols-1;       // Exclude the optimality variable
   int col_start, row_start1, row_start2;
   int start, stop;

   //fprintf(stderr, "\n**rhslpPtr->nrows = %d**\n", rhslpPtr->nrows);
   //fprintf(stderr, "\n**rhslpPtr->ncols = %d**\n", rhslpPtr->ncols);

   //******* Set object coefs ************//
   // delta column
   rhslpPtr->obj[0] = -1;

   // sigma_0 and sigma_i's
   for (i = 1; i < n1+2; i++)
       rhslpPtr->obj[i] = 1;

   // tau_0's and tau_1's
   for (i = n1+2; i < ncols; i++) // Without t_00 and t_01 cols
       rhslpPtr->obj[i] = 0;

   //*********** Set the ctype array continuous vars ********//
   for (i = 0; i < ncols; i++)    // Without t_00 and t_01 cols
       rhslpPtr->ctype[i] = 'C';

   //*********** Set the sense array to >= ***************//
   for (i = 0; i < nrows; i++) {      // Without t_00 and t_01 cols
       if (i == 2)
           rhslpPtr->sense[i] = 'E';
       else
           rhslpPtr->sense[i] = 'G';
   }
   //******* Set rhs coefs ************//
   for (i = 0; i < nrows; i++){
       if (i == 2)
           rhslpPtr->rhs[i] = 1;
        else
           rhslpPtr->rhs[i] = 0;
   }

   //*********** Set the bounds on each variable ********//
   // delta variable
   rhslpPtr->lb[0] = -1*CPX_INFBOUND;
   rhslpPtr->ub[0] =  CPX_INFBOUND;

   // sigma_0 and sigma_i's
   for (i = 1; i < n1+2; i++) {
       rhslpPtr->lb[i] = -1*CPX_INFBOUND;
       rhslpPtr->ub[i] =  CPX_INFBOUND;
   }

   // tau_0's and tau_1's
   for (i = n1+2; i < ncols; i++) { // Without t_00 and t_01 cols
       rhslpPtr->lb[i] = 0;
       rhslpPtr->ub[i] = CPX_INFBOUND;
   }

   //************** Set the constraint matrix coefs in sparse format ***********//


   // Set the column count for cols delta to sigma_i's
   // and initialize to zero for cols tau_0 to tau_1's.
   for (i = 0; i < 2+n1; i++) {
       rhslpPtr->cmatcnt[i] = 2;
       //fprintf(stderr, "rhslpPtr->cmatcnt[%d] = %d\n", i, rhslpPtr->cmatcnt[i]);
   }
   // Initialize to zero for cols tau_0 to tau_1's.
   for (i = n1+2; i < ncols; i++) {
       rhslpPtr->cmatcnt[i] = 0;   // Initialize to 0
       //fprintf(stderr, "rhslpPtr->cmatcnt[%d] = %d \n", i, rhslpPtr->cmatcnt[i]);
   }


   // delta column
   rhslpPtr->cmatbeg[0] = 0;
   rhslpPtr->cmatval[0] = -1;
   rhslpPtr->cmatval[1] = -1;
   rhslpPtr->cmatind[0] = 3+n1;
   rhslpPtr->cmatind[1] = nrows-1;

   // sigma_0 column
   rhslpPtr->cmatbeg[1] = 2;
   rhslpPtr->cmatval[2] = 1;
   rhslpPtr->cmatval[3] = 1;
   rhslpPtr->cmatind[2] = 0;
   rhslpPtr->cmatind[3] = 1;

   // Count the number of nonzeros (delta and sigma cols)
   rhslpPtr->nzcnt = 4;

   // sigma_i's columns
   row_start1 = 3;
   row_start2 = 4+n1;
   for (i = 2; i < n1+2; i++) {
   	rhslpPtr->cmatbeg[i] = rhslpPtr->nzcnt;
        rhslpPtr->cmatval[rhslpPtr->nzcnt] = 1;
        rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start1;
        rhslpPtr->nzcnt++;
        rhslpPtr->cmatval[rhslpPtr->nzcnt] = 1;
        rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start2;
        rhslpPtr->nzcnt++;
        row_start1++;
        row_start2++;
   }
   // tau_o's columns: add the matrix A coefs and b coefs
   //cols start at c = n1+2 to c < n1+2+m1
   //rols start at r = 3 to r < 4+n1

   // Loop over A nonzeros rows  and extract each row data
   row_start1 = 3;
   for (i = 0; i < masterprobPtr->nrows; i++) {
       if (i == masterprobPtr->nrows-1 ) {
           start = masterprobPtr->rmatbeg_A[i];
           stop  = masterprobPtr->nzcnt_A;
       } else {
           start = masterprobPtr->rmatbeg_A[i];
           stop  = masterprobPtr->rmatbeg_A[i+1];
       }
       rhslpPtr->cmatbeg[i+n1+2] = rhslpPtr->nzcnt;
       for (j = start; j < stop; j++){
            rhslpPtr->cmatval[rhslpPtr->nzcnt] = -1*masterprobPtr->rmatval_A[j];
            rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start1 + masterprobPtr->rmatind_A[j];
            rhslpPtr->nzcnt++;      	// count nonzeros
            rhslpPtr->cmatcnt[i+n1+2]++; 	// count nonzeros in this column
       } // End inner for loop

       // Add the rhs b[i] nonzero value to this column
       if (masterprobPtr->rhs[i] != 0) {
            rhslpPtr->cmatval[rhslpPtr->nzcnt] = masterprobPtr->rhs[i];
            rhslpPtr->cmatind[rhslpPtr->nzcnt] = 3+n1;
            rhslpPtr->nzcnt++;
            rhslpPtr->cmatcnt[i+n1+2]++; 	// count nonzeros in this column
        }
   } // End inner for loop

   // tau_1's columns: add the matrix A coefs and b coefs
   //cols start at c = 2+n1+m1 to c < 2+n1+2m1
   //rols start at r = 4+n1 to r < 5+2n1

   // Loop over A nonzeros rows  and extract each row data
   row_start1 = 4+n1;
   for (i = 0; i < masterprobPtr->nrows; i++) {
       if (i == masterprobPtr->nrows-1 ) {
           start = masterprobPtr->rmatbeg_A[i];
           stop  = masterprobPtr->nzcnt_A;
       } else {
           start = masterprobPtr->rmatbeg_A[i];
           stop  = masterprobPtr->rmatbeg_A[i+1];
       }
       rhslpPtr->cmatbeg[i+2+n1+m1] = rhslpPtr->nzcnt;
       for (j = start; j < stop; j++){
            rhslpPtr->cmatval[rhslpPtr->nzcnt] = -1*masterprobPtr->rmatval_A[j];
            rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start1 + masterprobPtr->rmatind_A[j];
            rhslpPtr->nzcnt++;      	// count nonzeros
            rhslpPtr->cmatcnt[i+2+n1+m1]++; 	// count nonzeros in this column
       } // End inner for loop

       // Add the rhs b[i] nonzero value to this column
       if (masterprobPtr->rhs[i] != 0) {
            rhslpPtr->cmatval[rhslpPtr->nzcnt] = masterprobPtr->rhs[i];
            rhslpPtr->cmatind[rhslpPtr->nzcnt] = 4+2*n1;
            rhslpPtr->nzcnt++;
            rhslpPtr->cmatcnt[i+2+n1+m1]++; 	// count nonzeros in this column
        }
   } // End inner for loop

   // ****** Ends here****

   /***************
   // Loop over A nonzeros and extract each row data
   row_start1 = 3;
   for (i = 0; i < m1; i++) {
       //fprintf(stderr, "\nrhslpPtr->cmatcnt[i+n1+2] = %d \n\n\n", rhslpPtr->cmatcnt[i+n1+2]);
       col_start = 1;
       for (j = 0; j < masterprobPtr->nzcnt_A; j++){
           if(masterprobPtr->rmatind_A[j] == i) {
                rhslpPtr->cmatval[rhslpPtr->nzcnt] = -1*masterprobPtr->rmatval_A[j];
                fprintf(stderr, "col = %d \n val[%d] = %f\n", i, rhslpPtr->nzcnt, rhslpPtr->cmatval[rhslpPtr->nzcnt]);
                rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start1 + masterprobPtr->rmatind_A[j];
                if (col_start == 1) { 		// set column beg
                     rhslpPtr->cmatbeg[i+n1+2] = rhslpPtr->nzcnt;
                     col_start = 0;          	// reset col beg
                }
                rhslpPtr->nzcnt++;      	// count nonzeros
                rhslpPtr->cmatcnt[i+n1+2]++; 	// count nonzeros in this column
           }
        }
        // Add the rhs b[i] nonzero value to this column
        if (masterprobPtr->rhs[i] != 0) {
              rhslpPtr->cmatval[rhslpPtr->nzcnt] = masterprobPtr->rhs[i];
              rhslpPtr->cmatind[rhslpPtr->nzcnt] = 3+n1;
              rhslpPtr->nzcnt++;
              rhslpPtr->cmatcnt[i+n1+2]++; 	// count nonzeros in this column
         }

   } // End outer for loop

   // tau_1's columns: add the matrix A coefs and b coefs
   //cols start at c = 2+n1+m1 to c < 2+n1+2m1
   //rols start at r = 4+n1 to r < 5+2n1

   // Loop over A nonzeros and extract each row data
   row_start1 = 4+n1;
   for (i = 0; i < m1; i++) {
       fprintf(stderr, "\nrhslpPtr->cmatcnt[i+2+n1+m1] = %d \n\n\n",rhslpPtr->cmatcnt[i+2+n1+m1]);
       col_start = 1;
       for (j = 0; j < masterprobPtr->nzcnt_A; j++){
           if(masterprobPtr->rmatind_A[j] == i) {
                rhslpPtr->cmatval[rhslpPtr->nzcnt] = -1*masterprobPtr->rmatval_A[j];
                fprintf(stdout, "A row = %d \n val[%d] = %f\n", i, rhslpPtr->nzcnt, rhslpPtr->cmatval[rhslpPtr->nzcnt]);
                rhslpPtr->cmatind[rhslpPtr->nzcnt] = row_start1 + masterprobPtr->rmatind_A[j];
                if (col_start == 1) {    // set column beg
                     rhslpPtr->cmatbeg[i+2+n1+m1] = rhslpPtr->nzcnt;
                     col_start = 0; // reset col beg
                }
                rhslpPtr->nzcnt++;  // count nonzeros
                rhslpPtr->cmatcnt[i+2+n1+m1]++; 	// count nonzeros in this column
           }
        }
        // Add the rhs b[i] nonzero value to this column
        if (masterprobPtr->rhs[i] != 0) {
              rhslpPtr->cmatval[rhslpPtr->nzcnt] = masterprobPtr->rhs[i];
              rhslpPtr->cmatind[rhslpPtr->nzcnt] = 4+2*n1;
              rhslpPtr->nzcnt++;
              rhslpPtr->cmatcnt[i+2+n1+m1]++; 	// count nonzeros in this column
         }

   } // End outer for loop

   ******/

 } //******************************* loadRHSLPdata ************************//

void
updateConstrMatrixWk(subproblem_t *subprobPtr, double *solnPi, int disjVar, int kIter)
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
 {
   int i, j;
   int rowstart;
   int newrow;

   newrow = subprobPtr->nrows_total;
   subprobPtr->cmatcnt_W[newrow] = 0;

   rowstart = subprobPtr->nzcnt_W;

   for (i = 0; i < subprobPtr->ncols; i++) {
       j = subprobPtr->nzcnt_W;
       if (solnPi[i] != 0) { // store nonzeros only
      	   subprobPtr->cmatbeg_W[newrow] = rowstart;
      	   subprobPtr->cmatval_W[j] = solnPi[i];
      	   subprobPtr->cmatind_W[j] = i;
      	   subprobPtr->nzcnt_W++;
      	   subprobPtr->cmatcnt_W[newrow]++;

   	}
   } // End for loop

   // Update the disj var index array for this pi
   //printf("\n\nMAKING DISJ VAR UPDATE\n\n");
   //printf("\n\nkIter = %d \n\n", kIter);
   subprobPtr->disjVarWindex[kIter] = disjVar;

   //printf("\n\n B: disj_var = %d\n\n", disjVar);
   //    for (i = 0; i < 10; i++)
   //       printf("B: subprobPtr->disjVarWindex[%d] = %d\n", i, subprobPtr->disjVarWindex[i]);


   // Update the  W matrix and lp nrows
   subprobPtr->nrows++;
   subprobPtr->nrows_total++;
   subprobPtr->indices_rowmip[subprobPtr->nrows_mip] = subprobPtr->nrows_mip;
   subprobPtr->nrows_mip++;

 } // ************ End updateConstrMatrixWk *********************//

int
addNewRowToSubProbWmat(CPXENVptr env, CPXLPptr lp_sub, double *solnPi, int len)
/**
 * This function appends the pi vector to the W^k matrix in the k-th iteration
 * Wk is the subproblem constraint matrix in scenario k and is in row sparse matrix format
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP problem object
 * @param solnPi array of pi values
 * @param len length of the solnPi array
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
   int i, nzcnt;
   int status;
	 
   double *rhs = NULL;
   char   *sense = NULL;
   int    *rmatbeg = NULL;
   int    *rmatind = NULL;
   double *rmatval= NULL ;

   // Allocate memory
   rhs     = (double*)malloc(sizeof(double));
   sense   = (char*)malloc(sizeof(char));
   rmatbeg = (int*)malloc(2*sizeof(int));
   rmatind = (int*)malloc(len*sizeof(int));
   rmatval = (double*)malloc(len*sizeof(double));

   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addNewRowToLP() \n");
     fprintf (stderr, " Failed to allocate memory to add row vars\n");
     return 1;
   }

   rhs[0] = 0;
   sense[0] = 'G';
   rmatbeg[0] = 0;

   nzcnt = 0;
   for (i = 0; i < len; i++) {
       if (solnPi[i] != 0) { // store nonzeros only
      	   rmatval[nzcnt] = solnPi[i];
      	   rmatind[nzcnt] = i;
      	   nzcnt++;
   	}
   } // End for loop
   rmatbeg[1] = nzcnt;

   status = CPXaddrows(env, lp_sub, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);

   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, "Function addNewRowToLP() \n");
      fprintf (stderr, " Failed to add new pi row to subproblem lp object\n");
      fprintf (stderr, " Error code = %d\n", status);
   }





   // Free up the memory allocated to the arrays
   if ( rhs != NULL ){
      free (rhs);
      rhs = NULL;
   }
   if ( sense != NULL ){
      free (sense);
      sense = NULL;
   }
   if ( rmatbeg != NULL ){
      free (rmatbeg);
      rmatbeg = NULL;
   }
   if ( rmatind != NULL ){
      free (rmatind);
      rmatind = NULL;
   }

	 // SMH: MemLeakFix (May 18, 2015)
	 if (rmatval != NULL)
	 {
		 free(rmatval);
		 rmatval = NULL;
	 }
	 
	 return status;

 } // ************ End addNewRowToSubProbWmat *********************//


double
scalProd(double *array1, double *array2, int len)

/**
 * This function computes the scalar product of two vectors given by the two arrays
 * @param array1 first array
 * @param array2 second array
 * @param len the length of the two arrays
 * @return returns the scala product of the two arrays
 */
{
  int i;
  double prod = 0;

  for (i = 0; i < len; i++) {
      //printf("solnC3Lambda_11[%d]*stochdataPtr->rhs[%d] = %6.6f * %6.6f\n",
      //        i, i, array1[i], array2[i]);
      prod += array1[i]*array2[i];
  }
  return prod;
} // ****** End scalProd *******//


double
scalProd2(double *array1, double *array2, int start, int len)

/**
 * This function computes the scalar product of two vectors of different sizes
 * given by the two arrays.
 * @param array1 first array
 * @param array2 second array
 * @param start the begining index for array1
 * @param len the length of array2, which starts at index 0
 * @return returns the scala product of the two arrays
 */
{
  int i;
  double prod = 0;

  for (i = 0; i < len; i++) {
      prod += array1[start]*array2[i];
      start++;
  }
  return prod;
} // ****** End scalProd2 *******//


double
getNu01(subproblem_t *subprobPtr, double *lambda_1, double *rhs,
        double lambda_2, double integral_bd, int disj_var)
/**
 * This function computes the nu_0 or nu_1 scalars
 * Must pass the -ve of lambda_02 for computing nu_0
 * @param array1 first array
 * @param array2 second array
 * @param len the length of the two arrays
 * @return returns the scala product of the two arrays
 */
{
  int i;
  int cnt;
  double nu;
  int newrows;

  // Use rhs due to fixed W matrix without explicit binary constraints
  nu = scalProd(lambda_1, rhs, subprobPtr->nrows_T) + lambda_2*integral_bd;

  // Add constribution due to binary constraints
  for (i = subprobPtr->nrows_T; i < subprobPtr->nrowsW; i++)
       nu += lambda_1[i]*-1;

  // Use rhs due to fixed W^k matrix rows for the current disj_var
  // Remember subprobPtr->nrows has been updated!
  newrows = subprobPtr->nrows_total - subprobPtr->nrowsW - 1;

  // count the number of rows add due to current disj var
  cnt = subprobPtr->nrowsW;
 // printf("New rows contribution \n");

  for (i = 0; i < newrows; i++) {
      if (subprobPtr->disjVarWindex[i] == disj_var) {
          nu += lambda_1[cnt]*rhs[subprobPtr->nrows_T+i];
         // printf("solnC3Lambda_11[%d]*stochdataPtr->rhs[%d] = %6.6f * %6.6f\n",
         //         subprobPtr->nrowsW+i, subprobPtr->nrowsW+i,
          //        lambda_1[subprobPtr->nrowsW+i], rhs[subprobPtr->nrowsW+i]);
          cnt++;
      }
  }
  //printf("\n \n");

  return nu;

} // ****** getNu01 *******//


void
getGammaH(double *lambda_h, stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
          double *gamma_h, int scenario, int random_T, int disj_var)
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
 * @param disj_var disjunctive variable
 */
{

  int i, j;
  int curr_row;
  int curr_col;
  int ncols;
  int nrows;
  int start, stop;
  int cnt;

  ncols = stochdataPtr->ncols;

  //count new rows added to T^k due to current disj_var
  // Remember subprobPtr->nrows has been updated!
  //newrows = subprobPtr->nrows - subprobPtr->nrowsW - 1;
  //cnt = 0;
  //for (i = 0; i < newrows; i++) {
  //    if (subprobPtr->disjVarWindex[i] == disj_var) {
  //        cnt++;
  //    }
  //}

  // fprintf(stdout, "Computing Gammas: \n");
   //Initialize result array with zeros
   for (i = 0; i < ncols; i++)
        gamma_h[i] = 0;

   if (random_T == 0) {  //Constant technolody matrix T
        cnt = 0;
       // Compute gamma_h coefs col by col due to the T matrix
       for (i = 0; i < ncols; i++) {
           if (i == ncols-1) {
               start = subprobPtr->cmatbeg_T[i];
               stop  = subprobPtr->nzcnt_T;
              // fprintf(stdout, "1. ncols = %d: start %d: stop %d \n", ncols, start, stop);
           } else {
               start = subprobPtr->cmatbeg_T[i];
               stop  = subprobPtr->cmatbeg_T[i+1];
               //fprintf(stdout, "2. ncols = %d: start %d: stop %d \n", ncols, start, stop);
           }
          // fprintf(stdout, "ncols = %d: start %d: stop %d \n", ncols, start, stop);
           for (j = start; j < stop; j++){
                curr_row = subprobPtr->cmatind_T[j];
                gamma_h[i] += subprobPtr->cmatval_T[j]*lambda_h[curr_row];
                //fprintf(stdout, "\nT: valT[%d] = %f, lambda[%d] = %f \n", j,
                //        subprobPtr->cmatval_T[j], curr_row, lambda_h[curr_row]);

           }
          //fprintf(stdout, "gamma_h[%d] = %f\n", i, gamma_h[i]);

        } // end out for loop

   } else { // Random technology matrix T(w)

        // Compute gamma_h coefs col by col due to the T(w) matrix
        for (i = 0; i < ncols; i++) {
            if (i == ncols-1) {
               start = stochdataPtr->cmatbeg_T[scenario][i];
               stop  = stochdataPtr->cnzcnt_T[scenario];
            } else {
               start = stochdataPtr->cmatbeg_T[scenario][i];
               stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
            }
            for (j = start; j < stop; j++){
               curr_row = stochdataPtr->cmatind_T[scenario][j];
               gamma_h[i] += stochdataPtr->cmatval_T[scenario][j]*lambda_h[curr_row];
            }
        } // end outer for loop

   } // End outer if\else statement

   // Continue computing gamma_h coefs row by row due to the newly added rows
   // to T(w) matrix
   if (subprobPtr->nrows_total-1 > subprobPtr->nrows_T) { // New rows have been added
      // fprintf(stdout, "New %d rows added!!\n", stochdataPtr->rnrows );
       nrows = stochdataPtr->rnrows;
       //cnt = 0;
       cnt = subprobPtr->nrows_T;
       ///////////////////
       //printf("\n\n disj_var = %d\n\n", disj_var);
       //for (i = 0; i < 10; i++)
       //   printf("subprobPtr->disjVarWindex[%d] = %d\n", i, subprobPtr->disjVarWindex[i]);

          //exit(0);
       //////////////////


       for (i = 0; i < nrows; i++) {
          if (subprobPtr->disjVarWindex[i] == disj_var) { // row due to curr disj_var
              if (i == nrows-1) {
                    start = stochdataPtr->rmatbeg_T[scenario][i];
                    stop  = stochdataPtr->rnzcnt_T[scenario];
              } else {
                    start = stochdataPtr->rmatbeg_T[scenario][i];
                    stop  = stochdataPtr->rmatbeg_T[scenario][i+1];
              }
              for (j = start; j < stop; j++){
                    curr_col = stochdataPtr->rmatind_T[scenario][j];
                    if (curr_col > stochdataPtr->ncols){
                        fprintf(stderr, "getGammaH():\n");
                        fprintf(stderr, "curr_col too large! curr_col = %d\n", curr_col);
                        fprintf(stderr, "curr_col too large! curr_col = %d\n", curr_col);
                        exit(1);
                    }
                           curr_row = subprobPtr->nrows_T + i;
                    //curr_row = subprobPtr->nrows_T + cnt;
                    gamma_h[curr_col] += stochdataPtr->rmatval_T[scenario][j]*lambda_h[cnt];
                    //gamma_h[curr_col] += stochdataPtr->rmatval_T[scenario][j]*lambda_h[curr_row];
                    //fprintf(stdout, "\nTr: stochdataPtr->rmatval_T[%d][%d] = %f, lambda[%d] = %f \n", scenario,
                    //    j, stochdataPtr->rmatval_T[scenario][j], curr_row, lambda_h[curr_row]);
                    //fprintf(stdout, "gamma_h[%d] = %f\n", curr_col, gamma_h[curr_col]);
                    //cnt++;
             }
             cnt++;
         } //End outer if

     } // end outer for loop
  } // End if statement



} // ******************** getGammaH ************************//



int
addNewColsToRHSlp(CPXENVptr env, CPXLPptr lp_rhs, double *gamma_0, double gamma_1[],
                  double nu_0, double nu_1, int len)
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
 {
   int i, nzcnt;
   int status;
   int row;
   double    *obj     = NULL;
   int    *cmatbeg = NULL;
   int    *cmatind = NULL;
   double *cmatval = NULL ;
   char **colnames = NULL;

   // Allocate memory
   obj = (double*)malloc(2*sizeof(double));
   cmatbeg = (int*)malloc(2*sizeof(int));
   cmatind = (int*)malloc((2*len+5)*sizeof(int));
   cmatval = (double*)malloc((2*len+5)*sizeof(double));
   colnames = (char**)malloc(2*sizeof(char*));

   if (cmatbeg == NULL || cmatind == NULL || cmatval == NULL || colnames == NULL) {
       fprintf (stderr, "\nd2algFuncs: \n");
       fprintf (stderr, "Function addNewColsToRHSlp \n");
       fprintf (stderr, " Failed to allocate memory to add col variables\n");
       return 1;
   }

   colnames[0] = (char*)malloc(NAMELEN*sizeof(char));
   colnames[1] = (char*)malloc(NAMELEN*sizeof(char));

   strcpy(colnames[0], "tau00");
   strcpy(colnames[1], "tau01");

   //fprintf (stderr, "colnames[0] = %s\n", colnames[0]);
  // fprintf (stderr, "colnames[1] = %s\n", colnames[1]);

   cmatbeg[0] = 0;
   nzcnt = 0;

   // Add obj coefs
   obj[0] = 0;
   obj[1] = 0;

   // *** Add nonzeros for the tau_00 column ***//
    // Add the -1 coef
    cmatval[nzcnt] = -1;
    cmatind[nzcnt] = 0;
    nzcnt++;

    // Add the 1 coef
    cmatval[nzcnt] = 1;
    cmatind[nzcnt] = 2;
    nzcnt++;


    //fprintf(stderr, "len = %d\n", len);
    // Add the gamma_0 coefs
    row = 3;
    for (i = 0; i < len; i++) {
       if (gamma_0[i] != 0) { // store nonzeros only
      	   cmatval[nzcnt] = -1*gamma_0[i];
      	   cmatind[nzcnt] = row;
      	   nzcnt++;
   	}
   	row++;
    } // End for loop
    // Add the nu_0 coef if nonzero
    if (nu_0 != 0) {
       cmatval[nzcnt] = nu_0;
       cmatind[nzcnt] = row;
       nzcnt++;
    }

   // Add nonzeros for the tau_01 column
    cmatbeg[1] = nzcnt;

    // Add the -1 coef
    cmatval[nzcnt] = -1;
    cmatind[nzcnt] = 1;
    nzcnt++;

    #ifdef DEBUB_FUNC
      fprintf(stderr, "*****WOW 2 gamma_1**: ");
  	           fprintf(stderr, "[");
  	           for (i = 0; i < len; i++)
  	                fprintf(stderr, "%f ", gamma_1[i]);
  	           fprintf(stderr, "]\n");
    #endif

    // Add the 1 coef
    cmatval[nzcnt] = 1;
    cmatind[nzcnt] = 2;
    nzcnt++;
 #ifdef DEBUB_FUNC
    fprintf(stderr, "*****WOW 3 gamma_1**: ");
  	           fprintf(stderr, "[");
  	           for (i = 0; i < len; i++)
  	                fprintf(stderr, "%f ", gamma_1[i]);
  	           fprintf(stderr, "]\n");

    fprintf(stderr, "[");
    for (i = 0; i < len; i++){
  	fprintf(stderr, "%f ", gamma_1[i]);
    }
    fprintf(stderr, "]\n");

  #endif

    // Add the gamma_1 coefs
    row = 4 + len;
    for (i = 0; i < len; i++) {
       if (gamma_1[i] != 0) { // store nonzeros only
      	   cmatval[nzcnt] = -1*gamma_1[i];
      	   cmatind[nzcnt] = row;
      	   nzcnt++;
      	   //fprintf(stderr, "row = %d: gamma_1[%d] = %f \n", row, i, gamma_1[i]);
   	}
   	row++;
    }
    //
    // Add the nu_1 coef if nonzero
    if (nu_1 != 0) {
       cmatval[nzcnt] = nu_1;
       cmatind[nzcnt] = row;
       nzcnt++;
    }

    //
    //fprintf(stderr, "Adding cols to RHS lp\n");

   // Add the cols to the RHS lp now
   status = CPXaddcols(env, lp_rhs, 2, nzcnt, obj, cmatbeg, cmatind, cmatval, NULL,
                       NULL, colnames);
   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, " Function addNewColsToRHSlp \n");
      fprintf (stderr, " Failed to add new tau_0h cols to RHS lp object\n");
      fprintf (stderr, " Error code = %d\n", status);
   }

   //fprintf(stderr, "Done adding cols to RHS lp\n");

   // Free up the memory allocated to the arrays
   if ( cmatbeg != NULL ){
      free (cmatbeg);
      cmatbeg = NULL;
   }
   if ( cmatind != NULL ){
      free (cmatind);
      cmatind = NULL;
   }

   return status;

 } // ************ End addNewColsToRHSlp() *********************//


void
updateRHSrAndMatT(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,
                  double nu, double *gamma, int scenario, double *solnX)
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
{
   int i;
   int curr_row = stochdataPtr->nrows + stochdataPtr->rnrows; // Add number of original/initial
   int temp;
                                                      // rows in subproblem plus the newly added

    //printf("... stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);

   //fprintf(stdout, "curr_row = %d \n", curr_row);
   // Append delta/sigma_0 coef to r(w)
   stochdataPtr->rhs[scenario][curr_row] = nu;
   //fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", scenario, curr_row, stochdataPtr->rhs[scenario][curr_row]);

   curr_row = stochdataPtr->rnrows; // Newly appended rows, initially zero

   //fprintf(stdout, "updateRHSrAndMatT(): \n");
   //fprintf(stdout, "rnzcnt_T[%d] = %d \n", scenario, stochdataPtr->rnzcnt_T[scenario] );
   // Append new row to T(w): new row = [sigma_i]/sigma_0
   temp = stochdataPtr->rnzcnt_T[scenario];
   stochdataPtr->rmatbeg_T[scenario][curr_row] = temp;
   //stochdataPtr->rmatbeg_T[scenario][curr_row] = stochdataPtr->rnzcnt_T[scenario];
   //stochdataPtr->rmatcnt_T[scenario][curr_row] = 0;

   //printf(".... stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);

   //fprintf(stdout, "** rnzcnt_T[%d] = %d \n",  scenario, stochdataPtr->rnzcnt_T[scenario] );
   //fprintf(stdout, "** rmatbeg_T[%d] = %d \n", scenario, stochdataPtr->rmatbeg_T[scenario][curr_row]);
   //fprintf(stdout, "** stochdataPtr->ncols = %d\n", stochdataPtr->ncols);
   for (i = 0; i < stochdataPtr->ncols; i++){
       //fprintf(stdout, "\n i = %d\n", i);
       ////////////
       //if (scenario == 81){
       //          printf("* stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);
       //}
       ///////////

       //if (sigma_i[i]/sigma_0 > NONZERO_LB || sigma_i[i]/sigma_0 < -NONZERO_LB) {
      //if (sigma_i[i] != 0) {
         //fprintf(stdout, "** sigma_i[%d] = %f\n",  i, sigma_i[i]);
         //fprintf(stdout, "** solnX[%d] = %f\n", i, solnX[i]);
               temp = stochdataPtr->rnzcnt_T[scenario];
                 //stochdataPtr->rmatval_T[scenario][stochdataPtr->rnzcnt_T[scenario]] = sigma_i[i]/sigma_0;
               stochdataPtr->rmatval_T[scenario][temp] = gamma[i];

               //printf("** stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);

                 //fprintf(stdout, "\n\nstochdataPtr->rmatval_T[%d][%d] = %f\n\n",
                 //fprintf(stdout, "\n\nstochdataPtr->rmatval_T[%d][%d] = %f\n\n",
                 //scenario,  stochdataPtr->rnzcnt_T[scenario],
                 //stochdataPtr->rmatval_T[scenario][stochdataPtr->rnzcnt_T[scenario]]);

                 //stochdataPtr->rmatind_T[scenario][stochdataPtr->rnzcnt_T[scenario]] = i;
               stochdataPtr->rmatind_T[scenario][temp] = i;

               //printf("*** stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);

               //fprintf(stdout, "> stochdataPtr->rmatind_T[%d][%d] = %d\n",
               //         scenario, stochdataPtr->rnzcnt_T[scenario], stochdataPtr->rmatind_T[scenario][stochdataPtr->rnzcnt_T[scenario]]);
                 //stochdataPtr->rmatcnt_T[scenario][curr_row]++;
                 //stochdataPtr->rnzcnt_T[scenario]++;
             //stochdataPtr->rmatcnt_T[scenario][curr_row] += 1;
               //printf("**** stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);
               stochdataPtr->rnzcnt_T[scenario] += 1;
               //printf("***** stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);
        //}
   } // End for loop

   // Compute rhs rho(w) for this scenario's newly appended row
   curr_row = stochdataPtr->nrows + stochdataPtr->rnrows; // newly appended row index
     //printf("stochdataPtr->nrows = %d \n", stochdataPtr->nrows);
     //printf("stochdataPtr->rnrows = %d \n", stochdataPtr->rnrows);
   subprobPtr->rhsRho[scenario][curr_row] = nu;
      //fprintf(stdout, "1. subprobPtr->rhsRho[%d][%d] = %f\n", scenario, curr_row, subprobPtr->rhsRho[scenario][curr_row]);
   for (i = 0; i < stochdataPtr->ncols; i++){
        //printf("i = %d ", i);
        subprobPtr->rhsRho[scenario][curr_row] -= gamma[i]*solnX[i];
        //fprintf(stdout, "subprobPtr->rhsRho[%d][%d] = %f\n", scenario, curr_row, subprobPtr->rhsRho[scenario][curr_row]);
   }
   // Update row indices array for rhs update for the scenario lp/mip
   subprobPtr->indices_row[curr_row] = curr_row;



  // fprintf(stdout, "Done updating T and RHS\n");

} //*************** End updateRHSrAndMatT **********************************//


int
addNewRowToMaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr,
                  double *solnX, double expObjval, double SUBPROB_LB)
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
 {
   int j, nzcnt;
   int status;
   double *rhs = NULL;
   char   *sense = NULL;
   int    *rmatbeg = NULL;
   int    *rmatind = NULL;
   double *rmatval= NULL ;
   int modeS_k;      // Number of nonzeros in the master prob solution
   double q_s;	     // subprob expected obj - lower bound

   int ncols = masterprobPtr->ncols;

   // Allocate memory
   rhs     = (double*)malloc(sizeof(double));
   sense   = (char*)malloc(sizeof(char));
   rmatbeg = (int*)malloc(sizeof(int));
   rmatind = (int*)malloc(ncols*sizeof(int));
   rmatval = (double*)malloc(ncols*sizeof(double));

   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addNewRowToLP() \n");
     fprintf (stderr, " Failed to allocate memory to add row vars\n");
     return 1;
   }

   //*********** DERIVE OPTIMALITY CUT **************//

   q_s = expObjval;

   sense[0] = 'G';
   rmatbeg[0] = 0;

   nzcnt = 0;

   // compute |S_k| modeS_k: total number of x soln components equal to zero
   modeS_k = 0;
   for (j = 0; j < ncols-1; j++) {
       if (solnX[j] >= 1 - INT_PRECISION) {
           modeS_k += 1;               // count nonzeros
           rmatval[nzcnt] = -q_s; // store nonzeros only
           rmatind[nzcnt] = j;
           nzcnt++;
       } else {
           rmatval[nzcnt] = q_s; // store nonzeros only
           rmatind[nzcnt] = j;
           nzcnt++;
       }
   }
   // Add the theta or opt column coef
   rmatval[nzcnt] = 1; // store nonzeros only
   rmatind[nzcnt] = ncols-1;
   nzcnt++;

  // fprintf(stderr, "L L - C U T:\n");
   //fprintf(stderr, "modeS_k = %d \n", modeS_k);
   //for (j = 0; j < ncols; j++) {
     //  fprintf(stderr, "rmatval[%d] = %f	", j, rmatval[j]);
     //  fprintf(stderr, "rmatind[%d] = %d\n", j, rmatind[j]);
   //}

   // Set the rhs
   rhs[0] = -q_s*(modeS_k - 1);

   status = CPXaddrows(env, lp_master, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);

   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, "Function addNewRowToLP() \n");
      fprintf (stderr, " Failed to add new pi row to subproblem lp object\n");
      fprintf (stderr, " Error code = %d\n", status);
   }

   // Free up the memory allocated to the arrays
   if ( rhs != NULL ){
      free (rhs);
      rhs = NULL;
   }
   if ( sense != NULL ){
      free (sense);
      sense = NULL;
   }
   if ( rmatbeg != NULL ){
      free (rmatbeg);
      rmatbeg = NULL;
   }
   if ( rmatind != NULL ){
      free (rmatind);
      rmatind = NULL;
   }

   return status;

 } // ************ addNewRowToMaster *********************//


int
addBendersCutToMaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr, double SUBPROB_LB)
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
 {
   //fprintf(stdout, "\n\n\n 1. addBendersCutToMaster: \n");
   int j, nzcnt;
   int status;
   double *rhs;
   char   *sense;
   int    *rmatbeg;
   int    *rmatind;
   double *rmatval;
   int ncols;

   ncols = masterprobPtr->ncols;

   //fprintf(stdout, "\n\n\n 2. addBendersCutToMaster: \n");
   //fprintf(stdout, "\n masterprobPtr->ncols: %d \n", ncols);

   // Allocate memory
   rhs     = (double*)malloc(sizeof(double));
   //fprintf(stdout, "\n\n\n 2a. addBendersCutToMaster: \n");
   sense   = (char*)malloc(sizeof(char));
   //fprintf(stdout, "\n\n\n 2b. addBendersCutToMaster: \n");
   rmatbeg = (int*)malloc(sizeof(int));
   //fprintf(stdout, "\n\n\n 2c. addBendersCutToMaster: \n");
   rmatind = (int*)malloc((ncols+1)*sizeof(int));
   //fprintf(stdout, "\n\n\n 2d. addBendersCutToMaster: \n");
   rmatval = (double*)malloc((ncols+1)*sizeof(double));
   //fprintf(stdout, "\n\n\n 2e. addBendersCutToMaster: \n");
   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addNewRowToLP() \n");
     fprintf (stderr, " Failed to allocate memory to add row vars\n");
     return 1;
   }
   //fprintf(stdout, "\n\n\n 3. addBendersCutToMaster: \n");

   //*********** DERIVE OPTIMALITY CUT **************//
   sense[0] = 'G';
   rmatbeg[0] = 0;

   nzcnt = 0;

   //fprintf(stdout, "\n\n\n 4. addBendersCutToMaster: \n");

   // store the x coefs
   for (j = 0; j < ncols-1; j++) {
       if (masterprobPtr->cutCoefs[j] > NONZERO_LB || masterprobPtr->cutCoefs[j] < -NONZERO_LB) {
           rmatval[nzcnt] = masterprobPtr->cutCoefs[j];   // store nonzeros only
           rmatind[nzcnt] = j;
           nzcnt++;
       }
   }
   // Add the theta or opt column coef
   rmatval[nzcnt] = 1; // store nonzeros only
   rmatind[nzcnt] = ncols-1;
   nzcnt++;

   //for (j = 0; j < ncols; j++) {
   //    fprintf(stderr, "rmatval[%d] = %f	", j, rmatval[j]);
   //    fprintf(stderr, "rmatind[%d] = %d\n", j, rmatind[j]);
   //}

   // Set the rhs: translate it
   rhs[0] = masterprobPtr->rhsCoef+SUBPROB_LB;

   //for (j = 0; j < ncols-1; j++) {
   //   fprintf(stdout, "masterprobPtr->cutCoefs[%d] = %f	\n", j, masterprobPtr->cutCoefs[j]);
   //
   //}
  //fprintf(stdout, "masterprobPtr->rhsCoef = %f\n", masterprobPtr->rhsCoef);
 // fprintf(stdout, "rhsCoef = %f\n", rhs[0]);

   status = CPXaddrows(env, lp_master, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);

   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, "Function addNewRowToLP() \n");
      fprintf (stderr, " Failed to add new pi row to subproblem lp object\n");
      fprintf (stderr, " Error code = %d\n", status);
   }

   // Free up the memory allocated to the arrays
   if ( rhs != NULL ){
      free (rhs);
      rhs = NULL;
   }
   if ( sense != NULL ){
      free (sense);
      sense = NULL;
   }
   if ( rmatbeg != NULL ){
      free (rmatbeg);
      rmatbeg = NULL;
   }
   if ( rmatind != NULL ){
      free (rmatind);
      rmatind = NULL;
   }
   if ( rmatval != NULL ){
      free (rmatval);
      rmatind = NULL;
   }

   //exit(1);

   return status;

 } // ************ addNewRowToMaster *********************//

 int
dropScenSolnFromC3obj(solnMatrix_t *solnPtr, double *c3objcoefs,
                      double* scenProb, int disj_scen, double * c3secondobj)
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
 * @return status: zero if scenario is dropped, otherwise returns
 * a nonzero if all scenarios have been dropped
 */
 {
   int i, j;
   int start, stop;
   int col, row;
   int ncols;
   double sum;


   ncols = solnPtr->ncols+1+solnPtr->nrows;

  // fprintf(stderr, "\nCurrent C^3 Objective Coefs: \n");
  // for (i = 0; i < ncols+1; i++) {
  //        fprintf(stderr," coef[%d] = %6.6g\n", i, c3objcoefs[i]);
  // }
   //fprintf(stderr,"\n");

   if (solnPtr->ncondScens == 1){ // No more scenarios to drop
       return 1;
   }
   if (solnPtr->scenToDrop == disj_scen){ // don't drop off disjunctive scenario
       solnPtr->scenToDrop++;
   }

   // Compute re-weighted conditional scenario probabilities
   sum = 0;
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++) {
       sum += solnPtr->condScenProbs[i];
   }
   if (solnPtr->disj_scen_index < solnPtr->scenToDrop){
      j = solnPtr->disj_scen_index;
      sum += solnPtr->condScenProbs[j];
   }
  // fprintf(stderr,"sum = %f\n", sum);
   // Compute the conditional scenario probabilities
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++) {
       solnPtr->condScenProbs[i] = solnPtr->condScenProbs[i]/sum;
       //fprintf(stderr," solnPtr->condScenProbs[%d] = %6.6g\n", i, solnPtr->condScenProbs[i]);
   }
   if (solnPtr->disj_scen_index < solnPtr->scenToDrop){
      j = solnPtr->disj_scen_index;
      solnPtr->condScenProbs[j] = solnPtr->condScenProbs[j]/sum;
      //fprintf(stderr," solnPtr->condScenProbs[%d] = %6.6g\n", j, solnPtr->condScenProbs[j]);
   }

   // Re-set the pi obj coefs to zeros
   ncols = solnPtr->ncols+1+solnPtr->nrows;
   for (j = 0; j < solnPtr->ncols + 1; j++)
       c3objcoefs[j] = 0;

   //Set the pi obj coefs
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++){
       row = solnPtr->condScens[i];
       if (row == solnPtr->nrows-1 ) {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->nzcnt_S;
       } else {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->rmatbeg_S[row+1];
       }
       for (j = start; j < stop; j++) {
           col = solnPtr->rmatind_S[j];
           c3objcoefs[col] += solnPtr->condScenProbs[i]*solnPtr->rmatval_S[j];
       }
   } // end outer for loop

   if (solnPtr->disj_scen_index < solnPtr->scenToDrop){
       row = solnPtr->condScens[solnPtr->disj_scen_index];
       if (row == solnPtr->nrows-1 ) {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->nzcnt_S;
       } else {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->rmatbeg_S[row+1];
       }
       for (j = start; j < stop; j++) {
           col = solnPtr->rmatind_S[j];
           c3objcoefs[col] += solnPtr->condScenProbs[solnPtr->disj_scen_index]*solnPtr->rmatval_S[j];
       }
    }

   // Set the pi_0 obj coefs as new conditional probs
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++){
	   j = solnPtr->condScens[i];
        c3secondobj[j] = solnPtr->condScenProbs[i];
   }
   if (solnPtr->disj_scen_index < solnPtr->scenToDrop){
       j = solnPtr->condScens[solnPtr->disj_scen_index];
       c3secondobj[j] = solnPtr->condScenProbs[solnPtr->disj_scen_index];
   }


   // Set the pi_0 obj coefs as new conditional probs
   //start = solnPtr->ncols + scenToDrop+1;
   //c3objcoefs[start] = 0; // Set the dropped scenario coef = 0
   //start += 1;
   //scenario = scenToDrop+1;;
   //for (i = start; i < ncols+1; i++) {
   //    c3objcoefs[i] = -scenCondProb[scenario];
   //    //c3objcoefs[i] = -scenProb[scenario];
   //    scenario++;
   //}

  // fprintf(stderr, "\nC^3 Objective Coefs After Dropping off Scenario %d soln coefs: \n",
  //         solnPtr->scenToDrop);
   ncols = solnPtr->ncols+1+solnPtr->nrows;
   //for (i = 0; i < ncols; i++) {
   //       fprintf(stderr," coef[%d] = %6.6g\n", i, c3objcoefs[i]);
   //}

   // Updates
   solnPtr->ncondScens--; // reduce the number of scenarios in computing the cond. expctn
   solnPtr->scenToDrop++; // Scenario to drop next

   return 0;

 } // ***************dropScenSolnFromC3obj ***********************//



/****
int
dropScenSolnFromC3obj(solnMatrix_t *solnPtr, double *c3objcoefs,
                      double* scenProb, int disj_scen)
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
 * @return status: zero if scenario is dropped, otherwise returns
 * a nonzero if all scenarios have been dropped

 {
   int i, j;
   int start, stop;
   int col, row;
   int ncols;
   double sum;


   ncols = solnPtr->ncols+1+solnPtr->nrows;

  // fprintf(stderr, "\nCurrent C^3 Objective Coefs: \n");
  // for (i = 0; i < ncols+1; i++) {
  //        fprintf(stderr," coef[%d] = %6.6g\n", i, c3objcoefs[i]);
  // }
   //fprintf(stderr,"\n");

   if (solnPtr->ncondScens == 1){ // No more scenarios to drop
       return 1;
   }
   if (solnPtr->scenToDrop == disj_scen){ // don't drop off disjunctive scenario
       solnPtr->scenToDrop++;
   }

   // Compute re-weighted conditional scenario probabilities
   sum = 0;
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++) {
       sum += solnPtr->condScenProbs[i];
   }
   if (disj_scen < solnPtr->scenToDrop){
      j = solnPtr->disj_scen_index;
      sum += solnPtr->condScenProbs[j];
   }
  // fprintf(stderr,"sum = %f\n", sum);
   // Compute the conditional scenario probabilities
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++) {
       solnPtr->condScenProbs[i] = solnPtr->condScenProbs[i]/sum;
       //fprintf(stderr," solnPtr->condScenProbs[%d] = %6.6g\n", i, solnPtr->condScenProbs[i]);
   }
   if (disj_scen < solnPtr->scenToDrop){
      j = solnPtr->disj_scen_index;
      solnPtr->condScenProbs[j] = solnPtr->condScenProbs[j]/sum;
      //fprintf(stderr," solnPtr->condScenProbs[%d] = %6.6g\n", j, solnPtr->condScenProbs[j]);
   }

   // Re-set the pi obj coefs to zeros
   ncols = solnPtr->ncols+1+solnPtr->nrows;
   for (j = 0; j < ncols; j++)
       c3objcoefs[j] = 0;

   //Set the pi obj coefs
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++){
       row = solnPtr->condScens[i];
       if (row == solnPtr->nrows-1 ) {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->nzcnt_S;
       } else {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->rmatbeg_S[row+1];
       }
       for (j = start; j < stop; j++) {
           col = solnPtr->rmatind_S[j];
           c3objcoefs[col] += solnPtr->condScenProbs[i]*solnPtr->rmatval_S[j];
       }
   } // end outer for loop

   if (disj_scen < solnPtr->scenToDrop){
       row = solnPtr->condScens[solnPtr->disj_scen_index];
       if (row == solnPtr->nrows-1 ) {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->nzcnt_S;
       } else {
           start = solnPtr->rmatbeg_S[row];
           stop  = solnPtr->rmatbeg_S[row+1];
       }
       for (j = start; j < stop; j++) {
           col = solnPtr->rmatind_S[j];
           c3objcoefs[col] += solnPtr->condScenProbs[solnPtr->disj_scen_index]*solnPtr->rmatval_S[j];
       }
    }

   // Set the pi_0 obj coefs as new conditional probs
   for (i = solnPtr->scenToDrop+1; i < solnPtr->ncondScens; i++){
        j = solnPtr->ncols + 1 + solnPtr->condScens[i];
        c3objcoefs[j] = -solnPtr->condScenProbs[i];
   }
   if (disj_scen < solnPtr->scenToDrop){
       j = solnPtr->ncols + 1 + solnPtr->condScens[solnPtr->disj_scen_index];
       c3objcoefs[j] = -solnPtr->condScenProbs[solnPtr->disj_scen_index];
   }


   // Set the pi_0 obj coefs as new conditional probs
   //start = solnPtr->ncols + scenToDrop+1;
   //c3objcoefs[start] = 0; // Set the dropped scenario coef = 0
   //start += 1;
   //scenario = scenToDrop+1;;
   //for (i = start; i < ncols+1; i++) {
   //    c3objcoefs[i] = -scenCondProb[scenario];
   //    //c3objcoefs[i] = -scenProb[scenario];
   //    scenario++;
   //}

  // fprintf(stderr, "\nC^3 Objective Coefs After Dropping off Scenario %d soln coefs: \n",
  //         solnPtr->scenToDrop);
   ncols = solnPtr->ncols+1+solnPtr->nrows;
   //for (i = 0; i < ncols; i++) {
   //       fprintf(stderr," coef[%d] = %6.6g\n", i, c3objcoefs[i]);
   //}

   // Updates
   solnPtr->ncondScens--; // reduce the number of scenarios in computing the cond. expctn
   solnPtr->scenToDrop++; // Scenario to drop next

   return 0;

 }
 / ***************dropScenSolnFromC3obj ***********************

***/

int
freeC3lpModelAndData(CPXENVptr env, CPXLPptr lp_c3, c3lpProb_t *c3lpPtr)
/**
 * This frees up the C^3 lp CPLEX object and data arrays
 * @param env  a pointer to the CPLEX environment object
 * @param lp_c3  a pointer to the CPLEX LP problem object
 * @param c3lpPtr  a pointer to the C^3 data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
   int status;

   /* Free up the problem as allocated by CPXcreateprob*/
   if ( lp_c3 != NULL ) {
       status = CPXfreeprob (env, &lp_c3);
       if ( status ) {
            fprintf (stderr, "CPXfreeprob lp_c3 failed, error code %d.\n", status);
       }
    }

   // Free up the memory allocated to the arrays
   if ( c3lpPtr->indices != NULL ){
      free (c3lpPtr->indices);
      c3lpPtr->indices = NULL;
   }

   if ( c3lpPtr->obj != NULL ){
      free (c3lpPtr->obj);
      c3lpPtr->obj = NULL;
   }

   if ( c3lpPtr->ctype != NULL ){
      free (c3lpPtr->ctype);
      c3lpPtr->ctype = NULL;
   }

   if ( c3lpPtr->sense != NULL ){
      free (c3lpPtr->sense);
      c3lpPtr->sense = NULL;
   }

   if ( c3lpPtr->rhs != NULL ){
      free (c3lpPtr->rhs);
      c3lpPtr->rhs = NULL;
   }

   if ( c3lpPtr->lb != NULL ){
      free (c3lpPtr->lb);
      c3lpPtr->lb = NULL;
   }

   if ( c3lpPtr->ub != NULL ){
      free (c3lpPtr->ub);
      c3lpPtr->ub = NULL;
   }

   if ( c3lpPtr->cmatbeg != NULL ){
      free (c3lpPtr->cmatbeg);
      c3lpPtr->cmatbeg = NULL;
   }

   if ( c3lpPtr->cmatcnt != NULL ){
      free (c3lpPtr->cmatcnt);
      c3lpPtr->cmatcnt = NULL;
   }

   if ( c3lpPtr->cmatind != NULL ){
      free (c3lpPtr->cmatind );
      c3lpPtr->cmatind  = NULL;
   }

   if ( c3lpPtr->cmatval != NULL ){
      free (c3lpPtr->cmatval );
      c3lpPtr->cmatval  = NULL;
   }

   return status;

 } // ************ End freeC3lpModelAndData *********************//


void
computeBendersCutCoefs(stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                       int curr_nrows, double *duals, double *cutCoefs, double *rhsCoef,
                       int random_T, double total_redcosts, FILE *fpout)
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
 {
   int i, j;
   int curr_row;
   int curr_col;
   int nrows;
   int ncols = stochdataPtr->ncols;
   int start, stop;
   double *coefs = NULL;
   double curr_rhs = 0;
   double prod;

   // Allocate memory
   coefs = (double*)malloc(ncols*sizeof(double));

   //Initialization
   if (scenario == 0){
      *rhsCoef = 0;
      for (i = 0; i < ncols; i++)
           cutCoefs[i] = 0;
   }
   for (i = 0; i < ncols; i++)
       coefs[i] = 0;

   #ifdef DEUBG_FUNC
      printf("\n\n\n\n Computing Benders cut for SCENARIO %d: \n\n", scenario);
      for (i = 0; i < curr_nrows; i++){
          fprintf(stdout, "duals[%d] = %f * stochdataPtr->rhs[%d][%d] = %f \n",
                       i, duals[i], scenario, i, stochdataPtr->rhs[scenario][i]);
      }
    #endif

   // Accumulate rhs coef
  // *rhsCoef += scalProd(duals, stochdataPtr->rhs[scenario], curr_nrows)*stochdataPtr->scenProb[scenario];
   prod = scalProd(duals, stochdataPtr->rhs[scenario], curr_nrows);


   //fprintf(stdout, "prod = %f \n", prod);

   //fprintf(stdout, "total_redcosts = %f \n",total_redcosts);

   prod += total_redcosts;

   //fprintf(stdout, "prod+total_redcosts = %f \n",prod);

   *rhsCoef += prod*stochdataPtr->scenProb[scenario];

   //fprintf(stdout, "prod = %f \n",prod*stochdataPtr->scenProb[scenario]);

   //fprintf(stdout, "*rhsCoef += %f \n",*rhsCoef);

   //curr_rhs = (scalProd(duals, stochdataPtr->rhs[scenario], curr_nrows) + total_redcosts)*stochdataPtr->scenProb[scenario];

   //*rhsCoef += curr_rhs;

   if (random_T == 0) {  // Constant technology matrix T
      // Compute rho(w) col by col : T matrix is column by column sparse matrix
      #ifdef DEBUG_FUNC
            for (i = 0; i < ncols; i++)
                 printf("subprobPtr->cmatbeg_T[%d] = %d  \n", i, subprobPtr->cmatbeg_T[i]);
      #endif

      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->nzcnt_T;
          } else {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->cmatbeg_T[i+1];
          }
          //fprintf(stdout, "\n\n");
          for (j = start; j < stop; j++){
             //fprintf(stdout, "start = %d stop = %d\n", start, stop);
             curr_row = subprobPtr->cmatind_T[j];
             coefs[j] += duals[curr_row]*subprobPtr->cmatval_T[j];
             #ifdef DEBUG_FUNC
               fprintf(stdout, "coefs[%d] = %f \n", j,coefs[j]);
               fprintf(stdout, "col = %d	cmatval_T[%d] = %f: \n", i, j,
                    subprobPtr->cmatval_T[j]);

               fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %6.6f\n", scenario, curr_row,
                              stochdataPtr->rhs[scenario][curr_row]);
             #endif
          }
      }

   //exit(0);

   } else { //Random Technology matrix T(w)

      // Compute rho(w) col by col : T(w) matrix is column by column sparse matrix
      // This is the T(w) read from the stoch file
      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
          } else {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cnzcnt_T[scenario];
          }
          for (j = start; j < stop; j++){
             curr_row = stochdataPtr->cmatind_T[scenario][j];
             coefs[j] += duals[curr_row]*stochdataPtr->cmatval_T[scenario][j];
             fprintf(stdout, "col = %d	val[%d] = %f: \n", i, j, stochdataPtr->cmatval_T[scenario][j]);
          }
      }

   } // End if/else statement


   // Compute rho(w) row by row for the newly added rows
   // The T(w) matrix for the newly added rows is in row-by-row sparse matrix format
   //fprintf(stdout, "\n\nHHHHHHHHHHHHHHHHH nrows > nrows_T: %d > %d \n\n\n", subprobPtr->nrows, subprobPtr->nrows_T);
   if (subprobPtr->nrows > subprobPtr->nrows_T) { // New rows have been added
         //fprintf(stdout, "\n Number of rows added = %d \n", stochdataPtr->rnrows );
         // Continue computing rho(w) row by row due to the T(w) matrix, which is
         // in sparse format
         nrows = stochdataPtr->rnrows;
         for (i = 0; i < nrows; i++) {
             if (i == nrows-1) {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rnzcnt_T[scenario];
              } else {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rmatbeg_T[scenario][i+1];
              }

              curr_row = subprobPtr->nrows_T + i;

              for (j = start; j < stop; j++){
                  curr_col = stochdataPtr->rmatind_T[scenario][j];

                  #ifdef DEBUG_FUNC
                       fprintf(stdout, "duals[%d] = %f * stochdataPtr->rmatval_T[%d][%d] = %f \n",
                          curr_row, duals[curr_row], scenario, j,
                          stochdataPtr->rmatval_T[scenario][j]);
                  #endif

                  coefs[curr_col] += duals[curr_row]*stochdataPtr->rmatval_T[scenario][j];

                  #ifdef DEBUG_FUNC
                      fprintf(stdout, "coefs[%d] = %f \n", curr_col,coefs[curr_col]);
                      fprintf(stdout, "\n stochdataPtr->cmatval_T[%d][%d] = %d \n", scenario, j,
                              stochdataPtr->cmatval_T[scenario][j]);
                  #endif

              }
         } // end outer for loop

    } //End if statement

    // Set the weighted coefs
    for (i = 0; i < ncols; i++){
        cutCoefs[i] += coefs[i]*stochdataPtr->scenProb[scenario];
    }

     //****** Print the cut coefs ******/
     #ifdef DEBUG_FUNC
        fprintf(stdout, "computeBendersCutCoefs (): \n");
        fprintf(stdout, "\nCut coefs at scenario %d \n", scenario);
        fprintf(stdout, "rhsCoef = %f \n", *rhsCoef);
  	for (i = 0; i < ncols; i++)
  	    fprintf(stdout, "coefs[%d] = %f \n", i, coefs[i]);
  	fprintf(stdout, "\nWeighted: \n");
  	for (i = 0; i < ncols; i++)
  	    fprintf(stdout, "cutCoefs[%d] = %f \n", i, cutCoefs[i]);

  	fprintf(stdout, " End computeBendersCutCoefs: \n\n");
    #endif

    if ( coefs != NULL ){
      free (coefs);
      coefs = NULL;
   }


 } //******************************* computeBendersCutCoefs *************************************//



int
addBinaryConstrsSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr)
/**
 * This function appends expliticty binary constraints to the
 * W matrix in the subprobl CPLEX lp object
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param subprobPtr a pointer to the subproblem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
   int i, nzcnt;
   int status;
   char lpname[20];

   double *rhs = NULL;
   char   *sense = NULL;
   int    *rmatbeg = NULL;
   int    *rmatind = NULL;
   double *rmatval= NULL;

   int ncols = subprobPtr->ncols;
   //printf("**subprobPtr->ncols = %d\n", subprobPtr->ncols);
   strcpy(lpname, "sub.lp");

   // Allocate memory
   rhs = (double*)malloc(ncols*sizeof(double));
   sense = (char*)malloc(ncols*sizeof(char));
   rmatbeg = (int*)malloc(ncols*sizeof(int));
   rmatind = (int*)malloc(ncols*sizeof(int));
   rmatval = (double*)malloc(ncols*sizeof(double));

   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addNewRowToLP() \n");
     fprintf (stderr, " Failed to allocate memory to add binary constraint row variables\n");
     fprintf (stderr, " to subproblem lp. \n");
     return 1;
   }

   //****** Add explicity binary constraints to the W matrix for the D^2-Algorithm ******//
   nzcnt = 0;
   for (i = 0; i < ncols; i++){
      if (subprobPtr->ctype[i] == 'B'){
	printf("**subprobPtr->ctype[%d] = %c\n", i, subprobPtr->ctype[i]);
    	rmatbeg[nzcnt] = nzcnt;
        rmatind[nzcnt] = i;
        rmatval[nzcnt] = -1;
        rhs[nzcnt] = -1;
        sense[nzcnt] = 'G';
        nzcnt++;
     } // End if

   } // End for loop

   // Update number of rows in the subproblem lp
   //printf("nzcnt = %d\n", nzcnt);

   status = CPXaddrows(env, lp_sub, 0, nzcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);

	 
	 /* Semih is debugging */
	 
	 CPXwriteprob(env, lp_sub, "hello.lp", "lp");
	 
	 
   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, "Function addBinaryConstrsSubProbLP \n");
      fprintf (stderr, " Failed to allocate memory to add binary constraint row variables\n");
      fprintf (stderr, " to subproblem lp. \n");
   }
   #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_sub, lpname, NULL);
  	if ( status ) {
   		fprintf (stderr, "addBinaryConstrsSubProbLP:\n");
      		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
      		return(status);
 	}
  #endif
   // Free up the memory allocated to the arrays
   if ( rhs != NULL ){
      free (rhs);
      rhs = NULL;
   }
   if ( sense != NULL ){
      free (sense);
      sense = NULL;
   }
   if ( rmatbeg != NULL ){
      free (rmatbeg);
      rmatbeg = NULL;
   }
   if ( rmatind != NULL ){
      free (rmatind);
      rmatind = NULL;
   }

   return status;

 } // ************ End addBinaryConstrsSubProbLP *********************//


int
addFeasColToSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr)
/**
 * This function adds a column to subprobl CPLEX lp object to enable
 * complete resourse as required by the D^2 algorithm
 * This column is penalized in the object function
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param subprobPtr a pointer to the subproblem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
   int i, nzcnt;
   int status;

   double *obj = NULL;
   int    *cmatbeg = NULL;
   int    *cmatind = NULL;
   double *cmatval= NULL;
   double *lb = NULL;
   double *ub = NULL;
   char **colname = NULL;
   int cur_numcols;

   nzcnt = subprobPtr->nrows_mip;

   // Allocate memory
   obj = (double*)malloc(sizeof(double));
   cmatbeg = (int*)malloc(sizeof(int));
   cmatind = (int*)malloc(nzcnt*sizeof(int));
   cmatval = (double*)malloc(nzcnt*sizeof(double));
   lb = (double*)malloc(sizeof(double));
   ub = (double*)malloc(sizeof(double));
   colname = (char**)malloc(sizeof(char*));

   if ( obj == NULL || cmatbeg == NULL || cmatind == NULL ||
        cmatval == NULL || lb == NULL || ub == NULL ||
        colname == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addFeasColToSubProbLP: \n");
     fprintf (stderr, " Failed to allocate memory to feasibility column variables\n");
     fprintf (stderr, " for subproblem lp. \n");
     return 1;
   }
   colname[0] = (char*)malloc(NAMELEN*sizeof(char*));
   if (colname == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addFeasColToSubProbLP: \n");
     fprintf (stderr, " Failed to allocate memory to feasibility colname variable\n");
     fprintf (stderr, " for subproblem lp. \n");
     return 1;
   }

   // Add a feasibility column to subproblem lp data to enable complete resourse
   obj[0] = PENALTY; // High penalty cost: see config.h

   //Add feasibility col name to sub prob data
   strcpy(colname[0], "slack");

   cur_numcols = CPXgetnumcols (env, lp_sub);

   strcpy(subprobPtr->colnames[cur_numcols], "slack");

   cmatbeg[0] = 0;

   //nzcnt = subprobPtr->nrows_mip;

   for (i = 0; i < nzcnt; i++) {
   	cmatind[i] = i;
   	cmatval[i] = 1;
   }

   lb[0] = 0  ;
   ub[0] = CPX_INFBOUND ;

//   status = CPXaddcols(env, lp_sub, 1, nzcnt, obj, cmatbeg, cmatind, cmatval,
//                       lb, ub, colname);
	 status = 0;
   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addFeasColToSubProbLP: \n");
     fprintf (stderr, " Failed to add feasibility column to subproblem lp. \n");
     return status;
   }

   // Free up the memory allocated to the arrays
   if ( obj != NULL ){
      free (obj);
      obj = NULL;
   }
   if ( cmatbeg != NULL ){
      free (cmatbeg);
      cmatbeg = NULL;
   }
   if ( cmatind != NULL ){
      free (cmatind);
      cmatind = NULL;
   }
   if ( lb != NULL ){
      free (lb);
      lb = NULL;
   }
   if ( ub != NULL ){
      free (ub);
      ub = NULL;
   }
   if ( colname != NULL ){
	  free (colname);
      colname = NULL;
   }
	 
	 // SMH: MemLeakFix (May 18, 2015)
	 if (cmatval != NULL) {
		 free(cmatval);
		 cmatval = NULL;
	 }

   return status;

 } // ************ addFeasColToSubProbLP *********************//


int
addOptColToMasterProbLP(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr)
/**
 * This function adds an optimality (theta) column to master problem CPLEX lp object
 * for adding optimality cuts as required by the D^2 algorithm
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr a pointer to the master problem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
  //******************************************************************************
 // ******************** Add theta (opt) variable to master problem ***************
{
 double obj[1];
 char **colname_theta;
 obj[0] = 1;
 int i;
 int status;
 char mastername_lp[NAMELEN];

 strcpy( mastername_lp, "mastermip.lp");

 colname_theta = (char**)malloc(sizeof(char*));
 if ( colname_theta == NULL ) {
     fprintf (stderr, "\n addOptColToMasterProbLP:\n");
     fprintf (stderr, " Failed to allocate memory to master col name ""theta"".\n");
     return 1;
  }

 colname_theta[0] = (char*)malloc(NAMELEN*sizeof(char*));
 if ( colname_theta[0] == NULL ) {
     fprintf (stderr, "\n addOptColToMasterProbLP:\n");
     fprintf (stderr, " Failed to allocate memory to master problem optimality col name ""theta"".\n");
     return 1;
  }
 strcpy(colname_theta[0], "theta");

 status = CPXnewcols(env, lp_master, 1, obj, NULL, NULL, NULL, colname_theta);
 if ( status ) {
     fprintf (stderr, "\n addOptColToMasterProbLP():\n");
     fprintf (stderr, " Failed to add new theta column (var) to master LP.\n");
     return(status);
  }


 i = masterprobPtr->ncols;
 strcpy(masterprobPtr->colnames[i], "theta");
 masterprobPtr->ncols +=1;

 //********************** Write prob to file in lp format **************************
 #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_master, mastername_lp, NULL);
  	if ( status ) {
   		fprintf (stderr, " addOptColToMasterProbLP:\n");
      		fprintf (stderr, "Failure to write lp_master to file, error %d.\n", status);
      		return(status);
 	}
  #endif

  // Free memory
   if ( colname_theta != NULL ){
      free (colname_theta);
      colname_theta = NULL;
   }

   return (0);

} // End addOptColToMasterProbLP

double
getTotalDualsDueToReducedCosts(subproblem_t *subprobPtr, double *redcosts, double *solnY)
/**
 * This function sums the total duals for binary constraints extracted as
 * the negative of the reduced costs
 * @param subprobPtr a pointer to the subproblem data structure
 * @param redcost a pointer to the subproblem reduced costs
 * @param solnY a pointer to the subproblem solution
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
    int i;
   double sum = 0;

   //printf("\n\ngetTotalDualsDueToReducedCosts: \n");

    for (i = 0; i < subprobPtr->ncols; i++){
        if (subprobPtr->ctype[i] == 'B' && solnY[i] >= 1-INT_PRECISION){ // soln[i] = 1, tight
            sum += -1*redcosts[i]*-1;
            //printf("\n\nRedcosts[%d] = %0.6f\n", i, redcosts[i]);
            //printf("sum          = %f",  sum);

        }

    }
    return sum;
    //return 0;

 } // end getTotalDualsDueToReducedCosts



int
DEPloadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,
                      int ncols_master, int nrows_core, int ncols_core, FILE *fpout)
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
 {

    int status;
    int surplus; // To check for sufficiency of the aloocated array size
    int i, j;

    #ifdef DEBUG_FUNC
        fprintf (stderr, "loadMasterProblemData(...):\n");
        fprintf(stderr, "masterprobPtr->nrows = %d\n masterprobPtr->ncols = %d\n", masterprobPtr->nrows, masterprobPtr->ncols);
    #endif

    //***** Create master problem lp by deleting stage_2 cols and rows from core file lp  *******
    status = CPXdelcols(env, lp_core, ncols_master, ncols_core-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to delete stage-2 cols from lp_core, error %d.\n", status);
       return status;
    }

    /* Delete stage-2 rows (constraints) */
    status = CPXdelrows(env, lp_core, nrows_master, nrows_core-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to delete stage-2 rows from lp_core rows, error %d.\n", status);
       return status;
    }

    // Write master to file in lp format
    #ifdef DEBUG_FUNC
        fprintf (stdout, "\nloadMasterProblemData():\n");
   	status = CPXwriteprob(env, lp_core, "master.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write lp_core problem lp to file, error %d.\n", status);
      		return status;
   	}
    #endif


    //************* Get and store the ctype (var type) array for the master problem lp ********
    status = CPXgetctype(env, lp_core, masterprobPtr->ctype, 0, ncols_master-1);
    if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to get ctype from lp_core for the master lp, error %d.\n", status);
       return (status);
    }
    #ifdef DEBUG_FUNC
        fprintf(fpout,"loadMasterProblemData(...): \n Master problem ctype: \n");
        for (j = 0; j < ncols_master; j++)
            fprintf(fpout, "col %d: %c\n", j, masterprobPtr->ctype[j]);
    #endif

   //************* Extract the obj from the lp_core from the core file ****************
   status = CPXgetobj(env, lp_core, masterprobPtr->obj, 0, ncols_master-1);
   if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to write get obj from lp_core for the master lp, error %d.\n", status);
       return(status);
   }

   //******************* Extract the A matrix from lp_core from the core file ****************
   // in row sparse format: more convenient for the formation of the RHS LP later
   status = CPXgetrows(env, lp_core, &masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatind_A,
                      masterprobPtr->rmatval_A, masterprobPtr->rmatspace_A, &surplus, 0, nrows_master-1);
   if ( status ) {
       fprintf (stderr, "loadMasterProblemData(...):\n");
       fprintf (stderr, "Failure to get rows from lp_core for A matrix, error %d.\n", status);
       fprintf (stderr, "Surplus value is: %d.\n", surplus);
       return(status);
   }


  //****** Get the sense vector for the range of contraints for the subproblem ******/
  status = CPXgetsense(env, lp_core, masterprobPtr->sense, 0, nrows_master-1);
  if ( status ) {
     fprintf (stderr, "loadMasterProblemData(...):\n");
     fprintf (stderr, "Failure to get the senses from lp_core problem, error %d.\n", status);
     return(status);
  }

  //****** Get the rhs (b) vector for the master problem *******************************
  status = CPXgetrhs(env, lp_core, masterprobPtr->rhs, 0, nrows_master-1);
  if ( status ) {
     fprintf (stderr, "loadMasterProblemData(...):\n");
     fprintf (stderr, "Failure to get the rhs vector (b) from lp_core problem, error %d.\n", status);
     return(status);
  }


  //**************************** Print the A matrix, senses, and rhs ********************
  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadMasterProblemData(...):\n");
     fprintf(stdout, "masterprobPtr->rmatspace_A = %d\n", masterprobPtr->rmatspace_A);
     fprintf(stdout, "masterprobPtr->nzcnt_A = %d\n", masterprobPtr->nzcnt_A);
     fprintf(stdout, "The A matrix is: \n");
     printSparseMatrix(nrows_master, masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatcnt_A, masterprobPtr->rmatind_A,
                       masterprobPtr->rmatval_A, stdout);

     fprintf(stdout, "master problem senses > \n");
     for (j = 0; j < masterprobPtr->nrows; j++)
           fprintf(stdout, "Row %d: %c\n", j, masterprobPtr->sense[j]);

     fprintf(stdout, "master problem rhs (b) > \n");
     for (j = 0; j < masterprobPtr->nrows; j++)
           fprintf(stdout, "Row %d: %f\n", j, masterprobPtr->rhs[j]);
  #endif

  return (0); // successful return

 } //************************** End DEPloadMasterProblemData function *****************************


int
DEPloadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
                   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
                   FILE *fpout)
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
 {

    int status;
    int surplus; // To check for sufficiency of the aloocated array size
    int i, j;
    char lpname[20];
    int curr_row;
    char *ub = NULL;
    int *indices = NULL;
    int numscens = stochdataPtr->nscens;

    strcpy(lpname, "subproblp.lp");;


   //Create subproblem lp by deleting stage_1 cols and rows from core file lp

   //***************** Delete stage-1 rows (constraints) *********************************
   if (nrows_master > 0) {
          status = CPXdelrows(env, lp_sub, 0, nrows_master-1);
          if ( status ) {
                fprintf(stderr,"\nFunction loadSubProblemData(...): \n");
                fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
                return(status);
           }
   }
   //**********Extract the T matrix from the core file and store it in T *****************
   status = CPXgetcols(env, lp_sub, &subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatind_T, subprobPtr->cmatval_T,
                       subprobPtr->cmatspace_T, &surplus, 0, ncols_master-1);
   if ( status ) {
      fprintf (stderr, "Failure to extraxt T from lp_sub from T(w) matrix, error %d.\n", status);
      fprintf (stderr, "Surplus value is: %d.\n", surplus);
      return(status);
   }

  //******************** Print the T matrix ******************************************
  #ifdef DEBUG_FUNC
     fprintf(stderr, "\nFunction loadSubProblemData(...): \n");
     fprintf(stderr, "subprobPtr->nzcnt_T = %d\n", subprobPtr->nzcnt_T);
     fprintf(stderr, "*******The T matrix is********: \n");
     printSparseMatrix(ncols_master, subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
                       subprobPtr->cmatbeg_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, stderr);
  #endif


  //********************** Delete stage-1 columns (vars) **********************************
  status = CPXdelcols(env, lp_sub, 0, ncols_master-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "Failure to delete lp_sub cols, error %d.\n", status);
     return(status);
  }

  //************* Get and store the ctype (var type) array from the subproblem lp ********
  status = CPXgetctype(env, lp_sub, subprobPtr->ctype, 0, ncols_sub-1);
  if ( status ) {
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf (stderr, "Failure to get ctype from lp_sub, error %d.\n", status);
     return (status);
  }
  #ifdef DEBUG_FUNC
     fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
     fprintf(stdout,"Subproblem ctype: \n");
     for (j = 0; j < ncols_sub; j++)
         fprintf(stdout, "col %d: %c\n", j, subprobPtr->ctype[j]);
  #endif

   // Make updates
   nrows_sub          = CPXgetnumrows (env, lp_sub);
   subprobPtr->nrows  = nrows_sub;
   subprobPtr->nrowsW = nrows_sub;
   subprobPtr->nrows  = nrows_sub;
   subprobPtr->nrows_T = nrows_sub;
   ncols_sub          = CPXgetnumcols (env, lp_sub);
   subprobPtr->ncols  = ncols_sub;

   stochdataPtr->nrows  = nrows_sub;	// Initialize: VERY IMPORTANT!!!!!!

  //********************** Write prob to file in lp format **************************
 #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_sub, lpname, NULL);
  	if ( status ) {
   		fprintf (stderr, "loadSubProblemData:\n");
      		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
      		return(status);
 	}
  #endif


  //************* Extract the obj from the lp_sub ****************
  status = CPXgetobj(env, lp_sub, subprobPtr->obj, 0, ncols_sub-1);
  if ( status ) {
      fprintf (stderr, "loadSubProblemData:\n");
      fprintf (stderr, "Failure to get obj from lp_sub, error %d.\n", status);
      return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem obj\n");
       for (j = 0; j < subprobPtr->ncols; j++) {
           fprintf(stdout, "subprobPtr->obj[%d]: %f\n", j, subprobPtr->obj[j]);
       }
  #endif


  // Get the rhs vector for the range of contraints for the subproblem  //
  // and initialize each scenario rhs with these values                 //
  status = CPXgetrhs(env, lp_sub, subprobPtr->rhs, 0, nrows_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get the rhs from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem rhs\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->rhs[%d]: %f\n", j, subprobPtr->rhs[j]);
       }
  #endif

  // Get the sense vector for the range of contraints for the subproblem       //
  status = CPXgetsense(env, lp_sub, subprobPtr->sense, 0, subprobPtr->nrows-1);
  if ( status ) {
     fprintf (stderr, "loadSubProblemData:\n");
     fprintf (stderr, "Failure to get the senses from the subproblem lp, error %d.\n", status);
     return(status);
  }
  #ifdef DEBUG_FUNC
       fprintf (stdout, "loadSubProblemData:\n");
       fprintf(stdout, "Subproblem senses\n");
       for (j = 0; j < subprobPtr->nrows; j++) {
           fprintf(stdout, "subprobPtr->sense[%d:] %c\n", j, subprobPtr->sense[j]);
       }
  #endif

  // Initialize each scenario rhs with these values: Need to rewrite this section for efficiency !
  // Don't know the number of scenarios yet!!

  // Initialize the scenario subproblem rhs(w)
  for (i = 0; i < numscens; i++){
      for (j = 0; j < subprobPtr->nrows; j++)
          stochdataPtr->rhs[i][j] = subprobPtr->rhs[j];
  }

  #ifdef DEBUG_FUNC
     fprintf (stdout, "loadSubProblemData:\n");
     for (i = 0; i < 3; i++){
        fprintf(stdout, "Subproblem rhs for scenario %d\n", i);
        for (j = 0; j < subprobPtr->nrows; j++) {
            fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", i, j, stochdataPtr->rhs[i][j]);
         }
     }
  #endif

  //Get the lower bounds on the subproblem cols (vars)
  status = CPXgetlb(env, lp_sub, subprobPtr->lb, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get lower bounds lb from lp_sub, error %d.\n", status);
     return (status);
  }
  //Get the upper bounds on the subproblem cols (vars)
  status = CPXgetub(env, lp_sub, subprobPtr->ub, 0, ncols_sub-1);
  if ( status ) {
     fprintf (stderr, "Failure to get upper bounds ub from lp_sub, error %d.\n", status);
     return (status);
  }


  #ifdef DEBUG_FUNC
    fprintf (stdout, "loadSubProblemData:\n");
    fprintf(stdout, "\n\n\n  Subproblem lp column ub and lb:\n");
    for(i = 0; i < ncols_sub; i++){
       fprintf(stdout, "lb[%d] = %f	ub[%d] = %f \n", i,subprobPtr->lb[i],i,subprobPtr->ub[i]);
    }
  #endif



 // Initialize the indices array corresponding to the constraints //
 // for which the rhs coefs are to be changed //
 //  for (j = 0; j < subprobPtr->nrows; j++) {
 //    subprobPtr->indices_row[j] = j;
 //  }


 // Extract the W matrix from the lp_sub in ROW SPARSE MATRIX FORMAT
  status = CPXgetrows(env, lp_sub, &subprobPtr->nzcnt_W, subprobPtr->cmatbeg_W, subprobPtr->cmatind_W,
                      subprobPtr->cmatval_W, subprobPtr->cmatspace_W, &surplus, 0, nrows_sub-1);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "Failure to write get rows from lp_sub for W matrix, error %d.\n", status);
      fprintf (stderr, "Surplus value is: %d.\n", surplus);
      return(status);
  }

  // Print the row sparse matrix W
  #ifdef DEBUG_MAIN
     fprintf (stderr, "\nd2algMain():\n");
     fprintf(stderr, "subprobPtr->cmatspace_W = %d\n", subprobPtr->cmatspace_W);
     fprintf(stderr, "subprobPtr->nzcnt_W = %d\n", subprobPtr->nzcnt_W);
     fprintf(stderr, "The W matrix in row sparse format is: \n");
     printMatrix(nrows_sub, subprobPtr->nzcnt_W, subprobPtr->cmatbeg_W, subprobPtr->cmatind_W,
                       subprobPtr->cmatval_W, stderr);
  #endif



  return (0); // Successful return

 } //************************** End DEPloadSubProblemData function *****************************


int
DEPsetupDEP(CPXENVptr env, CPXLPptr lp_dep, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
            double random_T, int nrows_master, int ncols_master, char** depcolnames, FILE *fpout)
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
 {

    int status = 0;
    int i, j, scenario;
    int cnt;
    double *obj;
    char **colname;
    int *rowlist;
    int *collist;
    double *vallist;
    int start, stop;
    int size;
    int index;

    size = subprobPtr->nrows*ncols_master;    // max # of nonzero in T matrix

    // Allocate memory
    obj = (double*)malloc(subprobPtr->ncols*sizeof(double));
    colname = (char**)malloc(subprobPtr->ncols*sizeof(char*));
    rowlist = (int*)malloc(size*sizeof(int));
    collist = (int*)malloc(size*sizeof(int));
    vallist = (double*)malloc(size*sizeof(double));

   // size    = stochdataPtr->nscens*subprobPtr->ncols;

   // depcolnames = (char**)malloc(size*sizeof(char*));

    if (obj     == NULL || colname == NULL || rowlist == NULL ||
        collist == NULL || vallist == NULL){
        fprintf(stderr, "DEPsetupDEP(): \n");
        fprintf(stderr, "Failed to allocate memory to obj and colname arrays.\n");
        return 1;
    }
    for(i = 0; i < subprobPtr->ncols; i++) {
    	colname[i] = (char*)malloc((NAMELEN)*sizeof(char));
    	if (colname[i] == NULL) {
        	fprintf(stderr, "DEPsetupDEP(): \n");
       		fprintf(stderr, "Failed to allocate memory to colname[%d] array.\n", i);
        	return 1;
    	}
    }
    //for(i = 0; i < size; i++) {
    //	depcolnames[i] = (char*)malloc((NAMELEN)*sizeof(char));
    //	if (depcolnames[i] == NULL) {
    //    	fprintf(stderr, "DEPsetupDEP(): \n");
    //   		fprintf(stderr, "Failed to allocate memory to depcolname[%d] array.\n", i);
    //    	return 1;
    //	}
    // }


    // To add the W matrix for each scenario (soon)
    // shift the cols by first stage ncols
    for (i = 0; i < subprobPtr->nzcnt_W; i++)
            subprobPtr->cmatind_W[i] += ncols_master;
    // To add the T matrix for each scenario (soon)
    // shift the rows by first stage nrows
    for (i = 0; i < subprobPtr->nzcnt_T; i++)
            subprobPtr->cmatind_T[i] += nrows_master;
   index = 0;
   for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
       //fprintf(stdout, "scenario = %d\n", scenario);
       // Setup this scenario objective
       for (i = 0; i < subprobPtr->ncols; i++)
           obj[i] = subprobPtr->obj[i]*stochdataPtr->scenProb[scenario];

       #ifdef DEBUG_FUNC
           for (i = 0; i < subprobPtr->ncols; i++)
                fprintf(stdout, "obj[%d] = %f\n", i, obj[i]);
       #endif

       // Setup this scenario column (var) names
       for (i = 0; i < subprobPtr->ncols; i++){
           strcpy(colname[i], subprobPtr->colnames[i]);
           strcpy(depcolnames[index], subprobPtr->colnames[i]);
           strcat(colname[i], stochdataPtr->scenName[scenario]);
           strcat(depcolnames[index], stochdataPtr->scenName[scenario]);
           //fprintf(stdout,"depcolnames[%d] = %s\n", index, depcolnames[index]);
           index++;
       }
       //printf("index = %d\n", index);

       #ifdef DEBUG_FUNC
           for (i = 0; i < subprobPtr->ncols; i++) {
                fprintf(stdout, "colname[%d] = %s\n", i, colname[i]);
           }
       #endif

       /////////////////////////////////////////////////
       // Add new columns to dep for current scenario //
       /////////////////////////////////////////////////
       status = CPXnewcols(env, lp_dep, subprobPtr->ncols, obj, subprobPtr->lb, subprobPtr->ub,
                           subprobPtr->ctype, colname);
       if (status) {
          fprintf(stderr, "Failed to add new columns to D.E.P lp for scenario %d\n", scenario);
          return status;
       }

       /////////////////////////////////////////////////
       //       Update W col and T row indices        //
       /////////////////////////////////////////////////
       if (scenario != 0) {
           for (i = 0; i < subprobPtr->nzcnt_W; i++) {
               subprobPtr->cmatind_W[i] += subprobPtr->ncols;
           }
           for (i = 0; i < subprobPtr->nzcnt_T; i++) {
               subprobPtr->cmatind_T[i] += subprobPtr->nrows;
           }
        }
        status = CPXaddrows(env, lp_dep, 0, subprobPtr->nrows, subprobPtr->nzcnt_W,
                            stochdataPtr->rhs[scenario], subprobPtr->sense,
                            subprobPtr->cmatbeg_W, subprobPtr->cmatind_W,
                            subprobPtr->cmatval_W, colname, NULL);
  	if ( status ) {
      		fprintf (stderr, "\n DEPsetupDEP():\n");
      		fprintf (stderr, "Failure to add rows to DEP lp for scenario %d matrix, error %d.\n",
      		         scenario, status);
      		return(status);
        }
        //////////////////////////////////////////////////
        //      Add T matrix coefs for this scenario    //
        //////////////////////////////////////////////////
        cnt = 0;
        if (random_T == 0) {  //Constant technolody matrix T
            for (i = 0; i < ncols_master; i++) {
               if (i == ncols_master-1) {
               	   start = subprobPtr->cmatbeg_T[i];
                   stop  = subprobPtr->nzcnt_T;
               } else {
                   start = subprobPtr->cmatbeg_T[i];
                   stop  = subprobPtr->cmatbeg_T[i+1];
               }
              for (j = start; j < stop; j++){
                  rowlist[cnt] = subprobPtr->cmatind_T[j];
                  collist[cnt] = i;
                  vallist[cnt] = subprobPtr->cmatval_T[j];
                  cnt++;
              }
           } // end out for loop
       } else { // Random technology matrix T(w)
        for (i = 0; i < ncols_master; i++) {
            if (i == ncols_master-1) {
               start = stochdataPtr->cmatbeg_T[scenario][i];
               stop  = stochdataPtr->cnzcnt_T[scenario];
            } else {
               start = stochdataPtr->cmatbeg_T[scenario][i];
               stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
            }
            for (j = start; j < stop; j++){
                  rowlist[cnt] = stochdataPtr->cmatind_T[scenario][j];
                  collist[cnt] = i;
                  vallist[cnt] = stochdataPtr->cmatval_T[scenario][j];
                  cnt++;
              }
        } // end outer for loop

      } // End outer if\else statement

      /////////////////////////////////////////////////
      //             Change DEP lp coefs             //
      /////////////////////////////////////////////////
      if (random_T == 0)   // constant technolody matrix T
            status = CPXchgcoeflist(env, lp_dep, subprobPtr->nzcnt_T, rowlist, collist, vallist);
       else  // random technolody matrix T
            status = CPXchgcoeflist(env, lp_dep, stochdataPtr->cnzcnt_T[scenario], rowlist, collist, vallist);

      if ( status ) {
      	   fprintf (stderr, "\n DEPsetupDEP():\n");
      	   fprintf (stderr, "Failure to change coef list (T) to DEP lp for scenario %d matrix, error %d.\n",
      		    scenario, status);
      	   return(status);
      }

   } // for scenario loop


    // Free arrays
    if ( **colname != NULL ) {
      free (colname);
      **colname = NULL;
    }
    if ( obj != NULL ) {
      free (obj);
      obj = NULL;
    }

    return status;

}

// Added Dece 26, 2003

int
getnumscenarios(char *filename)
/**
 * Reads the STOCH file data and determines the total number of scenarios
 * @param filename  STOCH file name
 * @return 0 number of scenarios
 */
 {

   FILE *stoch;
   char field1[NAMELEN], field2[NAMELEN];
   char buffer[LENGTH];

   int numscenarios = 0;  // counters

   // Open the TIME file
   stoch = fopen(filename,"r" );
   if(stoch == NULL) {
   	fprintf(stderr, "\ngetnumscenarios(...): \n");
     	fprintf(stderr, " Could not open the STOCH file %s for reading!\n", filename);
     	return(1);
    }


   //************************* Read the TIME file ************************************
   // Read first line: e.g. STOCH	example
   if (fgets (buffer, LENGTH, stoch) != NULL) {
        sscanf(buffer, "%s", field1);
        if (strcmp(field1, "STOCH") != 0) {
           fprintf(stderr, "\ngetnumscenarios(...): \n");
           fprintf(stderr, " The first line of the STOCH file %s must start with ""STOCH"" \n!", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
   }
   // Read the second line: e.g SCENARIOS	DISCRETE
   if (fgets (buffer, LENGTH, stoch) != NULL){
        sscanf(buffer, "%s%s", field1, field2);
       if ( strcmp(field1, "SCENARIOS") != 0){
       	   fprintf(stderr, "\nloadStochFile(...): \n");
           fprintf(stderr, " The second line of file %s must start with ""SCENARIOS""!\n", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
        if ( strcmp(field2, "DISCRETE") != 0){
       	   fprintf(stderr, "\ngetnumscenarios(...): \n");
           fprintf(stderr, " The second word in the second line of file %s must be ""DISCRETE""!\n", filename);
           fprintf(stderr, "Exiting...\n");
           return(1);
        }
   }

   // Read the rest of the file and count scenarios
   while (fgets (buffer, LENGTH, stoch) != NULL) {
      //fprintf(stderr, "BUFFER: %s\n", buffer);
      sscanf(buffer, "%s", field1);
      if (strcmp(field1, "SC") == 0)  // Read in new scenario data
      	    numscenarios++; 	// count scenario
   }

   // close the STOCH file
   fclose(stoch);

   return numscenarios;

} //************************** End getnumscenarios() function *****************************


double
getrandobjsubproblowerbound(stochfile_info *stochdataPtr, subproblem_t *subprobPtr)
/**
 * Computes and returns the value of the lower bound on the second stage problem.
 * This is required by the D^2 algorithm to keep the subproblem obj positive
 * for all scenario subproblems
 * @param stochdataPtr pointer to the subproblem lp stoch data structure
 * @param subprobPtr   pointer to the subproblem lp data structure
 * @return the lower bound SUBPROB_LB_L on second stage problem
 */
{
 int i, j, k;
 double lb_scen;
 double lb = 0;

  for (i = 0; i < stochdataPtr->nscens; i++) {
      lb_scen = 0;
      // Do non-random obj parts first
      for (j = 0; j < stochdataPtr->obj_index[0]; j++) {
          //fprintf(stdout, "subprobPtr->ctype[%d] =  %c: subprobPtr->obj[%d] = %f\n", j, subprobPtr->ctype[j], j, subprobPtr->obj[j]);
          if (subprobPtr->ctype[j] == 'B' || subprobPtr->ctype[j] == 'I') {
              if (subprobPtr->obj[j] < 0)
   	          lb_scen += subprobPtr->obj[j];
          }
      }
      if (stochdataPtr->obj_cnt < stochdataPtr->ncols-1) {
          for (j = stochdataPtr->obj_index[stochdataPtr->obj_cnt]; j < stochdataPtr->ncols-1; j++) { // don't include artificial variable
              //fprintf(stdout, "subprobPtr->ctype[%d] =  %c: subprobPtr->obj[%d] = %f\n", j, subprobPtr->ctype[j], j, subprobPtr->obj[j]);
              if (subprobPtr->ctype[j] == 'B' || subprobPtr->ctype[j] == 'I') {
                  if (subprobPtr->obj[j] < 0)
   	               lb_scen += subprobPtr->obj[j];
              }
          }
      }
      // Do random obj parts now
      for (j = 0; j < stochdataPtr->obj_cnt; j++) { // don't include artificial variable
          k = stochdataPtr->obj_index[j];
          //fprintf(stdout, "subprobPtr->ctype[%d] =  %c: stochdataPtr->obj[%d][%d] = %f\n",
          //                 k, subprobPtr->ctype[k], k, i, stochdataPtr->obj[i][k]);
          if (subprobPtr->ctype[k] == 'B' || subprobPtr->ctype[k] == 'I') {
              if (stochdataPtr->obj[i][k] < 0)
   	          lb_scen += stochdataPtr->obj[i][k];
          }
      }
      //fprintf(stdout, "lb_scen %d = %f\n", i+1, lb_scen);
      lb +=  stochdataPtr->scenProb[i]*lb_scen;

  } // End outer for loop

     			//else { // if ctype[i] is 'G' need to opt the subproblems!!
               		// Will take care of this case later!!
      			//}
  return -1*lb;

} // ******************* getrandobjsubproblowerbound() ************************** //


int
resetconstraints(CPXENVptr env, CPXLPptr lp)
/**
 * Resets <= to <= constraints in the CORE FILE lp. This is required by the D2 Algorithm.
 * This is done row by row by multiplying both sides of the constraint by -1.
 * The scenario RHS must also be multiplied by -1 in the function loadStochFile function
 * @param env CPLEX environment
 * @param lp pointer to the lp model
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */
 {
   int i, j;
   int cur_numcols;
   int cur_numrows;
   int status =  0;
   int nzcnt;
   int  *rmatbeg;
   int  *rmatind;
   double *rmatval;
   char *sense;

   int surplus;

   int *rowindex;
   char *newsense;
   double *rhs;

   int *rows;

   cur_numcols = CPXgetnumcols (env, lp);
   cur_numrows = CPXgetnumrows (env, lp);

   rmatbeg = (int *) malloc (cur_numrows*sizeof(int));
   rmatind = (int *) malloc (cur_numcols*sizeof(int));
   rmatval = (double *) malloc (cur_numcols*sizeof(double));
   rows    = (int *) malloc (cur_numcols*sizeof(int));
   sense   = (char *) malloc (cur_numrows*sizeof(char));

   if (rmatbeg == NULL || rmatind == NULL || rmatval == NULL || rows == NULL ||
       sense == NULL) {
        fprintf(stderr, "\nresetconstraints: ");
        fprintf(stderr, "\nMemory allocation failure for arrays.");
        return 1;
   }
   rowindex   = (int *) malloc (sizeof(int));
   rhs        = (double *) malloc (sizeof(double));
   newsense   = (char *) malloc (sizeof(char));
   if (rowindex == NULL || rhs == NULL || newsense == NULL) {
        fprintf(stderr, "\nresetconstraints: ");
        fprintf(stderr, "\n**Memory allocation failure for arrays.");
        return 1;
   }

   newsense[0] = 'G';

   //printf("\n  cur_numrows = %d\n", cur_numrows);
   //printf("   cur_numcols = %d\n", cur_numcols);

   // Get sense
   status = CPXgetsense (env, lp, sense, 0, cur_numrows-1);
   if (status) {
       fprintf(stderr, "\nresetconstraints: ");
       fprintf(stderr, "\nCPXgetsense failure at row %d. status = %d\n", i, status);
       return status;
   }

   // Loop thru all the contraints
   for (i = 0; i < cur_numrows; i++) {

       //printf("   cur_row = %d\n", i);
       // Change sense
       //printf("   sense[%d] = %c\n", i, sense[i]);
        if (sense[i] == 'L') {

            rowindex[0] = i;
            status =  CPXchgsense(env, lp, 1, rowindex, newsense);
            if (status) {
                 fprintf(stderr, "\nresetconstraints: ");
                 fprintf(stderr, "\nCPXchgsense failure at row %d. status = %d\n", i, status);
                 return status;
            }

            // Change RHS
            status = CPXgetrhs (env, lp, rhs, i, i);
            if (status) {
                 fprintf(stderr, "\nresetconstraints: ");
                 fprintf(stderr, "\nCPXgetrhs failure at row %d. status = %d\n", i, status);
                 return status;
            }
            rhs[0] *= -1;
            status = CPXchgrhs (env,lp,1,rowindex, rhs);
            if (status) {
                 fprintf(stderr, "\nresetconstraints: ");
                 fprintf(stderr, "\nCPXchgrhs failure at row %d. status = %d\n", i, status);
                 return status;
            }

            // Change constraint
            status =  CPXgetrows (env,lp, &nzcnt, rmatbeg, rmatind, rmatval, cur_numcols,
                             &surplus, i, i);
            if (status) {
                  fprintf(stderr, "\nresetconstraints: ");
                  fprintf(stderr, "\nCPXgetrows failure at row %d. status = %d\n", i, status);
                 return status;
            }
            for (j = 0; j < nzcnt; j++){
                  rmatval[j] *= -1;
                  rows[j] = i;
            }

            status = CPXchgcoeflist (env, lp, nzcnt, rows, rmatind, rmatval);
            if (status) {
                  fprintf(stderr, "\nresetconstraints: ");
                  fprintf(stderr, "\nCPXchgcoeflist failure at row %d. status = %d\n", i, status);
                 return status;
            }

       } // End if


   } // End for loop



   #ifdef DEBUG_FUNC
        //Write sub prob to file in lp format
   	status = CPXwriteprob(env, lp, "prob.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "resetconstraints:\n");
      		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		return status;
   	}
    #endif


   free(rmatbeg);
   free(rmatind);
   free(rmatval);
   free(rows);
   free(sense);

   free(rowindex);
   free(rhs);
   free(newsense);


  // printf("\n Done resetting constraints...\n");

   return 0;

 } // resetconstraints




/************************************************************************
 *			ORDERED LIST FUNCS				*
 *									*
 * This part contains function for an implementation of an ordered     	*
 * linked list. Ordering is by nondecreasing objective value of the 	*
 * list node struct variable.						*
 *  									*
 *   Headers: orderedlist.h						*
 * 									*
 *  Author: Lewis Ntaimo						*
 *  Date:   December 21, 2003						*
 ************************************************************************/

 /**
  * This function initialize the ordered list
  * @param list pointer to the lost
  */
int initialise(ORDEREDLIST *list)
{
	list->front = (NODE *) malloc(sizeof(NODE));
	if (list->front == NULL) {
		fprintf(stderr, "\nCannot initialize list! \n");
		return 0;
	}
	list->front->next = NULL;
	list->size = 0;
	return 1;
}

 /**
  * This function adds a node to the ordered list
  * @param list pointer to the lost
  * @param newnode pointer to the new node to add to list
  * @return returns nonzero on success, zero on failure
  */
  int add(ORDEREDLIST *list, NODE *newnode)
  {
        NODE *temp;

        if (isempty(list)) { /* Add at front of node */
              #ifdef DEBUG
                   printf("Inserting first node\n");
              #endif
           newnode->next = NULL;
           list->front->next = newnode;
        } else {
             if (newnode->objval <= list->front->next->objval) { /* insert in front list */
                  #ifdef DEBUG
                        printf("Inserting new node at front of list \n");
                  #endif
                 newnode->next = list->front->next;
                 list->front->next = newnode;
             } else { /* insert inside list */
                temp = list->front->next;
                do {

                   if (temp->next == NULL) { /* Append at end of list */
                       #ifdef DEBUG
                            printf("Inserting new node at end of list \n");
                       #endif
                       temp->next = newnode;
                       newnode->next = NULL;
                       break;
                   } else if (newnode->objval <= temp->next->objval) {
                       #ifdef DEBUG
                            printf("Inserting new node inside of list \n");
                       #endif
                      newnode->next = temp->next;
                      temp->next = newnode;
                      break;
                   }
                   temp = temp->next;

                } while (temp != NULL);

             } /* end inner if/else */

        } /* end outer if/else */

        list->size++;
        return 1;
  }

  /**
   * This function removes the node at the front of ordered node list
   * @param list pointer to the lost
   */
  void removefrontnode(ORDEREDLIST *list)
  {
        NODE *oldnode;
        double thisobjval;
        oldnode = list->front->next;
        thisobjval = oldnode->objval;
        /* Check if removing the last node from the list */
        if (list->front->next->next == NULL)
            list->front->next = NULL;
        else
			list->front->next = list->front->next->next;

        /////// 3/1/04 ///////

     /*   if (oldnode->varind != NULL)
             free(oldnode->varind);
        if (oldnode->varval != NULL)
             free(oldnode->varval);
        if (oldnode->varbnds != NULL)
             free(oldnode->varbnds);
       */
        /////////////////////
        if (oldnode != NULL)
		{
			// SMH: MemLeakFix (May 18, 2015)
			free(oldnode->varind);
			free(oldnode->varval);
			free(oldnode->varbnds);
			free(oldnode->gamma);
			
			free(oldnode);
		}

        list->size--;
  }

  /**
   * This function checks if ordered list is empty
   * @param list pointer to the lost
   * @return returns nonzero if list is empty, otherwise returns zero
   */
   int isempty(ORDEREDLIST *list)
    {
        return (list->front->next == NULL);
    }

  /**
   * This function creates a node and sets its member values
   * @param objval pointer to the lost
   * @param nodedepth pointer to the new node to add to list
   * @param max_numnodes maximum number of nodes to explore
   * @param rmatind pointer to array of branching variable indices
   * @param rmatval pointer to array of branching variable matrix coefs
   * @param varbnds pointer to array of branching variable rhs bounds
   * @param nu the nu value for this node
   * @param numcols_master number of first stage variables
   * @param gamma array containing gamma values for this node
   * @return returns a pointer to the newly created node
   */
   NODE *createnode(double objval, int nodedepth, int max_numnodes, int *varind,
                    double *varval, double * varbnds, int ncols_A, double nu, double *gamma)
   {
        int j;
        int size;
        NODE *newnode = NULL;
        size = max_numnodes;

        // Allocate memory
        //printf("Allocating memory to node vars...\n");
        //printf("size = %d\n", size);

        newnode = (NODE *)malloc(sizeof(NODE));
        newnode->varind  = (int*)malloc(size*sizeof(int));
        newnode->varval  = (double *)malloc(size*sizeof(double));
        newnode->varbnds = (double *)malloc(size*sizeof(double));

        newnode->gamma = (double *)malloc(ncols_A*sizeof(double));

        if (newnode == NULL || varind == NULL || varval == NULL || varbnds == NULL || gamma == NULL) {
             fprintf(stderr, "makenode(): Failure to allocate memory to new node struct. Bailing out...");
             return 0;
        }
        newnode->objval     = objval;
        newnode->node_depth = nodedepth;
        for (j = 0; j < nodedepth; j++) {
             newnode->varind[j] = varind[j] ;
             newnode->varval[j] = varval[j];
             newnode->varbnds[j] = varbnds[j];
        }
        // Copy nu and gamma
        newnode->nu = nu;
        for (j = 0; j < ncols_A; j++)
             newnode->gamma[j] = gamma[j] ;

        newnode->next = NULL;

        #ifdef DEBUG_FUNC
             printf("newnode->node_depth = %3d\t", newnode->node_depth);
             printf("copied arrays...\n");
             for (j = 0; j < nodedepth; j++) {
                 printf("newnode->varind[%d] = %3d\t", j, newnode->varind[j]);
                 printf("newnode->varval[%d] = %3.3f\t", j, newnode->varval[j]);
                 printf("newnode->varbnds[%d] = %3.3f\n", j, newnode->varbnds[j]);
             }
             printf("Done copying arrays...\n");
        #endif
        #ifdef DEBUG_FUNC
             printf("\n\nNEW NODE: id = %d\n", newnode->node_depth);
             printf("newnode->nu = %f \n", newnode->nu);
             for (j = 0; j < ncols_A; j++)
                 printf("newnode->gamma[%d] = %3f\n", j, newnode->gamma[j]);
        #endif

        return newnode;
 } // End createnode

/*********************END ORDERED LIST FUNCS **************************/



/************************************************************************
 *			BRANCH-AND-BOUND FUNCS				*
 *									*
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
BBgetMaxFractionIndex(double *solnY, int numcols, int *index, double *value, char* ctype)
{
  int j;
  int status = 0;
  double lval, uval;
  double curr_max;
  double max = 0;

  *index = -1;
  *value = -1.0;

  for (j = 0; j < numcols; j++)
  {
	   if (ctype[j] != 'C') {
	      lval = solnY[j] - floor(solnY[j]);
	      uval = ceil(solnY[j]) - solnY[j];
              if (solnY[j] > INT_PRECISION && solnY[j] < 1-INT_PRECISION){ // Fractional component
		  curr_max = minimum(lval, uval);
		  if (curr_max > max){
		      max = curr_max;
		      *value = solnY[j]; // Fractional soln
		      *index = j;     // Fractional index
		      status = 1;
	           }
             } // end outer if
           } // End outer if
   }

   return status;

} //****************** BBgetMaxFractionIndex ***************************

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
BBaddbranchconstrs(CPXENVptr env, CPXLPptr lp_sub, int num_constrs, int *varind, double *varval,
                   double *varbnds)
 {
   int i;
   int nzcnt;
   int status;

   char   *sense   = NULL;
   int    *rmatbeg = NULL;
   int    *rmatind = NULL;
   double *rmatval = NULL;
   double *rhs     = NULL;

   // Allocate memory
   rhs     = (double*)malloc(num_constrs*sizeof(double));
   sense   = (char*)malloc(num_constrs*sizeof(char));
   rmatbeg = (int*)malloc(num_constrs*sizeof(int));
   rmatind = (int*)malloc(num_constrs*sizeof(int));
   rmatval = (double*)malloc(num_constrs*sizeof(double));

   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "BBaddbranchconstrs: \n");
     fprintf (stderr, " Failed to allocate memory to constraint row variables\n");
     fprintf (stderr, " to subproblem lp. \n");
     return 1;
   }
   nzcnt = 0;


   // Add branching constraints to the subproblem lp
   for (i = 0; i < num_constrs; i++){
    	rmatbeg[i] = i;
        rmatind[i] = varind[i];
        rmatval[i] = varval[i];
        rhs[i]     = varbnds[i];
        sense[i] = 'G';
        nzcnt++;

   } // End for loop

   status = CPXaddrows(env, lp_sub, 0, num_constrs, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);
   if (status) {
        fprintf (stderr, "\nBBaddbranchconstrs:\n");
        fprintf (stderr, " Failed to add branching rows to subproblem lp. Bailing out...\n");
   }

   #ifdef DEBUG_FUNC
  	status = CPXwriteprob(env, lp_sub, "nodelp.lp", NULL);
  	if ( status ) {
   		fprintf (stderr, "\nBBaddbranchconstrs:\n");
      		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
      		return(status);
 	}
  #endif

  // Free up the memory allocated to the arrays

  free (sense);
  free (rmatbeg);
  free (rmatind);
  free (rmatval);
  free (rhs);

  return status;

 } // ************ End BBaddbranchconstrs() *********************//


/**
 * This function performs a branch-and-bound (B&B) alg that uses the best  node
 * strategy. The size of the B&B tree is based on the maximum # of nodes to
 * explore passed to this function.
 * @param env CPLEX environment pointer
 * @param lp_sub subproblem lp pointer
 * @param subprobPtr   pointer to the subproblem lp data structure
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param
 * @return status zero for success and non-zero otherwise
 */
int
BBdobranchbound(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, stochfile_info *stochdataPtr,
                int max_numnodes, int scenario, double *bb_nu, double **bb_gamma, int *num_nodesexplored,
                int *num_termnodes, double *cputime_sublp, double *bestbound, double *solnY_sub, int *fract_sub, int random_T)
{
   int i, j;
   //double bestbound;
   double objval;
   double *solnY;
   double *duals;
   double *redcosts;		// Reduced costs
   int size;
   int cur_numcols;
   int cur_numrows;
   int orig_numrows;
   int status = 0;
   int solstat;

    ORDEREDLIST list;		// TB&B list
    NODE *newnode; 		// A node to hold nodal information

   int    frac_index;		// Initial node subproblem max fractional index
   double frac_value;		// Initial node subproblem max fractional value
   int    lfrac_index;		// Left branch subproblem max fractional index
   double lfrac_value;		// Left branch subproblem max fractional value
   int    rfrac_index;		// Right branch subproblem max fractional index
   double rfrac_value;		// Right branch subproblem max fractional value

   int    fractional;		// flag = 0 if integralm soln found, nonzero otherwise
   int node_id;			// Nodes are indexed from 1 starting with the initial node
   int node_depth;		// The depth of a node in the TB&B Tree
   int cur_node_depth;		// The depth of a currently removed node in the TB&B Tree
   int *varind;                 // Store branching variable index as double
   double *varval;              // Store branching variable matrix coef
   double *varbnds;             // Store branching variable rhs bound

   double *rhs;		   // subproblem rhs

   int numcols_master;
   double nu;
   double *gamma;
   double total_redcosts;
   double cputime_start;
   double cputime_stop;
   int fathomed;		// Flag: 1 if node from list is fathomed by bound, 0 otherwise
   int start;			// Flag: 1 if initially starting TB&B search, 0 otherwise;
   int newnode_created = 0;        // Flag: 1 if node has been created, will need to free it

    cur_numcols = CPXgetnumcols (env, lp_sub);
    cur_numrows = CPXgetnumrows (env, lp_sub);
    orig_numrows =  CPXgetnumrows (env, lp_sub);
    numcols_master = stochdataPtr->ncols;

   ////////////////////////////////////////////////////////////
   //                  Memory Allocations                    //
   //                   		    		     //
   ////////////////////////////////////////////////////////////
   solnY    = (double *)malloc(cur_numcols*sizeof(double));
   redcosts = (double *)malloc(cur_numcols*sizeof(double));
   duals = (double *)malloc((cur_numrows+max_numnodes)*sizeof(double));
   if (solnY == NULL || duals == NULL || redcosts == NULL) {
       fprintf(stderr, "\nBBdobranchbound:");
       fprintf(stderr, "\nFailure to allocate memory to solnY. Bailing out...");
       return 1;
    }
    size = max_numnodes;//2 + 1;
    varind  = (int *)malloc(size*sizeof(int));
    varval  = (double *)malloc(size*sizeof(double));
    varbnds = (double *)malloc(size*sizeof(double));
    if (varind == NULL || varval == NULL || varbnds == NULL) {
         fprintf(stderr, "\nBBdobranchbound:");
         fprintf(stderr, "\nFailure to allocate memory to node arrays. Bailing out...");
         return 1;
    }
    rhs = (double *)malloc((cur_numcols+max_numnodes)*sizeof(double));
    //gamma = (double *)malloc(cur_numcols*sizeof(double));
    gamma = (double *)malloc(numcols_master*sizeof(double));
    if (rhs == NULL || gamma == NULL) {
         fprintf(stderr, "\nBBdobranchbound:");
         fprintf(stderr, "\nFailure to allocate memory to gamma arrays. Bailing out...");
         return 1;
    }



    ////////////////////////////////////////////////////////////
    //                         Start TB&B                     //
    //                   		    		      //
    ////////////////////////////////////////////////////////////
    //fprintf (stdout, "\nStarting B&B for scenario %d: \n", scenario);
    //status = CPXwriteprob(env, lp_sub, "BBsub.lp", NULL);
    //if ( status ) {
    //	 fprintf (stderr, "BBdobranchbound():\n");
    //   	 fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
    //  	 fprintf (stderr, "Error code: %d\n", status);
    //  	 return status;
    //}

    // Initializations
    *bestbound = CPX_INFBOUND;
    node_id = 1;		// Start counting the nodes explored
    //fprintf (stdout, "\nInitializing list\n");
    if (!initialise(&list))
    {
       fprintf(stderr, "\nBBdobranchbound():\n");
       fprintf(stderr, "Unable to initialise the TB&B ordered list\n");
       return 1;
    }

    // Solve initial node LP
    //fprintf(stdout, "\nOptimizing node 0 subproblem scenario index %d  \n", scenario);
    cputime_start = clock();
    status = CPXdualopt (env, lp_sub);
    //status = CPXlpopt (env, lp_sub);
    cputime_stop = clock();
    *cputime_sublp += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
    //fprintf(stdout, "\nDone Re-ptimizing node 0 subproblem scenario index %d  \n", scenario);

    if ( status ) {
        fprintf(stderr, "\nBBdobranchbound():\n");
        fprintf (stderr, "Failed to optimize scenario %d subproblem.\n", scenario);
        fprintf(stderr, "Error %d", status);
        return status;
    }
    //solstat = CPXgetstat (env, lp_sub);
    //fprintf (stdout, "Solution status %d.\n", solstat);
    status  = CPXgetobjval (env, lp_sub, &objval);
    if ( status ) {
         fprintf(stderr, "\nBBdobranchbound():\n");
         fprintf (stderr,"Failed to obtain objective value for subproblem scenario index %d .\n", scenario);
         fprintf(stderr, "Error code: %d", status);
         return status;
     }
     //fprintf (stdout, "objval %f.\n", objval);
     status = CPXgetx (env, lp_sub, solnY, 0, cur_numcols-1);
     if ( status ) {
          fprintf(stderr, "\nBBdobranchbound():\n");
          fprintf (stderr, "Failed to obtain subproblem LP solution.\n");
          fprintf (stderr, "status: %d\n", status);
          return status;
     }
     status = CPXgetdj (env, lp_sub, redcosts, 0, cur_numcols-1);
     if ( status ) {
          fprintf(stderr, "\nBBdobranchbound():\n");
          fprintf (stderr, "Failed to obtain subproblem reduced costs.\n");
          fprintf (stderr, "status: %d\n", status);
          return status;
      }
      status = CPXgetpi (env, lp_sub, duals, 0, cur_numrows-1);
      if ( status ) {
          fprintf(stderr, "\nBBdobranchbound():\n");
          fprintf (stderr, "Failed to obtain subproblem dual solution.\n");
          fprintf (stderr, "status: %d\n", status);
          return status;
      }

     //status = CPXgetrhs (env, lp_sub, rhs, 0, cur_numrows-1);
     // if ( status ) {
     //    fprintf(stderr, "\nBBdobranchbound():\n");
     //     fprintf (stderr, "Failed to obtain subproblem rhs.\n");
     //     fprintf (stderr, "status: %d\n", status);
     //     return status;
     //}

  /*   //////////////////////////////////////////
      fprintf(stdout, "\nPrimal Solution: Root Node\n");
      for (j = 0; j < cur_numcols; j++) {
           fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], solnY[j]);
      }
       fprintf(stdout, "\nDual Solution:\n");
      for (j = 0; j < cur_numrows; j++) {
           fprintf (stdout, " Row %d: = %6.6f\n", j, duals[j]);
      }
      fprintf(stdout, "\nReduced Costs:\n");
      for (j = 0; j < cur_numcols-1; j++) {
           fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], redcosts[j]);
      }

      fprintf(stdout, "\nrhs:\n");
      for (j = 0; j < cur_numrows; j++) {
           fprintf (stdout, " Row %d: = %6.6f\n", j, rhs[j]);
      }
      fprintf(stdout, "\nRHS:\n");
      for (j = 0; j < cur_numrows; j++) {
           fprintf (stdout, " Row %d: = %6.6f\n", j, stochdataPtr->rhs[scenario][j]);
      }
    */  /////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////
      //                      Compute nu and gamma                                //
      //////////////////////////////////////////////////////////////////////////////

      total_redcosts = getTotalDualsDueToReducedCosts(subprobPtr, redcosts, solnY);
      BBcomputegammanu(stochdataPtr, subprobPtr, scenario, cur_numrows, duals,
                             gamma, &nu, random_T, total_redcosts);
      *num_nodesexplored += 1; // count this node as exlored

      //printf("\n\n nu = %3.3g\n", nu);
      //for (j = 0; j < numcols_master; j++) {
      //    printf(" gamma[%d] = %3.3g\n", j, gamma[j]);
      //}
      //printf("\n");

      // Store the nu and gammas
      //bb_nu[*num_nodesexplored] = nu;
      //for (j = 0; j < numcols_master; j++) {
      //    bb_gamma[*num_nodesexplored][j] = gamma[j];
      ///}

      //////////////////////////////////////////////////////////////////////////////


     // Get max fractional index and value
     fractional = BBgetMaxFractionIndex(solnY, cur_numcols-1, &frac_index, &frac_value, subprobPtr->ctype);
     //fprintf(stdout, "\n Fractional status = %d \n", status);

     ////////////////////////////////////////////////////////////
     //                   Explore the TB&B                     //
     //                   		    		       //
     ////////////////////////////////////////////////////////////

     if (fractional) {
          //fprintf(stdout, "\nSolution is fractional\n");
          //fprintf(stdout, "frac_index = %d \t frac_value = %f \n", frac_index, frac_value);

          // Create 2 branches
          node_depth = 1;
          //fprintf(stdout, "list.size = %d\n", list.size);

    	 start = 1;


         do {



         	//////////////// Special Case ////////
               	if (max_numnodes == 1){
               	      //fprintf(stdout, "\n **node_id = %d \n", node_id);
               	      //fprintf(stdout, "\n **Nodes explored equals max_numnodes = %d \n", max_numnodes);
               	      //fprintf(stdout, "\n\ncomputing nu and gamma...\n");
               	      // get nu
                      // get gamma
                      // Store the nu and gammas
                       bb_nu[*num_termnodes] = nu;
                       for (j = 0; j < numcols_master; j++) {
                           bb_gamma[*num_termnodes][j] = gamma[j];
                       }
                       *num_termnodes += 1; // count this node as a terminal node
               	       break;
               	 }
               	//////////////////////////////////////
              fathomed = 0;
              ////////////////////////////////////////
              //	 Left branch -y_i >= 0      //
              ////////////////////////////////////////
              if (start == 1) { // Initial iteration
                   varind[0] = frac_index;

                   // Left branch -y_i >= 0
                   varval[0]  = -1.0;
                   varbnds[0] = floor(frac_value);
                   //fprintf(stdout, "\nLeft branch\n");
                   //fprintf(stdout, "varind[0] = %d \t varval[0] = %f \t varbnds[0] = %f\n",
                   //                 varind[0], varval[0], varbnds[0]);
                   status = BBaddbranchconstrs(env, lp_sub, 1, varind, varval, varbnds);
                   if ( status ) {
                        fprintf(stderr, "\nBBdobranchbound():\n");
                        fprintf (stderr, "Failed to add branching constraints to subproblem CPLEX LP object.\n");
                        fprintf (stderr, "status: %d\n", status);
                        return status;
                   }

               } else { // Get node from from of list
                     objval = list.front->next->objval;
                     //fprintf(stdout, "list.size = %d\n", list.size);
                     //fprintf(stdout, "list.front->next->objval= %f\n", list.front->next->objval);

                     if (objval >= *bestbound) { //  Fathom by bound
                          fathomed = 1;
                          //fprintf(stdout, "\n\nFATHOMING NODE %d BY BOUND!\n", list.front->next->node_id);
                          //fprintf(stdout, "\nobjval >= *bestbound:%f >= %f\n", objval, *bestbound);
                          //fprintf(stdout, "\nGetting nu and gamma...\n");
                          // Get and store nu and gamma
                          bb_nu[*num_termnodes] = list.front->next->nu;
                          for (j = 0; j < numcols_master; j++) {
                               bb_gamma[*num_termnodes][j] = list.front->next->gamma[j];
                          }
                          *num_termnodes += 1; // count this node as a terminal node

                          //fprintf(stdout, "\n\n *************REMOVING NODE TO LIST!\n");
                          removefrontnode(&list); // Delete the node from the list
                     } else { // Need to explore this node
                          j = list.front->next->node_depth;
                          // fprintf(stdout, "list.front->next->node_depth = %d\n", list.front->next->node_depth);
                          // Copy node information
                          cur_node_depth = list.front->next->node_depth;
                          frac_index     = list.front->next->frac_index;

                          frac_value     = list.front->next->frac_value;
                          //fprintf(stdout, "list.front->next->frac_index= %d\n", list.front->next->frac_index);
                          //fprintf(stdout, "list.front->next->frac_value= %f\n", list.front->next->frac_value);
                          for (i = 0; i < j; i++) {
                             // varind[i]  = list.front->next->varind[i];
                              varval[i]  = list.front->next->varval[i];
                              varbnds[i] = list.front->next->varbnds[i];
                              varind[i]  = list.front->next->varind[i];
                              //fprintf(stdout, "list.front->next->varind[i] = %d\n", list.front->next->varind[i]);
                              //fprintf(stdout, "list.front->next->varval[i]   = %f\n", list.front->next->varval[i]);
                              //fprintf(stdout, "list.front->next->varbnds[i]  = %f\n\n", list.front->next->varbnds[i]);
                          }
                          //fprintf(stdout, "\nNode depth from list = %d\n", j);
                          //fprintf(stdout, "list.front->next->frac_index = %d\n", list.front->next->frac_index);
                          //fprintf(stdout, "list.front->next->frac_value = %f\n", list.front->next->frac_value);
                          //fprintf(stdout, "list.front->next->node_id    = %d\n", list.front->next->node_id);

                          varind[j]  = frac_index;
                          varval[j]  = -1.0;
                          varbnds[j] = floor(frac_value);;
                          //fprintf(stdout, "\nLeft branch\n");
                          //for (i = 0; i <= j; i++) {
                          //     fprintf(stdout, "varind[%d] = %d \t varval[%d] = %f \t varbnds[%d] = %f\n",
                          //                i, varind[i], i, varval[i], i, varbnds[i]);
                          //}

                          status = BBaddbranchconstrs(env, lp_sub, j+1, varind, varval, varbnds);
                          if ( status ) {
                               fprintf(stderr, "\nBBdobranchbound():\n");
                               fprintf (stderr, "Failed to add branching constraints to subproblem CPLEX LP object.\n");
                               fprintf (stderr, "status: %d\n", status);
                               return status;
                           }
                           //fprintf(stdout, "\nDone adding branching constrs to sub lp\n");
                       } // end if(objval >= *bestbound)
               } // End if (start == 1)/else



               ////////////////////////////////////////////////////
               if (!(fathomed)) { // Explore this node
                    // Solve this node LP
                    //fprintf(stdout, "\nOptimizing node %d subproblem scenario index %d  \n", node_id, scenario);
    		    cputime_start = clock();
                    status = CPXdualopt (env, lp_sub);
                    //status = CPXlpopt (env, lp_sub);
                    cputime_stop = clock();
                    *cputime_sublp += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
                    //fprintf(stdout, "\nDone Re-ptimizing node %d subproblem scenario index %d  \n", node_id, scenario);
                    if ( status ) {
                        fprintf(stderr, "\nBBdobranchbound():\n");
                        fprintf (stderr, "Failed to optimize scenario %d subproblem.\n", scenario);
                        fprintf(stderr, "Error %d", status);
                        return status;
                     }
                     //solstat = CPXgetstat (env, lp_sub);
                     //fprintf (stdout, "Solution status %d.\n", solstat);
                     status  = CPXgetobjval (env, lp_sub, &objval);
                     if ( status ) {
                          fprintf(stderr, "\nBBdobranchbound():\n");
                          fprintf (stderr,"Failed to obtain objective value for node subproblem scenario index %d .\n", scenario);
                          fprintf(stderr, "Error code: %d", status);
                          return status;
                      }
                      //fprintf (stdout, "objval %f.\n", objval);

                      status = CPXgetx (env, lp_sub, solnY, 0, cur_numcols-1);
                      if ( status ) {
                           fprintf(stderr, "\nBBdobranchbound():\n");
                           fprintf (stderr, "Failed to obtain subproblem primal solution.\n");
                           fprintf (stderr, "status: %d\n", status);
                           return status;
                       }
                       status = CPXgetdj (env, lp_sub, redcosts, 0, cur_numcols-1);
               		if ( status ) {
                    	     fprintf(stderr, "\nBBdobranchbound():\n");
                    	     fprintf (stderr, "Failed to obtain subproblem reduced costs.\n");
                    	     fprintf (stderr, "status: %d\n", status);
                    	     return status;
                	}
                	cur_numrows = CPXgetnumrows (env, lp_sub);
                	status = CPXgetpi (env, lp_sub, duals, 0, cur_numrows-1);
                	if ( status ) {
                   	     fprintf(stderr, "\nBBdobranchbound():\n");
                   	     fprintf (stderr, "Failed to obtain subproblem dual solution.\n");
                   	     fprintf (stderr, "status: %d\n", status);
                   	     return status;
                	}

                	status = CPXgetrhs (env, lp_sub, rhs, 0, cur_numrows-1);
      			if ( status ) {
          		      fprintf(stderr, "\nBBdobranchbound():\n");
          		      fprintf (stderr, "Failed to obtain subproblem rhs.\n");
          		      fprintf (stderr, "status: %d\n", status);
          		      return status;
      			}


      			// Delete the branching rows from subproblem lb object
                	status = CPXdelrows (env, lp_sub, orig_numrows, cur_numrows-1);
                	if ( status ) {
                   	     fprintf(stderr, "\nBBdobranchbound():\n");
                   	     fprintf (stderr, "Failed to delele subproblem rows.\n");
                   	     fprintf (stderr, "status: %d\n", status);
                   	     return status;
                	}

      		/*	/////////////////////////////////////////////
                	fprintf(stdout, "\nPrimal Solution: Left Node\n");
                	for (j = 0; j < cur_numcols; j++) {
                             fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], solnY[j]);
                	}
                	//fprintf(stdout, "\nDual Solution:\n");
                	for (j = 0; j < cur_numrows; j++) {
                             fprintf (stdout, " Row %d: = %6.6f\n", j, duals[j]);
                	}
                	//fprintf(stdout, "\nReduced Costs:\n");
                	for (j = 0; j < cur_numcols-1; j++) {
                             fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], redcosts[j]);
                	}

      			//fprintf(stdout, "\nrhs:\n");
      			for (j = 0; j < cur_numrows; j++) {
           		      fprintf (stdout, " Row %d: = %6.6f\n", j, rhs[j]);
      			}
      			fprintf(stdout, "\nOriginal subproblem RHS:\n");
      			for (j = 0; j < orig_numrows; j++) {
           		      fprintf (stdout, " stochdataPtr[%d][%d]: = %6.6f\n", scenario, j, stochdataPtr->rhs[scenario][j]);
      			}

      		*/	////////////////////////////////////////////////////////////



      			//////////////////////////////////////////////////////////////////////////////
      			//                       Compute nu and gamma				    //
      			//////////////////////////////////////////////////////////////////////////////
      			total_redcosts = getTotalDualsDueToReducedCosts(subprobPtr, redcosts, solnY);
      			//fprintf(stdout, "\nCalling BBcomputegammanu:\n");
                        BBcomputegammanu(stochdataPtr, subprobPtr, scenario, orig_numrows, duals,
                                         gamma, &nu, random_T, total_redcosts);
                        *num_nodesexplored += 1; // count this node as exlored

                        //printf("\n Before adding branching duals: nu = %3.3g\n", nu);
                        for (j = orig_numrows; j < cur_numrows; j++)
                            nu += rhs[j]*duals[j];

                        //printf("\n\n nu = %3.3g\n", nu);
                        //for (j = 0; j < numcols_master; j++) {
                        ////}
                        //printf("\n");

                        // Store nu and gamma values
      			//bb_nu[*num_nodesexplored] = nu;
      			//for (j = 0; j < numcols_master; j++) {
          		//     bb_gamma[*num_nodesexplored][j] = gamma[j];
      			//}

      			//////////////////////////////////////////////////////////////////////////////
      			//Write sub prob to file in lp format
                        #ifdef DEBUG_FUNC
                              fprintf (stdout, "orig_numrows  = %d\n", orig_numrows);
                              fprintf (stdout, "cur_numrows   = %d\n", cur_numrows);
                              fprintf (stdout, "cur_numrows-1 = %d\n", cur_numrows-1);
   	                      status = CPXwriteprob(env, lp_sub, "SUB.lp", NULL);
   	                      if ( status ) {
   		                  fprintf (stderr, "BBdo():\n");
      		                  fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		                  exit(1);
   	                       }
                        #endif

                	// Get max fractional index and value
                	fractional = BBgetMaxFractionIndex(solnY, cur_numcols-1, &lfrac_index, &lfrac_value,
                                                   subprobPtr->ctype);
                	//fprintf(stdout, "\n Fractional status = %d \n", fractional);

                        node_id++;
                        ////////////////
               		if (node_id == max_numnodes){
               		        //fprintf(stdout, "\n **node_id = %d \n", node_id);
               		        //fprintf(stdout, "\n **Nodes explored equals max_numnodes = %d \n", max_numnodes);
               		        if (!(fractional)) {
               		             //fprintf(stdout, "\nINTEGRAL SOLUTION\n");
                      		     if (objval < *bestbound) {
                            		  *bestbound = objval;
                            		  /// store incumbent solution
                                          *fract_sub = 0;    // set flag to false
                                          for (j = 0; j < cur_numcols; j++) {
                                              solnY_sub[j] = solnY[j];
                                          }
                      		      }
                       		      //fprintf(stdout, "\n\nSolution is integral at node %d!\n", node_id);
                       		      //fprintf(stdout, "\n\n*bestbound = %f\n", *bestbound);
               		        }
               		        // Get and store nu and gamma
                                bb_nu[*num_termnodes] = nu;
                                for (j = 0; j < numcols_master; j++) {
                                     bb_gamma[*num_termnodes][j] = gamma[j];
                                }
                                *num_termnodes += 1; // count this node as a terminal node

               		        ////// NEED TO CHECK THIS AGAIN //////// 3/8/04
                       		if (start != 1) {
                       		      //fprintf(stdout, "\n\n *************REMOVING NODE FROM LIST!\n");
                       		      removefrontnode(&list); // Delete the node from the list
                       		}
               		        break;
               		}
               		 /////////////
                        if (fractional) { // Check if we can fathom this node by bound
                             //fprintf(stdout, "\nSolution is fractional\n");
                             //fprintf(stdout, "lfrac_index = %d \t lfrac_value = %f \n", lfrac_index, lfrac_value);
                             if (objval >= *bestbound) { //  Fathom by bound
                                 //fprintf(stdout, "\n\nFATHOMING NODE BY BOUND!\n");
                                 //fprintf(stdout, "\n\nobjval >= *bestbound:%f >= %f\n", objval, *bestbound);
                                 // Get and store nu and gamma
                                 bb_nu[*num_termnodes] = nu;
                                 for (j = 0; j < numcols_master; j++) {
                                      bb_gamma[*num_termnodes][j] = gamma[j];
                                 }
                                 *num_termnodes += 1; // count this node as a terminal node

                                 ////// NEED TO CHECK THIS AGAIN /////// 3/8/04
                                 if (start != 1) {
                                     // Delete old node from the front of the list
                                     //fprintf(stdout, "\n\n *************REMOVING NODE TO LIST!\n");
                                     removefrontnode(&list);
                                 }
                              } else { // Add problem to list
                                     //fprintf(stdout, "\n\nADDING LEFT NODE TO LIST at start = %d\n", start);
                                     if (start == 1) { // Initial iteration
                                         newnode = createnode(objval, node_depth, max_numnodes, varind, varval, varbnds,
                                                              numcols_master, nu, gamma);
                                     } else {
                                         //node_depth = list.front->next->node_depth + 1;
                                         newnode = createnode(objval, cur_node_depth+1, max_numnodes,
                                                              varind, varval, varbnds, numcols_master, nu, gamma);
                                         //fprintf(stderr, "\n *****REMOVING NODE FROM FRONT OF LIST******:\n");

                                         // Delete old node from the front of the list
                                         removefrontnode(&list);
                                     }
                                     newnode_created = 1;   // Flag node has been created, will need to free it

                                     newnode->frac_index = lfrac_index;
                                     newnode->frac_value = lfrac_value;
                                     newnode->node_id    = node_id;

                                     // Add new node now
                                     if (!add(&list, newnode))
                                     {
                                          fprintf(stderr, "\nBBdobranchbound():\n");
                                          fprintf (stderr, "Unable to allocate space for the new node");
                                          return 1;
                                      }

                                } // End inner if/else
                          } else { // Fathom node: Integral solution.

                      	         //fprintf(stdout, "\nFATHOM NODE BY INTEGRAL SOLUTION\n");
                      	         if (objval < *bestbound) {
                                      *bestbound = objval;
                                      /// store incumbent solution
                                      *fract_sub = 0;    // set flag to false
                                      for (j = 0; j < cur_numcols; j++) {
                                           solnY_sub[j] = solnY[j];
                                      }
                                  }
                                  // Get and store nu and gamma
                                  bb_nu[*num_termnodes] = nu;
                                  for (j = 0; j < numcols_master; j++) {
                                       bb_gamma[*num_termnodes][j] = gamma[j];
                                  }
                                  *num_termnodes += 1; // count this node as a terminal node

                                  //fprintf(stdout, "\n\nSolution is integral at node %d!\n", node_id);
                                  //fprintf(stdout, "\n\n*bestbound = %f\n", *bestbound);
                                  if (start != 1 && !isempty(&list)) {
                                     //fprintf(stderr, "\n *****REMOVING NODE FROM FRONT OF LIST******:\n");

                                     // Delete old node from the front of the list
                                     removefrontnode(&list);
                                   }
                          } // end if(fractional)/else



              		//////////////////////////////////
                        //    Right branch y_i >= 1     //
                        //////////////////////////////////
                	if (start == 1) { // Initial iteration
                   	      	varind[0] = frac_index;

                   	     	// Right branch y_i >= 1
                             	varval[0]  = 1.0;
                             	varbnds[0] = ceil(frac_value);
                             	//fprintf(stdout, "\n*Right branch\n");
                             	//fprintf(stdout, "varind[0] = %d \t varval[0] = %f \t varbnds[0] = %f\n",
                                //    		 varind[0], varval[0], varbnds[0]);
                   		status = BBaddbranchconstrs(env, lp_sub, 1, varind, varval, varbnds);
                   		if ( status ) {
                        		fprintf(stderr, "\nBBdobranchbound():\n");
                        		fprintf (stderr, "Failed to add branching constraints to subproblem CPLEX LP object.\n");
                        		fprintf (stderr, "status: %d\n", status);
                        		return status;
                   		}
                   		//fprintf(stdout, "\nEnd Right branch var setttings\n");
               		} else { // Use node from front of list
                          	//j = list.front->next->node_depth;
                          	j = cur_node_depth;
                          	//fprintf(stdout, "Node depth from list = %d\n", j);
                          	varind[j]  = frac_index;
                          	varval[j]  = 1.0;
                          	varbnds[j] = ceil(frac_value);
                          	//fprintf(stdout, "\nRight branch\n");
                          	//for (i = 0; i <= j; i++) {
                                //      fprintf(stdout, "varind[%d] = %d \t varval[%d] = %f \t varbnds[%d] = %f\n",
                                //           i, varind[i], i, varval[i], i, varbnds[i]);
                                 //}
                          	status = BBaddbranchconstrs(env, lp_sub, j+1, varind, varval, varbnds);
                          	if ( status ) {
                               		fprintf(stderr, "\nBBdobranchbound():\n");
                               		fprintf (stderr, "Failed to add branching constraints to subproblem CPLEX LP object.\n");
                               		fprintf (stderr, "status: %d\n", status);
                               		return status;
                           	}

                        } // End if (start == 1)/else

                	// Solve this node LP
                	//fprintf(stdout, "\nOptimizing node %d subproblem scenario index %d  \n", node_id, scenario);
    			cputime_start = clock();
    			status = CPXdualopt (env, lp_sub);
    			//status = CPXlpopt (env, lp_sub);
    			cputime_stop = clock();
    			*cputime_sublp += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
    			//fprintf(stdout, "\nDone Re-ptimizing node %d subproblem scenario index %d  \n", node_id, scenario);
                	if ( status ) {
                   		fprintf(stderr, "\nBBdobranchbound():\n");
                   		fprintf (stderr, "Failed to optimize scenario %d subproblem.\n", scenario);
                   		fprintf(stderr, "Error %d", status);
                   		return status;
                	}
                	//solstat = CPXgetstat (env, lp_sub);
                	//fprintf (stdout, "Solution status %d.\n", solstat);
                	status  = CPXgetobjval (env, lp_sub, &objval);
                	if ( status ) {
                   		fprintf(stderr, "\nBBdobranchbound():\n");
                   		fprintf (stderr,"Failed to obtain objective value for node subproblem scenario index %d .\n", scenario);
                   		fprintf(stderr, "Error code: %d", status);
                   		return status;
                	}
                	//fprintf (stdout, "objval %f.\n", objval);

                	status = CPXgetx (env, lp_sub, solnY, 0, cur_numcols-1);
                	if ( status ) {
                   		fprintf(stderr, "\nBBdobranchbound():\n");
                   		fprintf (stderr, "Failed to obtain subproblem primal solution.\n");
                   		fprintf (stderr, "status: %d\n", status);
                   		return status;
                 	}
                 	status = CPXgetdj (env, lp_sub, redcosts, 0, cur_numcols-1);
                 	if ( status ) {
                      		fprintf(stderr, "\nBBdobranchbound():\n");
                      		fprintf (stderr, "Failed to obtain subproblem reduced costs.\n");
                      		fprintf (stderr, "status: %d\n", status);
                      		return status;
                 	}
                 	cur_numrows = CPXgetnumrows (env, lp_sub);
                 	status = CPXgetpi (env, lp_sub, duals, 0, cur_numrows-1);
                 	if ( status ) {
                   		fprintf(stderr, "\nBBdobranchbound():\n");
                   		fprintf (stderr, "Failed to obtain subproblem dual solution.\n");
                   		fprintf (stderr, "status: %d\n", status);
                   		return status;
                	}
                	status = CPXgetrhs (env, lp_sub, rhs, 0, cur_numrows-1);
      			if ( status ) {
          			fprintf(stderr, "\nBBdobranchbound():\n");
          			fprintf (stderr, "Failed to obtain subproblem rhs.\n");
          			fprintf (stderr, "status: %d\n", status);
          			return status;
      			}
                /*	//////////////////////////////////////////////////
                	fprintf(stdout, "\nPrimal Solution: Right Node\n");
                	for (j = 0; j < cur_numcols; j++) {
                    		fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], solnY[j]);
                	}
                	fprintf(stdout, "\nDual Solution:\n");
                	for (j = 0; j < cur_numrows; j++) {
                    		fprintf (stdout, " Row %d: = %6.6f\n", j, duals[j]);
                	}
                	fprintf(stdout, "\nReduced Costs:\n");
                	for (j = 0; j < cur_numcols-1; j++) {
                     		fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], redcosts[j]);
                	}

      			fprintf(stdout, "\nRHS:\n");
      			for (j = 0; j < cur_numrows; j++) {
           			fprintf (stdout, " Row %d: = %6.6f\n", j, rhs[j]);
      			}
      			fprintf(stdout, "\nOriginal subproblem RHS:\n");
      			for (j = 0; j < orig_numrows; j++) {
           		      fprintf (stdout, " stochdataPtr[%d][%d]: = %6.6f\n", scenario, j, stochdataPtr->rhs[scenario][j]);
      			}
      		*/
      			//////////////////////////////////////////////////////////////////////////////
      			//                       Compute nu and gamma				    //
      			//////////////////////////////////////////////////////////////////////////////

      			fprintf(stdout, "\norig_numrows: %d\n", orig_numrows);
      			total_redcosts = getTotalDualsDueToReducedCosts(subprobPtr, redcosts, solnY);
				   if (total_redcosts <= 1e-6)
					   printf("hello");
      			printf("\n\n total_redcosts = %3.4f\n", total_redcosts);

                        BBcomputegammanu(stochdataPtr, subprobPtr, scenario, orig_numrows, duals,
                                         gamma, &nu, random_T, total_redcosts);
                        *num_nodesexplored += 1; // count this node as exlored

                        printf("\n\n Before: nu = %3.4f\n\n", nu);
                        for (j = orig_numrows; j < cur_numrows; j++)
                            nu += rhs[j]*duals[j];

                        //printf("\n\n nu = %3.4f\n", nu);
                        //for (j = 0; j < numcols_master; j++) {
                        //       printf(" gamma[%d] = %3.6f\n", j, gamma[j]);
                        //}
                        //printf("\n\n");

                        // Store nu and gamma values
      			//bb_nu[*num_nodesexplored] = nu;
      			//for (j = 0; j < numcols_master; j++) {
          		//     bb_gamma[*num_nodesexplored][j] = gamma[j];
      			//}

      			//////////////////////////////////////////////////////////////////////////////


                	// Delete the rows
                	status = CPXdelrows (env, lp_sub, orig_numrows, cur_numrows-1);
                	if ( status ) {
                   		fprintf(stderr, "\nBBdobranchbound():\n");
                   		fprintf (stderr, "Failed to delele subproblem rows for right node lp.\n");
                   		fprintf (stderr, "status: %d\n", status);
                   		return status;
                	}
                	// Get max fractional index and value
                	fractional = BBgetMaxFractionIndex(solnY, cur_numcols-1, &rfrac_index, &rfrac_value,
                                                   	   subprobPtr->ctype);
               		//fprintf(stdout, "\n Fractional status = %d \n", fractional);


               		 node_id++;
               		 ////////////////
               		  if (node_id == max_numnodes){
               		        //fprintf(stdout, "\n **node_id = %d \n", node_id);
               		        //fprintf(stdout, "\n ***Nodes explored equals max_numnodes = %d \n", max_numnodes);
               		        if (!(fractional)) {
               		             //fprintf(stdout, "\nINTEGRAL SOLUTION\n");
                      		     if (objval < *bestbound) {
                            		  *bestbound = objval;
                            		  /// store incumbent solution
                                          *fract_sub = 0;    // set flag to false
                                          for (j = 0; j < cur_numcols; j++) {
                                              solnY_sub[j] = solnY[j];
                                          }
                      		      }
                       		      //fprintf(stdout, "\n\nSolution is integral at node %d!\n", node_id);
                       		      //fprintf(stdout, "\n\n*bestbound = %f\n", *bestbound);
               		        }
               		        // Store nu and gamma
                                bb_nu[*num_termnodes] = nu;
                                for (j = 0; j < numcols_master; j++) {
                                    bb_gamma[*num_termnodes][j] = gamma[j];
                                }
                                *num_termnodes += 1; // count this node as a terminal node

               		        break;
               		  }
               		 /////////////

               		 //printf("\n\n *nu = %3.4f\n", nu);
                         //for (j = 0; j < numcols_master; j++) {
                         //     printf(" *gamma[%d] = %3.6f\n", j, gamma[j]);
                         //}
                         //printf("\n\n");


                	if (fractional) { // Check if we can fathom this node by bound
                     		//fprintf(stdout, "\nSolution is fractional\n");
                     		//fprintf(stdout, "rfrac_index = %d \t rfrac_value = %f \n", rfrac_index, rfrac_value);
                     		//nu = -999;
                     		//for (j = 0; j < numcols_master; j++) {
                             	//	gamma[j] = 0.0;
                     		//}
                     		if (objval >= *bestbound) { //  Fathom by bound

                          		//fprintf(stdout, "\n\nFATHOMING NODE BY BOUND!\n");
                          		//fprintf(stdout, "\n\nobjval >= *bestbound:%f >= %f\n", objval, *bestbound);
                          		// Store nu and gamma
                                        bb_nu[*num_termnodes] = nu;
                                        for (j = 0; j < numcols_master; j++) {
                                              bb_gamma[*num_termnodes][j] = gamma[j];
                                        }
                                        *num_termnodes += 1; // count this node as a terminal node
                     		} else { // Add problem to list
                          		//fprintf(stdout, "\n\nADDING RIGHT NODE TO LIST!\n");

                          		//printf("\n\n **nu = %3.4f\n", nu);
                                        //for (j = 0; j < numcols_master; j++) {
                                        //     printf(" **gamma[%d] = %3.6f\n", j, gamma[j]);
                                        //}
                                        //printf("\n\n");
                          		if (start == 1) { // Initial iteration
                                             newnode = createnode(objval, node_depth, max_numnodes, varind, varval,
                                                          varbnds, numcols_master, nu, gamma);
                                         } else {
                                             //node_depth = list.front->next->node_depth + 1;
                                             newnode = createnode(objval, cur_node_depth+1, max_numnodes,
                                                              varind, varval,  varbnds,numcols_master, nu, gamma);
                                         }
                                        // printf("\n\n DONE ADDING RIGHT NODE TO LIST\n");
                         		newnode->frac_index = rfrac_index;
                          		newnode->frac_value = rfrac_value;
                          		newnode->node_id    = node_id;
                          		if (!add(&list, newnode))
                          		{
                              			fprintf(stderr, "\nBBdobranchbound():\n");
                              			fprintf (stderr, "Unable to allocate space for the new node");
                             			return 1;
                          		}
                          		newnode_created = 1;   // Flag node has been created, will need to free it
                     		 } // End if(objval >= *bestbound)/else
                 	} else { // Fathom node: Integral solution.

                      		//fprintf(stdout, "\nFATHOM NODE BY INTEGRAL SOLUTION\n");
                      		if (objval < *bestbound) {
                                     *bestbound = objval;
                            	     /// store incumbent solution
                                     *fract_sub = 0;    // set flag to false
                                     for (j = 0; j < cur_numcols; j++) {
                                          solnY_sub[j] = solnY[j];
                                     }
                      		 }
                       		//fprintf(stdout, "\n\n Solution is integral at node %d!\n", node_id);
                       		//fprintf(stdout, "  *bestbound = %f\n", *bestbound);
                       		// Get and store nu and gamma
                                bb_nu[*num_termnodes] = nu;
                                for (j = 0; j < numcols_master; j++) {
                                      bb_gamma[*num_termnodes][j] = gamma[j];
                                }
                                *num_termnodes += 1; // count this node as a terminal node

                 	} // end if (fractional)/else

             } // end if (!(fathomed))
             ////////////////////////////////////////////////////



             fprintf(stdout, "\n\tAbout to explore node from list");
             fprintf(stdout, "\n\tlist.size = %d\n", list.size);
             fprintf(stdout, "\tNodes explored: %d\n", node_id);
             if (list.size != 0)
                  fprintf(stdout, "\tlist.front->next->objval: %f\n", list.front->next->objval);

             start = 2;      // Done first iteration

        } while (node_id < max_numnodes && list.size != 0);



     } else { // Integral solution: done! Fathom initial node
          //fprintf(stdout, "\n\nSolution is integral at initial node!\n");

          // Get the nu and gammas
          bb_nu[*num_termnodes] = nu;
          for (j = 0; j < numcols_master; j++) {
              bb_gamma[*num_termnodes][j] = gamma[j];
          }
          *num_termnodes += 1; // count this node as a terminal node


          *bestbound = objval;
          /// store incumbent solution
          *fract_sub = 0;    // set flag to false
          for (j = 0; j < cur_numcols; j++) {
               solnY_sub[j] = solnY[j];
          }

          //printf("\n\nnode_id = %d\n", node_id);
          //printf(" nu = %3.3g\n", nu);
          //for (j = 0; j < numcols_master; j++) {
          //    printf(" gamma[%d] = %3.6f\n", j, gamma[j]);
          //}


          //fprintf(stdout, "\n\nSolution is integral at initial node!\n");
          //fprintf(stdout, "\n\n*bestbound = %f\n", *bestbound);



     } // if (fractional)/else



     // Process and delete active nodes on list
     //fprintf(stdout, "\n\n %d active node(s) on list\n", list.size);
     while (!isempty(&list))
     {
            //fprintf(stdout, "\n\n Processing node %d at front of list\n", list.front->next->node_id);
            //printf("\n\nlist size: %d \t", list.size);
            //printf("list.front->next->objval: %f\n", list.front->next->objval);
            //printf("list.front->next->node_depth: %d\n", list.front->next->node_depth);

            //fprintf(stdout, "\n\nGetting nu and gamma...\n");
            // Get and store nu and gamma
            bb_nu[*num_termnodes] = list.front->next->nu;
            for (j = 0; j < numcols_master; j++) {
                 bb_gamma[*num_termnodes][j] = list.front->next->gamma[j];
            }
            *num_termnodes += 1; // count this node as a terminal node

            //fprintf(stdout, "\n\n *************REMOVING NODE TO LIST!\n");
            removefrontnode(&list);
     }
     //fprintf(stdout, "\n\n *bestbound = %f\n", *bestbound);
     //fprintf(stdout, "\n\n num_termnodes = %d\n", *num_termnodes);



     //Write sub prob to file in lp format
      #ifdef DEBUG_FUNC
   	     status = CPXwriteprob(env, lp_sub, "bbsub.lp", NULL);
   	     if ( status ) {
   		 fprintf (stderr, "BBdobranchbound:\n");
      		 fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		 return 1;
   	     }
      #endif

     // free memory
     //printf("start = %d\n", start);
     //if (newnode_created == 1 && newnode != NULL){
        //printf("newnode != NULL\n");
       // exit(1);
       // free(newnode);
     //}

     //fprintf(stdout, "\n free soln arrays\n");
     free(solnY);
     free(duals);
     free(redcosts);
     free(rhs);
     free(gamma);
     if (varind != NULL)
         free(varind);
     if (varval != NULL)
         free(varval);
     if (varbnds != NULL)
         free(varbnds);

     //printf("\n num_nodesexplored = %d\n", *num_nodesexplored);
     //fprintf(stdout, "\n Exiting BBdobranchbound()\n");

     return status;
 } // End BBdobranchbound



int
BBmalloclpstruct(bblpProb_t *bblpPtr, int nrows_A, int ncols_A, int num_nodes)
/**
 * Allocates memory to B&B reverse polar lp data structure variables
 * @param bblpPtr  pointer to data structure
 * @param nrows_A integer indicating number of stage 1 rows (constraints)
 * @param ncols_A integer indicating number of stage 1 columns (vars)
 * @param num_nodes integer indicating number of explored nodes in the TB&B tree
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
 {
    // initializations
    int status;
    int j;
    int ncols, nrows;

     bblpPtr->ncols  = 2 + ncols_A + num_nodes*(1 + nrows_A);
     bblpPtr->nrows  = 1 + num_nodes*(2 + ncols_A);

    nrows = bblpPtr->nrows;
    ncols = bblpPtr->ncols;
    bblpPtr->cmatspace  = nrows*ncols;  // assumed max number of nonzeros in W


    //printf("mallocBBlpstruct: nrows = %d\n", nrows);
    //printf("mallocBBlpstruct: ncols = %d\n", ncols);

    if (bblpPtr->obj      != NULL || bblpPtr->ctype   != NULL ||
        bblpPtr->sense    != NULL || bblpPtr->rhs     != NULL ||
        bblpPtr->cmatbeg  != NULL || bblpPtr->cmatind != NULL ||
        bblpPtr->cmatind  != NULL || bblpPtr->cmatval != NULL ||
        bblpPtr->lb       != NULL || bblpPtr->ub      != NULL)
    {
           printf("mallocBBlpstruct: bblpPtr->obj != NULL\n", ncols);
           exit(1);
    }

    bblpPtr->obj   = (double*)malloc(ncols*sizeof(double));
    bblpPtr->ctype = (char*)malloc(ncols*sizeof(char));
    bblpPtr->sense = (char*)malloc(nrows*sizeof(char));
    bblpPtr->rhs   = (double*)malloc(nrows*sizeof(double));
    bblpPtr->lb    = (double*)malloc(ncols*sizeof(double));
    bblpPtr->ub    = (double*)malloc(ncols*sizeof(double));

    bblpPtr->cmatbeg = (int*)malloc(ncols*sizeof(int));
    bblpPtr->cmatcnt = (int*)malloc(ncols*sizeof(int));
    bblpPtr->cmatind = (int*)malloc(bblpPtr->cmatspace*sizeof(int));
    bblpPtr->cmatval = (double*)malloc(bblpPtr->cmatspace*sizeof(double));


    if (bblpPtr->obj      == NULL || bblpPtr->ctype   == NULL ||
        bblpPtr->sense    == NULL || bblpPtr->rhs     == NULL ||
        bblpPtr->cmatbeg  == NULL || bblpPtr->cmatind == NULL ||
        bblpPtr->cmatind  == NULL || bblpPtr->cmatval == NULL ||
        bblpPtr->lb       == NULL || bblpPtr->ub      == NULL)
    {
        fprintf(stderr, "Function memAllocC3DualLPproblemStruct(...): \n");
        fprintf(stderr, "Failure to allocate memory to subproblem data arrays\n");
        fprintf(stderr, "Exiting...\n");
        return(1);
    }

    return(0); // successfull return

} //************************** End memAllocRHSLPproblemStruct function ********************




void
BBloadreversepolarLP(bblpProb_t *bblpPtr, masterproblem_t *masterprobPtr, int num_nodesexplored,
                     double *nu, double **gamma, double *solnX, double eta_k)
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
 {
   int i, j, k, n;
   int nrows = bblpPtr->nrows;
   int ncols = bblpPtr->ncols;
   int numrows_A    = masterprobPtr->nrows;
   int numcols_A    = masterprobPtr->ncols-1;       // Exclude the optimality variable
   int col_start, row_start1, row_start2, row_start3;
   int start, stop;

   //fprintf(stdout, "\nloadbbreversepolarLP: nrows_A = %d\n", numrows_A);
   //fprintf(stdout, "\nloadbbreversepolarLP: ncols_A = %d\n", numcols_A);
   //for (i = 0; i < masterprobPtr->ncols; i++)
   //   fprintf(stdout, "solnX[%d] = %f\n", i,solnX[i]);

   //******* Set object coefs ************//
   // delta column
   bblpPtr->obj[0] = -1;
   // sigma_0 column
   bblpPtr->obj[1] = eta_k;
   //sigma_i's
   for (i = 2; i < numcols_A+2; i++)
       bblpPtr->obj[i] = solnX[i-2];

   // tau_0's and tau_1's
   for (i = numcols_A+2; i < ncols; i++)
       bblpPtr->obj[i] = 0;

   //*********** Set the ctype array continuous vars ********//
   for (i = 0; i < ncols; i++)
       bblpPtr->ctype[i] = 'C';

   //*********** Set the sense array to >= ***************//
   for (i = 0; i < nrows; i++)
       bblpPtr->sense[i] = 'G';
   bblpPtr->sense[num_nodesexplored] = 'E';  // Normalization contraint

   //******* Set RHS coefs ************//
   for (i = 0; i < nrows; i++)
       bblpPtr->rhs[i] = 0;
   bblpPtr->rhs[num_nodesexplored] = 1;

   //*********** Set the bounds on each variable ********//
   // delta variable
   bblpPtr->lb[0] = -1*CPX_INFBOUND;
   bblpPtr->ub[0] =  CPX_INFBOUND;

   // sigma_0 and sigma_i's
   for (i = 1; i < numcols_A+2; i++) {
       bblpPtr->lb[i] = -1*CPX_INFBOUND;
       bblpPtr->ub[i] =  CPX_INFBOUND;
   }

   // tau_0's and tau_1's
   for (i = numcols_A+2; i < ncols; i++) {
	   bblpPtr->lb[i] = 0;
       bblpPtr->ub[i] = CPX_INFBOUND;
   }

   //************** Set the constraint matrix coefs in sparse format ***********//

   // Set the column count for cols delta
    bblpPtr->cmatcnt[0] = num_nodesexplored;

   // Set the column count for cols sigma_o
    bblpPtr->cmatcnt[1] = num_nodesexplored;

   // Set the column count for cols sigma_i's
   for (i = 2; i < numcols_A+2; i++) {
       bblpPtr->cmatcnt[i] = num_nodesexplored;
   }
   // Initialize to zero for cols tau_0's to tau_i's.
   for (i = numcols_A+2; i < ncols; i++) {
       bblpPtr->cmatcnt[i] = 0;
   }

   // delta column
   bblpPtr->nzcnt = 0;
   row_start1 = nrows - num_nodesexplored;
   bblpPtr->cmatbeg[0] = 0;
   for (i = 0; i < num_nodesexplored; i++) {
       bblpPtr->cmatval[i] = -1;
       bblpPtr->cmatind[i] = row_start1;
       bblpPtr->nzcnt++;
       row_start1++;
   }

   // sigma_0 column
   bblpPtr->cmatbeg[1] = num_nodesexplored;
   row_start1 = 0;
   for (i = num_nodesexplored; i < 2*num_nodesexplored; i++) {
	   bblpPtr->cmatval[i] = 1;
           bblpPtr->cmatind[i] = row_start1;
           bblpPtr->nzcnt++;
           row_start1++;
   }



   // sigma_i's columns
   row_start1 = num_nodesexplored+1; //skip normalization constraint
   for (i = 2; i < numcols_A+2; i++) {
        bblpPtr->cmatbeg[i] = bblpPtr->nzcnt;
        row_start2 = row_start1;
        for (j = 0; j < num_nodesexplored; j++) {
            bblpPtr->cmatval[bblpPtr->nzcnt] = 1;
            bblpPtr->cmatind[bblpPtr->nzcnt] = row_start2;
            bblpPtr->nzcnt++;
            row_start2 += numcols_A;
        }
        row_start1++;
   }
   //printSparseMatrix(numcols_A+2, bblpPtr->nzcnt, bblpPtr->cmatbeg, bblpPtr->cmatcnt, bblpPtr->cmatind,
   //                    bblpPtr->cmatval, stdout);

   // tau_0's columns
   row_start1 = 0;
   row_start2 = num_nodesexplored+1; // skip normalization constraint
   row_start3 = nrows - num_nodesexplored;
   k = 0;// termnodes count
   for (i = numcols_A+2; i < numcols_A+2+num_nodesexplored; i++) {
          // sigma_0 constraint
	   bblpPtr->cmatbeg[i] = bblpPtr->nzcnt;
	   bblpPtr->cmatval[bblpPtr->nzcnt] = -1;
	   bblpPtr->cmatind[bblpPtr->nzcnt] = row_start1;
	   bblpPtr->cmatcnt[i]++;
	   bblpPtr->nzcnt++;
	   row_start1++;

	   // normalization constraint
	   bblpPtr->cmatval[bblpPtr->nzcnt] = 1;
	   bblpPtr->cmatind[bblpPtr->nzcnt] = num_nodesexplored;
	   bblpPtr->cmatcnt[i]++;
	   bblpPtr->nzcnt++;

	   // sigma_i constraints: store nonzeros only
	   for (j = 0; j < numcols_A; j++){
               if (gamma[k][j] > NONZERO_LB || gamma[k][j] < -NONZERO_LB) {
	           bblpPtr->cmatval[bblpPtr->nzcnt] = -1*gamma[k][j];
	           bblpPtr->cmatind[bblpPtr->nzcnt] = row_start2;
	           bblpPtr->cmatcnt[i]++;
	           bblpPtr->nzcnt++;
	       }
	       row_start2++;
	   }

	   // delta contraint
	   if (nu[k] > NONZERO_LB || nu[k] < -NONZERO_LB){
		   bblpPtr->cmatval[bblpPtr->nzcnt] = nu[k];
	   	   bblpPtr->cmatind[bblpPtr->nzcnt] = row_start3;
	   	   bblpPtr->cmatcnt[i]++;
	   	   bblpPtr->nzcnt++;
	   }
	   row_start3++;
	   k++;
    }
   // tau_i's columns: add the matrix A coefs and b coefs
   //cols start at c = numcols_A+2 to c+num_nodesexplored
   //rols start at r = 3 to r < 4+numcols_A

   // Loop over A nonzeros rows  and extract each row data
   row_start2 = nrows - num_nodesexplored;
   k = 2 + numcols_A + num_nodesexplored; // column index
   for (n = 0; n < num_nodesexplored; n++) {  // Loop for each terminal node in the TB&B Tree
       row_start1 = num_nodesexplored + 1 + (n)*numcols_A;    // skip normalization constraint
       for (i = 0; i < numrows_A; i++) {
           if (i == numrows_A-1 ) {
               start = masterprobPtr->rmatbeg_A[i];
               stop  = masterprobPtr->nzcnt_A;
           } else {
               start = masterprobPtr->rmatbeg_A[i];
               stop  = masterprobPtr->rmatbeg_A[i+1];
           }
           bblpPtr->cmatbeg[k] = bblpPtr->nzcnt;
           for (j = start; j < stop; j++){
               bblpPtr->cmatval[bblpPtr->nzcnt] = -1*masterprobPtr->rmatval_A[j];
               bblpPtr->cmatind[bblpPtr->nzcnt] = row_start1;
               bblpPtr->nzcnt++;      	// count nonzeros
               bblpPtr->cmatcnt[k]++; 	// count nonzeros in this column
               row_start1++;               // increment row count
            } // End inner for loop

            // Add the rhs b[i] nonzero value to this column
            if (masterprobPtr->rhs[i] > NONZERO_LB || masterprobPtr->rhs[i] < -NONZERO_LB){
                 bblpPtr->cmatval[bblpPtr->nzcnt] = masterprobPtr->rhs[i];
                 bblpPtr->cmatind[bblpPtr->nzcnt] = row_start2;
                 bblpPtr->nzcnt++;
                 bblpPtr->cmatcnt[k]++; 	// count nonzeros in this column
            }
            k++; // next column
        } // End for loop
        row_start2++;
  } // End outer n for loop

 } //******************************* loadbbreversepolarLP ************************//


int
BBfreelpmodel(CPXENVptr env, CPXLPptr lp_bb, bblpProb_t *bblpPtr)
/**
 * This frees up the BB lp CPLEX object and data arrays
 * @param env  a pointer to the CPLEX environment object
 * @param lp_bb  a pointer to the CPLEX LP problem object
 * @param bblpPtr  a pointer to the BB data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
 {
   int status = 0;

   /* Free up the problem as allocated by CPXcreateprob*/
   /*if ( lp_bb != NULL ) {
       status = CPXfreeprob (env, &lp_bb);
       if ( status ) {
           fprintf (stderr, "CPXfreeprob lp_bb failed, error code %d.\n", status);
       }
   }*/

   // Free up the memory allocated to the arrays
   if ( bblpPtr->obj != NULL ){
      free (bblpPtr->obj);
      bblpPtr->obj = NULL;
   }

   if ( bblpPtr->ctype != NULL ){
      free (bblpPtr->ctype);
      bblpPtr->ctype = NULL;
   }

   if ( bblpPtr->sense != NULL ){
      free (bblpPtr->sense);
      bblpPtr->sense = NULL;
   }

   if ( bblpPtr->rhs != NULL ){
      free (bblpPtr->rhs);
      bblpPtr->rhs = NULL;
   }

   if ( bblpPtr->lb != NULL ){
      free (bblpPtr->lb);
      bblpPtr->lb = NULL;
   }

   if ( bblpPtr->ub != NULL ){
      free (bblpPtr->ub);
      bblpPtr->ub = NULL;
   }

   if ( bblpPtr->cmatbeg != NULL ){
      free (bblpPtr->cmatbeg);
      bblpPtr->cmatbeg = NULL;
   }

   if ( bblpPtr->cmatcnt != NULL ){
      free (bblpPtr->cmatcnt);
      bblpPtr->cmatcnt = NULL;
   }

   if ( bblpPtr->cmatind != NULL ){
      free (bblpPtr->cmatind );
      bblpPtr->cmatind  = NULL;
   }

   if ( bblpPtr->cmatval != NULL ){
      free (bblpPtr->cmatval );
      bblpPtr->cmatval  = NULL;
   }

   return status;

 } //



void
BBcomputegammanu (stochfile_info *stochdataPtr, subproblem_t *subprobPtr,int scenario,
                  int curr_nrows, double *duals, double *cutCoefs, double *rhsCoef,
                  int random_T, double total_redcosts)
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
 {
   int i, j;
   int curr_row;
   int curr_col;
   int nrows;
   int ncols = stochdataPtr->ncols;
   int start, stop;
   double *coefs = NULL;
   double curr_rhs = 0;
   double prod;

   // Allocate memory
   coefs = (double*)malloc(ncols*sizeof(double));
   if (coefs == NULL){
      fprintf(stderr, "BBcomputegammanu: Cannot allocate memory to coefs: exiting\n");
      exit(1);
   }
   //Initialization
   if (scenario == 0){
      *rhsCoef = 0;
      for (i = 0; i < ncols; i++)
           cutCoefs[i] = 0;
   }
   for (i = 0; i < ncols; i++)
       coefs[i] = 0;

   #ifdef DEUBG_FUNC
      printf("\n\n\n\n Computing Benders cut for SCENARIO %d: \n\n", scenario);
      for (i = 0; i < curr_nrows; i++){
          fprintf(stdout, "duals[%d] = %f * stochdataPtr->rhs[%d][%d] = %f \n",
                       i, duals[i], scenario, i, stochdataPtr->rhs[scenario][i]);
      }
    #endif

   // Accumulate rhs coef
   prod = scalProd(duals, stochdataPtr->rhs[scenario], curr_nrows);

   //fprintf(stdout, "prod = %f \n", prod);

   *rhsCoef = prod + total_redcosts;

   //fprintf(stdout, "*rhsCoef += %f \n",*rhsCoef);



   if (random_T == 0) {  // Constant technology matrix T
      // Compute rho(w) col by col : T matrix is column by column sparse matrix
      #ifdef DEBUG_FUNC
            for (i = 0; i < ncols; i++)
                 printf("subprobPtr->cmatbeg_T[%d] = %d  \n", i, subprobPtr->cmatbeg_T[i]);
      #endif

      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->nzcnt_T;
          } else {
             start = subprobPtr->cmatbeg_T[i];
             stop  = subprobPtr->cmatbeg_T[i+1];
          }
          //fprintf(stdout, "\n\n");
          for (j = start; j < stop; j++){
             //fprintf(stdout, "start = %d stop = %d\n", start, stop);
             curr_row = subprobPtr->cmatind_T[j];
             coefs[j] += duals[curr_row]*subprobPtr->cmatval_T[j];
             #ifdef DEBUG_FUNC
               fprintf(stdout, "coefs[%d] = %f \n", j,coefs[j]);
               fprintf(stdout, "col = %d	cmatval_T[%d] = %f: \n", i, j,
                    subprobPtr->cmatval_T[j]);

               fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %6.6f\n", scenario, curr_row,
                              stochdataPtr->rhs[scenario][curr_row]);
             #endif
          }
      }

   //exit(0);

   } else { //Random Technology matrix T(w)

      // Compute rho(w) col by col : T(w) matrix is column by column sparse matrix
      // This is the T(w) read from the stoch file
      for (i = 0; i < ncols; i++){
          if (i == ncols-1) {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cmatbeg_T[scenario][i+1];
          } else {
             start = stochdataPtr->cmatbeg_T[scenario][i];
             stop  = stochdataPtr->cnzcnt_T[scenario];
          }
          for (j = start; j < stop; j++){
             curr_row = stochdataPtr->cmatind_T[scenario][j];
             coefs[j] += duals[curr_row]*stochdataPtr->cmatval_T[scenario][j];
             fprintf(stdout, "col = %d	val[%d] = %f: \n", i, j, stochdataPtr->cmatval_T[scenario][j]);
          }
      }

   } // End if/else statement


   // Compute rho(w) row by row for the newly added rows
   // The T(w) matrix for the newly added rows is in row-by-row sparse matrix format
   //fprintf(stdout, "\n\nHHHHHHHHHHHHHHHHH nrows > nrows_T: %d > %d \n\n\n", subprobPtr->nrows, subprobPtr->nrows_T);
   if (subprobPtr->nrows > subprobPtr->nrows_T) { // New rows have been added
         //fprintf(stdout, "\n Number of rows added = %d \n", stochdataPtr->rnrows );
         // Continue computing rho(w) row by row due to the T(w) matrix, which is
         // in sparse format
         nrows = stochdataPtr->rnrows;
         for (i = 0; i < nrows; i++) {
             if (i == nrows-1) {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rnzcnt_T[scenario];
              } else {
                  start = stochdataPtr->rmatbeg_T[scenario][i];
                  stop  = stochdataPtr->rmatbeg_T[scenario][i+1];
              }

              curr_row = subprobPtr->nrows_T + i;

              for (j = start; j < stop; j++) {
                  curr_col = stochdataPtr->rmatind_T[scenario][j];
                  #ifdef DEBUG_FUNC
                       fprintf(stdout, "duals[%d] = %f * stochdataPtr->rmatval_T[%d][%d] = %f \n",
                          curr_row, duals[curr_row], scenario, j,
                          stochdataPtr->rmatval_T[scenario][j]);
                  #endif

                  coefs[curr_col] += duals[curr_row]*stochdataPtr->rmatval_T[scenario][j];

                  #ifdef DEBUG_FUNC
                      fprintf(stdout, "coefs[%d] = %f \n", curr_col,coefs[curr_col]);
                      fprintf(stdout, "\n *stochdataPtr->rmatval_T[%d][%d] = %f \n", scenario, j,
                              stochdataPtr->rmatval_T[scenario][j]);
                  #endif

              }
         } // end outer for loop

    } //End if statement

    // Set the weighted coefs
    for (i = 0; i < ncols; i++){
        cutCoefs[i] = coefs[i];
    }

     //****** Print the cut coefs ******/
     #ifdef DEBUG_FUNC
        fprintf(stdout, "computeBendersCutCoefs (): \n");
        fprintf(stdout, "\nCut coefs at scenario %d \n", scenario);
        fprintf(stdout, "rhsCoef = %f \n", *rhsCoef);
  	for (i = 0; i < ncols; i++)
  	    fprintf(stdout, "coefs[%d] = %f \n", i, coefs[i]);
  	fprintf(stdout, "\nWeighted: \n");
  	for (i = 0; i < ncols; i++)
  	    fprintf(stdout, "cutCoefs[%d] = %f \n", i, cutCoefs[i]);

  	fprintf(stdout, " End computeBendersCutCoefs: \n\n");
     #endif

    if ( coefs != NULL ){
      free (coefs);
      coefs = NULL;
   }


 } //******************************* BBcomputegammanu *************************************//


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
 BBaddD2optcuttomaster(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr, double rhscoef,
                      double *coefs)
 {
   int j, nzcnt;
   int status;
   double *rhs;
   char   *sense;
   int    *rmatbeg;
   int    *rmatind;
   double *rmatval;
   int ncols;

   ncols = masterprobPtr->ncols;

   //fprintf(stdout, "\n\n\n 2. addBendersCutToMaster: \n");
   //fprintf(stdout, "\n masterprobPtr->ncols: %d \n", ncols);

   // Allocate memory
   rhs     = (double*)malloc(sizeof(double));
   sense   = (char*)malloc(sizeof(char));
   rmatbeg = (int*)malloc(sizeof(int));
   rmatind = (int*)malloc((ncols+1)*sizeof(int));
   rmatval = (double*)malloc((ncols+1)*sizeof(double));
   if ( rhs == NULL || sense == NULL || rmatbeg == NULL ||
        rmatind == NULL || rmatval == NULL) {
     fprintf (stderr, "\nd2algFuncs:\n");
     fprintf (stderr, "Function addNewRowToLP() \n");
     fprintf (stderr, " Failed to allocate memory to add row vars\n");
     return 1;
   }
   //fprintf(stdout, "\n\n\n 3. addBendersCutToMaster: \n");

   //*********** DERIVE OPTIMALITY CUT **************//
   sense[0] = 'G';
   rmatbeg[0] = 0;

   nzcnt = 0;

   //fprintf(stdout, "\n\n\n 4. addBendersCutToMaster: \n");

   // store the x coefs
   for (j = 0; j < ncols-1; j++) {
       if (coefs[j] > NONZERO_LB || coefs[j] < -NONZERO_LB) {
           rmatval[nzcnt] = coefs[j];   // store nonzeros only
           rmatind[nzcnt] = j;
           nzcnt++;
       }
   }
   // Add the theta or opt column coef
   rmatval[nzcnt] = 1; // store nonzeros only
   rmatind[nzcnt] = ncols-1;
   nzcnt++;

   //for (j = 0; j < ncols; j++) {
   //    fprintf(stderr, "rmatval[%d] = %f	", j, rmatval[j]);
   //    fprintf(stderr, "rmatind[%d] = %d\n", j, rmatind[j]);
   //}

   // Set the rhs
   rhs[0] = rhscoef;

   status = CPXaddrows(env, lp_master, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
                       NULL, NULL);

   if (status) {
      fprintf (stderr, "\nd2algFuncs:\n");
      fprintf (stderr, "Function addNewRowToLP() \n");
      fprintf (stderr, " Failed to add new pi row to subproblem lp object\n");
      fprintf (stderr, " Error code = %d\n", status);
   }

   // Free up the memory allocated to the arrays
   if ( rhs != NULL ){
      free (rhs);
      rhs = NULL;
   }
   if ( sense != NULL ){
      free (sense);
      sense = NULL;
   }
   if ( rmatbeg != NULL ){
      free (rmatbeg);
      rmatbeg = NULL;
   }
   if ( rmatind != NULL ){
      free (rmatind);
      rmatind = NULL;
   }
   if ( rmatval != NULL ){
      free (rmatval);
      rmatind = NULL;
   }

   return status;

 } // ************ BBaddD2optcuttomaster *********************//


/*********************END BRANCH-AND-BOUND FUNCS **************************/
///////////////////////////////////////////////////////



