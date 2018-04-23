 /***********************************************************************************************
 *  Program d2bacMain.c implements the Disjuctive Decomposition with 					*
 *  Branch And Cut algorithm for 2-stage SMIP problems.													*
 *																																								*
 *         Author: 			Yang Yuan																										*
 *         Date : 			Sep 10, 2006																									*
 *	   		Finished :  	Oct 30, 2006																									*
 *         Revised :  	Feb 20, 2007																									*
 *																																								*
 *			Based on the preivous work by Lewis Ntaimo															*
 ************************************************************************************************/




#include <ilcplex/cplex.h>  /* Bring in the CPLEX function declarations */

/* Bring in the declarations for the string and character functions and malloc
   and math functions */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> 		/* Computation time functions */
#include "d2bacFuncs.h"  /* Prototypes for SMPS reader functions    */
#include "d2bac.h"  		/* Function prototypes                     */
#include "d2bacConfig.h"        /* Constant params that can be configured  */

//#include "orderedlist.h"        /* header file for ordered linked list c file  */

//#define DEBUG_MAIN  // Debugging mode for main  if uncommented
//#define DISPLAY_CUTCOEFS //if uncommented displays common cut coefs

// Initialize the structures for storing input file data
stochfile_info    struct_s;             // STOCH file data structure
masterproblem_t   master_prob_t ;       // Master problem lp data structure
subproblem_t      sub_prob_t ;          // Subproblem lp data  structure
solnMatrix_t      soln_matrix_t;		// Subproblems solution matrix data structure
c3lpProb_t        c3lp_prob_t;	        // Pointer to C3-LP problem lp data  structure
rhslpProb_t       rhsLP_prob_t;	        // Pointer to RHS-LP problem lp data  structure

bblpProb_t       bbLP_prob_t;	        // Pointer to B&B reverse polar-LP problem lp data  structure
/*********************************** START MAIN FUNCTION **************************/

#ifndef  CPX_PROTOTYPE_MIN
int
main (int argc, char *argv[])
#else
int
main (argc, argv)
int   argc;
char  *argv[];
#endif
{
   //****** Output file variables *******
   FILE *fpout;      														// Pointer to DEBUG output file
   FILE *fpSolnOut;  												// Pointer to optimal solution output file
   char debugoutfname[NAMELEN];
   //char solnFileName[NAMELEN];
   char corename_lp[NAMELEN];
   char mastername_lp[NAMELEN];
   char subprobname_lp[NAMELEN];
   char subprobname_mip[NAMELEN];
   char c3probname_lp[NAMELEN];
   char rhsprobname_lp[NAMELEN];
   char bbprobname_lp[NAMELEN];

   stochfile_info *stochdataPtr;	// Structure to hold STOCH file info
   stochdataPtr = &struct_s;		// set its address

   masterproblem_t *masterprobPtr;	// Structure to hold master problem data
   masterprobPtr = &master_prob_t;	// set its address

   subproblem_t *subprobPtr;		// Structure to hold subproblem data
   subprobPtr = &sub_prob_t;		// set its address

   solnMatrix_t *solnPtr;	// Second stage solution sparse matrix struct
   solnPtr = &soln_matrix_t;	// set is address

   c3lpProb_t *c3lpPtr;		// Structure to hold C3-dual LP data
   c3lpPtr = &c3lp_prob_t;	// set its address

   rhslpProb_t *rhslpPtr;		// Structure to hold RHS LP data
   rhslpPtr = &rhsLP_prob_t;	// set its address

   bblpProb_t *bblpPtr;		// Structure to hold B&B reverse polar LP data
   bblpPtr = &bbLP_prob_t;	// set its address

   //****** CORE file LP model info storage variables  *****

   int  surplus;
   int  surplus_row;
   int ncols_master;
   int nrows_master;
   int nrows_sub;
   int ncols_sub;
   int nrows_submip;
   int nrows_core;	// number of rows in CORE file lp
   int ncols_core;	// number of cols in CORE file lp
   int nscens;		// number of subproblem scenarios
   int random_T = 0;  	// A nonzero value indicates that the T(w) is random
    			// while a zero value indicates that T(w) = T for all scenarios
   int zero = 0;
   int random_obj = 0;  // A nonzero value indicates that the T(w) is random
    			// while a zero value indicates that T(w) = T for all scenarios
   int disj_var_floor;  // Floor of the second stage disj var solution
   int disj_var_ceil;   // Ceiling of the second stage disj var solution

 /* Declare and allocate space for the variables and arrays where we will
      store the optimization results including the status, objective value,
      and variable values. */

   int      solstat;
   double   objval;
   double   objval_m; 		// obj value of the master problem
   double   *solnX     = NULL;
   double   *solnXprev     = NULL;
   double   *solnY     = NULL;
   double   *redcosts     = NULL;   // reduced costs
   double   **scen_rhs_rho   = NULL;
   double   *solnC3pi        = NULL;
   double   *solnC3pi_0      = NULL;
   double   *C3secstageQ=NULL     ;
   double   *C3subgrad01 = NULL;
   double   *C3subgrad11 = NULL;
   double   *C3subgrad02_12 = NULL;
    int  *C3addrowind         = NULL;
   double   *C3addrowval         = NULL;
   double   *C3secobj=NULL;
   double   *solnC3Lambda_01 = NULL;
   double   *solnC3Lambda_11 = NULL;
   double   *solnC3Lambda_02_12 = NULL; // solnC3Lambda_02_12[0] = lambda_02;
   																					// solnC3Lambda_02_12[1] = lambda_12;
   int cntmark1, cntmark2;		// counter markers
   int solnCutOff; // Indicates whether scenario fractional soln is cut off or not.
   int scenDropped; // Indicates which scenario soln is drop from the c^3 obj
   int fracSolnCutOff; // Nonezero indicates that the fractional soln is cut off

   CPXENVptr     env = NULL;
   CPXLPptr      lp_master = NULL;	// pointer to lp from core file
   CPXLPptr      lp_sub_clone = NULL;   // a clone of lp from core file
   CPXLPptr      lp_sub = NULL;		// pointer to subproblem lp to be created
   CPXLPptr      lp_submip = NULL;	// pointer to subproblem mip to be created
   //CPXLPptr      lp_c3 = NULL;	// pointer to C^3 lp to be created later
   CPXLPptr      lp_rhs = NULL;		// pointer to RHS lp to be created later

   CPXLPptr lp_bb = NULL;		// pointer to subproblem lp to be created

   int           status;
   int           cur_numcols;
   int           cur_numrows;

   // ***** D2 algorithm optimization variables ***** //

  // Define the lower bound on the second stage scenario subproblems
  // This value is require to make sure that the obj value of all
  // scenario subproblems remains positive as required by the D2 algorithm
  double SUBPROB_LB;

  // Bound on the obj value of the master problem
  double BOUND_V_k = CPX_INFBOUND;
  double C3UpperBnd;
  double C3LowerBnd;

  // Flag for displaying log info to the screen
  int DISPLAY_LOG_INFO = 1;

  int i, j, k, scenario; 			// Counters
  int temp;
   int disj_var = -1;			// Disjunction variable
   int disj_scen;			// Disjunction scenario
   int disj_ind;			// Disjunction variable index in the soln matrix
   int fractional;			// Nonzero value indicates that fractional soln found
   double expObjval;			// Averaged solution for all scenario subproblems
   double expmipObjval;			// Average obj values for all scenario mips
   double lpObjval;			// Averaged solution for all scenario lps
   double totalRedCost;
   					// For testing purposes
   double *aveSubprobSoln;		// Array to store subproblems averaged solution
   double nu_0;				// Coefs for the RHS LP
   double nu_1;
   double nu;
   double nu_max;
   int    nu_nega;
   double gamma_max;
   double gamma_temp;
   double *gamma_index;
   double *gamma;
   double *gamma_0;
   double *gamma_1;
   double **beta;

   double solnRHSdelta;
   double solnRHSsigma_0;
   double   *solnRHSlp = NULL;
   double   *solnRHSsigma_i = NULL;
   int   *delstat = NULL;
   int nodecount;
   int nodeint;
   int nodeleftcnt;
   double curr_lb_V;

   //Computation time variables
   time_t wallclock_start;
   time_t wallclock_stop;
   clock_t cputime_start;
   clock_t cputime_stop;
   double total_wctime;
   double cputime_master = 0;
   double cputime_sublp  = 0;
   double cputime_submip = 0;
   double cputime_c3     = 0;
   double cputime_rhs    = 0;
   double total_cputime;

   int nd2cuts = 0;
   int n_iterations = 0;
   int no_nd2cuts_iters = 0;
   double C3epsino = 1E-1;
   int startToDropScens;
   int fractSoln;

   double prev_lb;
   int MIP_SOLVE;
   int mycount;
   double sumcount1, sumcount2, sumcount;
   double C3addrhs;
   int C3_index;
   int C3_iteration;
   char C3_sense = 'G';
   int nmipsolves = 0;
   int diffMasterSoln;
   int prevDisjVar;
   double curr_percent_gap = 100.0;
   int ADD_LL = 0;
   int ADD_CUT = 0;
   int addCutsToMip = 0;
   int check_obj_chg = 0;
   int writeNextLP;
   double space;

   char soln_filename[NAMELEN];
   char time_filename[NAMELEN];
   char core_filename[NAMELEN];
   char stoch_filename[NAMELEN];

   char rowname_start[NAMELEN];  // First row name in subproblem
   char colname_start[NAMELEN];  // First column name in subproblem

   ///////////////////////
   //   TB&B Variables  //
   ///////////////////////
   int num_nodesexplored; 	// Number of of expored nodes in the TB&T for a given scenario subproblem
   int max_numnodes;	        // Maximum number of nodes to explore in the TB&B tree
   int num_termnodes;           // Number of terminal nodes in the TB&B tree (# of disjunctions)
   double *bb_nu;	// Array to contain nu values for each terminal node in the TB&B tree
   double **bb_gamma;	// 2D Array to contain gamma vectors for each terminal node in the TB&B tree
   double eta_k;	// Current iteration k eta value (expected value of second stage)

   double soln_bbdelta;     // delta value
   double  soln_bbsigma_0;  //  sigma_0 value
   double *soln_bbsigma_i;  //  sigma_i values
   double *soln_bblp;       // TB&B Optimality cut LP (Reverse Polar LP)
   double d2optcut_rhs;     // D^2 "Optimality cut" RHS coef
   double *d2optcut_coefs;  // D^2 "Optimality cut" LHS coefs
   double cputime_bblp  = 0; // The TB&B LP
   int cnt = 0;
   double PERC_START_TBB = 0;   // Start TB&B if percent gap is less than this value
   double ADD_D2CUTS = 0;   // Add D2 cuts if this is set to 1
   int my_max_nodes = 1;
   int fract_sub;
   double bestbound;
   int lp_intsoln = 0;
   char soln_fileseqname[10];

   int cur_numnzs;
   double matrixdensity;
   double prev_percent_gap = 101.0;

   // Start wall clock
   wallclock_start = clock();

   // Initialize CPU time to zero: will only accumulate time for
   // calls to cplex solve routines
   total_cputime  = 0;

   //***** Check the command line arguments *****
   if ( argc != 2 ) {
      usage (argv[0]);
   }

/*
   printf("\nEnter solution output file name sequence number or name> ");
   scanf("%s", soln_fileseqname);
   do {
       printf("\nEnter a nonzero maximum number of nodes to explore in the TB&B Tree> ");
       scanf("%d", &max_numnodes);
       if (max_numnodes <= 0) {
           printf("\nPlease enter a nonzero number!!\n");
       }
   } while (max_numnodes <= 0);

   //////////////////////////////////////////////////////////////////////////////////////////////
   //					Get user options			               //
   //////////////////////////////////////////////////////////////////////////////////////////////
      printf("\n < Set Algorithm Option A. >\n");
      printf("\n  1. Enter 0 for NO Laporte and Louveaux optimality cut in the D2 Algorithm.");
      printf("\n  2. Enter 1 to add Laporte and Louveaux optimality cut during last iteration of");
      printf("\n     the D2 Algorithm.");
      printf("\n  Option A:");
      scanf("%d", &ADD_LL);


      printf("\n < Set Algorithm Option B. >\n");
      printf("\n  1. Enter 0 NOT to add D2 cuts to subproblem MIP for upper bounding.");
      printf("\n  2. Enter 1 to add D2 cuts to subproblem MIP for upper bounding.");
      printf("\n  Option B:");
      scanf("%d", &addCutsToMip);
*/
	ADD_CUT = 1;
	ADD_LL  = 1;
	addCutsToMip = 0;
	check_obj_chg = 0;
	DISPLAY_LOG_INFO = 0;
	PERC_START_TBB = 10;
	max_numnodes = 3;

/*
      printf("\n < Set Algorithm Option C. >\n");
      printf("\n  1. Enter 0 NOT to start MIP solves immediately master obj stays constant.");
      printf("\n  2. Enter 1 to start MIP solves immediately master obj stays constant.");
      printf("\n  Option C:");
      scanf("%d", &check_obj_chg);

      printf("\n < Set Algorithm Option D. >\n");
      printf("\n  1. Enter 0 for NOT displaying algorithm log information.");
      printf("\n  2. Enter 1 to display algorithm log information.");
      printf("\n  Option D:");
      scanf("%d", &DISPLAY_LOG_INFO);

      if (max_numnodes > 1) {
          do {
             printf("\n < Set Algorithm Option E. >\n");
             printf("\n  1. Enter percent gap value [0-100] below which TB&B is initiated>");
             printf("\n  Option E:");
             scanf("%lf", &PERC_START_TBB);
             if (PERC_START_TBB > 100 || PERC_START_TBB < 0) {
                 printf("\n  1. Percent gap value MUST be in the range [0-100]!");
                 printf("\n  1. Re-enter value in the range [0-100]!");
             }
          } while (PERC_START_TBB > 100 || PERC_START_TBB < 0);
      }

*/

      //printf("\n  PERC_START_TBB:  %f\n", PERC_START_TBB);

      strcpy(time_filename,  argv[1]); strcat(time_filename, ".tim");
      strcpy(core_filename,  argv[1]); strcat(core_filename, ".cor");
      strcpy(stoch_filename, argv[1]); strcat(stoch_filename, ".sto");

      // Create output file name as: "d2bac<inputfilename>.out"
      strcpy(soln_filename,  "d2bac");
      strcat(soln_filename, argv[1]);
      strcat(soln_filename, "_1");
      strcat(soln_filename, ".out");

      //printf("\n  time_filename:  %s", time_filename);
      //printf("\n  core_filename:  %s", core_filename);
      //printf("\n  stoch_filename: %s", stoch_filename);
      //printf("\n  soln_filename:  %s\n", soln_filename);

   //////////////////////////////////////////////////////////////////////////////////////////////



  //************* Initialize the CPLEX environment ******************
   env = CPXopenCPLEX (&status);

   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */

   if ( env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      goto TERMINATE;
   }

   //************ Turn on output to the screen **************
   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr,
               "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }

	//************ SMH: set the number of threads to 1 **************
//	status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
//	if ( status ) {
//		fprintf (stderr, "\nd2algMain():\n");
//		fprintf (stderr,
//				 "Failure to set the number of threads, error %d.\n", status);
//		goto TERMINATE;
//	}

  //************* Copy file and lp names ********************
    strcpy(debugoutfname,  "debug.out");
    strcpy(corename_lp,  "core.lp");
    strcpy(mastername_lp,  "master.lp");
    strcpy(subprobname_lp, "subprob.lp");
    strcpy(subprobname_mip, "subprobmip.lp");
    strcpy(c3probname_lp, "c3prob.lp");
    strcpy(rhsprobname_lp, "rhsprob.lp");
    strcpy(bbprobname_lp, "bb.lp");

  // Open debug output file
   fpout = fopen(debugoutfname,"w" );
   if(fpout == NULL) {
       fprintf (stderr, "\nd2algMain():\n");
       fprintf(stderr, "Could not open default DEBUG output file %s for writing!\n", debugoutfname);
       fprintf(stderr, "Terminating...\n");
       return(0);
    }

   // Open solution output file
   fpSolnOut = fopen(soln_filename,"w" );
   if(fpSolnOut == NULL) {
       fprintf (stderr, "\nd2algMain():\n");
       fprintf(stderr, "Could not open solution output file %s for writing!\n", soln_filename);
       fprintf(stderr, "Terminating...\n");
       return(0);
    }

  /////////////////////////////////////////////////////////////////////////
  //    			Load CORE file LP   			 //
  //      								 //
  /////////////////////////////////////////////////////////////////////////

  //****** Load the CORE file into an MIP model for creating the master ******
  lp_master = CPXcreateprob (env, &status, "master_lp");
  if ( lp_master == NULL ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failed to create LP.\n");
      goto TERMINATE;
  }

  //* Now read the file, and copy the data into the created lp
  status = CPXreadcopyprob (env, lp_master, core_filename, NULL);
  if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
     goto TERMINATE;
  }
    //****** Load the CORE file into an MIP model for creating subproblem mip ******
  lp_submip = CPXcreateprob (env, &status, "subprobmip");
  if ( lp_submip == NULL ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failed to create CPX subproblem lp_submip LP.\n");
      goto TERMINATE;
  }

  //* Now read the file, and copy the data into the created lp
  status = CPXreadcopyprob (env, lp_submip, core_filename, NULL);
  if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
     fprintf (stderr, " for the lp_submip LP.\n");
     goto TERMINATE;
  }

   //****** Load the CORE file into an MIP model for creating subproblem lp ******
  lp_sub = CPXcreateprob (env, &status, "subprob");
  if ( lp_sub == NULL ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failed to create CPX subproblem lp_sub LP.\n");
      goto TERMINATE;
  }

  //* Now read the file, and copy the data into the created lp
  status = CPXreadcopyprob (env, lp_sub, core_filename, NULL);
  if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
     fprintf (stderr, " for the lp_sub LP.\n");
     goto TERMINATE;
  }

	CPXwriteprob(env, lp_submip, "subproblem", "lp");

 /* lp_sub_clone = CPXcreateprob (env, &status, "subprob_clone");
  if ( lp_sub_clone == NULL ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failed to create CPX subprob_clone LP.\n");
      goto TERMINATE;
  }
*/
  ////////////////////////////////////////////////////////////
  //     Load TIME file data for splitting CORE lp into     //
  //           master and subproblem LP and MIP             //
  ////////////////////////////////////////////////////////////

  if (DISPLAY_LOG_INFO) {
       fprintf (stdout, "\n Loading TIME file info...\n");
   }
  //************** Load TIME file **************
  status = loadTimeFile(rowname_start, colname_start, time_filename, fpout);
  if (status) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf(stderr, "Failed to read TIME file %s! \n", time_filename);
      fprintf(stderr, "Exiting...\n");
      goto TERMINATE;
  }

  #ifdef DEBUG_MAIN
       fprintf(stderr, " rowname_start = %s\n", rowname_start);
       fprintf(stderr, " colname_start = %s\n", colname_start);
  #endif

 //**************** Get the number of cols and rows in the master and subproblem *********
    nrows_core = CPXgetnumrows (env, lp_master);
    ncols_core = CPXgetnumcols (env, lp_master);
    status = CPXgetcolindex (env, lp_master, colname_start, &ncols_master);
    if ( status ) {
   	fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "Failure to get column index from CORE LP, error %d.\n", status);
      	goto TERMINATE;
    }
    fprintf(stdout, " ncols_master = %d\n", ncols_master);

    status = CPXgetrowindex (env, lp_master, rowname_start, &nrows_master);
    if ( status ) {
   	fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "Failure to get row index from CORE LP, error %d.\n", status);
      	goto TERMINATE;
    }
    nrows_sub = nrows_core - nrows_master;
    ncols_sub = ncols_core - ncols_master;

   //#ifdef DEBUG_MAIN
       fprintf (stdout, "\nd2algMain():\n");
       fprintf(stdout, " nrows_master = %d\n", nrows_master);
       fprintf(stdout, " ncols_master = %d\n", ncols_master);
       fprintf(stdout, " nrows_sub = %d\n", nrows_sub);
       fprintf(stdout, " ncols_sub = %d\n", ncols_sub);
       fprintf(stdout, " nrows_core = %d\n", nrows_core);
       fprintf(stdout, " ncols_core = %d\n", ncols_core);
   //#endif

  //  exit(2);

   ////////////////////////////////////////////////////////////
   //     Get number of scenarios from the STOCH files       //
   //                   		    		     //
   ////////////////////////////////////////////////////////////
   stochdataPtr->nscens = getnumscenarios(stoch_filename);
   //printf("\nstochdataPtr->nscens = %d\n", stochdataPtr->nscens);


  ////////////////////////////////////////////////////////////
  //     Allocate memory for problem structures             //
  //                   		    			    //
  ////////////////////////////////////////////////////////////

  //****************** Allocate memory to STOCH file data structure ********************
  if (DISPLAY_LOG_INFO) {
       fprintf (stdout, "\n Allocating memory to STOCH data struct...\n");
   }
  status = memAllocStochFileStruct(&struct_s, nrows_sub+ncols_sub, ncols_master, ncols_sub);
  if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Failed to allocate memory to STOCH file data structures.\n");
     goto TERMINATE;
  }


  //**************** Allocate memory to master problem data structure *******************
  if (DISPLAY_LOG_INFO) {
       fprintf (stdout, "\n Allocating memory to Master data struct...\n");
   }
  status = memAllocMasterProblemStructs(masterprobPtr, nrows_master, ncols_master);
  if ( status ) {
  	fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "d2algMain(): Failed to intialize master problem data structures!\n");
      	goto TERMINATE;
   }

//**************** Allocate memory to sub problem data structure ***********************
  if (DISPLAY_LOG_INFO) {
       fprintf (stdout, "\n Allocating memory to subproblem struct...\n");
   }
  status = memAllocSubProblemStruct(subprobPtr, nrows_sub, ncols_sub,  ncols_master, stochdataPtr->nscens);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "d2algMain(): Failed to intialize subproblem data structures!\n");
      goto TERMINATE;
   }

   ///////////////////////////////////////////////////////////
   //        Memory  allocations for TB&B Variables         //
   ///////////////////////////////////////////////////////////
   bb_nu   = (double *) malloc (max_numnodes*sizeof(double));
   bb_gamma = (double **) malloc (max_numnodes*sizeof(double *));
   if ( bb_nu == NULL || bb_gamma == NULL) {
        fprintf (stderr, "Memory Allocation failure for bb_nu and bb_gamma arrays. Terminating...\n");
        goto TERMINATE;
   }
   for (i = 0; i < max_numnodes; i++){
        bb_gamma[i] = (double *) malloc ((masterprobPtr->ncols)*sizeof(double));
        if ( bb_gamma[i] == NULL) {
             fprintf (stderr, "Memory Allocation failure for bb_gamma[%d] array.\n", i);
             fprintf(stderr, "Terminating...\n");
             goto TERMINATE;
        }
    }
    soln_bblp       = (double *) malloc ((ncols_master+3)*sizeof(double));
    soln_bbsigma_i  = (double *) malloc (ncols_master*sizeof(double ));
    d2optcut_coefs  = (double *) malloc (ncols_master*sizeof(double ));
    if ( soln_bblp == NULL || soln_bbsigma_i == NULL || d2optcut_coefs == NULL) {
        fprintf (stderr, "Memory Allocation failure for soln_bbsigma_i arrays. Terminating...\n");
        goto TERMINATE;
    }

    ////////////////////////////////////////////////
    //          Allocate memory to BB LP          //
    ////////////////////////////////////////////////


    status = BBmalloclpstruct(bblpPtr, masterprobPtr->nrows, masterprobPtr->ncols, max_numnodes);
    if ( status ) {
        fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "Failure to allocate memory to BB LP struct\n");
      	goto TERMINATE;
    }


   /////////////////////////////////////////////////////////////
   //     Set up master and subproblem MIPs and subproblem LP //
   //                   		    			     //
   /////////////////////////////////////////////////////////////
   //*******Extract master problem lp data and store into master lp structure*******
   status = loadMasterProblemData(env, lp_master, masterprobPtr, nrows_master,
                              ncols_master, nrows_core, ncols_core, fpout);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "d2algMain(): Failed to load master problem data into lp data structures.\n");
      goto TERMINATE;
   }

   status = addOptColToMasterProbLP(env, lp_master, masterprobPtr);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "d2algMain(): Failed to add optimality column to master lp.\n");
      goto TERMINATE;
   }

 //*******Set up subproblem mip data  *******
  status = setupSubProbMip(env, lp_submip, subprobPtr, nrows_master,
                              ncols_master, nrows_sub, ncols_sub, fpout);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "d2algMain(): Failed to load subproblem data into lp data structures.\n");
      goto TERMINATE;
   }

   //printf("ncols_sub    = %d \n", CPXgetnumcols (env, lp_submip));
   //printf("nrows_submip = %d \n", CPXgetnumrows (env, lp_submip));
   nrows_submip = CPXgetnumrows (env, lp_submip);




    //printf("LOADING SUB \n");
   //printf("1. subprobPtr->ncols= %d \n", subprobPtr->ncols);
   //******* Set up subproblem lp data and extract subproblem lp data, ******************
   //******* store into subproblem lp structure and created subproblem ******************
   //******* CPX lp obj lp_sub. Also add new feasibility column to      ******************
   //******* CPX lp obj lp_sub as required by the D^2 algorithm and     ******************
   //******* add explicity binary constraints for all binary variables  ******************
   status = loadSubProblemData(env, lp_sub, subprobPtr, nrows_master, ncols_master,
                               nrows_sub, ncols_sub, stochdataPtr, fpout);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "d2algMain(): Failed to load subproblem data into lp data structures.\n");
      goto TERMINATE;
   }
   //printf("DONE LOADING SUB \n");

   // Make updates: Feasibility column has been added to subproblem  lp
   // as well as explicity binary constraints as required in generating the
   // D^2 cut.
   nrows_sub = CPXgetnumrows (env, lp_sub);
   ncols_sub = CPXgetnumcols (env, lp_sub);

   #ifdef DEBUG_MAIN
        fprintf (stderr, "\nd2algMain():\n");
        printf("nrows_sub    = %d \n", CPXgetnumrows (env, lp_sub));
        printf("ncols_sub    = %d \n", CPXgetnumcols (env, lp_sub));
        printf("stochdataPtr->nrows = %d \n", stochdataPtr->nrows);

        //Write sub prob to file in lp format
   	status = CPXwriteprob(env, lp_sub, subprobname_lp, NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	}
   #endif



  ////////////////////////////////////////////////////////////
  //     Get master and subproblem column and row names     //
  //                   		    			    //
  ////////////////////////////////////////////////////////////

  //****** Get master column names *****
  status = CPXgetcolname(env, lp_master, masterprobPtr->colnames, masterprobPtr->colnamestore,
                         ncols_master*FIELDLEN, &surplus, 0, ncols_master-1);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failure to get master lp colnames, error %d.\n", status);
      fprintf (stderr, " surplus: %d.\n", surplus);
      goto TERMINATE;
  }
  status = CPXgetcolname(env, lp_sub, subprobPtr->colnames, subprobPtr->colnamestore, ncols_sub*FIELDLEN,
                         &surplus, 0, ncols_sub-1);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failure to get subproblem lp colnames, error %d.\n", status);
      fprintf (stderr, " surplus: %d.\n", surplus);
      goto TERMINATE;
  }
	
  status = CPXgetrowname(env, lp_submip, subprobPtr->rownames, subprobPtr->rownamestore, nrows_submip*FIELDLEN,
                         &surplus, 0, nrows_submip-1);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, " Failure to get subproblem lp rownames, error %d.\n", status);
      fprintf (stderr, " surplus: %d.\n", surplus);
      goto TERMINATE;
  }

  #ifdef DEBUG_MAIN
       fprintf (stdout, "\nd2algMain():\n");
       fprintf(stdout, " \n Master column names:\n");
       for (i = 0; i < ncols_master; i++)
            fprintf(stdout, "   masterprobPtr->colnames[%d] = %s\n", i, masterprobPtr->colnames[i]);
       fprintf(stdout, " \n Subproblem column names:\n");
       for (i = 0; i < ncols_sub; i++)
            fprintf(stdout, "   subprobPtr->colnames[%d] = %s\n", i, subprobPtr->colnames[i]);
       fprintf(stdout, " \n Subproblem row names:\n");
       for (i = 0; i < nrows_submip; i++)
            fprintf(stdout, "   subprobPtr->rownames[%d] = %s\n", i, subprobPtr->rownames[i]);
   #endif


   ///////////////////////////////////////////////////////////////////////////
   //     Extract the W matrix from the lp_sub in ROW SPARSE MATRIX FORMAT  //
   //     This is convenient for constructing the C^3 LP.   	            //
   ///////////////////////////////////////////////////////////////////////////
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
       printf("subprobPtr->nrows_mip = %d \n", subprobPtr->nrows_mip);
       printf("subprobPtr->nrows = %d \n", subprobPtr->nrows);
       printf("subprobPtr->nrows_total = %d \n", subprobPtr->nrows_total);
   #endif

   ////////////// Added June 5, 2003 ///////////////////////////////////////////
   //    Delete the explicity binary constraints added to the subproblem      //
   //									      //
   /////////////////////////////////////////////////////////////////////////////
	

//	status = CPXwriteprob(env, lp_sub, "before.lp", "lp");
	status = CPXdelrows(env, lp_sub, subprobPtr->nrows_mip, subprobPtr->nrows_total-1);
//	status = CPXwriteprob(env, lp_sub, "after.lp", "lp");

   if ( status ) {
       fprintf(stderr,"\nMain(...): \n");
       fprintf (stderr,"Failure to delete suproblem lp binary rows, error %d.\n", status);
       return(status);
   }
   // update subproblem rows
   subprobPtr->nrows = subprobPtr->nrows_mip;

   //////////////////////////////////////////////////////////////
   //      Load STOCH file data into the struct stochdataPtr   //
   //        	                                               //
   //////////////////////////////////////////////////////////////
   if (DISPLAY_LOG_INFO) {
        printf("loading StochFile ...\n");
   }
   status = loadStochFile(stochdataPtr, subprobPtr, stoch_filename,
                         &random_T, &random_obj, masterprobPtr->ncols, stdout);
   if (status){
      fprintf (stderr, "\nd2algMain():\n");
      fprintf(stderr, "Failed to read STOCH file %s \n", stoch_filename);
      fprintf(stderr, "Exiting...\n");
      goto TERMINATE;
   }

   //exit(1);

   #ifdef DEBUG_MAIN
       printf("\n\n Done loadStochFile: random_obj = %d \t stochdataPtr->obj_cnt = %d \n",
           random_obj, stochdataPtr->obj_cnt);
       if (random_obj) {
             printf("\n RANDOM OBJECT COEFS:");
             for (i = 0 ; i < stochdataPtr->nscens; i ++ ) {
                   printf("\n");
                   for (j = 0 ; j < stochdataPtr->obj_cnt; j ++ ) {
                        printf("obj_index[%d] = %d \t obj[%d][%d] = %f\n", j,  stochdataPtr->obj_index[j],
                                i, j, stochdataPtr->obj[i][j]);
                   }
             }
        }
   #endif



   #ifdef DEBUG_MAIN
      if (random_T == 1) { // Random techology matrix present

           for (i = 0; i < stochdataPtr->nscens; i++) {
               fprintf(stderr, "\nTechnology Matrix T(%d): \n", i);
               fprintf(stderr, "cmatbeg: %d %d \n", stochdataPtr->cmatbeg_T[i][0],
                                             stochdataPtr->cmatbeg_T[i][1]);
               fprintf(stderr, "\nTechnology Matrix T(%d): \n", i);
  	       printSparseMatrix(masterprobPtr->ncols, stochdataPtr->cnzcnt_T[i], stochdataPtr->cmatbeg_T[i],
                          stochdataPtr->cmatcnt_T[i], stochdataPtr->cmatind_T[i], stochdataPtr->cmatval_T[i], stderr);
            }

  	    fprintf (stderr, "\nd2algMain():\n");
            fprintf(stderr, "random_T = %d \n", random_T);
            fprintf(stderr, "if random_T = 0, then T(w) = T (constant),");
            fprintf(stderr, "otherwise it is random. \n");
       }
   #endif

   //printf("stochdataPtr->nscens = %d \n", stochdataPtr->nscens);



 ///////////////////////////////////////////////////////////////////////////
 //    Allocate memory to the T matrix to be added by the D2 Algorithm    //
 //    This is a row by row sparse matrix .				  //
 ///////////////////////////////////////////////////////////////////////////

 if (DISPLAY_LOG_INFO) {
       fprintf (stdout, "\n Allocating memory to TIME struct...\n");
  }





  //************ Allocate memory to sub problem soln sparse matrix data structure *****
  // Do not store the feas var soln: usually zero anyway, otherwise subproblem is infeasible.

  ///printf("subprobPtr->ncols= %d \n", subprobPtr->ncols);
  status = memAllocSolnMatrix(solnPtr, stochdataPtr->nscens, subprobPtr->ncols-1,
                              masterprobPtr->ncols);
  if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Failed to allocate memory to subproblem solns data structure.\n");
     fprintf (stderr, " Status: %d.\n", status);
     goto TERMINATE;
  }



 //********************* Allocate memory to RHS lp data structure **********************
 status = memAllocRHSLPproblemStruct(rhslpPtr, masterprobPtr->nrows, masterprobPtr->ncols-1);
 if ( status ) {
     fprintf (stderr, "\nd2algMain():\n");
     fprintf (stderr, " Failed to allocate memory to RHS dual lp data structure.\n");
     goto TERMINATE;
  }


 //************************** load data into RHS lp data structure **********************
 loadRHSLPdata(rhslpPtr, masterprobPtr, fpout);

 #ifdef DEBUG_MAIN
    fprintf (stderr, "\n\n*********d2algMain()************:\n");
    fprintf(stderr, "\nRHS LP Constraint Matrix: \n");
    fprintf(stderr, "rows = %d\n cols = %d\n", rhslpPtr->nrows, rhslpPtr->ncols);
    printSparseMatrix(rhslpPtr->ncols, rhslpPtr->nzcnt, rhslpPtr->cmatbeg, rhslpPtr->cmatcnt,
                      rhslpPtr->cmatind, rhslpPtr->cmatval, stderr);
 #endif

 //******************** Create RHS problem CPLEX lp model *********************************
 lp_rhs = CPXcreateprob(env, &status, rhsprobname_lp);
 if ( status ) {
   	fprintf (stderr, "\nd2algMain():\n");
      	fprintf (stderr, "Failure to create CPLEX lp object lp_rhs, error %d.\n", status);
    	goto TERMINATE;
 }
 status = CPXcopylp(env, lp_rhs, rhslpPtr->ncols, rhslpPtr->nrows,
                          CPX_MIN, rhslpPtr->obj, rhslpPtr->rhs, rhslpPtr->sense,
                          rhslpPtr->cmatbeg, rhslpPtr->cmatcnt, rhslpPtr->cmatind,
                          rhslpPtr->cmatval, rhslpPtr->lb, rhslpPtr->ub, NULL);
 if ( status ) {
   	fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "Failure to copy lp_rhs data into CPLEX lp object, error %d.\n", status);
    	goto TERMINATE;
 }
 // Write sub prob to file in lp format
   #ifdef DEBUG_MAIN
        fprintf (fpout, "\nd2algMain():\n");
   	status = CPXwriteprob(env, lp_rhs, rhsprobname_lp, NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write rhs problem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	}
   #endif


  // ********* Reset the r(w) and T(w) or T matrix for those constraints whose ************
  // ********* sense have been changed from 'L' to 'G'				***********
 // resetSubprobRhsAndT(stochdataPtr, subprobPtr, solnX, random_T, fpout);


//************************** Allocate space for subproblem average solution *******************
   aveSubprobSoln = (double *) malloc (3*subprobPtr->ncols*sizeof(double));
   if ( aveSubprobSoln == NULL ) {
      fprintf (stderr, "No memory for subproblem average solution values.\n");
      goto TERMINATE;
   }






   if (random_obj)
       SUBPROB_LB = getrandobjsubproblowerbound(stochdataPtr, subprobPtr);
   else
       SUBPROB_LB = getSubProbLowerBound(subprobPtr);


  //SUBPROB_LB += 4000;

  #ifdef DEBUG_MAIN
     fprintf (fpout, "SUBPROB_LB = %f \n", SUBPROB_LB);
     fprintf (stdout, "\nd2algMain():\n");
     fprintf (stdout, "SUBPROB_LB = %f \n", SUBPROB_LB);
   #endif



   // goto Terminate

  //* Allocate space for master solution *
   solnX = (double *) malloc (masterprobPtr->ncols*sizeof(double));
   if ( solnX == NULL ) {
      fprintf (stderr, "No memory for master solution values.\n");
      goto TERMINATE;
   }
   //* Allocate space for master solution copy *
   solnXprev = (double *) malloc (masterprobPtr->ncols*sizeof(double));
   if ( solnXprev == NULL ) {
      fprintf (stderr, "No memory for master solution copy values.\n");
      goto TERMINATE;
   }

   //* Allocate space for solution *
    solnY = (double *) malloc (subprobPtr->ncols*sizeof(double));
    redcosts = (double *) malloc (subprobPtr->ncols*sizeof(double));
    if ( solnY == NULL || redcosts == NULL) {
          fprintf (stderr, "No memory for subproblem solution and reduced cost arrays.\n");
          goto TERMINATE;
     }


//**************** Optimize the master problem and obtain solution ********************
     //#ifdef DEBUG_MAIN
   	status = CPXwriteprob(env, lp_master, mastername_lp, NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write subproblem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	}
   //#endif
   //printf("Opt Master \n");

   cputime_start = clock();
   status = CPXmipopt (env, lp_master);
   cputime_stop = clock();
   cputime_master += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;

   //printf("End Opt Master \n");

   if ( status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      goto TERMINATE;
   }

   solstat = CPXgetstat (env, lp_master);
   //printf ("Solution status %d.\n", solstat);
   //printf("End sol stat for Opt Master \n");
   status  = CPXgetmipobjval (env, lp_master, &objval_m);

   if ( status ) {
      fprintf (stderr,"Failed to obtain objective value. Error code %d\n", status);
      goto TERMINATE;
   }


   fprintf (fpSolnOut, "Master Obj value = %.10g\n", objval_m);

   //* The size of the problem should be obtained by asking CPLEX what
   //   the actual size is. cur_numcols stores the current number
   //   of columns.

   cur_numcols = CPXgetnumcols (env, lp_master);

   status = CPXgetmipx (env, lp_master, solnX, 0, cur_numcols-1);
   if ( status ) {
      fprintf (stderr, "Failed to obtain master solution.\n");
      goto TERMINATE;
   }
   status = CPXgetmipx (env, lp_master, solnXprev, 0, cur_numcols-1);
   if ( status ) {
      fprintf (stderr, "Failed to obtain master solution copy.\n");
      goto TERMINATE;
   }

   //* Write out the solution *
   fprintf(fpSolnOut, "\nMaster Problem Solution:\n");
   fprintf(stdout, "\nMaster Problem Solution:\n");
   for (j = 0; j < cur_numcols; j++) {
      //printf ( " Column %s:  Value = %17.10g\n", masterprobPtr->colnames[j], solnX[j]);
      //printf ( " Column %d:  Value = %4.4f\n", j, solnX[j]);
      fprintf (fpSolnOut, " %s: = %6.6g\n", masterprobPtr->colnames[j], solnX[j]);
      fprintf (stdout, " %s: = %6.6g\n", masterprobPtr->colnames[j], solnX[j]);
   }

   ////////////////////////////////////////////////
   //	Change subproblem mip to LP		 //
   ////////////////////////////////////////////////

   status = CPXchgprobtype(env, lp_sub, CPXPROB_LP);
   if ( status ) {
      fprintf (stderr, "Failed to change subproblem mip into LP.\n");
      goto TERMINATE;
   }
   // Write sub prob to file in lp format
   //#ifdef DEBUG_MAIN
   	status = CPXwriteprob(env, lp_sub, subprobname_lp, NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write subproblem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	}
   //#endif

   ///////////////////////////////////////////////





/////

//solnX[0] = 0;
//solnX[1] = 1;
//solnX[2] = 0;
//solnX[3] = 0;
//solnX[4] = 0;
//solnX[5] = 1;
//solnX[6] = 0;
//solnX[7] = 1;
//solnX[8] = 0;
//solnX[9] = 0;
/////

   cur_numrows = CPXgetnumrows (env, lp_submip);
   cur_numcols = CPXgetnumcols (env, lp_submip);
   cur_numnzs  = CPXgetnumnz (env, lp_submip);

   matrixdensity = (double)cur_numrows;
   matrixdensity *= cur_numcols;
   matrixdensity = (double)cur_numnzs/matrixdensity;

   fprintf(fpSolnOut, "\n\n cur_numrows =  %d\n\n", cur_numrows);
   fprintf(fpSolnOut, "\n\n cur_numcols =  %d\n\n", cur_numcols);
   fprintf(fpSolnOut, "\n\n cur_numnzs  =  %d\n\n", cur_numnzs);

   fprintf(fpSolnOut, "\n\nConstraint Matrix Density = %f\n\n", matrixdensity);






fprintf (fpout, "SUBPROB_LB = %f \n", SUBPROB_LB);
fprintf (stdout, "SUBPROB_LB = %f \n", SUBPROB_LB);



// ************************************ START ALGORITHM *********************************

prev_lb = -99999;
// objval_m = 0; // This is wrong!!!

MIP_SOLVE = 0;
mycount = 0;
n_iterations = 0;



do {

  fprintf(fpSolnOut, "prev_lb = %f\n", prev_lb);
  fprintf(fpSolnOut, "objval_m = %f\n", objval_m);


  if (check_obj_chg) { // Check if master obj has not changed
      if ((objval_m - prev_lb) < 0.001){
            MIP_SOLVE = 1;
            if (DISPLAY_LOG_INFO) {
                   fprintf (stdout, "\n Initializing scenario subproblem MIP solves...\n");
            }
      }
   }

    fprintf (fpSolnOut, "ITERATION %d\n", n_iterations);
    if (DISPLAY_LOG_INFO) {
        fprintf (stdout, "\n\nITERATION %d\n", n_iterations);
        fprintf(stdout, " [LB, UB]: [%f, %f]  ", objval_m, BOUND_V_k);
        fprintf (stdout,"gap:     %f%%\n", curr_percent_gap);
        fprintf (stdout, " Optimizing scenario subproblem lps...\n");
        fprintf (stdout, " Current accumulated subproblem lps CPU time: %f seconds\n", cputime_sublp);
    }



 //**************************************************************************************
 //*************** Solve subproblems for all the scenarios ******************************


    expObjval = 0;   // Initialize subprob expected obj value
    lpObjval = 0;

  ///////////////////////////////////////////////////////////
  //               SOLVE LPs FOR ALL SCENARIOS             //
  ///////////////////////////////////////////////////////////
  if (MIP_SOLVE == 0) {

      lp_intsoln = 0;
      for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
          // Reset the subproblem RHS based on current stage one X solution
          // The rhs is h(w) = r(w) - T(w)x = subprobPtr->rhsRho[scenario]
	  createsubprobRHS(stochdataPtr, subprobPtr, scenario, solnX, random_T, fpout);
      }
      //printf("2. stochdataPtr->rmatind_T[18][0] = %d\n", stochdataPtr->rmatind_T[18][0]);
      for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
         // printf(" \nSCENARIO: %d\n\n", scenario);
          // Reset the subproblem RHS based on current stage one X solution
          // The rhs is h(w) = r(w) - T(w)x = subprobPtr->rhsRho[scenario]

          //createsubprobRHS(stochdataPtr, subprobPtr, scenario, solnX, random_T, fpout);

          #ifdef DEBUG_MAIN
              fprintf (stdout, "d2algMain():\n");
              for(i=0; i < subprobPtr->nrows;i++){
                    fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", scenario, i,
                           stochdataPtr->rhs[scenario][i]);
               }
           #endif

           #ifdef DEBUG_MAIN
           fprintf (stdout, "d2algMain():\n");
           for(i=0; i < subprobPtr->nrows;i++){
                fprintf(stdout, "index[%d] = %d   subprobPtr->rhsRho[%d][%d] = %f\n",i,
                          subprobPtr->indices_row[i], scenario, i,
                        subprobPtr->rhsRho[scenario][i]);
            }
           #endif

	//for (i = 0; i<subprobPtr->nrows; i++) {
	//	printf ("index[%d] = %d, subprobPtr->rhsRho[%d][%d] = %f\n", i, subprobPtr->indices_row[i], scenario, i,
//subprobPtr->rhsRho[scenario][i]);
	//}
          // Set the RHS for this scenario sub problem
          status = CPXchgrhs(env, lp_sub, subprobPtr->nrows, subprobPtr->indices_row,
                          subprobPtr->rhsRho[scenario]);
          if ( status ) {
              fprintf (stderr, "Failed to set the rhs for the index scenario %d subproblem.\n", scenario);
              fprintf(stderr, "Error: %d\n", status);
              goto TERMINATE;
          }

          // Set the obj for this scenario sub problem if random obj
          if (random_obj) {
        //      fprintf(stdout, " Setting obj for scenario %d  \n", scenario);
        //      fprintf(stdout, " stochdataPtr->obj_cnt %d  \n", stochdataPtr->obj_cnt);
              status = CPXchgobj(env, lp_sub, stochdataPtr->obj_cnt, stochdataPtr->obj_index, stochdataPtr->obj[scenario]);
              if ( status ) {
                  fprintf (stderr,"Failed to change objective value for scenario subproblem index %d .\n", scenario);
                  fprintf(stderr, "Error code: %d", status);
                  goto TERMINATE;
              }
          }
              // Write sub prob to file in lp format
              #ifdef DEBUG_MAIN
   	         status = CPXwriteprob(env, lp_sub, "SP.lp", NULL);
   	         if ( status ) {
   		       fprintf (stderr, "d2algMain():\n");
      		       fprintf (stderr, "Failure to write subproblem lp to file, error %d.\n", status);
      		       goto TERMINATE;
   	          }
               #endif
          //}
          #ifdef DEBUG_MAIN
            fprintf(stdout, "Optimizing subproblem scenario index %d  \n", scenario);
          #endif

          cputime_start = clock();
          //fprintf(stdout, "Start optimizing subproblem scenario index %d  \n", scenario);
          ////if (scenario == 1)
          ////    exit(1);
          //status = CPXdualopt (env, lp_sub);
          status = CPXdualopt (env, lp_sub);
          //fprintf(stdout, "Done optimizing subproblem scenario index %d  \n", scenario);
          if ( status ) {
               fprintf (stderr, "Failed to optimize scenario %d subproblem.\n", scenario);
               fprintf(stderr, "Error %d\n", status);
               goto TERMINATE;
           }
           //status = CPXlpopt (env, lp_sub);
           //status = CPXprimopt (env, lp_sub);
           cputime_stop = clock();
           cputime_sublp += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
           ///////////////////
           //fprintf(stdout, "CPU for scenario %d  lp: %f\n", scenario,(double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC);

            status  = CPXgetobjval (env, lp_sub, &objval);
	   if ( status ) {
               fprintf (stderr,"Failed to obtain objective value for subproblem scenario index %d .\n", scenario);
               fprintf(stderr, "Error code: %d", status);
               goto TERMINATE;
           }

           //fprintf (fpSolnOut, "LP objective value   : %.10g\n", objval);
           // Update subprob expected obj value
           expObjval += objval*stochdataPtr->scenProb[scenario];

           cur_numcols = CPXgetnumcols (env, lp_sub);

           status = CPXgetx (env, lp_sub, solnY, 0, cur_numcols-1);
           if ( status ) {
              fprintf (stderr, "Failed to obtain subproblem LP solution.\n");
              goto TERMINATE;
           }
           status = CPXgetdj (env, lp_sub, redcosts, 0, cur_numcols-1);
           if ( status ) {
              fprintf (stderr, "Failed to obtain subproblem LP reduced costs.\n");
              goto TERMINATE;
           }
           // Store the solution
           storeSubProbSoln(solnPtr, solnY, scenario);

           // Get the dual solution
           cur_numrows = CPXgetnumrows (env, lp_sub);
           //fprintf (stdout, "Getting the duals....\n");
           status = CPXgetpi (env, lp_sub, subprobPtr->duals , 0, cur_numrows-1);
           if ( status ) {
                fprintf (stderr, "Failed to obtain subproblem dual solution.\n");
                fprintf(stderr, "Error code %d \n", status);
                goto TERMINATE;
            }

           #ifdef DEBUG_MAIN
           fprintf (stdout, "d2algMain():\n");
               fprintf(fpSolnOut, "\nSubproblem Problem Dual Solution:\n");
               for (j = 0; j < cur_numfrows; j++)
                    fprintf (fpSolnOut, "row%d: = %6.6g\n", j+1, subprobPtr->duals[j]);
            #endif

            //  Compute Benders' cut coefs
            //fprintf(fpout, "\nComputing Benders cut after optimization:\n");
            totalRedCost = getTotalDualsDueToReducedCosts(subprobPtr, redcosts, solnY);
             computeBendersCutCoefs(stochdataPtr, subprobPtr, scenario, cur_numrows,
                          subprobPtr->duals, masterprobPtr->cutCoefs, &masterprobPtr->rhsCoef,
                          random_T, totalRedCost, fpout);
       } // End scenario for loop
      //**************** Done optimizing all scenario subproblems **************************




     //exit(3);






      #ifdef DEBUG_MAIN
         fprintf (fpout, "\nd2algMain():\n");
         fprintf (stdout, "Expected subproblem obj = %f.\n", expObjval);
         fprintf(fpout, "\nSolution for all scenario subproblems:\n");
         printSparseMatrix(solnPtr->nrows, solnPtr->nzcnt_S,
                           solnPtr->rmatbeg_S, solnPtr->rmatcnt_S,
                           solnPtr->rmatind_S, solnPtr->rmatval_S, fpout);
      #endif

     //fprintf (stdout, "Checking solution if fractional....\n");
     //********** Check if solution solution is fractional and get the disjunction **********
     //********** variable and disjunction scenario *****************************************

     //prevDisjVar = disj_var;

     //fprintf (stdout, "\n\n Before disj_var: %d\n\n", disj_var);
     //fractional = getDisjunctionVar(solnPtr, &disj_var, &disj_scen, &disj_ind, subprobPtr->ctype);

     fractional = getMaxDisjunctionVar(solnPtr, &disj_var, &disj_scen, &disj_ind, subprobPtr->ctype);
     //fprintf (stdout, "\n\n Real disj_var: %d\n\n", disj_var);

  }

  else {

  fractional = 0;

  } // end if (MIP_SOLVE == 0)
  ///////////////////////////////////////////////////////////




	if (fractional ) { // Fractional solution, more work to do!!
       //fprintf(stdout, "\n Fractional solution: need to set up C3 LP\n");
   	#ifdef DEBUG_MAIN
             fprintf (stdout, "\nd2algMain():\n");
             fprintf (stdout, "FRACTIONAL solution found.\n");
             fprintf(stdout, " disj_scen = %d, disj_var = %d\n", disj_scen, disj_var);
             fprintf(stdout, " disj_ind = %d\n", disj_ind);
        #endif

        //*************Get averaged solution of all the scenario subproblems **************
        //*********This is the objective coef array for the C^3-LP coming up next *********

        getCondAverageSubProblemsSoln(solnPtr, stochdataPtr->scenProb, disj_var, disj_scen, aveSubprobSoln);

        //getConditionC3objCoefs(solnPtr, stochdataPtr->scenProb, aveSubprobSoln,
        //                      disj_scen);

        #ifdef DEBUG_MAIN
             fprintf (stdout, "\nd2algMain():\n");
             fprintf (stdout, "Average subproblem solution is > \n");
             for (i = 0; i < solnPtr->ncols; i++)
                   fprintf(stdout, "aveSubprobSoln[%d] = %f\n", i, aveSubprobSoln[i] );
        #endif


        //fprintf(stdout, "Allocating memory to C^3 LP\n");

        //********************* Allocate memory to C3 Dual lp data structure **********************
        status = memAllocC3LPproblemStruct(c3lpPtr, subprobPtr->nrows_total, subprobPtr->ncols,
                                        stochdataPtr->nscens);
        if ( status ) {
            fprintf (stderr, "\nd2algMain():\n");
            fprintf (stderr, " Failed to allocate memory to c3 dual lp data structure.\n");
            goto TERMINATE;
         }
         //fprintf(stderr, "Done allocating memory to C^3 LP\n");


        //*************************** Load C^3 lp data **********************************
        disj_var_floor = floor(solnPtr->rmatval_S[disj_ind]);
        disj_var_ceil =  ceil(solnPtr->rmatval_S[disj_ind]);


         loadC3LP(c3lpPtr, subprobPtr, stochdataPtr, aveSubprobSoln, solnX, disj_var,
                     disj_scen, disj_var_floor, disj_var_ceil, random_T, fpout);

//	printf ("\n c3lpPtr->nrowsWk = %d, in main", c3lpPtr->nrowsWk);

		//allocate memory for C3 LP information

        solnC3Lambda_01 = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
        solnC3Lambda_11 = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
        solnC3Lambda_02_12 = (double *) malloc (3*sizeof(double));
        solnC3pi = (double *) malloc (subprobPtr->ncols*sizeof(double));
        C3subgrad01   = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
		C3subgrad11   = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
		C3subgrad02_12= (double *) malloc (2*sizeof(double));
        C3addrowind = (int *) malloc ((2*c3lpPtr->nrowsWk+3)*sizeof(int));
        C3addrowval = (double *) malloc ((2*c3lpPtr->nrowsWk+3)*sizeof(double));
        C3secstageQ = (double *) malloc (stochdataPtr->nscens * sizeof(double));
        C3secobj        = (double *) malloc (stochdataPtr->nscens * sizeof(double));

                if (solnC3Lambda_01 == NULL || solnC3Lambda_11 == NULL ||
                    solnC3Lambda_02_12 == NULL || solnC3pi == NULL || C3subgrad01 == NULL || C3subgrad11==NULL || C3subgrad02_12 == NULL
                    || C3addrowind == NULL || C3addrowval == NULL || C3secstageQ == NULL || C3secobj == NULL)
          {
                fprintf (stderr, "No memory for C^3 solution values.\n");
              goto TERMINATE;
         }


        #ifdef DEBUG_MAIN
               	fprintf (stdout, "\n\n*********d2algMain()************:\n");
               	fprintf(stdout, "**disj_var = %d \n", disj_var);
               	fprintf(stdout, "disj_var_floor = %d \n", disj_var_floor);
               	fprintf(stdout, "disj_var_ceil = %d\n", disj_var_ceil);
               	printf("c3lpPtr->cmatspace = %d", c3lpPtr->cmatspace);
    	 	fprintf(stdout, "\nC^3 LP Constraint Matrix: \n");
    		fprintf(stdout, "rows = %d\n cols = %d\n", c3lpPtr->nrows, c3lpPtr->ncols);
    		printMatrix(c3lpPtr->ncols, c3lpPtr->nzcnt, c3lpPtr->cmatbeg,
                            c3lpPtr->cmatind, c3lpPtr->cmatval, stdout);
 	#endif

        //fprintf(stdout, "DONE PRINTING C^3 LP Constraint Matrix\n");



 	//******************** Create C^3 problem CPLEX lp model *********************************

 	CPXLPptr lp_c3 = NULL;		// pointer to subproblem lp to be created
 	lp_c3 = CPXcreateprob(env, &status, c3probname_lp);
        if ( status ) {
   	      fprintf (stderr, "\nd2algMain():\n");
      	      fprintf (stderr, "Failure to create CPLEX lp object lp_c3, error %d.\n", status);
    	      goto TERMINATE;
        }
        #ifdef DEBUG_MAIN
           fprintf(stdout, "\n AFTER CPXcreateprob\n");
           fprintf(stdout, " \n**c3lpPtr->nrows = %d**\n", c3lpPtr->nrows);
           fprintf(stdout, " **c3lpPtr->ncols = %d**\n\n", c3lpPtr->ncols);
        #endif
        status = CPXcopylp(env, lp_c3, c3lpPtr->ncols, c3lpPtr->nrows, CPX_MIN, c3lpPtr->obj,
        		   c3lpPtr->rhs, c3lpPtr->sense, c3lpPtr->cmatbeg, c3lpPtr->cmatcnt,
        		   c3lpPtr->cmatind, c3lpPtr->cmatval, c3lpPtr->lb, c3lpPtr->ub, NULL);

        //fprintf(stdout, "\n AFTER CPXcopylp\n");
        if ( status ) {
   	     fprintf (stderr, "d2algMain():\n");
      	     fprintf (stderr, "Failure to copy lp_c3 data into CPLEX lp object, error %d.\n", status);
    	     goto TERMINATE;
         }



     //    status = CPXlpopt (env, lp_c3);
        // Write sub prob to file in lp format
        #ifdef DEBUG_MAIN
             fprintf (fpout, "\nd2algMain():\n");
   	     status = CPXwriteprob(env, lp_c3, c3probname_lp, NULL);
   	    if ( status ) {
   		 fprintf (stderr, "d2algMain():\n");
      		 fprintf (stderr, "Failure to write c3 problem lp to file, error %d.\n", status);
      		 goto TERMINATE;
   	     }
         #endif

       //if (n_iterations == 3)
       //    goto TERMINATE;


        //fprintf(stdout, "\n AFTER CPXwriteprob\n");

       // Optimize the C^3 LP until fractional solution is cut off
       // then obtain solution: Pi_s, pi_0's, and lambda's.
       // If fractional soln is not cut off drop scenario solns from
       // the obj coef one at a time.

      //scenDropped = 0; // Counter for dropped off c^3 scenario obj coefs
      startToDropScens = 0;
      fracSolnCutOff = 1;

      if (DISPLAY_LOG_INFO) {
          fprintf (stdout, "\n Optimizing the Common-Cut-Coefficients lp...\n");
          fprintf (stdout, " Current accumulated Common-Cut-Coefficients lp CPU time: %f seconds\n", cputime_c3);
      }


		for ( i = 0; i < stochdataPtr->nscens; i++ ) {
			C3secobj[i] = 0;
		}

		for (i = 0; i<solnPtr->ncondScens; i++) {
			j = solnPtr->condScens[i];
			C3secobj[j] = solnPtr->condScenProbs[i];
		}


      do {
//		  printf ("\n getting into the C3lp loop! \n");
          solnCutOff = 0; // Assume that fractional soln is cut off


          // Write sub prob to file in lp format
          #ifdef DEBUG_MAIN
              fprintf (stdout, "\nd2algMain():\n");
              fprintf (stdout, "\nWriting c3 lp format \n");
   	          status = CPXwriteprob(env, lp_c3, c3probname_lp, NULL);
   	          if ( status ) {
   	              fprintf (stderr, "d2algMain():\n");
                  fprintf (stderr, "Failure to write c3 problem lp to file, error %d.\n", status);
                  goto TERMINATE;
   	          }
          #endif
		  C3UpperBnd = CPX_INFBOUND;
		  C3LowerBnd = -CPX_INFBOUND;

          cputime_start = clock();
          status = CPXlpopt (env, lp_c3);
          cputime_stop = clock();

		  cputime_c3 += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
          if ( status ) {
              fprintf (stderr, "Failed to optimize C^3 LP! Error code %d .\n", status);
              goto TERMINATE;
          }
          status  = CPXgetobjval (env, lp_c3, &C3LowerBnd);
          if ( status ) {
			  fprintf (stderr,"Failed to obtain objective value for C^3 lp .\n");
              fprintf (stderr,"Error code: %d\n", status);
              goto TERMINATE;
          }



		 do {
	//		printf ("\n getting into the innermost loop, and start algorithms! \n");

	//		cputime_c3 += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;

	//		printf ("Number of Columns is %d ", CPXgetnumcols(env, lp_c3));
	//		printf ("Number of rows is %d \n", CPXgetnumrows(env, lp_c3));

			status = CPXgetx (env, lp_c3, solnC3pi, 0, subprobPtr->ncols-1);
			if ( status ) {
				fprintf (stderr, "Failed to obtain C^3 LP pi's solution.\n");
				fprintf (stderr, "Error code: %d\n", status);
				goto TERMINATE;
			}

			cntmark1 = subprobPtr->ncols + c3lpPtr->nrowsWk;

			status = CPXgetx(env, lp_c3, solnC3Lambda_01, subprobPtr->ncols, cntmark1-1);
		    if ( status ) {
			      fprintf (stderr, "failed to obtain C^3 LP lambda_01's solution.\n");
				  fprintf (stderr, "Error code: %d\n", status);
				  goto TERMINATE;
			}

			cntmark2 = cntmark1 + c3lpPtr->nrowsWk;

			status = CPXgetx(env, lp_c3, solnC3Lambda_11, cntmark1, cntmark2-1);
		    if ( status ) {
			      fprintf (stderr, "Failed to obtain C^3 LP lambda_11's solution.\n");
				  fprintf (stderr, "Error code: %d\n", status);
				  goto TERMINATE;
			}


		    status = CPXgetx (env, lp_c3, solnC3Lambda_02_12, cntmark2, cntmark2+2);
         	   if ( status ) {
               		fprintf (stderr, "Failed to obtain C^3 LP lambda_02_02's solution.\n");
      		        fprintf (stderr, "Error code: %d\n", status);
              	 	goto TERMINATE;
            	}


			for (i = 0; i<c3lpPtr->nrowsWk; i++)
			{
				C3subgrad01[i] = 0;
          		        C3subgrad11[i] = 0;
			}

			C3subgrad02_12[0] = 0;
			C3subgrad02_12[1] = 0;

			C3addrhs = 0;


			C3_iteration = 0;

			for (i = 0; i < solnPtr->ncondScens; i++)
			{
				k = solnPtr->condScens[i];

				sumcount1 = 0;
				sumcount2 = 0;

           		        for (j = 0; j<2*c3lpPtr->nrowsWk+3; j++)
                                {
                                        C3addrowind[j] = 0;
                                        C3addrowval[j] = 0;
                                }


				for (j = 0; j < c3lpPtr->nrowsWk; j++)
				{
					if (j < subprobPtr->nrows_T) {
						if (subprobPtr->rhsRho[k][j] > NONZERO_LB || subprobPtr->rhsRho[k][j] < -NONZERO_LB) {
							sumcount1 = sumcount1 - solnC3Lambda_01[j]*subprobPtr->rhsRho[k][j];
							sumcount2 = sumcount2 - solnC3Lambda_11[j]*subprobPtr->rhsRho[k][j];
						}
					} else if (j < subprobPtr->nrowsW) {
						sumcount1 = sumcount1 - solnC3Lambda_01[j]*(-1);
						sumcount2 = sumcount2 - solnC3Lambda_11[j]*(-1);
					}
         			           else {
						C3_index = j - subprobPtr->nrowsW + subprobPtr->nrows_T;
						if (subprobPtr->rhsRho[k][C3_index] > NONZERO_LB || subprobPtr->rhsRho[k][C3_index] < -NONZERO_LB) {
							sumcount1 = sumcount1 - solnC3Lambda_01[j]*subprobPtr->rhsRho[k][C3_index];
							sumcount2 = sumcount2 - solnC3Lambda_11[j]*subprobPtr->rhsRho[k][C3_index];
						}
					}
				}

				sumcount1 = sumcount1 + solnC3Lambda_02_12[0]*disj_var_floor;
				sumcount2 = sumcount2 - solnC3Lambda_02_12[1]*disj_var_ceil;


				C3secstageQ[i] = maximum(sumcount1, sumcount2);

				mycount = 0;
				if (C3secstageQ[i] > 1+NONZERO_LB)
				{
					sumcount = -1;
					if (sumcount1 >= sumcount2) {
						for (j=0; j<c3lpPtr->nrowsWk; j++) {
							if (j < subprobPtr->nrows_T) {
								if (subprobPtr->rhsRho[k][j]>NONZERO_LB || subprobPtr->rhsRho[k][j]<-NONZERO_LB)
								{
									C3addrowind[mycount] =  j + subprobPtr->ncols;
									C3addrowval[mycount] =  subprobPtr->rhsRho[k][j];
									mycount++;
								}
							} else if (j < subprobPtr->nrowsW) {
								C3addrowind[mycount] = j + subprobPtr->ncols;
								C3addrowval[mycount] = -1;
								mycount++;
							} else {
								C3_index = j - subprobPtr->nrowsW + subprobPtr->nrows_T;
								 if (subprobPtr->rhsRho[k][C3_index] > NONZERO_LB || subprobPtr->rhsRho[k][C3_index] < -NONZERO_LB) {
									C3addrowind[mycount] = j + subprobPtr->ncols;
									C3addrowval[mycount] = subprobPtr->rhsRho[k][C3_index];
									mycount++;
								}
							}
						}
						C3addrowind[mycount] = subprobPtr->ncols + 2*c3lpPtr->nrowsWk;
						C3addrowval[mycount] = -disj_var_floor;
						mycount++;
					}
					else if (sumcount2 > sumcount1) {
                                                for (j=0; j<c3lpPtr->nrowsWk; j++) {
                                                        if (j < subprobPtr->nrows_T) {
                                                                if (subprobPtr->rhsRho[k][j]>NONZERO_LB || subprobPtr->rhsRho[k][j]<-NONZERO_LB)
                                                                {
                                                                        C3addrowind[mycount] =  j + c3lpPtr->nrowsWk + subprobPtr->ncols;
                                                                        C3addrowval[mycount] =  subprobPtr->rhsRho[i][j];
                                                                        mycount++;
                                                                }
                                                        } else if (j < subprobPtr->nrowsW) {
                                                                C3addrowind[mycount] = j + c3lpPtr->nrowsWk + subprobPtr->ncols;
                                                                C3addrowval[mycount] = -1;
                                                                mycount++;
                                                        } else {
                                                                C3_index = j - subprobPtr->nrowsW + subprobPtr->nrows_T;
                                                                 if (subprobPtr->rhsRho[k][C3_index] > NONZERO_LB || subprobPtr->rhsRho[k][C3_index] < -NONZERO_LB) {
                                                                        C3addrowind[mycount] = j + subprobPtr->ncols + c3lpPtr->nrowsWk;
                                                                        C3addrowval[mycount] = subprobPtr->rhsRho[k][C3_index];
                                                                        mycount++;
                                                                }
                                                        }
                                                }

						C3addrowind[mycount] = subprobPtr->ncols + 2*c3lpPtr->nrowsWk + 1;
						C3addrowval[mycount] = disj_var_ceil;
						mycount++;
					}

					status = CPXaddrows(env, lp_c3, 0, 1, mycount, &sumcount, &C3_sense, &zero, C3addrowind, C3addrowval, NULL, NULL);
				    if ( status ) {
						fprintf (stderr, "\nd2algMain():\n");
						fprintf (stderr, " Main: Failed to add cut row to C3 master lp.\n");
						fprintf (stderr, "Error code: %d\n", status);
						goto TERMINATE;
					}
				}
				else if (C3secstageQ[i] <= -1-NONZERO_LB) {
					C3secstageQ[i] = -1;
					C3_iteration++;
				}
				else
					C3_iteration++;




				C3addrhs = C3addrhs + C3secobj[k]*C3secstageQ[i];

                for (j = 0; j < c3lpPtr->nrowsWk;j++)
                {
//                     printf ("i = %d, k = %d \n", i, j);
                     if (sumcount1 >= sumcount2) {
						sumcount = C3subgrad01[j];
                        if (j < subprobPtr->nrows_T)
                           sumcount = sumcount - C3secobj[k]*subprobPtr->rhsRho[k][j];
                        else if (j < subprobPtr->nrowsW)
                           sumcount = sumcount - C3secobj[k]*(-1);
                        else {
                           C3_index = j - subprobPtr->nrowsW + subprobPtr->nrows_T;
                           sumcount = sumcount - C3secobj[k]*subprobPtr->rhsRho[k][C3_index];
                        }
						C3subgrad01[j] = sumcount;
                      }
					 else {
                        sumcount = C3subgrad11[j];
                        if (j< subprobPtr->nrows_T)
                            sumcount = sumcount - C3secobj[k]*subprobPtr->rhsRho[k][j];
                        else if (j<subprobPtr->nrowsW)
                            sumcount = sumcount-C3secobj[k]*(-1);
                        else {
                            C3_index = j - subprobPtr->nrowsW + subprobPtr->nrows_T;
                            sumcount = sumcount - C3secobj[k]*subprobPtr->rhsRho[k][C3_index];
                        }
                        C3subgrad11[j] = sumcount;
                     }

				}

				if (sumcount1 >= sumcount2)
					C3subgrad02_12[0] = C3subgrad02_12[0] + C3secobj[k] * disj_var_floor;
				else
					C3subgrad02_12[1] = C3subgrad02_12[1] - C3secobj[k] * disj_var_ceil;

			}

			if (C3_iteration == solnPtr->ncondScens) {
				sumcount = 0;

				for (i = 0; i<2*c3lpPtr->nrowsWk+3; i++)
				{
					C3addrowind[i] = 0;
					C3addrowval[i] = 0;
				}

				for (i = 0; i < subprobPtr->ncols-1; i++) {
					sumcount = sumcount + c3lpPtr->obj[i]*solnC3pi[i];
					//fprintf(stdout, "c3lpPtr->obj[%d] = %f\n", i, c3lpPtr->obj[i]);
				}

	//			printf ("sumcount = %f, C3addrhs = %f, Theta = %f. \n", sumcount, C3addrhs, solnC3Lambda_02_12[2]);
				sumcount = sumcount + C3addrhs;
				mycount = 0;

				C3UpperBnd = minimum(C3UpperBnd, sumcount);

				for (i = 0; i<c3lpPtr->nrowsWk; i++)
				{
					C3addrhs = C3addrhs - C3subgrad01[i] * solnC3Lambda_01[i];
					C3addrhs = C3addrhs - C3subgrad11[i] * solnC3Lambda_11[i];
				}


				C3addrhs = C3addrhs - C3subgrad02_12[0]*solnC3Lambda_02_12[0] - C3subgrad02_12[1]*solnC3Lambda_02_12[1];
	//			printf ("f(lambda)-subgrad*lambda = %f.\n", C3addrhs);

				for (i = 0; i<2*c3lpPtr->nrowsWk+2; i++)
				{
					if (i<c3lpPtr->nrowsWk){
						if (C3subgrad01[i] >0.00001 || C3subgrad01[i] <-0.00001 ) {
							C3addrowind[mycount] = i + subprobPtr->ncols;
							C3addrowval[mycount] = -C3subgrad01[i];
							mycount++;
						}
					}
					else if (i<2*c3lpPtr->nrowsWk) {
						if (C3subgrad11[i-c3lpPtr->nrowsWk] >0.00001 || C3subgrad11[i-c3lpPtr->nrowsWk] < -0.00001 ) {
							C3addrowind[mycount] = i + subprobPtr->ncols;
							C3addrowval[mycount] = -C3subgrad11[i-c3lpPtr->nrowsWk];
							mycount++;
						}
					}
					else if (i == 2*c3lpPtr->nrowsWk){
						if (C3subgrad02_12[0] > 0.00001 || C3subgrad02_12[0] < -0.00001)  {
							C3addrowind[mycount] = i + subprobPtr->ncols;
							C3addrowval[mycount] = -C3subgrad02_12[0];
							mycount++;
						}
					}
					else if (i == 2*c3lpPtr->nrowsWk+1){
		                if (C3subgrad02_12[1] > 0.00001 || C3subgrad02_12[1] < -0.00001) {
		                    C3addrowind[mycount] = i + subprobPtr->ncols;
		                    C3addrowval[mycount] = -C3subgrad02_12[1];
		                    mycount++;
                       	}
					}
				}

				C3addrowind[mycount] = 2*c3lpPtr->nrowsWk + 2 + subprobPtr->ncols;
				C3addrowval[mycount] = 1;
				mycount++;

				status = CPXaddrows(env, lp_c3, 0, 1, mycount, &C3addrhs, &C3_sense, &zero, C3addrowind, C3addrowval, NULL, NULL);
			    if ( status ) {
					fprintf (stderr, "\nd2algMain():\n");
					fprintf (stderr, " Main: Failed to add cut row to C3 master lp.\n");
					fprintf (stderr, "Error code: %d\n", status);
					goto TERMINATE;
				}

				c3lpPtr->nrows++;
			}

            cputime_start = clock();
//			 status = CPXwriteprob(env, lp_c3, "c3_lp", "lp");
//			 if (status) {
//				 fprintf(stderr, "Failed to write c3 lp %d", status);
//			 }
			 
            status = CPXlpopt (env, lp_c3);
            cputime_stop = clock();
            cputime_c3 += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;

            if ( status ) {
                fprintf (stderr, "Failed to optimize C^3 LP! Error code %d .\n", status);
                goto TERMINATE;
            }
            status  = CPXgetobjval (env, lp_c3, &C3LowerBnd);
            if ( status ) {
                fprintf (stderr,"Failed to obtain objective value for C^3 lp .\n");
                fprintf (stderr,"Error code: %d\n", status);
                goto TERMINATE;
            }
			status = CPXwriteprob(env, lp_c3, "C3.lp", NULL);


        //                printf ("\n finishing one innermost loop! \n");

//			printf ("C3UpperBnd = %lf, C3LowerBnd = %lf \n", C3UpperBnd, C3LowerBnd);

		}  while ( C3UpperBnd - C3LowerBnd > 0.000001);


                    status = CPXwriteprob(env, lp_c3, "C3.lp", NULL);
                    if ( status ) {
                        fprintf (stderr, "d2algMain():\n");
                        fprintf (stderr, "Failure to write c3 problem lp to file, error %d.\n", status);
                        goto TERMINATE;
                     }

          //*** if C^3 LP objective value is greater than zero (min prob) ***
          //*** then current fractional solution is not cut off !!!!!!!!!!***
          //*** Use the conditional expected scenario soln as objective   ***
          if (C3UpperBnd > -NONZERO_LB) {


                    solnCutOff = 1;


                    //***** Drop of scenario******
                    // Set number of obj coefs to change
        	    temp = subprobPtr->ncols;

        	    // Reset c^3 obj coefs
        	    for (j = 0; j < temp; j++){
                        c3lpPtr->obj[j] = 0;
        	    }
				for (i =0; i<stochdataPtr->nscens; i++)
					C3secobj[i] = 0;

//				printf ("time for drop the scenario.\n");

			printf ("solnPtr->ncondScens = %d. \n", solnPtr->ncondScens);
                status = dropScenSolnFromC3obj(solnPtr, c3lpPtr->obj,
                                      stochdataPtr->scenProb, disj_scen, C3secobj);
                    if (status) {
                        fprintf (stdout, "***** FAILED TO CUT OFF FRACTIONAL SOLN ****\n");
                        fprintf (stdout, "DROPPED ALL SCENARIOS!!!\n");
                        solnCutOff = 1;
                        fracSolnCutOff = 0;
                        no_nd2cuts_iters++;
        	         //goto TERMINATE;
                         break;
                    }
                    //fprintf (stdout, "Dropped off scenario\n");

                    //exit(0);




             // fprintf (stdout, "Changing c3 lp obj coefs\n");

              // Change the objective coefs for the C3 lp
             temp = subprobPtr->ncols+2*c3lpPtr->nrowsWk+3;
             status = CPXchgobj(env, lp_c3, temp, c3lpPtr->indices, c3lpPtr->obj);
             if (status) {
        	 fprintf (stderr,"**Failed to change objective coefs for C^3 lp** .\n");
                 fprintf (stderr,"Error code: %d\n", status);
                 goto TERMINATE;
             }


        	// Write sub prob to file in lp format
                #ifdef DEBUG_MAIN
                    fprintf (fpout, "\nd2algMain():\n");
   	            status = CPXwriteprob(env, lp_c3, c3probname_lp, NULL);
   	            if ( status ) {
   		        fprintf (stderr, "d2algMain():\n");
      		        fprintf (stderr, "Failure to write c3 problem lp to file, error %d.\n", status);
      		        goto TERMINATE;
   	             }
                #endif

                //scenDropped = 2;


                //exit(1);

          } else { // Update number of cols in next C^3 lp
               //c3lpPtr->ncols += 2;
               solnCutOff = 0;
               startToDropScens = 0;

               if (DISPLAY_LOG_INFO) {
                   fprintf (stdout, "\n Common-Cut-Coefficients lp objective value %.10g\n", C3UpperBnd);
               }
          }

     } while (solnCutOff); // Re-optimize c^3 lp

    //******************************
    if (fracSolnCutOff ) {

        cur_numcols = CPXgetnumcols (env, lp_c3);

        // Allocate space for solution
//        solnC3pi = (double *) malloc (subprobPtr->ncols*sizeof(double));
        solnC3pi_0 = (double *) malloc (stochdataPtr->nscens*sizeof(double));
//        solnC3Lambda_01 = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
//        solnC3Lambda_11 = (double *) malloc (c3lpPtr->nrowsWk*sizeof(double));
//        solnC3Lambda_02_12 = (double *) malloc (2*sizeof(double));

        if ( solnC3pi_0      == NULL ) {
              fprintf (stderr, "No memory for C^3 solution values.\n");
              goto TERMINATE;
         }

		for (i=0; i<stochdataPtr->nscens; i++)
			solnC3pi_0[i] = 0;

		for (i=0; i< solnPtr->ncondScens; i++) {
			k = solnPtr->condScens[i];
			solnC3pi_0[k] = -C3secstageQ[i];
		}

			status = CPXgetx (env, lp_c3, solnC3pi, 0, subprobPtr->ncols-1);
			if ( status ) {
				fprintf (stderr, "Failed to obtain C^3 LP pi's solution.\n");
				fprintf (stderr, "Error code: %d\n", status);
				goto TERMINATE;
			}

			cntmark1 = subprobPtr->ncols + c3lpPtr->nrowsWk;

			status = CPXgetx(env, lp_c3, solnC3Lambda_01, subprobPtr->ncols, cntmark1-1);
		    if ( status ) {
			      fprintf (stderr, "Failed to obtain C^3 LP lambda_01's solution.\n");
				  fprintf (stderr, "Error code: %d\n", status);
				  goto TERMINATE;
			}

			cntmark2 = cntmark1 + c3lpPtr->nrowsWk;

			status = CPXgetx(env, lp_c3, solnC3Lambda_11, cntmark1, cntmark2-1);
		    if ( status ) {
			      fprintf (stderr, "Failed to obtain C^3 LP lambda_11's solution.\n");
				  fprintf (stderr, "Error code: %d\n", status);
				  goto TERMINATE;
			}

		    status = CPXgetx (env, lp_c3, solnC3Lambda_02_12, cntmark2, cntmark2+1);
            if ( status ) {
               fprintf (stderr, "Failed to obtain C^3 LP lambda_02_02's solution.\n");
               fprintf (stderr, "Error code: %d\n", status);
               goto TERMINATE;
            }



          //**** Updata the subproblem constaint matrix W^k by adding the pi row ****
          // Also update the row indices for subproblem lp and mip
          updateConstrMatrixWk(subprobPtr, solnC3pi, disj_var, nd2cuts);

          //Update number of d^2 cuts added thus far
          nd2cuts++;
          subprobPtr->nd2cuts = nd2cuts;

         // fprintf(stdout, "\n ****** \n nd2cuts = %d \n \n************\n", nd2cuts);

          // Print the row sparse matrix W
  	  #ifdef DEBUG_MAIN
                fprintf (stdout, "\nd2algMain():\n");
     		fprintf(stdout, "subprobPtr->cmatspace_W = %d\n", subprobPtr->cmatspace_W);
     		fprintf(stdout, "subprobPtr->nzcnt_W = %d\n", subprobPtr->nzcnt_W);
     		fprintf(stdout, "\n\nW^k matrix in row sparse format is: \n");
     		printSparseMatrix(subprobPtr->nrows_total, subprobPtr->nzcnt_W, subprobPtr->cmatbeg_W, subprobPtr->cmatcnt_W, subprobPtr->cmatind_W,
                       subprobPtr->cmatval_W, stdout);
  	   #endif

          //**** Also update the subproblem lp instance by adding the pi row to the lp ****
          status = addNewRowToSubProbWmat(env, lp_sub, solnC3pi, subprobPtr->ncols);
          if (status) {
      		fprintf (stderr, "\nd2algMain:\n");
      		fprintf (stderr, " Failed to add new pi row to subproblem lp\n");
      		fprintf (stderr, "Error code: %d\n", status);
      		goto TERMINATE;
           }

           //Write sub prob to file in lp format
   	   #ifdef DEBUG_MAIN
                  fprintf (stderr, "\nd2algMain():\n");
   	          status = CPXwriteprob(env, lp_sub, subprobname_lp, NULL);
   	          if ( status ) {
   		     fprintf (stderr, "d2algMain():\n");
      		     fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		     fprintf (stderr, "Error code: %d\n", status);
      		     goto TERMINATE;
   	           }
            #endif


            ////////////////////////////////

         if (addCutsToMip == 1) {
             //**** Also update the subproblem mip instance by adding the pi row to the lp ****
             status = addNewRowToSubProbWmat(env, lp_submip, solnC3pi, subprobPtr->ncols);
             if (status) {
      		 fprintf (stderr, "\nd2algMain:\n");
      		 fprintf (stderr, " Failed to add new pi row to subproblem mip\n");
      		 fprintf (stderr, "Error code: %d\n", status);
      		 goto TERMINATE;
              }
           }

         ////////////////////////////////


           //Write sub prob to file in lp format
   	   #ifdef DEBUG_MAIN
                  fprintf (stderr, "\nd2algMain():\n");
   	          status = CPXwriteprob(env, lp_submip, subprobname_mip, NULL);
   	          if ( status ) {
   		     fprintf (stderr, "d2algMain():\n");
      		     fprintf (stderr, "Failure to write sub problem mip to file, error %d.\n", status);
      		     fprintf (stderr, "Error code: %d\n", status);
      		     goto TERMINATE;
   	           }
            #endif



    } //************************** end if(fracSolnCutOff)



 //    printf("Freeing lp_c3 mode and data structure\n");

     status = freeC3lpModelAndData(env, lp_c3, c3lpPtr);
     if (status){
         fprintf(stderr, "d2algMain():");
         fprintf(stderr, "Could not free C^3 LP model. Exiting....");
     }
   //  printf("lp_c3 Freed\n");

    //*** GO TO SOLVE RHS LPs ***


 } // End if (fractional)
 else if (MIP_SOLVE == 0)
 {  // INTEGER SOLUTION FOUND BY LP_RELAXATION, COOL!!
    lp_intsoln = 1;
    //printf("\n\n Integer Solution Found without adding any d^2 cuts!!\n");
    no_nd2cuts_iters++;
    //   printf("no_nd2cuts_iters = %d \n\n", no_nd2cuts_iters);
    ////////////////////////////////////////////////////////////////
    //          ADD Benders-type optimality cut                   //
    ////////////////////////////////////////////////////////////////
    if (DISPLAY_LOG_INFO) {
        fprintf (stdout, "\n Done adding Benders Cut to master problem MIP...\n");
    }
    status = addBendersCutToMaster(env, lp_master, masterprobPtr, SUBPROB_LB);
    if (status) {
               fprintf (stderr, "\nd2algMain:\n");
      	       fprintf (stderr, " Failed to add benders cut to master lp\n");
      	       goto TERMINATE;
    }

 } // if else statement





    //************ SOLVE RHS LPs FOR ALL SCENARIOS IF THERE WAS A FRACTIONAL SOLUTION **************
    if (fractional && fracSolnCutOff) {

        if (DISPLAY_LOG_INFO) {
            fprintf (stdout, "\n Optimizing RHS lps for all scenarios...\n");
            fprintf (stdout, " Current accumulated RHS lps CPU time: %f seconds\n", cputime_rhs);
        }

        // Remember that number of rows in the subproblem have been updated
        int nrows_prev = subprobPtr->nrows-1;

        // Allocate space for the gamma_0, gamma and gamma_1 vectors for the RHS LP
        //fprintf(stdout, "ncols = %d\n", stochdataPtr->ncols);
             gamma_0 = (double *) malloc (stochdataPtr->ncols*sizeof(double));
             gamma_1 = (double *) malloc (stochdataPtr->ncols*sizeof(double));
			 gamma = (double *) malloc (stochdataPtr->ncols*sizeof(double));
			 gamma_index = (double *) malloc (stochdataPtr->ncols*sizeof(double));
			 beta        = (double **) malloc (stochdataPtr->ncols* sizeof(double *));
			 if ( gamma_index == NULL || beta == NULL) {
			    fprintf (stderr, "Memory Allocation failure for gamma_index and beta arrays. Terminating...\n");
			    goto TERMINATE;
			 }

			 for (i = 0; i < stochdataPtr->ncols; i++){
			     beta[i] = (double *) malloc (stochdataPtr->ncols*sizeof(double));
			     if ( beta[i] == NULL) {
			        fprintf (stderr, "Memory Allocation failure for beta[%d] array.\n", i);
					fprintf(stderr, "Terminating...\n");
					goto TERMINATE;
				 }
			 }

        for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
             // Compute nu_0 and nu_1 parameters
             //fprintf(stdout, " computing nu_0 \n");
             //nu_0 = scalProd(solnC3Lambda_01, stochdataPtr->rhs[scenario], nrows_prev)
             //       - solnC3Lambda_02_12[0]*disj_var_floor;
             nu_0 = getNu01(subprobPtr, solnC3Lambda_01, stochdataPtr->rhs[scenario],
                            -solnC3Lambda_02_12[0], disj_var_floor, disj_var);

             //fprintf(stdout, " computing nu_1 \n");
             //nu_1 = scalProd(solnC3Lambda_11, stochdataPtr->rhs[scenario], nrows_prev)
             //       + solnC3Lambda_02_12[1]*disj_var_ceil;
             nu_1 = getNu01(subprobPtr, solnC3Lambda_11, stochdataPtr->rhs[scenario],
                            solnC3Lambda_02_12[1], disj_var_ceil, disj_var);


             //fprintf(stdout, " \nnu_0 = %6.6f\n", nu_0);
             //fprintf(stdout, " nu_1 = %6.6f\n", nu_1);

             //fprintf(stdout, "T matrix is: \n");
             //printSparseMatrix(stochdataPtr->ncols, subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
             //                   subprobPtr->cmatcnt_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, stdout);

             //fprintf(stdout,"Computing gamma_o and gamma_1 \n");
             // Compute gamma_0 and gamma_1 vectors
             getGammaH(solnC3Lambda_01, stochdataPtr, subprobPtr, gamma_0, scenario, random_T, disj_var);
	     //fprintf(stdout,"Done computing gamma_o  \n");
	     //fprintf(stdout,"Computing gamma_1 and gamma_1 \n");
	     getGammaH(solnC3Lambda_11, stochdataPtr, subprobPtr, gamma_1, scenario, random_T, disj_var);
	     //fprintf(stdout,"Done computing gamma_1  \n");

   	     #ifdef DEBUG_MAIN
                   fprintf(stdout, "d2alMain(): \n");
                   fprintf(stdout, "gamma_0 for scenario %d = ", scenario);
                   fprintf(stdout, "[");
  	           for (i = 0; i < stochdataPtr->ncols; i++)
  	                fprintf(stdout, "%f ", gamma_0[i]);
  	           fprintf(stdout, "]\n");
  	           fprintf(stdout, "gamma_1 for scenario %d = ", scenario);
  	           fprintf(stdout, "[");
  	           for (i = 0; i < stochdataPtr->ncols; i++)
  	                fprintf(stdout, "%f ", gamma_1[i]);
  	           fprintf(stdout, "]\n");
             #endif

             //**** Add the tau_00 and tau_01 cols to the RHS ***
         #ifdef DEBUG_MAIN
             fprintf(stdout, "Adding new cols (tau_00 and tau_01) to RHS lp...\n");

             fprintf(stdout, "***gamma_1 for scenario %d = ", scenario);
  	           fprintf(fpSolnOut, "[");
  	           for (i = 0; i < stochdataPtr->ncols; i++)
  	                fprintf(fpSolnOut, "%f ", gamma_1[i]);
  	           fprintf(fpSolnOut, "]\n");
  	  #endif
//*********************************************** ADD CUT *******************************


           if (ADD_CUT == 1)
           {
			   nu = minimum (nu_0, nu_1);
			   for (i = 0; i < stochdataPtr->ncols; i++)
			   {
				   gamma[i] = minimum(gamma_0[i], gamma_1[i]);
			   }
           }
		   else if ( ADD_CUT == 2)
		   {
			   nu_nega = 2;

			   if (nu_0 <= 0 || nu_1<= 0)
			   {
				   nu_max = maximum(-nu_0+1, -nu_1+1);
				   nu_0 = nu_0 + nu_max;
				   nu_1 = nu_1 + nu_max;
				   nu_nega = 1;
			   }
			   else
			   {
				   nu_nega = 0;
				   nu_max = maximum(1/nu_0, 1/nu_1);
			   }

			   for (i = 0; i<stochdataPtr->ncols; i++)
			   {
				   gamma_0[i] = gamma_0[i]/nu_0;
				   gamma_1[i] = gamma_1[i]/nu_1;
				   if (gamma_0[i] > 0 || gamma_1[i] > 0)
					   gamma_index[i] = 1;
				   else if (gamma_0[i] < 0 && gamma_1[i] < 0)
					   gamma_index[i] = 2;
				   else
					   gamma_index[i] = 3;
			   }

			   // calculate beta[i][j]
			   for (i = 0; i<stochdataPtr->ncols; i++)
			   {
				   for (j = 0; j<stochdataPtr->ncols; j++)
				   {
					   if (gamma_index[i]==1 && gamma_index[j]==2)
					   {
						   if (gamma_0[i] > 0 && gamma_1[i] > 0)
							   beta[i][j] = minimum(-gamma_0[j]/gamma_0[i], -gamma_1[j]/gamma_1[i]);
						   else if(gamma_0[i] > 0 && gamma_1[i] <=0)
							   beta[i][j] = -gamma_0[j]/gamma_0[i];
						   else if(gamma_1[i] > 0 && gamma_0[i] <=0)
							   beta[i][j] = -gamma_1[j]/gamma_1[i];
					   }
					   else
						   beta[i][j] = 0;
				   }
			   }

			   // calculate pi0
			   for (j = 0; j<stochdataPtr->ncols; j++)
			   {
				   if (gamma_index[j] == 1)
					   gamma[j] = maximum(gamma_0[j], gamma_1[j]);
				   else if (gamma_index[j] == 2)
				   {
					   gamma_max = -9999999999;
					   for (i = 0; i<stochdataPtr->ncols; i++)
					   {
						   if (gamma_index[i] == 1)
						   {
							   gamma_temp =  (-gamma[i])*beta[i][j];
							   if (gamma_temp > gamma_max)
								   gamma_max = gamma_temp;
						   }
					   }
					   gamma[j] = gamma_max;
				   }
				   else
					   gamma[i] = 0;
			   }

			   if (nu_nega == 1)
				   nu = 1 - nu_max;
			   else if (nu_nega == 0)
			   {
				   nu = 1/nu_max;
				   for (i = 0; i<stochdataPtr->ncols; i++)
					   gamma[i] = gamma[i]/nu_max;
			   }
		   }
              updateRHSrAndMatT(&struct_s, subprobPtr, nu, gamma, scenario, solnX);

              //updateRHSrAndMatT(stochdataPtr, subprobPtr, solnRHSdelta, solnRHSsigma_0,
              //                  solnRHSsigma_i, scenario, solnX);

              //printf("4. stochdataPtr->rmatind_T[18][0] = %d\n\n\n", stochdataPtr->rmatind_T[18][0]);

                  //fprintf(stdout, " END   updateRHSrAndMatT \n" );
              #ifdef DEBUG_MAIN
                  for (j = 0; j < subprobPtr->nrows; j++)
                       fprintf(stdout, "2. r[%d][%d] = %f \n", scenario, j, stochdataPtr->rhs[scenario][j]);

                   fprintf(fpSolnOut, "updated r(w): \n");
                   temp = stochdataPtr->nrows + stochdataPtr->rnrows +1;
                   for (j = 0; j < subprobPtr->nrows; j++)
                       fprintf(stdout, "r[%d][%d] = %f \n", scenario, j, stochdataPtr->rhs[scenario][j]);
                    temp = stochdataPtr->rnrows +1;
                    fprintf(stdout, "\nUpdated Technology Matrix T(%d): \n", scenario);
  	            printSparseMatrix(temp, stochdataPtr->rnzcnt_T[scenario],
  	                              stochdataPtr->rmatbeg_T[scenario], stochdataPtr->rmatcnt_T[scenario],
  	                              stochdataPtr->rmatind_T[scenario], stochdataPtr->rmatval_T[scenario],
  	                              stdout);
              #endif

    /******************************* no RHSLP needed ****************************

              //*** Delete the tau_00 and tau_01 cols added to the lp_rhs lp ***
              delstat = (int*) malloc (cur_numcols*sizeof(int));
              for (i = 0; i < cur_numcols-2; i++)
                   delstat[i] = 0;       // Just to make sure no elem is equal to 1
              delstat[cur_numcols-2] = 1;
              delstat[cur_numcols-1] = 1;
              status = CPXdelsetcols(env, lp_rhs, delstat);

               if (status){
                  fprintf(stderr, "Failed to delete set of cols for RHS lp. status: %d", status);
                  goto TERMINATE;
               }

               // free memory
               free(delstat);
               delstat = NULL;

              // Write RHS lp to file in LP format
             #ifdef DEBUG_MAIN
   	          status = CPXwriteprob(env, lp_rhs, rhsprobname_lp, NULL);
   	          if ( status ) {
   		       fprintf (stderr, "d2algMain():\n");
      		       fprintf (stderr, "Failure to write RHS problem lp to file, error %d.\n", status);
      		       goto TERMINATE;
   	          }
             #endif
 ***********************************************************************************/

         } // End for scenario loop









         if (fracSolnCutOff) {
            //******* Update the number of rows in the updated T(w) matrix ********
            stochdataPtr->rnrows++;
          }


    } // End if (fractional && and solnCutOff)


        //************ OPTIMIZE THE RE-UPDATED SUBPROBLEM LPs ***********************


  //printf("fractional = %d\n", fractional);
 // printf("MIP_SOLVE = %d\n", MIP_SOLVE);
  //exit(2);

  if (fractional || MIP_SOLVE) {

        ///////////////////////////////////////////////////////////
        //	Solve subproblems for all the scenarios 	 //
        ///////////////////////////////////////////////////////////
        if (DISPLAY_LOG_INFO) {
            fprintf (stdout, "\n Re-optimizing scenario subproblems...\n");
            fprintf (stdout, " Current accumulated subproblem lps CPU time: %f seconds\n", cputime_sublp);
            fprintf (stdout, " Current accumulated subproblem MIPs CPU time: %f seconds\n", cputime_submip);
        }

        expObjval    = 0;    // Initialize subprob expected obj value
        expmipObjval = 0;    // Initialize mip subprob expected obj value
        fractional   = 0;     // Assume integral solutions will be found

       fprintf(fpSolnOut, " \n\n\n stochdataPtr->nscens =  %d 		\n", stochdataPtr->nscens);

       ////////////////////////////////////////////////////////////
       for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
	//	printf ("getting the loop for scenario subproblems.\n");
             // fprintf(fpSolnOut, " \n\n\n ////////////////////// SCENARIO %d :		\n", scenario);
            #ifdef DEBUG_MAIN
                fprintf (stdout, "d2algMain(): rhsRho \n");
                for(i=0; i < subprobPtr->nrows;i++){
                       fprintf(stdout, "index[%d] = %d   subprobPtr->rhsRho[%d][%d] = %f\n",i, subprobPtr->indices_row[i],
                       scenario, i, subprobPtr->rhsRho[scenario][i]);
                }
            #endif

	//	printf ("start to change rhs of the subproblem. ");
            // Set the RHS for this scenario sub problem
            status = CPXchgrhs(env, lp_sub, subprobPtr->nrows, subprobPtr->indices_row,
                          subprobPtr->rhsRho[scenario]);
           if ( status ) {
                 fprintf (stderr, "Failed to set the rhs for the index scenario %d subproblem.\n", scenario);
                 fprintf(stderr, "Error: %d\n", status);
                 goto TERMINATE;
           }
	//	printf ("finish to change the rhs of the subproblem.\n");

           // Set the obj for this scenario sub problem if random obj
           if (random_obj) {
                 status = CPXchgobj(env, lp_sub, stochdataPtr->obj_cnt, stochdataPtr->obj_index, stochdataPtr->obj[scenario]);
                if ( status ) {
                      fprintf (stderr,"Failed to change objective value for scenario subproblem index %d .\n", scenario);
                      fprintf(stderr, "Error code: %d", status);
                      goto TERMINATE;
                }


           }

           //Write sub prob to file in lp format
           #ifdef DEBUG_MAIN
                fprintf (fpout, "d2algMain():\n");
   	        status = CPXwriteprob(env, lp_sub, subprobname_lp, NULL);
   	        if ( status ) {
   		   fprintf (stderr, "d2algMain():\n");
      		   fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		   goto TERMINATE;
   	       }
            #endif


        /////////////////////////////////////////////////////////
        //              BRANCH-AND-BOUND HERE		        //
        /////////////////////////////////////////////////////////
        bestbound = 0.0;
        if (MIP_SOLVE == 0) {

              if (curr_percent_gap > PERC_START_TBB)
                   my_max_nodes = 1;
               else
                   my_max_nodes = max_numnodes;

               num_nodesexplored = 0; // Initialized
               num_termnodes     = 0; // Initialized
               fract_sub         = 1;  // Initialized flag to fractional solution
	//	printf ("start bbdobranchbound.");
               //fprintf(stdout, "\n Main: Calling BBdobranchbound() for scenario %d\n", scenario);
               status = BBdobranchbound(env, lp_sub, subprobPtr, stochdataPtr, my_max_nodes, scenario,
                                    bb_nu, bb_gamma, &num_nodesexplored, &num_termnodes,
                                    &cputime_sublp,  &bestbound, solnY, &fract_sub,random_T);
              if ( status ) {
   		   fprintf (stderr, "d2algMain():\n");
      		   fprintf (stderr, "Failure to do branch-and-bound, error %d.\n", status);
      		   goto TERMINATE;
   	      }
	//	printf ("finish bbdobranchbound.\n");

   	      // Check if this scenario yielded an integral solution
   	      //fprintf(stdout, " bestbound = %f \n", bestbound);
   	      //fprintf(stdout, " fract_sub = %d \n", fract_sub);
   	      if (fractional == 0) {
   	           if (fract_sub == 0) { // integral subproblem scenario
   	               // Store the solution
   	               expmipObjval += bestbound*stochdataPtr->scenProb[scenario];
                       storeSubProbSoln(solnPtr, solnY, scenario);
   	               //fprintf(stdout, "\n Subproblem Best INTEGRAL Solution for scenario %d \n", scenario);
   	               //for (j = 0; j < ncols_sub; j++) {
                       //   fprintf (stdout, " %s: = %6.6f\n", subprobPtr->colnames[j], solnY[j]);
                       //}
   	           } else {
   	               fractional = 1;  // No incumbent integral soln for all scenarios
   	               //fprintf(stdout, "\n Subproblem FRACTIONAL Solution for scenario %d \n", scenario);
   	           }
   	      }

              //Write sub prob to file in lp format
              #ifdef DEBUG_MAIN
                fprintf (fpout, "d2algMain():\n");
   	        status = CPXwriteprob(env, lp_sub, "BBsub.lp", NULL);
   	        if ( status ) {
   		   fprintf (stderr, "d2algMain():\n");
      		   fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		   goto TERMINATE;
   	         }
               #endif

            //printf("\n\n MAIN: num_nodesexplored = %d\n", num_nodesexplored);
            //printf(" MAIN: num_termnodes = %d\n", num_termnodes);

            // Translate the nu's by the lower bound on the expected value
            for (i = 0; i < num_termnodes; i++)
                  bb_nu[i] += SUBPROB_LB;

            //for (i = 0; i < num_termnodes; i++){
            //      printf("\n bb_nu[%d] = %f\n", i, bb_nu[i]);
            //      for (j = 0; j < masterprobPtr->ncols-1; j++){
            //           printf("\n bb_gamma[%d] = %f\n", j, bb_gamma[i][j]);
            //      }
            //}

            //////////////////////////////////////////////////////////////////////
            //         Compute the expected "D^2 optimality cut" coefs          //
            //									//
            //////////////////////////////////////////////////////////////////////

            if (scenario == 0) { // Re-initialize
                  d2optcut_rhs = 0;
                  for (i = 0; i < masterprobPtr->ncols-1; i++)
                      d2optcut_coefs[i] = 0;
            }
        if (num_termnodes == 1) { ////// No need to solver Reverse Polar LP

             d2optcut_rhs +=  stochdataPtr->scenProb[scenario]*bb_nu[0];
             for (i = 0; i < masterprobPtr->ncols-1; i++)
                  d2optcut_coefs[i] += stochdataPtr->scenProb[scenario]*bb_gamma[0][i];

          } else {

             ////////////////////////////////////////////////
            //          Allocate memory to BB LP          //
            ////////////////////////////////////////////////
   	     bblpPtr->ncols  = 2 + masterprobPtr->ncols - 1 + num_termnodes*(1 + masterprobPtr->nrows);
             bblpPtr->nrows  = 1 + num_termnodes*(2 + masterprobPtr->ncols - 1);
             //printf("\n bblpPtr->nrows = %d\n", bblpPtr->nrows);
             //printf("\n bblpPtr->ncols = %d\n", bblpPtr->ncols);

   	     ////////////////////////////////////////////////
             //          Load BB LP data                   //
             ////////////////////////////////////////////////


             eta_k = solnX[masterprobPtr->ncols-1];

             if (n_iterations == 0) { //(eta_k == 0){
                 //eta_k = 1;
                 eta_k = SUBPROB_LB;
              }

	//	printf ("start bbloadreversepolarLP. ");
   	     BBloadreversepolarLP(bblpPtr, masterprobPtr, num_termnodes,
                     bb_nu, bb_gamma, solnX, eta_k);
	//	printf ("finish bbloadreversepolarLP. \n");


   	     ////////////////////////////////////////////////
             //          Create BB CPLEX LP model          //
             ////////////////////////////////////////////////
 	     //CPXLPptr lp_bb = NULL;		// pointer to subproblem lp to be created

 	     lp_bb = CPXcreateprob(env, &status, bbprobname_lp);
             if ( status ) {
   	         fprintf (stderr, "\nd2algMain():\n");
      	         fprintf (stderr, "Failure to create CPLEX lp object lp_bb, error %d.\n", status);
    	         goto TERMINATE;
              }
              #ifdef DEBUG_MAIN
                   fprintf(stdout, "\n AFTER CPXcreateprob\n");
                   fprintf(stdout, " \n**bblpPtr->nrows = %d**\n", bblpPtr->nrows);
                   fprintf(stdout, " **bblpPtr->ncols = %d**\n\n", bblpPtr->ncols);
              #endif
              status = CPXcopylp(env, lp_bb, bblpPtr->ncols, bblpPtr->nrows, CPX_MIN, bblpPtr->obj,
        		   bblpPtr->rhs, bblpPtr->sense, bblpPtr->cmatbeg, bblpPtr->cmatcnt,
        		   bblpPtr->cmatind, bblpPtr->cmatval, bblpPtr->lb, bblpPtr->ub, NULL);

               //fprintf(stdout, "\n AFTER CPXcopylp\n");
              if ( status ) {
   	          fprintf (stderr, "d2algMain():\n");
      	          fprintf (stderr, "CPXcopylp failure: copy lp_bb data into CPLEX lp object failed, error %d.\n", status);
    	          goto TERMINATE;
              }



             // Write sub prob to file in lp format
             #ifdef DEBUG_MAIN
                  fprintf (fpout, "\nd2algMain():\n");
   	          status = CPXwriteprob(env, lp_bb, "bbrevpolar.lp", NULL);
   	          if ( status ) {
   		       fprintf (stderr, "d2algMain():\n");
      		       fprintf (stderr, "Failure to write bb problem lp to file, error %d.\n", status);
      		       goto TERMINATE;
   	           }

                    //printSparseMatrix(bblpPtr->ncols, bblpPtr->nzcnt, bblpPtr->cmatbeg,
                    //                  bblpPtr->cmatcnt, bblpPtr->cmatind, bblpPtr->cmatval, stdout);
             #endif

             //if (n_iterations == 0 && scenario == 1)
             //    exit(1);


              ///////////////////////////////////////////////
              //           OPTIMIZE BB OPT CUT LP          //
              ///////////////////////////////////////////////

              //fprintf(stdout, "Optimizing TB&B Optimaliy Cut LP...   \n");
              cputime_start = clock();
              //status = CPXdualopt (env, lp_bb);
              status = CPXprimopt (env, lp_bb);
              cputime_stop = clock();
              cputime_bblp += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;

              if ( status ) {
                 fprintf (stderr, "Failed to optimize TB&B Optimaliy Cut LP. Error code: %d\n", status);
                 fprintf (stderr, "Error code: %d\n", status);
                 goto TERMINATE;
              }

              status  = CPXgetobjval (env, lp_bb, &objval);
	      if ( status ) {
                 fprintf (stderr,"Failed to obtain objective value for RHS lp .\n");
                 fprintf(stderr, "Error code: %d", status);
                 goto TERMINATE;
              }
              //fprintf (stdout, "TB&B Opt Cut lp objective value %.10g\n", objval);
              //cur_numcols = CPXgetnumcols (env, lp_bb);

 	      // Get lp solution: delta, sigma_0, and sigma_i's only
              status = CPXgetx (env, lp_bb, soln_bblp, 0, masterprobPtr->ncols+1);
              if ( status ) {
                  fprintf (stderr, "Failed to obtain TB&B Opt Cut LP solution. Terminating...\n");
                  fprintf(stderr, "Error code: %d", status);
                  goto TERMINATE;
              }

             // Get the delta solution
             soln_bbdelta = soln_bblp[0];
             // Get the sigma_0 solution
             soln_bbsigma_0 = soln_bblp[1];

             // Get the sigma_i's solution
             for (j = 0; j < masterprobPtr->ncols-1; j++){
                   soln_bbsigma_i[j] = soln_bblp[j+2];
             }

             //* Write out the solution *
             //#ifdef DEBUG_MAIN
                   fprintf (stdout, "d2algMain():\n");
                   fprintf(stdout, "\nTB&B Optimality cut LP  solution for SCENARIO %d:\n", scenario);
                   fprintf (stdout," delta = %6.6g\n", soln_bbdelta);
                   fprintf (stdout," sigma_0 = %6.6g\n", soln_bbsigma_0);
                   for (j = 0; j < masterprobPtr->ncols-1; j++)
                        fprintf (stdout, " bbsigma_i[%d] = %6.6g\n",  j, soln_bbsigma_i[j]);
                   fprintf (stdout," rhs = %6.6g\n", soln_bbdelta/soln_bbsigma_0);
                   for (j = 0; j < masterprobPtr->ncols-1; j++)
                        fprintf (stdout, " lhs[%d] = %6.6g\n",  j, soln_bbsigma_i[j]/soln_bbsigma_0);
             //#endif

              //Compute the expected "D^2 optimality cut" coefs
              //if (scenario == 0) { // Re-initialize
              //    d2optcut_rhs = 0;
              //    for (i = 0; i < masterprobPtr->ncols-1; i++)
              //        d2optcut_coefs[i] = 0;
              //}
              d2optcut_rhs +=  stochdataPtr->scenProb[scenario]*soln_bbdelta/soln_bbsigma_0;
              for (i = 0; i < masterprobPtr->ncols-1; i++)
                  d2optcut_coefs[i] += stochdataPtr->scenProb[scenario]*soln_bbsigma_i[i]/soln_bbsigma_0;


             //BBfreelpmodel(env, lp_bb, bblpPtr);
             // Free up the problem as allocated by CPXcreateprob
             if ( lp_bb != NULL ) {
                    //fprintf(stdout, "FREE-ING lp__bb......................  \n");
                    status = CPXfreeprob (env, &lp_bb);
                    if ( status ) {
                         fprintf (stderr, "CPXfreeprob lp_bb failed, error code %d.\n", status);
                    }
                }
          } // End if(num_termnodes == 1)/else
            /////////////////////////////////////////////////



       }

       else  // MIP_SOLVE == 1

       {
          /////////////////////////////////////////////////////
          //                   SOLVE MIPS                    //
          /////////////////////////////////////////////////////

	  // Set the RHS for this scenario mip problem
	  if (nd2cuts > 0) { // Cuts added
	     i = subprobPtr->nrowsWmip;
	     for (j = subprobPtr->nrowsW; j < subprobPtr->nrows; j++) {
	            //fprintf(stdout, "subprobPtr->rhsRho[%d] = %f - ", i, subprobPtr->rhsRho[scenario][i]);
	            subprobPtr->rhsRho[scenario][i] = subprobPtr->rhsRho[scenario][j];
	            //fprintf(stdout, "subprobPtr->rhsRho[%d] = %f: subprobPtr->rhsRho[%d] = %f\n",
	            //               i,subprobPtr->rhsRho[scenario][i], j, subprobPtr->rhsRho[scenario][j]);
	            i++;
	      }
	  }
          #ifdef DEBUG_MAIN
              fprintf (stdout, "d2algMain(): MIP\n");
              fprintf (stdout, "subprobPtr->nrows_mip: %d\n", subprobPtr->nrows_mip);
              for(i=0; i < subprobPtr->nrows_mip;i++){
                   fprintf(stdout, "indexmip[%d] = %d   rhsmip[%d] = %f\n",i, subprobPtr->indices_rowmip[i], i,
                        subprobPtr->rhsRho[scenario][i]);
              }
          #endif

          if (addCutsToMip == 1) {
               status = CPXchgrhs(env, lp_submip, subprobPtr->nrows_mip, subprobPtr->indices_rowmip,
                                  subprobPtr->rhsRho[scenario]);
           } else {
               status = CPXchgrhs(env, lp_submip, subprobPtr->nrowsWmip, subprobPtr->indices_rowmip,
                                  subprobPtr->rhsRho[scenario]);
           }
           if ( status ) {
                 fprintf (stderr, "Failed to set the rhs for the index scenario %d subproblem mip.\n", scenario);
                 fprintf(stderr, "Error: %d\n", status);
                 goto TERMINATE;
           }

           //Write sub prob to file in lp format
 //          #ifdef DEBUG_MAIN
                fprintf (fpout, "d2algMain():\n");
   	        status = CPXwriteprob(env, lp_submip, subprobname_mip, NULL);
   	        if ( status ) {
   		   fprintf (stderr, "d2algMain():\n");
      		   fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		   goto TERMINATE;
   	       }

 //           #endif

          #ifdef DEBUG_MAIN
               fprintf(fpSolnOut, "Optimizing subproblem MIP scenario index %d  \n", scenario);
           #endif

		printf ("current MIP_SOLVE = %d. \n", MIP_SOLVE);

//   status = CPXsetintparam (env, CPX_PARAM_BNDSTRENIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr,
               "Failure to turn off the bound strengthing, error %d.\n", status);
      goto TERMINATE;
   }
		   

//   status = CPXsetintparam (env, CPX_PARAM_COEREDIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr,
               "Failure to turn off the coefficient reduce, error %d.\n", status);
      goto TERMINATE;
   }

         cputime_start = clock();
         status = CPXmipopt (env, lp_submip);
         cputime_stop = clock();
         cputime_submip += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;

		   status = CPXsetintparam(env, CPX_PARAM_EPGAP, 1e-4);
		   status = CPXsetintparam (env, CPX_PARAM_SCRIND, 0) ;
		   status = CPXsetintparam(env, CPX_PARAM_PREIND, 0);
		   

		   
         if ( status ) {
            fprintf (stderr, "Failed to optimize scenario %d subproblem.\n", scenario);
            fprintf(stderr, "Error %d", status);
            goto TERMINATE;
         }
         status  = CPXgetmipobjval (env, lp_submip, &objval);
//		    CPXwriteprob(env, lp_submip, "what.lp", "lp");
		   
         fprintf (fpSolnOut, "\nSCENARIO %d MIP obj value = %.10g", scenario, objval);
	 if ( status ) {
             fprintf (stderr,"Failed to obtain objective value for subproblem scenario index %d .\n", scenario);
		 CPXwriteprob(env, lp_submip, "what.lp", "lp");
		 
             fprintf(stderr, "Error code: %d", status);
             goto TERMINATE;
         }

         //fprintf (fpSolnOut, "Scenario %d Obj value = %.10g\n", scenario, objval);
         // Update subprob expected obj value
         expmipObjval += objval*stochdataPtr->scenProb[scenario];

         cur_numcols = CPXgetnumcols (env, lp_submip);

         status = CPXgetmipx (env, lp_submip, solnY, 0, cur_numcols-1);
         if ( status ) {
              fprintf (stderr, "Failed to obtain subproblem mip solution.\n");
              goto TERMINATE;
          }

         // * Write out the solution *
         #ifdef DEBUG_MAIN
             fprintf (stdout, "d2algMain():\n");
             fprintf(fpSolnOut, "\nScenario %d\n", scenario);
             fprintf(fpSolnOut, "\nSubproblem problem MIP Solution:\n");
             for (j = 0; j < cur_numcols; j++)
                   fprintf (fpSolnOut, " %s: = %6.6g\n", subprobPtr->colnames[j], solnY[j]);
         #endif

         // Store the solution
         storeSubProbSoln(solnPtr, solnY, scenario);
         // Count subproblem mip solves
            nmipsolves++;

      } // End if (MIP_SOLVE == 0)/else

      /////////////////////////////////////////////////////////////////////
   } // End scenario for loop

    ////////////////////////////////////////////////////////
    //              Set incumbent solution                //
    ////////////////////////////////////////////////////////
    if (fractional == 0) {
   	fprintf (fpSolnOut, "\nExpected Objective value %.10g\n\n", expmipObjval);
        expmipObjval += SUBPROB_LB;
        //fprintf (stdout, "\nExpected mip objective after translation value %.10g\n\n", expmipObjval);
        curr_lb_V = objval_m - solnX[masterprobPtr->ncols-1] + expmipObjval;
        if (curr_lb_V < BOUND_V_k){
     	     fprintf (fpSolnOut, "\nSetting incumbent soln\n");
             setIncumbent(solnX, objval_m, expmipObjval-SUBPROB_LB, masterprobPtr->ncols, solnPtr);
             BOUND_V_k = curr_lb_V;
        }
     	fprintf(fpSolnOut, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);
      } // end if fractional

      // Print "D2-BAC Optimalit"y Cut Coefs
      #ifdef DEBUG_MAIN
          if (max_numnodes > 1 && curr_percent_gap <= PERC_START_TBB) {
             fprintf (stdout," D2-BAC Optimality cut:\n");
             fprintf (stdout,"\n eta ");
             for (j = 0; j < masterprobPtr->ncols-1; j++)
                 fprintf (stdout, " + %3.6g x_%d", d2optcut_coefs[j], j);
             fprintf (stdout," >= %3.6g\n", d2optcut_rhs);
           }
      #endif


   }

   else  //if (fractional == 0) // INTEGER SOLUTION

   {
      ////////////////////////////////////////////////////////
      //     Update upper bound & set incumbent solution    //
      ////////////////////////////////////////////////////////

      fprintf(fpSolnOut, "INCUMBENT INTEGER SOLUTION FOUND FOR ALL SCENARIOS:)\n");
      /////
         expmipObjval = expObjval;
      ////
      // **** Translate subproblem expected objective value *******
      // This makes the subprob lower bound = zero as required by the D^2 algorithm
      // For the Laporte and Louveaux opt cut L = zero
      fprintf (fpSolnOut, "\nExpected Objective value %.10g\n\n", expObjval);
      expObjval += SUBPROB_LB;

      ///// Added June 3, 2003
         expmipObjval = expObjval;
     ////

      fprintf (fpSolnOut, "\nExpected Translated Objective value %.10g\n\n", expObjval);

      //******* Updata the lower bound on the master problem obj ********
      curr_lb_V = objval_m - solnX[masterprobPtr->ncols-1] + expObjval;
      if (curr_lb_V < BOUND_V_k){
            fprintf (fpSolnOut, "\nSetting incumbent soln\n");
            setIncumbent(solnX, objval_m, expObjval-SUBPROB_LB, masterprobPtr->ncols, solnPtr);
      }
      BOUND_V_k = minimum(curr_lb_V, BOUND_V_k );
      prev_lb = curr_lb_V;

      //fprintf (stdout, "\nTotal Objective value %.10g\n\n", curr_lb_V-SUBPROB_LB);
      //fprintf (fpSolnOut, "\nTotal Objective value %.10g\n\n", curr_lb_V-SUBPROB_LB);
      //fprintf(fpSolnOut, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);



  } // End if/else


    fprintf(fpSolnOut, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);
    //fprintf(stdout, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);

    if (DISPLAY_LOG_INFO) {
          fprintf (stdout, "\n Adding Benders Cut to master problem MIP...\n");
     }
     //fprintf (stdout, "\nSUBPROB_LB = %f\n", SUBPROB_LB);
     //fprintf(stdout, "masterprobPtr->rhsCoef = %f\n", masterprobPtr->rhsCoef);
     //for (j = 0; j < masterprobPtr->ncols-1; j++) {
      //  fprintf(stdout, "masterprobPtr->cutCoefs[%d] = %f	\n", j, masterprobPtr->cutCoefs[j]);
     //}

    ////////////////////////////////////////////////////////////////////////
    //              Add optimality cuts to master program                 //
    ////////////////////////////////////////////////////////////////////////
    if (MIP_SOLVE == 0 && lp_intsoln == 0) { // Add "D2 optimality cut"
          status = BBaddD2optcuttomaster(env, lp_master, masterprobPtr, d2optcut_rhs, d2optcut_coefs);
          if (status) {
              fprintf (stderr, "\nd2algMain:\n");
      	      fprintf (stderr, " Failed to add benders cut to master lp\n");
      	      goto TERMINATE;
          }
          if (DISPLAY_LOG_INFO) {
              fprintf (stdout, "\n Done ""D2 Optimality cut"" to master problem MIP...\n");
          }
     }
     //////////////////////////////////////////////////////////////////////////////////////


    //Write master prob to file in lp format
    #ifdef DEBUG_MAIN
        if (DISPLAY_LOG_INFO) {
              fprintf (stdout, "\n Writing master problem MIP...\n");
         }
         //fprintf (stdout, "\nd2algMain():\n");
         status = CPXwriteprob(env, lp_master, mastername_lp, NULL);
         if ( status ) {
              fprintf (stderr, "d2algMain():\n");
              fprintf (stderr, "Failure to write update master problem lp to file, error %d.\n", status);
              goto TERMINATE;
         }
         //exit(0);
     #endif

     ////////////////////////////////////////////////////////////////////////
     //              Add Laporte and Louveaux Cut to close the  gap        //
     //              between the UB and the LB 			           //
     ////////////////////////////////////////////////////////////////////////
     if (MIP_SOLVE && ADD_LL == 1 || fractional == 0 && ADD_LL == 1) {
          if (DISPLAY_LOG_INFO) {
                  fprintf (stdout, "\n Adding LL Cut to master problem MIP...\n");
          }
          // ADD LAPORTE AND LOUVEAUX OPTIMALITY CUT TO MASTER PROBLEM
          status = addNewRowToMaster(env, lp_master, masterprobPtr, solnX, expmipObjval, SUBPROB_LB);
          if (status) {
             fprintf (stderr, "\nd2algMain:\n");
      	     fprintf (stderr, " Failed to add new pi row to subproblem lp\n");
      	     goto TERMINATE;
          }
          MIP_SOLVE = 0;
     }
     /////////////////////////////////////////////////////////////

     //////////
      //if (n_iterations == 0)
          prev_lb = objval_m;

      ////////

     /////////////////////////////////////////////////////////////////////////////////////
      //		      Solve MASTER problem                     		         //
      /////////////////////////////////////////////////////////////////////////////////////
      if (DISPLAY_LOG_INFO) {
            fprintf (stdout, "\n Optimizing master problem MIP...\n");
            fprintf (stdout, " Current accumulated master problem MIP CPU time: %f seconds\n", cputime_master);
      }


     //*********** OPTIMIZE UPDATED MASTER PROBLEM ******************
//     printf("Re-optizing Master Problem...\n");
     //status = CPXmipopt (env, lp_master);
     cputime_start = clock();
     status = CPXmipopt (env, lp_master);
     cputime_stop = clock();
     cputime_master += (double)( cputime_stop - cputime_start )/CLOCKS_PER_SEC;
     if ( status ) {
        fprintf (stderr, "Failed to optimize MIP.\n");
        goto TERMINATE;
     }


//     printf ("finishing optimizing lp_master again. \n");
     status  = CPXgetmipobjval (env, lp_master, &objval_m);
     if ( status ) {
         fprintf (stderr,"##Failed to obtain master objective value.\n");
         fprintf (stderr,"Error Code: %d\n", status);
         goto TERMINATE;
     }
     fprintf (fpSolnOut, "\nMaster Obj value = %.10g\n", objval_m);
     curr_percent_gap = 100.0*(BOUND_V_k - objval_m)/BOUND_V_k;
     fprintf(fpSolnOut, "\n Percent gap = %0.3f%% \n", curr_percent_gap);
     fprintf(stdout, "\n ***Percent gap = %0.5f%% \n", curr_percent_gap);
     fprintf(fpSolnOut, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);
     fprintf(stdout, "[lb, ub] = [%f, %f]  \n", objval_m, BOUND_V_k);

     //***** The size of the problem should be obtained by asking CPLEX what
     //      the actual size is. cur_numcols stores the current number of columns****
     cur_numcols = CPXgetnumcols (env, lp_master);
     status = CPXgetmipx (env, lp_master, solnX, 0, cur_numcols-1);
     if ( status ) {
        fprintf (stderr, "Failed to obtain solution.\n");
        goto TERMINATE;
     }

     //* Write out the solution
     fprintf(fpSolnOut, "Master Problem Solution:\n");
   // fprintf(stdout, "Master Problem Solution:\n");
    // fprintf(stdout, " obj = %f:\n", objval_m);
     for (j = 0; j < cur_numcols; j++) {
          fprintf (fpSolnOut, " %s: = %6.6f\n", masterprobPtr->colnames[j], solnX[j]);
          //fprintf (stdout, " %s: = %6.6f\n", masterprobPtr->colnames[j], solnX[j]);
     }
     //fprintf(stdout, " Previous solution:\n");
     for (j = 0; j < cur_numcols; j++) {
          fprintf (fpSolnOut, " %s: = %6.6f\n", masterprobPtr->colnames[j], solnXprev[j]);
          //fprintf (stdout, " %s: = %6.6f\n", masterprobPtr->colnames[j], solnXprev[j]);
     }

     #ifdef DEBUG_MAIN
         fprintf (stdout, "\nd2algMain():\n");
         fprintf(stdout, "[BOUND_V_k, BOUND_v_k] = [%f, %f]  \n", BOUND_V_k, objval_m);
    #endif

     //*********** Check TERMINATION CONDITION ******************
     // count algorithm iterations
     n_iterations++;


     // Accumulate computational time
     total_cputime = cputime_master + cputime_sublp + cputime_submip +
                  cputime_c3 + cputime_rhs;

     // LOG TIME to File
     if (DISPLAY_LOG_INFO) {
       fprintf(stdout, "\nCurrent Total CPU Time at iteration %d is %0.2f\n", n_iterations-1, total_cputime);
     }
     fprintf(fpSolnOut, "\nCurrent Total CPU Time at iteration %d is %0.2f\n", n_iterations-1, total_cputime);

  /*   fprintf(fpSolnOut,"\n\n<Specific CPU TIMES>: \n");
     fprintf(fpSolnOut, " Master MIP : %0.2f seconds \n", cputime_master);
     fprintf(fpSolnOut, " Subprob LP : %0.2f seconds \n", cputime_sublp);
     fprintf(fpSolnOut, " Subprob MIP: %0.2f seconds \n", cputime_submip);
     fprintf(fpSolnOut, " C^3 LP     : %0.2f seconds \n", cputime_c3);
     fprintf(fpSolnOut, " RHS LP     : %0.2f seconds \n\n", cputime_rhs);
   */
     //exit(0);

     ///////////////////////////////////////////////////////
     // Check if current master solution is same as before//
     diffMasterSoln = isEqual(solnX, solnXprev, cur_numcols-1);
	
	for (int i=0; i<cur_numcols; i++)	printf("%f ", solnX[i]);
	
     if (diffMasterSoln == 0 && MIP_SOLVE) {
        fprintf(fpSolnOut, "Incumbent solution found by terminating due to same master soln\n");
        fprintf(fpSolnOut, "Being generated \n");
        curr_percent_gap = 100.0*(BOUND_V_k - objval_m)/BOUND_V_k;
        fprintf(fpSolnOut, "\n Percent gap = %0.3f%% \n", curr_percent_gap);
        break;
     }
     if (diffMasterSoln == 0) {
        //cnt++;
        //if (cnt == 2)
             MIP_SOLVE = 1;
        if (DISPLAY_LOG_INFO) {
            fprintf (stdout, "\n Initializing scenario subproblem MIP solves...\n");
        }
        fprintf (fpSolnOut, "\n Same master solution as in previous iteration...\n");
        fprintf (fpSolnOut, "\n Initializing scenario subproblem MIP solves...\n");
        curr_lb_V = objval_m - SUBPROB_LB;
        fprintf (fpSolnOut, "\n curr_lb_V = %g\n", curr_lb_V);
        //break;
     }


     if ((myabs(prev_percent_gap - curr_percent_gap) < PERCENT_GAP_THRES) &&  curr_percent_gap < PERCENT_GAP_LIM){
          cnt++;
          if (cnt == 2) {
               MIP_SOLVE = 1;
               if (DISPLAY_LOG_INFO) {
                     fprintf (stdout, "\n ->Gap not changing significantly...\n");
                     fprintf (stdout, "\n ->Initializing subproblem MIP solves in the next iterations...\n");
                }
          }

     } else {
         cnt = 0;
     }
     if (DISPLAY_LOG_INFO) {
         fprintf(stdout, "\ Prev_percent_gap - Curr_percent_gap at iteration %d is %g\n", n_iterations-1, prev_percent_gap - curr_percent_gap);
     }

     prev_percent_gap = curr_percent_gap;



     // Store current solution
     status = CPXgetmipx (env, lp_master, solnXprev, 0, cur_numcols-1);
     if ( status ) {
        fprintf (stderr, "Failed to obtain solution.\n");
        goto TERMINATE;
     }

     fprintf (fpSolnOut, "\n curr_percent_gap > PERCENT_GAP: %g > %g\n", curr_percent_gap, PERCENT_GAP);

     if (total_cputime > 3600)
	break;

 } while (curr_percent_gap > PERCENT_GAP) ; // ******** End while loop ******

	solnPtr->best_lb = objval_m;
        solnPtr->best_ub = BOUND_V_k;

  //Write sub prob to file in lp format
   #ifdef DEBUG_MAIN
        fprintf (stderr, "\nd2algMain():\n");
   	status = CPXwriteprob(env, lp_sub, subprobname_lp, NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	}
   #endif

   fprintf(fpout, "%d	%f	%f\n  ", n_iterations, objval_m, BOUND_V_k);

  //fprintf (stdout, "\n Optimal Solution Found!!\n");







  //for (i = 0; i < nd2cuts; i++)
   //    fprintf(stdout, "subprobPtr->disjVarWindex[%d] = %d \n", i, subprobPtr->disjVarWindex[i]);

  // stop wallclock time
  wallclock_stop = clock();

  // Set total cpu time
  total_cputime = cputime_master + cputime_sublp + cputime_submip +
                  cputime_c3 + cputime_rhs + cputime_bblp;

   total_wctime = (double)( wallclock_stop - wallclock_start )/CLOCKS_PER_SEC;

   printf( "\n Optimal solution found: Percent gap = %0.3f%% \n", curr_percent_gap);
   printf( "\n Total CPU Time: %0.2f seconds \n", total_cputime );
   printf( " Wallcklock    : %0.2f seconds \n", total_wctime );

   printf("\n<Specific CPU TIMES>: \n");
   printf( " Master MIP     : %0.2f seconds \n", cputime_master);
   printf( " Subprob LP     : %0.2f seconds \n", cputime_sublp);
   printf( " Subprob MIP    : %0.2f seconds \n", cputime_submip);
   printf( " C^3 LP         : %0.2f seconds \n", cputime_c3);
   printf( " RHS LP         : %0.2f seconds \n", cputime_rhs);
   printf( " TB&B OPT CUT LP: %0.2f seconds \n", cputime_bblp);
   printf(" Algorithm Iterations                  : %d \n", n_iterations);
   printf(" Algorithm Iterations with no D^2 cuts : %d \n", no_nd2cuts_iters);
   printf(" D^2 Cuts added                        : %d \n\n", nd2cuts);
   printf(" Number of subproblem MIP solves       : %d \n\n", nmipsolves);
   printf( "\n obj = %0.2f  \n", solnPtr->obj);
   // Print Optimal Solution to File
   fprintf(fpSolnOut, "\n Percent gap = %0.3f%% \n", curr_percent_gap);
   fprintf(fpSolnOut, "\nCPU Time: %0.2f seconds \n", total_cputime );
   fprintf(fpSolnOut, "Wallclock    : %0.2f seconds \n", total_wctime );

   fprintf(fpSolnOut,"\n\n<Specific CPU TIMES>: \n");
   fprintf(fpSolnOut, "Master MIP      : %0.2f seconds \n", cputime_master);
   fprintf(fpSolnOut, "Subprob LP      : %0.2f seconds \n", cputime_sublp);
   fprintf(fpSolnOut, "Subprob MIP     : %0.2f seconds \n", cputime_submip);
   fprintf(fpSolnOut, "C^3 LP          : %0.2f seconds \n", cputime_c3);
   fprintf(fpSolnOut, "RHS LP          : %0.2f seconds \n\n", cputime_rhs);
   fprintf(fpSolnOut, " TB&B OPT CUT LP: %0.2f seconds \n", cputime_bblp);
   fprintf(fpSolnOut, " Algorithm Iterations                  : %d \n", n_iterations);
   fprintf(fpSolnOut, " Algorithm Iterations with no D^2 cuts : %d \n", no_nd2cuts_iters);
   fprintf(fpSolnOut, " D^2 Cuts added                        : %d \n\n", nd2cuts);
   fprintf(fpSolnOut," Number of subproblem MIP solves        : %d \n\n", nmipsolves);

   printSoln(masterprobPtr->ncols, masterprobPtr->colnames, solnPtr,
             solnPtr->nrows, solnPtr->ncols,
             subprobPtr->colnames, fpSolnOut);


   /////////////////////////////////////////////////////////////
   //    PRINT SUBPROBLEM LP TO SEE THE D2 CUTS               //
   /////////////////////////////////////////////////////////////
   fprintf(fpout, "subprobPtr->nrows_mip = %d \n", subprobPtr->nrows_mip);
  fprintf(fpout, "subprobPtr->nrows = %d \n", subprobPtr->nrows);
  fprintf(fpout, "subprobPtr->nrows_total = %d \n", subprobPtr->nrows_total);

   nrows_sub = CPXgetnumrows (env, lp_sub);
   fprintf(fpout, "Before: nrows_sub = %d \n", nrows_sub);
   status = CPXdelrows (env, lp_sub, 0, nrows_sub-nd2cuts-1);
   if ( status ) {
   	fprintf (stderr, "d2algMain():\n");
      	fprintf (stderr, "Failure to delete subproblem rows, error %d.\n", status);
      	goto TERMINATE;
   }
   fprintf(fpout, "nrows_sub = %d \n", nrows_sub);
   #ifdef DEBUG_MAIN
   	status = CPXwriteprob(env, lp_sub, "SUBPROB.lp", NULL);
   	if ( status ) {
   		fprintf (stderr, "d2algMain():\n");
      		fprintf (stderr, "Failure to write subproblem lp to file, error %d.\n", status);
      		goto TERMINATE;
   	 }
    #endif

    ///////////////////////////////////////////////////////////////////////////
   //     Extract the W matrix from the lp_sub in ROW SPARSE MATRIX FORMAT  //
   //     This is convenient for constructing the C^3 LP.   	            //
   ///////////////////////////////////////////////////////////////////////////
   nrows_sub = CPXgetnumrows (env, lp_sub);
   ncols_sub = CPXgetnumcols (env, lp_sub);
   fprintf(fpout, "\n\n  nnrows_sub = %d \n", nrows_sub);
   fprintf(fpout, "  ncols_sub = %d \n", ncols_sub);

   status = CPXgetrows(env, lp_sub, &subprobPtr->nzcnt_W, subprobPtr->cmatbeg_W, subprobPtr->cmatind_W,
                      subprobPtr->cmatval_W, subprobPtr->cmatspace_W, &surplus, 0, nrows_sub-1);
  if ( status ) {
      fprintf (stderr, "\nd2algMain():\n");
      fprintf (stderr, "Failure to write get rows from lp_sub for W matrix, error %d.\n", status);
      fprintf (stderr, "Surplus value is: %d.\n", surplus);
      return(status);
  }
  fprintf(fpout, "  Density = %g \n\n", 1.0*subprobPtr->nzcnt_W/(ncols_sub*nrows_sub));

  // Print the row sparse matrix W
  //fprintf(fpout, "\n\nsubprobPtr->cmatspace_W = %d\n", subprobPtr->cmatspace_W);
  fprintf(fpout, "subprobPtr->nzcnt_W = %d\n", subprobPtr->nzcnt_W);
  fprintf(fpout, "The W matrix in row sparse format is: \n");
  printWMatrix(nrows_sub, subprobPtr->nzcnt_W, subprobPtr->cmatbeg_W, subprobPtr->cmatind_W,
                   subprobPtr->cmatval_W, fpout);



    fprintf(fpout, "\n\nSUBPROBLEM T(w): \n");
    /////////////////////////////////////////////////////////////
    //    PRINT SUBPROBLEM T(w) TO SEE THE D2 CUTS RHS         //
    /////////////////////////////////////////////////////////////

    for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
         fprintf(fpout, "\nSCENARIO %d RHS VECTOR r(w): \n", scenario);
         for (j = 0; j < stochdataPtr->nrows+nd2cuts; j++)
             fprintf(fpout, " r[%d][%d] = %f\n", scenario, j, stochdataPtr->rhs[scenario][j]);
         /////////////////////////////////////////////////////////////
         //    PRINT SUBPROBLEM LP RHS TO SEE THE D2 CUTS           //
         /////////////////////////////////////////////////////////////
         fprintf(fpout, "\nSCENARIO %d TECHNOLOGY MATRIX: \n", scenario);
         printMatrix(stochdataPtr->rnrows, stochdataPtr->rnzcnt_T[scenario], stochdataPtr->rmatbeg_T[scenario],
                stochdataPtr->rmatind_T[scenario], stochdataPtr->rmatval_T[scenario], fpout);
    }


   // close the output file
   fclose(fpout);
   fclose(fpSolnOut);


   //printf("\nTo TERMINATE...\n");



TERMINATE:

   #ifdef DEBUG_MAIN
     printf("\nDEBUG output written to file: %s\n", debugoutfname);
   #endif

   printf("\nOutput written to file: %s\n\n", soln_filename);

	// Free up the SMPS file data storage structures
	freeStochFileStruct(stochdataPtr, nrows_sub, ncols_master, ncols_sub);

	// Free up the masterproblem data structures
	freeMasterProblemStructs(masterprobPtr, ncols_master);

	// Free up the subproblem data structures
	freeSubProblemStruct(subprobPtr, stochdataPtr->nscens);
	
   // Free up the subproblem data structures
   if ( solnPtr != NULL ) {
      //free ((solnMatrix_t*)solnPtr);
      solnPtr = NULL;
   }


   // Free up the solution
   free_and_null ((char **) &solnX);
   free_and_null ((char **) &solnY);
   //if ( *solnX != NULL ) {
    //  free (solnX);
    // *solnX = NULL;
   //}
   //printf("To FREE solnY \n");
   //if ( *solnY != NULL ) {
   //   free (solnY);
   //   *solnY = NULL;
  // }

  // printf("To FREE lp_master \n");

   // Free up the problem as allocated by CPXcreateprob, if necessary
   if ( lp_master != NULL ) {
      //printf("Free master lp \n");
      status = CPXfreeprob(env, &lp_master);
      //printf("freed \n");
      if ( status ) {
         fprintf (stderr, "CPXfreeprob lp_master failed, error code %d.\n", status);
      }
   }

  //  Free up the problem as allocated by CPXcreateprob, if necessary
   if ( lp_sub != NULL ) {
      //printf("Free subproblem lp \n");
      status = CPXfreeprob(env, &lp_sub);
      //printf("freed \n");
      if ( status ) {
         fprintf (stderr, "CPXfreeprob lp_sub failed, error code %d.\n", status);
      }
   }

    //  Free up the problem as allocated by CPXcreateprob, if necessary
   if ( lp_submip != NULL ) {
      //printf("Free subproblem mip \n");
      status = CPXfreeprob(env, &lp_submip);
      //printf("freed \n");
      if ( status ) {
         fprintf (stderr, "CPXfreeprob lp_submip failed, error code %d.\n", status);
      }
   }




   printf("Freeing up the CPLEX environment...\n");

   // Free up the CPLEX environment, if necessary
   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      // Note that CPXcloseCPLEX produces no output,
      //   so the only way to see the cause of the error is to use
       //  CPXgeterrorstring.  For other CPLEX routines, the errors will
       //  be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.

      if ( status ) {
         char  errmsg[1024];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }


   //printf("TO RETURN\n");



   return (0);

}

//******************************* END MAIN FUNCTION **************************************

// This simple routine frees up the pointer *ptr, and sets *ptr to NULL

#ifndef  CPX_PROTOTYPE_MIN
static void
free_and_null (char **ptr)
#else
static void
free_and_null (ptr)
char  **ptr;
#endif
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */

#ifndef  CPX_PROTOTYPE_MIN
static void
usage (char *progname)
#else
static void
usage (progname)
char *progname;
#endif
{
   fprintf (stderr,"\nUsage: %s prefix\n", progname);
   fprintf (stderr,"   where prefix is the prefix name for the STOCH files:\n");
   fprintf (stderr,"         prefix.tim is the SMPS TIME file name \n");
   fprintf (stderr,"         prefix.cor is the SMPS CORE file name \n");
   fprintf (stderr,"         prefix.sto is the SMPS STOCH file name \n");
   fprintf (stderr," Exiting...\n");
   exit(1);
} /* END usage */

void
freeAndNull (double *ptr)
/**
 * Frees the memory allocated to a double array pointed to by ptr
 *
 */
{
   if ( ptr != NULL ) {
      free (ptr);
      ptr = NULL;
   }
} // ****** END free_and_null ********
