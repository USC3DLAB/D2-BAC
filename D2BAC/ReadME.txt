
# Author: Yang Yuan
This folder contains source code for an implementation of the D2-BAC algorithms.

To compile the program edit the makefile as follows:
  change the label : YANG_D2BAC=/directory to
                     <YOUR_LABEL>=<current directory>
  where, <YOUR_LABEL> is your choice of name and <current directory> is the directory
  where you have your source files. Make sure your change the label YANG_D2BAC to <YOUR_LABEL>
  throughout the makefile.

  you might have to change the compiler (CC and cc) directory and the cplex directory
  based on your own computer setting.

  To compile the program type: make d2bac
 
  To run the program type: d2bac corefile_prefix_name outputfilename
  
  A sample output for parameter setting:

/*******************************************************************************
hopi ] d2bac sslp_5_50_100  

Enter solution output file name sequence number or name> _3

Enter a nonzero maximum number of nodes to explore in the TB&B Tree> 3

 < Set Algorithm Option A. >

  1. Enter 0 for NO Laporte and Louveaux optimality cut in the D2 Algorithm.
  2. Enter 1 to add Laporte and Louveaux optimality cut during last iteration of
     the D2 Algorithm.
  Option A:1

 < Set Algorithm Option B. >

  1. Enter 0 NOT to add D2 cuts to subproblem MIP for upper bounding.
  2. Enter 1 to add D2 cuts to subproblem MIP for upper bounding.
  Option B:0

 < Set Algorithm Option C. >

  1. Enter 0 NOT to start MIP solves immediately master obj stays constant.
  2. Enter 1 to start MIP solves immediately master obj stays constant.
  Option C:0

 < Set Algorithm Option D. >

  1. Enter 0 for NOT displaying algorithm log information.
  2. Enter 1 to display algorithm log information.
  Option D:1

 < Set Algorithm Option E. >

  1. Enter percent gap value [0-100] below which TB&B is initiated>
  Option E:10

**************************************************************************/


