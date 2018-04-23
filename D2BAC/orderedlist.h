    #include <stdio.h> 
    #include <ctype.h>
    #include <stdlib.h>
    #include <malloc/malloc.h>
    
    #define MAXDEPTH 20
    
    
    /**
     * Node structure
     */
    typedef struct node
    {
        double objval;          /* Store some objective value */
        int    node_id;         /* The index number of this node in the B&B tree = # of vars branched on */
        int    node_depth; 	/* The depth of this node in the branch and bound tree */
        int frac_index;		/* Store next branching variable index */
        double frac_value;	/* Store next branching value */
        int    *varind;         /* Store branching variable indices */
        double *varval;         /* Store branching variable matrix coef */
        double *varbnds;        /* Store branching variable rhs bound */
        double nu;              /* Store nu coef for form reverse polar LP */
        double *gamma;          /* Store gamma values for form reverse polar LP*/
        struct node *next;      /* Pointer to the next node */
    } NODE;
    
    /**
     * Ordered list structure: list of nodes of type NODE
     */
    typedef struct list
    {  
        int size;           /* Number of nodes in the list */
        NODE *front;        /* Front of the list */
    } ORDEREDLIST;
    
    
    /**
     * This function initialize the ordered list
     * @param list pointer to the lost
     */
    int initialise(ORDEREDLIST *list);
    
    /**
     * This function adds a node to the ordered list
     * @param list pointer to the lost
     * @param newnode pointer to the new node to add to list
     * @return returns nonzero on success, zero on failure
     */
    int add(ORDEREDLIST *list, NODE *newnode);
    
    /**
     * This function removes the node at the front of ordered node list 
     * @param list pointer to the lost
     */
    void removefrontnode(ORDEREDLIST *list);
    
    /**
     * This function checks if ordered list is empty
     * @param list pointer to the lost
     * @return returns nonzero if list is empty, otherwise returns zero
     */
    int isempty(ORDEREDLIST *list);
    
    /**
     * This function creates a node and sets its member values
     * @param objval pointer to the lost
     * @param nodedepth pointer to the new node to add to list
     * @param rmatind pointer to array of branching variable indices
     * @param rmatval pointer to array of branching variable matrix coefs
     * @param rhsbnds pointer to array of branching variable rhs bounds
     * @return returns a pointer to the newly created node
     */
    NODE *createnode(double objval, int nodedepth, int max_numnodes, int *varind,
                    double *varval, double * varbnds, int ncols_A, double nu, double *gamma);
                 

