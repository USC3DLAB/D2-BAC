/** This program contains function for an implementation of an ordered linked list.
 *  Ordering is by nondecreasing objective value of the list node struct variable
 *  
 *   Headers: orderedlist.h
 * 
 *  Author: Lewis Ntaimo
 *  Date:   December 21, 2003
 **/
   
 #include "orderedlist.h"   
  
 /* Functions */
  
 /**
  * This function initialize the ordered list
  * @param list pointer to the lost
  */
  int initialise(ORDEREDLIST *list)
   {
        if ((list->front = (NODE *)malloc(sizeof(NODE))) == NULL)
        return 0;
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
        free(oldnode);
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
   * @param rmatind pointer to array of branching variable indices
   * @param rmatval pointer to array of branching variable matrix coefs
   * @param rhsbnds pointer to array of branching variable rhs bounds
   * @return returns a pointer to the newly created node
   */
   NODE *createnode2(double objval, int nodedepth, int *rmatind,
                    double *rmatval, double * rhsbnds)
   {
        int j;
        NODE * newnode;
        
        printf("Allocating memory to node vars...\n");
        newnode = (NODE *)malloc(sizeof(NODE));
        newnode->rmatind = (int *)malloc(nodedepth*sizeof(int));
        newnode->rmatval = (double *)malloc(nodedepth*sizeof(double));
        newnode->rhsbnds = (double *)malloc(nodedepth*sizeof(double));
        
        if (newnode == NULL || rmatind == NULL || rmatval == NULL || rhsbnds == NULL) {
             fprintf(stderr, "makenode(): Failure to allocate memory to new node struct. Bailing out...");
             return 0;
        }
        newnode->objval     = objval;  
        newnode->node_depth = nodedepth; 
        for (j = 0; j < nodedepth; j++) {
             newnode->rmatind[j] = rmatind[j] ;        
             newnode->rmatval[j] = rmatval[j];        
             newnode->rhsbnds[j] = rhsbnds[j]; 
        } 
        newnode->next = NULL;
           
        #ifdef DEBUG
             printf("newnode->node_depth = %3d\t", newnode->node_depth); 
             printf("copied arrays...\n");
             for (j = 0; j < nodedepth; j++) {
                 printf("newnode->rmatind[%d] = %3d\t", j, newnode->rmatind[j]); 
                 printf("newnode->rmatval[%d] = %3.3f\t", j, newnode->rmatval[j]); 
                 printf("newnode->rhsbnds[%d] = %3.3f\n", j, newnode->rhsbnds[j]); 
             }
             printf("Done copying arrays...\n"); 
        #endif

        return newnode;
 } // End createnode


/* End functions */

