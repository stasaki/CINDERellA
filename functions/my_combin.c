#include <stdio.h>
#include "mex.h"

/* return the combination number
 * all the combination under k (include 0)
 * last updated by S.TASKI 2012/8/8
 */


unsigned long long pascal(int n, int r, int dup);    /* Pascal*/


void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    int n = mxGetScalar(prhs[0]); /* The size of the set; for {1, 2, 3, 4} it's 4 */
    int k = mxGetScalar(prhs[1]); /* The size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
    /*int comb[16]; 
    /* comb[i] is the index of the i-th element in the
     * combination */
    int s=2;
    
    
    if (n==0){
    plhs[0] = mxCreateDoubleScalar(1);
    }else{
        if (n<k){
            k=n;
        }
    int i;
    unsigned long long pascal1=0;
    for (i=0; i<=k; i++ ) {
        pascal1 = pascal1+pascal(n,i,0);
    }
    
    plhs[0] = mxCreateDoubleScalar(pascal1);
    }
}

/* �?��せ数 nCr の計算[Pascalの3角形による方法] */
unsigned long long pascal(int n, int r, int dup)
{
    int nn;
    unsigned long long* a = (unsigned long long*)calloc(n,sizeof(unsigned long long));
    unsigned long long retv;
    
    nn = (dup) ? (n+r-1) : n;
    if ( nn-r < r ) r = nn-r;
    if ( r == 0 ) return 1;
    if ( r == 1 ) return nn;
    int i;
    for ( i=1; i<r; i++ ) a[i] = i+2;
    
    for ( i=3; i<=nn-r+1; i++ ) {
        a[0] = i;
        int j;
        for ( j=1; j<r; j++ ) a[j] += a[j-1];
    }
    retv = a[r-1];
    free(a);
    
    return retv;
}
