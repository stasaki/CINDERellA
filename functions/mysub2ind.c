#include "mex.h"
/* sub2ind
/* Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
/* Code covered by the 3-clause BSD License

/* nlhs = output number*/
/* plhs = output array consisting pointers*/
/* nrhs = input number */
/* prhs = input array consisting pointers*/

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
	{
        int i, nrows, ncols, maxN, indx;
        double *xp, *yp;
        
        maxN = 1;
        /*nrows = mxGetM(prhs[0]);*/
        ncols = mxGetN(prhs[0])-1;
        
               
        xp = mxGetPr(prhs[0]);
        yp = mxGetPr(prhs[1]);
        
        indx = *(yp);
        for (i = 0; i < ncols; i++){
            maxN = maxN * *(xp+i);
            indx += maxN * *(yp+i+1);
            indx -= maxN;
            /*mexPrintf("%1f\n",*(xp+i));*/
        }
		
        plhs[0]=mxCreateDoubleScalar(indx);
	}
    