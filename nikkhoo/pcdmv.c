#include "mex.h"
#include "math.h"
#define PI 3.141592653589793103245029

/*
 * pcdmv.c
 * Point Compound Dislocation Model fully vectorized
 *
 * The calling syntax is:
 *
 *	[ue,un,uv]=pcdmv(x,y,d,oX,oY,oZ,dV,a,b,nu)
 *
 * All input and output variables are scalars, vectors or N-D matrix of the
 * same number of elements, excepted nu which must be a scalar.
 *
 * This is a MEX-file for Matlab/Octave.
 *
 * Author: François Beauducel, IPGP/IRD
 * Reference: pcdm.m by F. Beauducel and A. Villié, after Nikkhoo et al. (2017)
 * Created: 2019-08-06
 * Updated: 2019-08-16
 */



/*
 * PTDdispSurf calculates surface displacements associated with a tensile 
 * point dislocation (PTD) in an elastic half-space (Okada, 1985).
 */

void PTDdispSurf(double x, double y, double d, double beta, double dip, 
		double dv, double nu, double* ue, double* un, double* uv) {
	
	double tmp, r;
	
	tmp = x * cos(beta) - y * sin(beta);
	y = x * sin(beta) + y * cos(beta);
	x = tmp;
	r = sqrt(x * x + y * y + d * d);
	tmp = 3 * pow(y*sin(dip) - d*cos(dip), 2)/pow(r, 5);

	dv /= 2*PI;
	*ue = dv * x * ( tmp - pow(sin(dip), 2) * (1 - 2*nu) * (
				 1/pow(r,3) - 1/(r * pow(r + d, 2))
				 + y * y * (3*r + d)/(pow(r, 3)*pow(r + d, 3)))
		       );
	*un = dv * y * ( tmp - pow(sin(dip), 2) * (1 - 2*nu) * (
				 1/(r * pow(r+d, 2)) 
				 - x * x * (3*r + d)/(pow(r, 3)*pow(r + d, 3)))
		       );
	*uv = dv * ( d  * tmp - pow(sin(dip), 2) * (1 - 2*nu) * ( 
			     1/(r * (r + d))
			     - x * x * (2*r + d)/(pow(r, 3)*pow(r + d, 2)))
		   );

	tmp =   *ue * cos(beta) + *un * sin(beta);
	*un = - *ue * sin(beta) + *un * cos(beta);
	*ue = tmp;

}

/*
 * pcdmv is the main function
 * NOTE: numel is the number or elements of input matrix so main loop uses i as
 * a linear index shift whatever the matrix dimensions are.
 */
void pcdm(double* x, double* y, double* d, double* ox, double* oy, double* oz,
		double* dv, double* a, double* b, double nu,
		double* ue, double* un, double* uv, size_t numel) {
	mwSize i;
	double dvx, dvy, dvz,
	       ax, ay, az, R1[3], R2[3], R3[3], norm, strike,
	       ue1, un1, uv1,
	       ue2, un2, uv2,
	       ue3, un3, uv3;

	for (i = 0; i < numel; i++) {
		/* recomputes 3 potencies from DVtot, A and B */
		dvz = *(dv + i) * *(a + i);
		dvy = ( *(dv + i) - dvz ) * *(b + i);
		dvx = ( *(dv + i) - dvz ) * ( 1 - *(b + i) );

		/* converts angles in radian */
		ax = *(ox + i) * PI / 180;
		ay = *(oy + i) * PI / 180;
		az = *(oz + i) * PI / 180;
		
		/* 3-D matrix of rotation */
		R1[0] =  cos(ay) * cos(az);
		R1[1] = -cos(ay) * sin(az);
		R1[2] =  sin(ay);
		R2[0] =  cos(az) * sin(ax) * sin(ay) + cos(ax) * sin(az);
		R2[1] = -sin(ax) * sin(ay) * sin(az) + cos(ax) * cos(az);
		R2[2] = -sin(ax) * cos(ay);
		R3[0] = -cos(ax) * sin(ay) * cos(az) + sin(ax) * sin(az);
		R3[1] =  cos(ax) * sin(ay) * sin(az) + sin(ax) * cos(az);
		R3[2] =  cos(ax) * cos(ay);

		/* calculates contribution of the first PTD */
		if (dvx != 0) {
			norm = sqrt(R1[0] * R1[0] + R1[1] * R1[1]);
			if ( norm != 0 ) strike = atan2(-R1[1]/norm, R1[0]/norm);
			else strike = 0;
			PTDdispSurf(*(x + i), *(y + i), *(d + i),
				strike - PI/2, acos(R1[2]), dvx, nu, &ue1, &un1, &uv1);
		} else {
			ue1 = 0;
			un1 = 0;
			uv1 = 0;
		}
		
		/* calculates contribution of the second PTD */
		if (dvy != 0) {
			norm = sqrt(R2[0] * R2[0] + R2[1] * R2[1]);
			if ( norm != 0 ) strike = atan2(-R2[1]/norm, R2[0]/norm);
			else strike = 0;
			PTDdispSurf(*(x + i), *(y + i), *(d + i),
				strike - PI/2, acos(R2[2]), dvy, nu, &ue2, &un2, &uv2);
		} else {
			ue2 = 0;
			un2 = 0;
			uv2 = 0;
		}
		
		/* calculates contribution of the third PTD */
		if (dvz != 0) {
			norm = sqrt(R3[0] * R3[0] + R3[1] * R3[1]);
			if ( norm != 0 ) strike = atan2(-R3[1]/norm, R3[0]/norm);
			else strike = 0;
			PTDdispSurf(*(x + i), *(y + i), *(d + i),
				strike - PI/2, acos(R3[2]), dvz, nu, &ue3, &un3, &uv3);
		} else {
			ue3 = 0;
			un3 = 0;
			uv3 = 0;
		}
		
		/* outputs */
		*(ue + i) = ue1 + ue2 + ue3;
		*(un + i) = un1 + un2 + un3;
		*(uv + i) = uv1 + uv2 + uv3;
	}
}

/*
 * the gateway function for MEX
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	double *x, *y, *d, *ox, *oy, *oz, *dv, *a, *b, *nu, *ue, *un, *uv;
	mwSize ndim;
	const mwSize *sz;
	size_t numel;
	int n;

	/* check for proper number of arguments */
	if (nrhs != 10)
		mexErrMsgIdAndTxt("MATLAB:pcdmv:invalidNumInputs", "Ten inputs required.");
	if (nlhs != 3)
		mexErrMsgIdAndTxt("MATLAB:pcdmv:invalidNumOutputs", "Three outputs required.");

	/*  get the dimensions of the matrix input x */
	ndim = mxGetNumberOfDimensions(prhs[0]);
	sz = mxGetDimensions(prhs[0]);
	numel = mxGetNumberOfElements(prhs[0]);

	/* check to make sure all input arguments are real, double matrix and same size as x */
	for (n = 0; n < nrhs; n++) {
		if (!mxIsDouble(prhs[n]) 
			|| (n < 9 && (mxGetNumberOfElements(prhs[n]) != numel ))
			|| (n == 9 && !mxIsScalar(prhs[n])))
			mexErrMsgIdAndTxt("MATLAB:pcdmv:fieldNotRealMatrix",
				"All input arguments must be real, double matrix with same size.");
	}

	/*  create pointers to each of the input matrices */
	x  = mxGetPr(prhs[0]);
	y  = mxGetPr(prhs[1]);
	d  = mxGetPr(prhs[2]);
	ox = mxGetPr(prhs[3]);
	oy = mxGetPr(prhs[4]);
	oz = mxGetPr(prhs[5]);
	dv = mxGetPr(prhs[6]);
	a  = mxGetPr(prhs[7]);
	b  = mxGetPr(prhs[8]);
	nu = mxGetPr(prhs[9]);

	/*  create the output matrices */
	plhs[0] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);

	/*  create pointers to a copy of the output matrices */
	ue = mxGetPr(plhs[0]);
	un = mxGetPr(plhs[1]);
	uv = mxGetPr(plhs[2]);

	/*  call the C subroutine */
	pcdm(x, y, d, ox, oy, oz, dv, a, b, *nu, ue, un, uv, numel);
}
