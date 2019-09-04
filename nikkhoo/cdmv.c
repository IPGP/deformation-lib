#include "mex.h"
#include "math.h"
#define PI 3.141592653589793103245029
#define EPS 1E-15

/*
 * cdmv.c
 * Compound Dislocation Model fully vectorized
 *
 * The calling syntax is:
 *
 *	[ue,un,uv,dv]=cdmv(x,y,d,oX,oY,oZ,aX,aY,aZ,opening,nu)
 *
 * All input and output variables are scalars, vectors or N-D matrix of the
 * same number of elements, excepted nu which must be a scalar.
 *
 * This is a MEX-file for Matlab/Octave.
 *
 * Author: François Beauducel, IPGP/IRD
 * Reference: cdmv.m by F. Beauducel and A. Villié, after Nikkhoo et al. (2017)
 * Created: 2019-08-31
 * Updated: 2019-09-03
 */

/*
 * AngDisDispSurf calculates the displacements associated with an angular
 * dislocation in a half-space.
 */
void AngDisDispSurf(double y1, double y2, double beta, double b1, double b2, double b3,
		double nu2, double a, double* v1, double* v2, double* v3) {

	double sinB, cosB, cotB, z1, z3, r, Fi,
		v1b1, v2b1, v3b1,
		v1b2, v2b2, v3b2,
		v1b3, v2b3, v3b3;

	sinB = sin(beta);
	cosB = cos(beta);
	cotB = cosB/sinB;
	z1 = y1*cosB + a*sinB;
	z3 = y1*sinB - a*cosB;
	r = sqrt(pow(y1,2) + pow(y2,2) + pow(a,2));

	Fi = 2*atan2(y2,(r + a)*cos(beta/2)/sin(beta/2) - y1); /* The Burgers function */

	v1b1 = b1*((1 - nu2*pow(cotB,2))*Fi 
			+ y2/(r + a)*(nu2*(cotB + y1/2/(r + a)) - y1/r)
			- y2*(r*sinB - y1)*cosB/r/(r - z3));
	v2b1 = b1*(nu2*((.5 + pow(cotB,2))*log(r+a) - cotB/sinB*log(r - z3))
			- 1/(r + a)*(nu2*(y1*cotB - a/2 - pow(y2,2)/2/(r + a)) + pow(y2,2)/r)
			+ pow(y2,2)*cosB/r/(r - z3));
	v3b1 = b1*(nu2*Fi*cotB + y2/(r + a)*(1 - nu2 + a/r) - y2*cosB/(r - z3)*(cosB + a/r));

	v1b2 = b2*(-nu2*((.5-pow(cotB,2))*log(r + a) + pow(cotB,2)*cosB*log(r - z3)) 
			- 1/(r + a)*(nu2*(y1*cotB + a/2 + pow(y1,2)/2/(r + a)) - pow(y1,2)/r)
			+ z1*(r*sinB - y1)/r/(r - z3));
	v2b2 = b2*((1 + nu2*pow(cotB,2))*Fi - y2/(r + a)*(nu2*(cotB + y1/2/(r + a)) - y1/r)
			- y2*z1/r/(r - z3));
	v3b2 = b2*(-nu2*cotB*(log(r + a) - cosB*log(r - z3)) 
			- y1/(r + a)*(1 - nu2 + a/r) + z1/(r - z3)*(cosB + a/r));

	v1b3 = b3*(y2*(r*sinB - y1)*sinB/r/(r - z3));
	v2b3 = b3*(-pow(y2,2)*sinB/r/(r - z3));
	v3b3 = b3*(Fi + y2*(r*cosB + a)*sinB/r/(r - z3));

	*v1 = (v1b1 + v1b2 + v1b3)/2/PI;
	*v2 = (v2b1 + v2b2 + v2b3)/2/PI;
	*v3 = (v3b1 + v3b2 + v3b3)/2/PI;

}


/*
 * AngSetupSurf calculates the displacements associated with an angular
 * dislocation pair on each side of an RD in a half-space.
 */

void AngSetupFSC(double x, double y, double bX, double bY, double bZ, 
		double* PA, double* PB, double nu2, double* ue, double* un, double* uv) {
	
	double A1, A2, norm, beta,
	       y1A, y2A, y1B, y2B, 
	       b1, b2, b3,
	       v1A, v2A, v3A,
	       v1B, v2B, v3B,
	       v1, v2, v3;

	beta = acos((PA[2]-PB[2])/sqrt(pow(PA[0]-PB[0],2) + pow(PA[1]-PB[1],2) + pow(PA[2]-PB[2],2)));

	if (fabs(beta)<EPS || fabs(PI - beta)<EPS) {
		*ue = 0;
		*un = 0;
		*uv = 0;
	} else {
		/* A is the sparse 3x3 transformation matrix (A33 = -1, other elements are 0) */
		norm = sqrt(pow(PB[0]-PA[0],2) + pow(PB[1]-PA[1],2));
		A1 = (PB[0]-PA[0])/norm; /* = A11 = -A22 */
		A2 = (PB[1]-PA[1])/norm; /* = A12 = A21 */

		/* Transform coordinates from EFCS to the first ADCS */
		y1A = A1*(x - PA[0]) + A2*(y - PA[1]);
		y2A = A2*(x - PA[0]) - A1*(y - PA[1]);

		/* Transform coordinates from EFCS to the second ADCS */
		y1B = y1A - (A1*(PB[0] - PA[0]) + A2*(PB[1] - PA[1]));
		y2B = y2A - (A2*(PB[0] - PA[0]) - A1*(PB[1] - PA[1]));

		/* Transform slip vector components from EFCS to ADCS */
		b1 = A1*bX + A2*bY;
		b2 = A2*bX - A1*bY;
		b3 = -bZ;

		/* artefact-free for the calculation points near the free surface */
		if (beta*y1A >= 0) {
			AngDisDispSurf(y1A, y2A, beta-PI, b1, b2, b3, nu2, -PA[2], &v1A, &v2A, &v3A);
			AngDisDispSurf(y1B, y2B, beta-PI, b1, b2, b3, nu2, -PB[2], &v1B, &v2B, &v3B);
		} else {
			AngDisDispSurf(y1A, y2A, beta, b1, b2, b3, nu2, -PA[2], &v1A, &v2A, &v3A);
			AngDisDispSurf(y1B, y2B, beta, b1, b2, b3, nu2, -PB[2], &v1B, &v2B, &v3B);
		}

		/* Calculate total displacements in ADCS */
		v1 = v1B - v1A;
		v2 = v2B - v2A;
		v3 = v3B - v3A;

		/* Transform total displacements from ADCS to EFCS */
		*ue = A1*v1 + A2*v2;
		*un = A2*v1 - A1*v2;
		*uv = -v3;
	}
	
}

/*
 * RDdispSurf calculates surface displacements associated with a rectangular
 * dislocation in an elastic half-space.
 */

void RDdispSurf(double x, double y, double P1[3], double P2[3], double P3[3], double P4[3], 
		double op, double nu2, double* ue, double* un, double* uv) {
	
	double norm,
	       bX, bY, bZ,
	       u1, v1, w1,
	       u2, v2, w2,
	       u3, v3, w3,
	       u4, v4, w4;

	bX = (P2[1]-P1[1])*(P4[2]-P1[2]) - (P2[2]-P1[2])*(P4[1]-P1[1]);
	bY = (P2[2]-P1[2])*(P4[0]-P1[0]) - (P2[0]-P1[0])*(P4[2]-P1[2]);
	bZ = (P2[0]-P1[0])*(P4[1]-P1[1]) - (P2[1]-P1[1])*(P4[0]-P1[0]);
	norm = sqrt(pow(bX,2) + pow(bY,2) + pow(bZ,2));
	bX *= op/norm;
	bY *= op/norm;
	bZ *= op/norm;

	AngSetupFSC(x, y, bX, bY, bZ, P1, P2, nu2, &u1, &v1, &w1); /* Side P1P2 */
	AngSetupFSC(x, y, bX, bY, bZ, P2, P3, nu2, &u2, &v2, &w2); /* Side P2P3 */
	AngSetupFSC(x, y, bX, bY, bZ, P3, P4, nu2, &u3, &v3, &w3); /* Side P3P4 */
	AngSetupFSC(x, y, bX, bY, bZ, P4, P1, nu2, &u4, &v4, &w4); /* Side P4P1 */

	*ue = u1 + u2 + u3 + u4;
	*un = v1 + v2 + v3 + v4;
	*uv = w1 + w2 + w3 + w4;

	/* printf("\nue = %g, un = %g, uv = %g\n\n",ue,un,uv); */
}

/*
 * cdmv is the main function
 * NOTE: numel is the number or elements of input matrix so main loop uses i as
 * a linear index shift whatever the matrix dimensions are.
 */
void cdm(double* x, double* y, double* d, double* ox, double* oy, double* oz,
		double* ax2, double* ay2, double* az2, double* op, double nu,
		double* ue, double* un, double* uv, double* dv, size_t numel) {
	mwSize i;
	double ax, ay, az,
	       oxd, oyd, ozd, norm, strike,
	       R11, R12, R13, R21, R22, R23, R31, R32, R33,
	       P1[3], P2[3], P3[3], P4[3],
	       Q1[3], Q2[3], Q3[3], Q4[3],
	       R1[3], R2[3], R3[3], R4[3],
	       ue1, un1, uv1,
	       ue2, un2, uv2,
	       ue3, un3, uv3;

	for (i = 0; i < numel; i++) {

		/* converts angles into radian */
		oxd = *(ox + i) * PI / 180;
		oyd = *(oy + i) * PI / 180;
		ozd = *(oz + i) * PI / 180;

		/* converts semi-axes to full axes */
		ax = *(ax2 + i)*2;
		ay = *(ay2 + i)*2;
		az = *(az2 + i)*2;

		/* converts nu to 1-2*nu */
		nu = 1 - 2*nu;
		
		/* 3-D matrix of rotation */
		R11 =  cos(oyd) * cos(ozd);
		R12 = -cos(oyd) * sin(ozd);
		R13 =  sin(oyd);
		R21 =  cos(ozd) * sin(oxd) * sin(oyd) + cos(oxd) * sin(ozd);
		R22 = -sin(oxd) * sin(oyd) * sin(ozd) + cos(oxd) * cos(ozd);
		R23 = -sin(oxd) * cos(oyd);
		R31 = -cos(oxd) * sin(oyd) * cos(ozd) + sin(oxd) * sin(ozd);
		R32 =  cos(oxd) * sin(oyd) * sin(ozd) + sin(oxd) * cos(ozd);
		R33 =  cos(oxd) * cos(oyd);

		/* coordinates for each RD summits */
		P1[0] = ay*R21/2 + az*R31/2;
		P1[1] = ay*R22/2 + az*R32/2;
		P1[2] = ay*R23/2 + az*R33/2 - *(d + i);
		P2[0] = P1[0] - ay*R21;
		P2[1] = P1[1] - ay*R22;
		P2[2] = P1[2] - ay*R23;
		P3[0] = P2[0] - az*R31;
		P3[1] = P2[1] - az*R32;
		P3[2] = P2[2] - az*R33;
		P4[0] = P1[0] - az*R31;
		P4[1] = P1[1] - az*R32;
		P4[2] = P1[2] - az*R33;

		Q1[0] = -ax*R11/2 + az*R31/2;
		Q1[1] = -ax*R12/2 + az*R32/2;
		Q1[2] = -ax*R13/2 + az*R33/2 - *(d + i);
		Q2[0] = Q1[0] + ax*R11;
		Q2[1] = Q1[1] + ax*R12;
		Q2[2] = Q1[2] + ax*R13;
		Q3[0] = Q2[0] - az*R31;
		Q3[1] = Q2[1] - az*R32;
		Q3[2] = Q2[2] - az*R33;
		Q4[0] = Q1[0] - az*R31;
		Q4[1] = Q1[1] - az*R32;
		Q4[2] = Q1[2] - az*R33;

		R1[0] = ax*R11/2 + ay*R21/2;
		R1[1] = ax*R12/2 + ay*R22/2;
		R1[2] = ax*R13/2 + ay*R23/2 - *(d + i);
		R2[0] = R1[0] - ax*R11;
		R2[1] = R1[1] - ax*R12;
		R2[2] = R1[2] - ax*R13;
		R3[0] = R2[0] - ay*R21;
		R3[1] = R2[1] - ay*R22;
		R3[2] = R2[2] - ay*R23;
		R4[0] = R1[0] - ay*R21;
		R4[1] = R1[1] - ay*R22;
		R4[2] = R1[2] - ay*R23;

		/* the CDM must be under the free surface */
		if ( P1[2]>0 || P2[2]>0 || P3[2]>0 || P4[2]>0 
		  || Q1[2]>0 || Q2[2]>0 || Q3[2]>0 || Q4[2]>0 
		  || Q1[2]>0 || Q2[2]>0 || Q3[2]>0 || Q4[2]>0 ) { 
			*(ue + i) = 0.0/0.0; /* NaN */
			*(un + i) = 0.0/0.0; /* NaN */
			*(uv + i) = 0.0/0.0; /* NaN */
		} else if (ax==0 && ay==0 && az==0) {
			*(ue + i) = 0;
			*(un + i) = 0;
			*(uv + i) = 0;
		} else if (ax==0 && ay>0 && az>0) {
			RDdispSurf(*(x + i), *(y + i), P1, P2, P3, P4, *(op + i), nu, (ue + i), (un + i), (uv + i));
		} else if (ax>0 && ay==0 && az>0) {
			RDdispSurf(*(x + i), *(y + i), Q1, Q2, Q3, Q4, *(op + i), nu, (ue + i), (un + i), (uv + i));
		} else if (ax>0 && ay>0 && az==0) {
			RDdispSurf(*(x + i), *(y + i), R1, R2, R3, R4, *(op + i), nu, (ue + i), (un + i), (uv + i));
		} else {
			RDdispSurf(*(x + i), *(y + i), P1, P2, P3, P4, *(op + i), nu, &ue1, &un1, &uv1);
			RDdispSurf(*(x + i), *(y + i), Q1, Q2, Q3, Q4, *(op + i), nu, &ue2, &un2, &uv2);
			RDdispSurf(*(x + i), *(y + i), R1, R2, R3, R4, *(op + i), nu, &ue3, &un3, &uv3);
			*(ue + i) = ue1 + ue2 + ue3;
			*(un + i) = un1 + un2 + un3;
			*(uv + i) = uv1 + uv2 + uv3;
		}
		
		/* calculates the CDM potency */
		*(dv + i) = *(op + i) * (ax*ay + ax*az + ay*az);
		
	}
}

/*
 * the gateway function for MEX
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	double *x, *y, *d, *ox, *oy, *oz, *ax, *ay, *az, *op, *nu, *ue, *un, *uv, *dv;
	mwSize ndim;
	const mwSize *sz;
	size_t numel;
	int n;

	/* check for proper number of arguments */
	if (nrhs != 11)
		mexErrMsgIdAndTxt("MATLAB:pcdmv:invalidNumInputs", "Eleven inputs required.");
	if (nlhs != 4)
		mexErrMsgIdAndTxt("MATLAB:pcdmv:invalidNumOutputs", "Four outputs required.");

	/*  get the dimensions of the matrix input x */
	ndim = mxGetNumberOfDimensions(prhs[0]);
	sz = mxGetDimensions(prhs[0]);
	numel = mxGetNumberOfElements(prhs[0]);

	/* check to make sure all input arguments are real, double matrix and same size as x */
	for (n = 0; n < nrhs; n++) {
		if (!mxIsDouble(prhs[n]) 
			|| (n < 10 && (mxGetNumberOfElements(prhs[n]) != numel ))
			|| (n == 10 && !mxIsScalar(prhs[n])))
			mexErrMsgIdAndTxt("MATLAB:pcdmv:fieldNotRealMatrix",
				"All input arguments must be real, double matrix with same size excepted NU that must be a scalar.");
	}

	/*  create pointers to each of the input matrices */
	x  = mxGetPr(prhs[0]);
	y  = mxGetPr(prhs[1]);
	d  = mxGetPr(prhs[2]);
	ox = mxGetPr(prhs[3]);
	oy = mxGetPr(prhs[4]);
	oz = mxGetPr(prhs[5]);
	ax = mxGetPr(prhs[6]);
	ay = mxGetPr(prhs[7]);
	az = mxGetPr(prhs[8]);
	op = mxGetPr(prhs[9]);
	nu = mxGetPr(prhs[10]);

	/*  create the output matrices */
	plhs[0] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
	plhs[3] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);

	/*  create pointers to a copy of the output matrices */
	ue = mxGetPr(plhs[0]);
	un = mxGetPr(plhs[1]);
	uv = mxGetPr(plhs[2]);
	dv = mxGetPr(plhs[3]);

	/*  call the C subroutine */
	cdm(x, y, d, ox, oy, oz, ax, ay, az, op, *nu, ue, un, uv, dv, numel);
}
