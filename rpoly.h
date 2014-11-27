/*      rpoly.cpp -- Jenkins-Traub real polynomial root finder.
 *
 *      (C) 2000, C. Bond.  All rights reserved.
 *
 *      Translation of TOMS493 from FORTRAN to C. This
 *      implementation of Jenkins-Traub partially adapts
 *      the original code to a C environment by restruction
 *      many of the 'goto' controls to better fit a block
 *      structured form. It also eliminates the global memory
 *      allocation in favor of local, dynamic memory management.
 *
 *      The calling conventions are slightly modified to return
 *      the number of roots found as the function value.
 *
 *      INPUT:
 *      op - double precision vector of coefficients in order of
 *              decreasing powers.
 *      degree - integer degree of polynomial
 *
 *      OUTPUT:
 *      zeror,zeroi - output double precision vectors of the
 *              real and imaginary parts of the zeros.
 *
 *      RETURN:
 *      returnval:   -1 if leading coefficient is zero, otherwise
 *                  number of roots found. 
 */


/*extern double *p_,*qp,*k,*qk,*svk;
extern double sr,si,u,v_,a_,b_,c,d_,a1,a2;
extern double a3,a6,a7,e_,f,g,h,szr,szi,lzr,lzi;
extern double eta,are,mre;
extern int n,nn,nmi,zerok;*/

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

int rpoly(double *op, int degree, double *zeror, double *zeroi);

//This is not part of the original rpoly. It was added by me.
//gives the roots of the cubic equation: a*x^2 + b*x + c = 0 (returns # of real roots found)
int quadratic_solver(double a, double b, double c, double *sol1, double *sol2);

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
}
#endif
