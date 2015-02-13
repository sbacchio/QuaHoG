#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <qhg_su3_mul.h>
#include <qhg_su3_linalg.h>

typedef struct {
  double re;
  double im;
} qcd_complex_16;

/* preprocessor macros, use carefully */
#define qcd_CONJ(x)    ( (qcd_complex_16) {x.re,-x.im} )
#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#define qcd_CMULR(x,y) ( (double) (x.re * y.re - x.im * y.im) )
#define qcd_CMULI(x,y) ( (double) (x.re * y.im + x.im * y.re) )
#define qcd_CADJOINTMUL(x,y)  ( (qcd_complex_16) {x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re } )

#define qcd_CADD(x,y)  ( (qcd_complex_16) {x.re+y.re, x.im+y.im})
#define qcd_CADDR(x,y) ( (double) (x.re+y.re))
#define qcd_CADDI(x,y) ( (double) (x.im+y.im))

#define qcd_CSUB(x,y)  ( (qcd_complex_16) {x.re-y.re, x.im-y.im})
#define qcd_CSUBR(x,y) ( (double) (x.re-y.re))
#define qcd_CSUBI(x,y) ( (double) (x.im-y.im))

#define qcd_CSCALE(x,a)  ( (qcd_complex_16) {x.re*(a), x.im*(a)})
#define qcd_CSCALER(x,a) ( (double) x.re*(a))
#define qcd_CSCALEI(x,a) ( (double) x.im*(a))

#define qcd_CDEV(x,y)  ( (qcd_complex_16) {(x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im), (x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) } )
#define qcd_CDEVR(x,y) ( (double)     ((x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im) ))
#define qcd_CDEVI(x,y) ( (double)     ((x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) ))

#define qcd_NORM(x)    ( (double) sqrt(x.re * x.re + x.im * x.im))
#define qcd_NORMSQUARED(x) ( (double) (x.re * x.re + x.im * x.im))
#define qcd_ARG(x)     ( (double) atan2(x.im,x.re))
#define qcd_CPOW(x,a)  ( (qcd_complex_16) {pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)), pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))})
#define qcd_CPOWR(x,a) ( (double) pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)))
#define qcd_CPOWI(x,a) ( (double) pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))) 


/*
  Projects the 3x3 matrices pointed to by u to SU(3). If M is an
  arbitary 3x3 matrix it can be factorised as the product of a
  hermitian and a unitary matrix, M = U H. We will calculate U as the
  projection of M to SU(3). With M^\dagger M = H^2 and diagonalising
  H^2 we calculate H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^\dagger where
  v is the 3x3 matrix whoes columns are the eigenvectors of H^2 and L
  are the eigenvalues. Thus: U = M H^-1

  [1] Phys. Lett. B307, 375-382, (1993)
*/

void
qhg_su3_project(_Complex double *u, int n)
{
  qcd_complex_16 H[3][3],detM,U[3][3],M[3][3];
  qcd_complex_16 b,D,v[3][3],vr[3][3];
  double a,ThirdRoot_18,ThirdRoot_12,ThirdRoot_2_3;
  double trace,e[3],de,cor[3];
  qcd_complex_16 w,ThirdRootOne[2];
  double sum;
  double norm,phase;

   /*
     Constants Used:
   */
   ThirdRootOne[0].re= 1;
   ThirdRootOne[0].im= sqrt(3.);

   ThirdRootOne[1].re= 1;
   ThirdRootOne[1].im=-sqrt(3.);

   ThirdRoot_12 =pow(12.,1./3.);
   ThirdRoot_18 =pow(18.,1./3.);
   ThirdRoot_2_3=pow((2./3.),1./3.);

   for(int i=0; i<n; i++) {
     for(int c1=0; c1<3; c1++)
       for(int c2=0; c2<3; c2++)
         M[c1][c2] = (qcd_complex_16){creal(u[i*NC*NC + c1*NC + c2]),
				      cimag(u[i*NC*NC + c1*NC + c2])};
     
      
     detM = qcd_CADD(qcd_CADD( qcd_CMUL(M[0][0],qcd_CMUL(M[1][1],M[2][2])),
			       qcd_CMUL(M[0][1],qcd_CMUL(M[1][2],M[2][0]))),
		     qcd_CMUL(M[0][2],qcd_CMUL(M[1][0],M[2][1])));
     
     detM = qcd_CSUB(detM,
		     qcd_CADD(qcd_CADD( qcd_CMUL(M[0][2],qcd_CMUL(M[1][1],M[2][0])),
					qcd_CMUL(M[0][0],qcd_CMUL(M[1][2],M[2][1]))),
			      qcd_CMUL(M[0][1],qcd_CMUL(M[1][0],M[2][2]))));
     phase = qcd_ARG(detM)/3.;
     
     H[0][0].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][0]);
     H[0][1].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][1]);
     H[0][2].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][2]);
     H[0][0].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][0]);
     H[0][1].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][1]);
     H[0][2].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][2]);
     
     H[1][0].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][0]);
     H[1][1].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][1]);
     H[1][2].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][2]);
     H[1][0].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][0]);
     H[1][1].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][1]);
     H[1][2].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][2]);

     H[2][0].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][0]);
     H[2][1].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][1]);
     H[2][2].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][2]);
     H[2][0].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][0]);
     H[2][1].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][1]);
     H[2][2].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][2]);

     /*
       Assure Hermiticity:
     */
     H[0][1].re = (H[0][1].re + H[1][0].re)/2.;
     H[0][1].im = (H[0][1].im - H[1][0].im)/2.;

     H[1][0] =  qcd_CONJ(H[0][1]);

     H[0][2].re = (H[0][2].re + H[2][0].re)/2.;
     H[0][2].im = (H[0][2].im - H[2][0].im)/2.;

     H[2][0] =  qcd_CONJ(H[0][2]);

     H[1][2].re = (H[1][2].re + H[2][1].re)/2.;
     H[1][2].im = (H[1][2].im - H[2][1].im)/2.;

     H[2][1] =  qcd_CONJ(H[1][2]);
     /*
       If H^2 is alread diagonal skip diagonalization and
       calculate U directly
     */
     sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);

     if(sum<=1e-08)
       {
         e[0]=1./sqrt(H[0][0].re);
         e[1]=1./sqrt(H[1][1].re);
         e[2]=1./sqrt(H[2][2].re);

         U[0][0] = (qcd_complex_16) { M[0][0].re*e[0], M[0][0].im*e[0] };
         U[0][1] = (qcd_complex_16) { M[0][1].re*e[0], M[0][1].im*e[0] };
         U[0][2] = (qcd_complex_16) { M[0][2].re*e[0], M[0][2].im*e[0] };

         U[1][0] = (qcd_complex_16) { M[1][0].re*e[1], M[1][0].im*e[1] };
         U[1][1] = (qcd_complex_16) { M[1][1].re*e[1], M[1][1].im*e[1] };
         U[1][2] = (qcd_complex_16) { M[1][2].re*e[1], M[1][2].im*e[1] };

         U[2][0] = (qcd_complex_16) { M[2][0].re*e[2], M[2][0].im*e[2] };
         U[2][1] = (qcd_complex_16) { M[2][1].re*e[2], M[2][1].im*e[2] };
         U[2][2] = (qcd_complex_16) { M[2][2].re*e[2], M[2][2].im*e[2] };

         for(int c1=0; c1<3; c1++)
	   for(int c2=0; c2<3; c2++)
	     u[i*NC*NC + c1*NC + c2] = U[c1][c2].re + _Complex_I*U[c1][c2].im; 
       }
     else
       {
         /*
           Make traceless to eliminate second order term in eigenvalue equation,
           i.e. eig^3 + a eig + b = 0, when H is traceless.
         */
         trace=(H[0][0].re+H[1][1].re+H[2][2].re)/3.;

         H[0][0].re-=trace;
         H[1][1].re-=trace;
         H[2][2].re-=trace;


         /*
           Solve for eigenvalues:
           e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
           e^3 + a*e + b = 0,

           a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
           b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

           D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

           e = D/(18^(1/3)) - ((2/3)^(1/3))/D
           e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
           e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
         */
         a = -(H[2][2].re*H[2][2].re - H[0][0].re*H[1][1].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));

         b.re  = - H[0][0].re*H[1][1].re*H[2][2].re + H[2][2].re*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
         b.im  =                                      H[2][2].re*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

         b.re +=   H[0][0].re*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
         b.im +=   H[0][0].re*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

         w.re=qcd_CPOWR(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
         w.im=qcd_CPOWI(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

         D=qcd_CPOW(((qcd_complex_16){-9.*b.re + w.re,-9.*b.im + w.im}),1./3.); 

         e[0] = D.re/(ThirdRoot_18) - qcd_CDEVR(((qcd_complex_16){a*ThirdRoot_2_3,0}),D);
         e[1] = a*qcd_CDEVR(ThirdRootOne[0],((qcd_complex_16){D.re*ThirdRoot_12,D.im*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
         e[2] = -e[0]-e[1];

         e[0]+= trace;
         e[1]+= trace;
         e[2]+= trace;

         H[0][0].re+=trace;
         H[1][1].re+=trace;
         H[2][2].re+=trace;

         /*
           Eigenvectors:
           v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
           v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
           v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
         */

         v[0][0].re = -(e[0]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[0][0].im = -(e[0]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[0][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[0][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[0][2].re =-e[0]*e[0] + e[0]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].re - H[0][0].re*H[1][1].re;
         v[0][2].im = 0.;

         v[1][0].re = -(e[1]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[1][0].im = -(e[1]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[1][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[1][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[1][2].re =-e[1]*e[1] + e[1]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].re - H[0][0].re*H[1][1].re;;
         v[1][2].im = 0.;

         /*
           Assure eigenvectors orthonormality:
           norm =  inner product v1.v1
           w    = (inner product v1.v2)/norm
           v2   = w*v1
         */

         norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
         w.re  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
         w.im  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
         w.re /= norm;
         w.im /= norm;

         v[1][0].re-= qcd_CMULR(w,v[0][0]);
         v[1][0].im-= qcd_CMULI(w,v[0][0]);

         v[1][1].re-= qcd_CMULR(w,v[0][1]);
         v[1][1].im-= qcd_CMULI(w,v[0][1]);

         v[1][2].re-= qcd_CMULR(w,v[0][2]);
         v[1][2].im-= qcd_CMULI(w,v[0][2]);

         norm=1./sqrt(norm);

         /*
           Normalize first and second eigenvector:
         */

         v[0][0].re*= norm;
         v[0][0].im*= norm;

         v[0][1].re*= norm;
         v[0][1].im*= norm;

         v[0][2].re*= norm;
         v[0][2].im*= norm;


         norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

         norm=1./sqrt(norm);

         v[1][0].re*= norm;
         v[1][0].im*= norm;

         v[1][1].re*= norm;
         v[1][1].im*= norm;

         v[1][2].re*= norm;
         v[1][2].im*= norm;

         /*
           v3 = v1 x v2
         */


         v[2][0].re =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
         v[2][0].im = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

         v[2][1].re = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
         v[2][1].im = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

         v[2][2].re =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
         v[2][2].im = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

         de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
         /*
	   cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
	   cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
	   cor[2] = - cor[0] - cor[1];
         */
         //to be compatible with Grenoble & Paris, don't apply corrections
         cor[0]=0;
         cor[1]=0;
         cor[2]=0;

         de = 1./sqrt(e[0]);
         b.re = de*cos(phase-cor[0]);
         b.im =-de*sin(phase-cor[0]);
         vr[0][0] = qcd_CMUL(b,v[0][0]);
         vr[0][1] = qcd_CMUL(b,v[0][1]);
         vr[0][2] = qcd_CMUL(b,v[0][2]);

         de = 1./sqrt(e[1]);
         b.re = de*cos(phase-cor[1]);
         b.im =-de*sin(phase-cor[1]);

         vr[1][0] = qcd_CMUL(b,v[1][0]);
         vr[1][1] = qcd_CMUL(b,v[1][1]);
         vr[1][2] = qcd_CMUL(b,v[1][2]);

         de = 1./sqrt(e[2]);
         b.re = de*cos(phase-cor[2]);
         b.im =-de*sin(phase-cor[2]);

         vr[2][0] = qcd_CMUL(b,v[2][0]);
         vr[2][1] = qcd_CMUL(b,v[2][1]);
         vr[2][2] = qcd_CMUL(b,v[2][2]);


         H[0][0].re= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].re= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].re= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

         H[0][0].im= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].im= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].im= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


         H[1][0].re= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].re= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].re= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

         H[1][0].im= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].im= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].im= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


         H[2][0].re= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].re= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].re= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

         H[2][0].im= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].im= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].im= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

         U[0][0].re= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
         U[0][1].re= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
         U[0][2].re= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

         U[0][0].im= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
         U[0][1].im= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
         U[0][2].im= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


         U[1][0].re= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
         U[1][1].re= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
         U[1][2].re= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

         U[1][0].im= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
         U[1][1].im= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
         U[1][2].im= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


         U[2][0].re= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
         U[2][1].re= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
         U[2][2].re= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

         U[2][0].im= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
         U[2][1].im= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
         U[2][2].im= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


         /*
           w    = inner product: col1.col2
           norm = inner product: col1.col1
         */

         norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
         w.re  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
         w.im  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
         w.re /= norm;
         w.im /= norm;


         U[0][1].re-=qcd_CMULR(w,U[0][0]);
         U[0][1].im-=qcd_CMULI(w,U[0][0]);

         U[1][1].re-=qcd_CMULR(w,U[1][0]);
         U[1][1].im-=qcd_CMULI(w,U[1][0]);

         U[2][1].re-=qcd_CMULR(w,U[2][0]);
         U[2][1].im-=qcd_CMULI(w,U[2][0]);

         norm = 1./sqrt(norm);

         U[0][0].re*= norm;
         U[0][0].im*= norm;
         U[1][0].re*= norm;
         U[1][0].im*= norm;
         U[2][0].re*= norm;
         U[2][0].im*= norm;

         norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
         norm = 1./sqrt(norm);

         U[0][1].re*= norm;
         U[0][1].im*= norm;
         U[1][1].re*= norm;
         U[1][1].im*= norm;
         U[2][1].re*= norm;
         U[2][1].im*= norm;

         /*
           col3 = col1 x col2
         */
         U[0][2].re =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
         U[0][2].im = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

         U[1][2].re = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
         U[1][2].im =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

         U[2][2].re =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
         U[2][2].im = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

         for(int c1=0; c1<3; c1++)
	   for(int c2=0; c2<3; c2++)
	     u[i*NC*NC + c1*NC + c2] = U[c1][c2].re + _Complex_I*U[c1][c2].im; 
       }
   }
   return;
}

/*
  Projects the 3x3 matrices pointed to by u to SU(3). If M is an
  arbitary 3x3 matrix it can be factorised as the product of a
  hermitian and a unitary matrix, M = U H. We will calculate U as the
  projection of M to SU(3). With M^\dagger M = H^2 and diagonalising
  H^2 we calculate H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^\dagger where
  v is the 3x3 matrix whoes columns are the eigenvectors of H^2 and L
  are the eigenvalues. Thus: U = M H^-1

  [1] Phys. Lett. B307, 375-382, (1993)
*/

void
_qhg_su3_project(_Complex double *u, int n)
{
  for(int i=0; i<n; i++) {
    qcd_complex_16 M[NC][NC];
    _Complex double *_M = &M[0][0].re;
    su3_linalg_ueqv(_M, &u[i*NC*NC]);
	
    qcd_complex_16 detM;
    _Complex double *_detM = &detM.re;
    *_detM = su3_linalg_det_u(_M);
    double phase = carg(*_detM)/3.;

    qcd_complex_16 H[NC][NC];
    _Complex double *_H = &H[0][0].re;
    
    su3_mul_du(_H, _M, _M);
     /*
       Assure Hermiticity:
     */

    _Complex double _Hd[NC*NC];
    su3_linalg_ueqvd(_Hd, _H);
    su3_linalg_upeqv(_H, _Hd);
    su3_linalg_au(0.5, _H);
     /*
       If H^2 is alread diagonal skip diagonalization and
       calculate U directly
     */
     double sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);
     
     if(sum<=1e-08)
       {
	 double e[NC];
	 qcd_complex_16 U[NC][NC];
         e[0]=1./sqrt(H[0][0].re);
         e[1]=1./sqrt(H[1][1].re);
         e[2]=1./sqrt(H[2][2].re);

         U[0][0] = (qcd_complex_16) { M[0][0].re*e[0], M[0][0].im*e[0] };
         U[0][1] = (qcd_complex_16) { M[0][1].re*e[0], M[0][1].im*e[0] };
         U[0][2] = (qcd_complex_16) { M[0][2].re*e[0], M[0][2].im*e[0] };

         U[1][0] = (qcd_complex_16) { M[1][0].re*e[1], M[1][0].im*e[1] };
         U[1][1] = (qcd_complex_16) { M[1][1].re*e[1], M[1][1].im*e[1] };
         U[1][2] = (qcd_complex_16) { M[1][2].re*e[1], M[1][2].im*e[1] };

         U[2][0] = (qcd_complex_16) { M[2][0].re*e[2], M[2][0].im*e[2] };
         U[2][1] = (qcd_complex_16) { M[2][1].re*e[2], M[2][1].im*e[2] };
         U[2][2] = (qcd_complex_16) { M[2][2].re*e[2], M[2][2].im*e[2] };

	 _Complex double *_U = &U[0][0].re;
	 su3_linalg_upeqv(&u[i*NC*NC], _U);
       }
     else
       {
	 double e[NC];
	 qcd_complex_16 U[NC][NC];
         /*
           Make traceless to eliminate second order term in eigenvalue equation,
           i.e. eig^3 + a eig + b = 0, when H is traceless.
         */

         double trace = su3_linalg_trace_u(_H)/3.;
	 _Complex double dtr[NC*NC];
	 su3_linalg_diag(dtr, (_Complex double[]){trace, trace, trace});
         /* _H[CC(0, 0)] -= dtr[CC(0, 0)]; */
	 /* _H[CC(1, 1)] -= dtr[CC(1, 1)]; */
	 /* _H[CC(2, 2)] -= dtr[CC(2, 2)]; */
	 
	 if(isnan(dtr[0])) {
	     /* for(int j=0; j<NC; j++) */
	     /*   for(int k=0; k<NC; k++) */
	     /* 	 printf("%+e%+e\n", creal(dtr[j*NC+k]), cimag(dtr[j*NC+k])); */
	     for(int j=0; j<NC; j++)
	       for(int k=0; k<NC; k++)
		 printf("%+e%+e\n", creal(_M[j*NC+k]), cimag(_M[j*NC+k]));
	     exit(0);
	 }
	 /* su3_linalg_umeqv(_H, dtr); */

	 /* for(int j=0; j<NC; j++) */
	 /*   for(int k=0; k<NC; k++) */
	 /*     printf("%+e%+e\n", creal(_H[j*NC+k]), cimag(_H[j*NC+k])); */
       /*
           Solve for eigenvalues:
           e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
           e^3 + a*e + b = 0,

           a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
           b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

           D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

           e = D/(18^(1/3)) - ((2/3)^(1/3))/D
           e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
           e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
         */
         double a = -(H[2][2].re*H[2][2].re - H[0][0].re*H[1][1].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));
	 qcd_complex_16 b;
         b.re  = - H[0][0].re*H[1][1].re*H[2][2].re + H[2][2].re*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
         b.im  =                                      H[2][2].re*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

         b.re +=   H[0][0].re*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
         b.im +=   H[0][0].re*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

	 qcd_complex_16 w;
         w.re=qcd_CPOWR(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
         w.im=qcd_CPOWI(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

         qcd_complex_16 D=qcd_CPOW(((qcd_complex_16){-9.*b.re + w.re,-9.*b.im + w.im}),1./3.); 


	 /*
	   Constants Used:
	 */
	 qcd_complex_16 ThirdRootOne[2];
	 double ThirdRoot_12;
	 double ThirdRoot_18;
	 double ThirdRoot_2_3;
	 ThirdRootOne[0].re= 1;
	 ThirdRootOne[0].im= sqrt(3.);
	 
	 ThirdRootOne[1].re= 1;
	 ThirdRootOne[1].im=-sqrt(3.);
	 
	 ThirdRoot_12 =pow(12.,1./3.);
	 ThirdRoot_18 =pow(18.,1./3.);
	 ThirdRoot_2_3=pow((2./3.),1./3.);
	 
         e[0] = D.re/(ThirdRoot_18) - qcd_CDEVR(((qcd_complex_16){a*ThirdRoot_2_3,0}),D);
         e[1] = a*qcd_CDEVR(ThirdRootOne[0],((qcd_complex_16){D.re*ThirdRoot_12,D.im*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
         e[2] = -e[0]-e[1];

         e[0]+= trace;
         e[1]+= trace;
         e[2]+= trace;

         H[0][0].re+=trace;
         H[1][1].re+=trace;
         H[2][2].re+=trace;

	 qcd_complex_16 v[NC][NC];
	 qcd_complex_16 vr[NC][NC];			      
         /*
           Eigenvectors:
           v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
           v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
           v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
         */

         v[0][0].re = -(e[0]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[0][0].im = -(e[0]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[0][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[0][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[0][2].re =-e[0]*e[0] + e[0]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].re - H[0][0].re*H[1][1].re;
         v[0][2].im = 0.;

         v[1][0].re = -(e[1]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[1][0].im = -(e[1]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[1][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[1][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[1][2].re =-e[1]*e[1] + e[1]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].re - H[0][0].re*H[1][1].re;;
         v[1][2].im = 0.;

         /*
           Assure eigenvectors orthonormality:
           norm =  inner product v1.v1
           w    = (inner product v1.v2)/norm
           v2   = w*v1
         */
	 double norm;
         norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
         w.re  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
         w.im  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
         w.re /= norm;
         w.im /= norm;

         v[1][0].re-= qcd_CMULR(w,v[0][0]);
         v[1][0].im-= qcd_CMULI(w,v[0][0]);

         v[1][1].re-= qcd_CMULR(w,v[0][1]);
         v[1][1].im-= qcd_CMULI(w,v[0][1]);

         v[1][2].re-= qcd_CMULR(w,v[0][2]);
         v[1][2].im-= qcd_CMULI(w,v[0][2]);

         norm=1./sqrt(norm);

         /*
           Normalize first and second eigenvector:
         */

         v[0][0].re*= norm;
         v[0][0].im*= norm;

         v[0][1].re*= norm;
         v[0][1].im*= norm;

         v[0][2].re*= norm;
         v[0][2].im*= norm;


         norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

         norm=1./sqrt(norm);

         v[1][0].re*= norm;
         v[1][0].im*= norm;

         v[1][1].re*= norm;
         v[1][1].im*= norm;

         v[1][2].re*= norm;
         v[1][2].im*= norm;

         /*
           v3 = v1 x v2
         */


         v[2][0].re =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
         v[2][0].im = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

         v[2][1].re = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
         v[2][1].im = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

         v[2][2].re =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
         v[2][2].im = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

         double de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
         /*
	   cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
	   cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
	   cor[2] = - cor[0] - cor[1];
         */
         //to be compatible with Grenoble & Paris, don't apply corrections
	 double cor[NC];
         cor[0]=0;
         cor[1]=0;
         cor[2]=0;

         de = 1./sqrt(e[0]);
         b.re = de*cos(phase-cor[0]);
         b.im =-de*sin(phase-cor[0]);
         vr[0][0] = qcd_CMUL(b,v[0][0]);
         vr[0][1] = qcd_CMUL(b,v[0][1]);
         vr[0][2] = qcd_CMUL(b,v[0][2]);

         de = 1./sqrt(e[1]);
         b.re = de*cos(phase-cor[1]);
         b.im =-de*sin(phase-cor[1]);

         vr[1][0] = qcd_CMUL(b,v[1][0]);
         vr[1][1] = qcd_CMUL(b,v[1][1]);
         vr[1][2] = qcd_CMUL(b,v[1][2]);

         de = 1./sqrt(e[2]);
         b.re = de*cos(phase-cor[2]);
         b.im =-de*sin(phase-cor[2]);

         vr[2][0] = qcd_CMUL(b,v[2][0]);
         vr[2][1] = qcd_CMUL(b,v[2][1]);
         vr[2][2] = qcd_CMUL(b,v[2][2]);


         H[0][0].re= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].re= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].re= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

         H[0][0].im= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].im= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].im= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


         H[1][0].re= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].re= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].re= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

         H[1][0].im= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].im= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].im= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


         H[2][0].re= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].re= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].re= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

         H[2][0].im= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].im= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].im= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

         U[0][0].re= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
         U[0][1].re= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
         U[0][2].re= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

         U[0][0].im= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
         U[0][1].im= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
         U[0][2].im= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


         U[1][0].re= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
         U[1][1].re= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
         U[1][2].re= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

         U[1][0].im= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
         U[1][1].im= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
         U[1][2].im= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


         U[2][0].re= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
         U[2][1].re= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
         U[2][2].re= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

         U[2][0].im= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
         U[2][1].im= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
         U[2][2].im= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


         /*
           w    = inner product: col1.col2
           norm = inner product: col1.col1
         */

         norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
         w.re  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
         w.im  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
         w.re /= norm;
         w.im /= norm;


         U[0][1].re-=qcd_CMULR(w,U[0][0]);
         U[0][1].im-=qcd_CMULI(w,U[0][0]);

         U[1][1].re-=qcd_CMULR(w,U[1][0]);
         U[1][1].im-=qcd_CMULI(w,U[1][0]);

         U[2][1].re-=qcd_CMULR(w,U[2][0]);
         U[2][1].im-=qcd_CMULI(w,U[2][0]);

         norm = 1./sqrt(norm);

         U[0][0].re*= norm;
         U[0][0].im*= norm;
         U[1][0].re*= norm;
         U[1][0].im*= norm;
         U[2][0].re*= norm;
         U[2][0].im*= norm;

         norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
         norm = 1./sqrt(norm);

         U[0][1].re*= norm;
         U[0][1].im*= norm;
         U[1][1].re*= norm;
         U[1][1].im*= norm;
         U[2][1].re*= norm;
         U[2][1].im*= norm;

         /*
           col3 = col1 x col2
         */
         U[0][2].re =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
         U[0][2].im = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

         U[1][2].re = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
         U[1][2].im =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

         U[2][2].re =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
         U[2][2].im = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

         for(int c1=0; c1<3; c1++)
	   for(int c2=0; c2<3; c2++)
	     u[i*NC*NC + c1*NC + c2] = U[c1][c2].re + _Complex_I*U[c1][c2].im; 
       }
   }
   return;
}

