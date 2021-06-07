#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include "cubature.h"
/*
 *    code for rate computations
 *    @author  GJ
 *    @version 0.2
 *
 */

#define SQR(x) ((x)*(x))
#define SGN(x) (double) ((x>0)-(x<0))
#define OOFP 0.079577471545947667884441881686
#define ZE2  1.644934066848226436472415166646

double mubar = 6.283185307179586476925286766559;

char E = 'U';

void check_ML(double,double,double,double,double);
double complex lam(double complex x, double complex y, double complex z) {
  double complex res = SQR(x) + SQR(y) + SQR(z) - 2.*( x*y + x*z + y*z );
  return res;
}
double complex n(int s, double complex x) {
  //if (creal(x)<0.) { return -1.-n(s,-x); }
  double complex e = cexp(-x);
  //if (creal(x)>10.) { return e*((double) s); }
  double complex denom ;
  if (s==1) { denom = - expm1(-creal(x)); }
  else { denom = 1.-((double) s)*e; }
  return e*((double) s)/denom;
}

double complex calN(int sA, int sB,
                    double complex xA, double complex xB) {

  double pA = SGN(creal(xA));
  double pB = SGN(creal(xB));

  if ((pA>0)&&(pB>0)) return + 1.+ n(sA,+xA) + n(sB,+xB);
  if ((pA<0)&&(pB>0)) return     - n(sA,-xA) + n(sB,+xB);
  if ((pA>0)&&(pB<0)) return     + n(sA,+xA) - n(sB,-xB);
  if ((pA<0)&&(pB<0)) return - 1.- n(sA,-xA) - n(sB,-xB);
}

double complex calF(double complex z, double *ab) {
  double a = ab[0], b = ab[1];
  double complex zz = z;// + I*1e-8;
  //printf("|z-a| = %g, |b| = %g\n", cabs(z-a), cabs(b) );
  double complex denom = csqrt( (zz-a)*(zz-a) - b*b  );
  return SGN(creal(z-a))/denom;
}

double calA(double al, double be, double ga) { // Eq.(C.16)

  double b2 = SQR(be),
         ag = 4.*al*ga,
         bg = be+2.*ga,
         D  = sqrt(fabs(b2-ag));

       if (b2>ag) { return log(fabs( (bg+D)/(bg-D) ))/D; }
  else if (b2<ag) { return 2.*atan(D/bg)/D; }
  else            { return 2./bg; }
}

double complex calG(double complex z, double ad, int nn, double *coeffs) {
  /*
   *  reproduce Eq. (C.14) :
   *
   *  ad = p_a*p_d
   *  nn = n
   *  coeffs = {c_0,c_1,...,c_n}
   *
   */
  if (nn<0) return 0.;
  double complex res = coeffs[nn], d[nn+1]; d[nn]=0.;
  for (int i=nn;i>0;i--) {
    res *= z; res += coeffs[i-1];   // construct R_n(z)
    d[i-1] =  - coeffs[i] + z*d[i]; // populate list: {d_0,d_1,...}
  }
  res *= .5*clog(cabs( (z+ad)/(z-ad) ))/ad; // = R_n(z)*log(...)/(2*p_a*p_d)
  for (int i=0;i<nn;i++) {
   if (i%2==0) { res += d[i]*pow(ad,i)/((double)i+1.); } // term ~ (p_a*p_d)^l
  }
  //printf("res = %g + I %g \n", creal(res), cimag(res));
  return res;
}

double complex calH(double complex z1, double complex z2,
                    double pa, double pd, double pe, double k, 
                    int nn, double *coeffs) {
  /*
   *  reproduce Eq. (C.4) :
   *
   *  nn = n
   *  coeffs = {a_0,a_1,...,a_n}
   *
   */
  double complex res = coeffs[nn];
  double b[nn+1]; b[nn]=0.;
  for (int i=nn;i>0;i--) {
    res *= z2; res += coeffs[i-1];           // construct Q_n(z2)
    b[i-1] =  creal( - coeffs[i] + z2*b[i]); // populate list: {b_0,b_1,...}
  }

  double k_pe = .5*( SQR(k)+SQR(pe)-SQR(pd) );
  double complex z12 = z1-z2;

  double al =    creal( SQR(z12) - SQR(pa*pe)   ),
         be = 2.*creal( z2*z12   + SQR(pa)*k_pe ),
         ga =    creal( SQR(z2)  - SQR(pa*k)    );

  //printf("calA = %.2f\n", calA(al,be,ga) );
  res = calA(al,be,ga)*res + calG(z1,pa*pd,nn,b); // or nn-1?

  return res;
}

//FILE *in;
FILE *out;
void rho_f_scan(double);   // k input
void virtual_check(double,double,double);   // omega,k,ma input
void Rate_fig2_scan(double,double,double,double); // M, ml, mQ, mS
void M_scan(double,double,double,double); // k, ml, mQ, mS

// integrating functions:
#include "Real.h"
#include "Virt.h"
#include "htl.h"



double rho_f(double o, double k, void *A_, void *B_, void *C_) {
  double K2 = SQR(o) - SQR(k);
  double res = 
    Rate_1_to_3(o,k,A_,B_,C_,_id)
   +Rate_2_to_2_sChan(o,k,A_,B_,C_,_id)
   +Rate_2_to_2_sChan(o,k,B_,A_,C_,_id)
   +Rate_2_to_2_sChan(o,k,C_,A_,B_,_id)
   +Rate_3_to_1_sChan(o,k,A_,B_,C_,_id)
   +Rate_3_to_1_sChan(o,k,B_,C_,A_,_id)
   +Rate_3_to_1_sChan(o,k,C_,A_,B_,_id);
  return ( res - pow(OOFP,3.)*K2/8.)*64.*M_PI;
}//*/

void rho_j_scan1();
void rho_j_scan2();
double rho_j_real(double M, double k, double m_l, double m_Q, double m_S) {

  double o = sqrt( SQR(M) + SQR(k) );

  //double m_l = .1, m_Q = .01, m_S = 1.;
  //double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};

  double res = 
    Rate_1_to_3(o,k,l_,Q_,S_,_j_1_3_s)
   +Rate_2_to_2_tChan(o,k,Q_,l_,S_,_j_2_2_t1)
   +Rate_2_to_2_tChan(o,k,l_,Q_,S_,_j_2_2_t2)
   +Rate_2_to_2_sChan(o,k,S_,Q_,l_,_j_2_2_s)
   +Rate_3_to_1_tChan(o,k,l_,S_,Q_,_j_3_1_t1)
   +Rate_3_to_1_tChan(o,k,Q_,S_,l_,_j_3_1_t2)
   +Rate_3_to_1_sChan(o,k,Q_,l_,S_,_j_3_1_s);
    ;
  return  res*32.*M_PI ;
}//*/
double rho_j_virt(double M, double k, double m_l, double m_Q, double m_S) {

  double o = sqrt( SQR(M) + SQR(k) );

  //double m_l = .1, m_Q = .01, m_S = 1.;
  //double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};
  double ccc[3];
  ccc[0] = 1.; ccc[1] = 0.; ccc[2] = 1.;
  E = 'K';

  // CHECK (1/2+n) -> 1/2
  //double temp1 = Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/SQR(M); //*/
  //double temp2 = Triangle_0_new(o,k,l_,Q_,S_,l_,S_,_j_virt); //*/
  //printf(" temp1 = %.8f \n", temp1);
  //printf(" temp2 = %.8f \n", temp2);
  //return 0.;

  double res = 
     Triangle_T(o,k,l_,Q_,S_,l_,S_,_j_virt) 
   + Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/SQR(M) //*/
    ;
  return  res*32.*M_PI ;
}//*/



/*--------------------------------------------------------------------*/




void test_C0();
int main () {
  //rho_j_virt(3.,1.,.1,.01,.4);
  //test_C0();
  //double c_[4] = {.1,.3,.001,.003};
  //double complex test = calG(1.+I*1e-5, .2, 3, c_);
  //printf("test = %g + I %g \n\n", test);
  E='U';
  //M_scan(2.,.1,.01,.4);
  //check_ML(3.,0.96875,.1,.01,1.);
  check_ML(.3,0.4,.1,.01,1.);
  //rho_j_scan1();

  //double m_l = .1, m_Q = .01, m_S = 1.;
  E='U';
  //Rate_fig2_scan(3.,.1,.01,1.); // l,Q,S
  //Rate_fig2_scan(3.,.1,.01,10.);
  //Rate_fig2_scan(.3,.1,.01,1.);
  E='K';
  //Rate_fig2_scan(3.,.1,.01,1.);
  //Rate_fig2_scan(3.,.1,.01,10.);
  //Rate_fig2_scan(.3,.1,.01,1.);//*/

  /*double aa[3] = {.1,0.,-1.};
  double bb[3] = {.01,0.,+1.};
  double cc[3] = {10.,0.,+1.};
  double dd[3] = {0.,0.,-1.};//*/

  //double ph = Phi(2.,1.1,aa,bb,0);
  //printf("1 -> 2 = %g\n\n", ph);

  //double ans1 = Rate_1_to_3(2.1,1.,aa,bb,cc);

  //double ans2s = Rate_2_to_2_sChan(3.,.01,aa,bb,cc,_id);
  //double ans2t = Rate_2_to_2_tChan(3.,.01,aa,bb,cc,_id);

  //double ans3s = Rate_3_to_1_sChan(2.1,1.1,aa,bb,cc,_id);
  //double ans3t = Rate_3_to_1_tChan(2.1,1.1,aa,bb,cc,_id);

  //printf("1 -> 3 = %g\n\n", ans1);
  //printf("2 -> 2 (s) = %g\n", ans2s);
  //printf("2 -> 2 (t) = %g\n\n", ans2t);
  //printf("3 -> 1 (s) = %g\n", ans3s);
  //printf("3 -> 1 (t) = %g\n", ans3t);//*/
  //
  /*
  double o = .12, ans;
  while (o<80.) {
    ans = rho_f(o,.1,aa,bb,cc);
    printf("omega = %g , rho_f = %g\n", o, ans);
    o*=1.3;
  }//*/

  //virtual_check(2.,1.1,.1);
  //virtual_check(2.,1.1,.01);
  //virtual_check(2.,1.1,.001);
  //virtual_check(2.,1.1,.0);

  //double ans1 = Rate_1_to_2(2.,1.1,aa,bb,cc,dd);
    //printf("omega = %g , rho_virt = %g\n", 2.1, ans1);

  //rho_f_scan(.0001);
  //rho_f_scan(.1);
  //rho_f_scan(1);
  //rho_f_scan(10);

  //
  /*double z = 2.;
  double ab[2] = {1.,-.2};
  double complex FF;
  while (ab[1]<1.9) {
    FF = calF(z,ab);
    printf("z=%.2f, a=%.2f, b=%.2f, result = %.2f + I %.2f \n", z, ab[0], ab[1], creal(FF), cimag(FF));
    ab[1]+=.1;
  }*/

  return 0;
}

/*--------------------------------------------------------------------*/

void virtual_check(double o, double k, double ma) {
  double aa[3] = {ma,0.,+1.};
  double bb[3] = {1.0,0.,+1.};
  double cc[3] = {0.0,0.,+1.};
  double dd[3] = {0.0,0.,+1.};

  int N_mb;
  double mb, mb_min, mb_max, step;
  double res;

  char *prefix=(char*)"out/virtual_test";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{ma=%g,o=%.1f,k=%.1f}.dat",ma,o,k);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  //fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_mb=200; 
  mb_min=1e-6;
  mb_max=1e-1;
  // don't change anything after that.

  mb = mb_min;
  step=pow(mb_max/mb_min,1./(double)(N_mb-1));


  printf(" Settings: o=%g, k=%g, with mb_min=%g, mb_max=%g\n",o,k,mb_min,mb_max); 
  double frac;

  for (int i=0; i<N_mb; i++) {
    frac = (double)i/(double)(N_mb-1);

    bb[0] = mb;
    res = Bubble_T(o,k,aa,bb,cc,dd);

    //printf(" mb = %.5e , [%2.2f%]\n", mb , 100.*frac); 
    fprintf( out, "%.8e    %.8e\n", mb, res );
    mb *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}


void rho_f_scan(double k) {
  double aa[3] = {0.0,0.,+1.};
  double bb[3] = {0.0,0.,+1.};
  double cc[3] = {0.0,0.,+1.};

  int N_o;
  double o, o_min, o_max, step;
  double res;

  char *prefix=(char*)"out/rho_f_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{k=%.1f}.dat",k);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  //fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_o=200; 
  o_min=fmax(k*1.01,.01);
  o_max=100.;
  // don't change anything after that.

  o = o_min;
  step=pow(o_max/o_min,1./(double)(N_o-1));


  printf(" Settings: k=%g, with o_min=%g, o_max=%g\n",k,o_min,o_max); 
  double frac;

  for (int i=0; i<N_o; i++) {
    frac = (double)i/(double)(N_o-1);

    res = rho_f(o,k,aa,bb,cc);
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" omega = %.5e , [%2.2f%]\n", o , 100.*frac); 
    fprintf( out, "%.8e    %.8e\n", o, res );
    o *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

#include <gsl/gsl_sf_bessel.h>
double k2av(double M) {
  // worth noting: k0=3*k for k0~26.
  return 3.*M*gsl_sf_bessel_Kn(3,M)/gsl_sf_bessel_Kn(2,M);
}

void rho_j_scan1() {
  int N_M;
  double M, M_min, M_max, step, o, k;
  double res_r, res_v;

  char *prefix=(char*)"out/rho_j_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"1.dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  //fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_M=20; 
  M_min=.1;
  M_max=100.;
  // don't change anything after that.

  M = M_min;
  step=pow(M_max/M_min,1./(double)(N_M-1));


  printf(" Settings: M_min=%g, M_max=%g\n",M_min,M_max); 
  double frac;

  for (int i=0; i<N_M; i++) {
    frac = (double)i/(double)(N_M-1);
    k = sqrt(k2av(M));
    //o = sqrt( M*M + k*k );

    res_r = rho_j_real(M,k,.001,.001,.01);
    res_v = rho_j_virt(M,k,.001,.001,.01);
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" M = %.5e , [%2.2f%]\n", M , 100.*frac); 
    fprintf( out, "%.8e    %.8e    %.8e\n", M, res_r, res_v );
    M *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

void rho_j_scan2() {
  int N_M;
  double mQ, M_min, M_max, step, o, k;
  double res_r, res_v;

  char *prefix=(char*)"out/rho_j_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"4.dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  //fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_M=30; 
  M_min=.0001;
  M_max=.01;
  // don't change anything after that.

  mQ = M_min;
  step=pow(M_max/M_min,1./(double)(N_M-1));
  double M=3.;
  k = sqrt(k2av(M));
  o = sqrt( SQR(M)+SQR(k) );

  printf(" Settings: M_min=%g, M_max=%g\n",M_min,M_max); 
  double frac;

  for (int i=0; i<N_M; i++) {
    frac = (double)i/(double)(N_M-1);

    res_r = rho_j_real(M,k,.001,mQ,.01); // l,Q,S
    res_v = rho_j_virt(M,k,.001,mQ,.01);
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" M = %.5e , [%2.2f%]\n", M , 100.*frac); 
    fprintf( out, "%.8e    %.8e    %.8e\n", mQ, res_r, res_v );
    mQ *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}



void Rate_fig2_scan(double M, double m_l, double m_Q, double m_S) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_Y = 2.*1e-2;
  double G12 = 2.*(SQR(g1)+3.*SQR(g2));

  //double m_l = .1, m_Q = .01, m_S = 1.;
  double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};

  double dm_S = 1e-7;
  double St_plus[3]  = {m_S+dm_S, mu_S, +1.};
  double St_minus[3] = {m_S-dm_S, mu_S, +1.};

  int N_k;
  double o, k, k_min, k_max, step;
  double res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1, tot;
  double res_virt;

  char *prefix=(char*)"out/rates_";
  char  suffix[50];
  char  filename[90];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"E=%c_{ml=%.2f,mQ=%.3f,mS=%.2f,M=%.2f}.dat",E,m_l,m_Q,m_S,M);
  strcat(filename,suffix);
  out=fopen(filename,"w");

  // Here are some parameters that can be changed:
  N_k=150; 
  k_min=0.001;
  k_max=15.;
  // don't change anything after that.

  k = k_min;
  step=(k_max-k_min)/((double) N_k-1);


  printf(" Settings: M=%g, with k_min=%g, k_max=%g\n",M,k_min,k_max); 
  double frac;
  double ccc[3];

  for (int i=0; i<N_k; i++) {
    frac = (double)i/(double)(N_k-1);
    o = sqrt( SQR(M) + SQR(k) );

    res_1_to_2  = 0.;
    res_1_to_2 += Phi(o,k,l_,S_,_448_);

    res_1_to_3  = 0.;
    res_1_to_3 += Rate_1_to_3(o,k,l_,Q_,S_,_B1_i );
    res_1_to_3 += Rate_1_to_3(o,k,S_,Q_,l_,_B1_ii);

    res_2_to_2  = 0.;
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,l_,S_,_45_i  );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,l_,Q_,S_,_45_ii );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,S_,l_,_45_iii);
    res_2_to_2 += Rate_2_to_2_tChan(o,k,S_,Q_,l_,_45_iii);

    res_2_to_2 += Rate_2_to_2_sChan(o,k,S_,Q_,l_,_B9_i);
    res_2_to_2 += Rate_2_to_2_sChan(o,k,l_,Q_,S_,_B9_ii);//*/

    res_3_to_1  = 0.;
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,l_,S_,_B27_i );
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,S_,l_,_B27_ii);

    res_3_to_1 += Rate_3_to_1_tChan(o,k,l_,S_,Q_,_B19_i );
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,S_,l_,_B19_ii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,S_,l_,Q_,_B19_iii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,l_,S_,_B19_iii);//*/

    res_virt  = 0.;
    res_virt += Bubble_T(o,k,l_,Q_,S_,l_,_virt_i);
    res_virt += Bubble_T(o,k,Q_,S_,l_,S_,_virt_ii);

    res_virt += (m_S)*( Bubble_T(o,k,Q_,S_,l_,St_plus ,_virt_iii)
                      - Bubble_T(o,k,Q_,S_,l_,St_minus,_virt_iii) )/(2.*dm_S);

    res_virt += (SQR(m_S)-SQR(M))*Triangle_T(o,k,l_,Q_,S_,l_,S_,_virt_iv);

    ccc[0] = 1.; ccc[1] = 0.; ccc[2] = 0.;
    res_virt += Bubble_0(o,k,l_,Q_,S_,l_,ccc);
    ccc[0] = 2.; ccc[1] = 1.; ccc[2] = 0.;
    res_virt += (SQR(m_S)-SQR(M))*Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 2.;
    res_virt += Bubble_0(o,k,Q_,S_,l_,S_,ccc);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 1.;
    res_virt += (m_S)*( Bubble_0(o,k,Q_,S_,l_,St_plus,ccc)
                      - Bubble_0(o,k,Q_,S_,l_,St_minus,ccc) )/(2.*dm_S);

    tot = res_1_to_3+res_2_to_2+res_3_to_1+res_virt;
    res_1_to_3 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_2_to_2 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_3_to_1 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_virt   *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    tot *= G12;
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" k = %.5e , [%2.2f%]\n", k , 100.*frac); 
    fprintf( out, "%.10e    %.10e    %.10e    %.10e    %.10e    %.10e    %.10e\n", 
                   k, res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1, res_virt, tot );
    k += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

void check_ML(double M, double k, double m_l, double m_Q, double m_S) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_B = 0., mu_Y = 2.*1e-2;

  //double m_l = .1, m_Q = .01, m_S = 1.;
  double mu_l = mu_L - .5*mu_Y, mu_Q = 0., mu_S = .5*mu_Y;
  //double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};

  double dm_S = 1e-7;
  double St_plus[3]  = {m_S+dm_S, mu_S, +1.};
  double St_minus[3] = {m_S-dm_S, mu_S, +1.};

  double res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1;
  double res_virt, res_virt_0;

  double G12 = 2.*(SQR(g1)+3.*SQR(g2));
  double o = sqrt( SQR(M) + SQR(k) );
  double res_i, res_ii, res_iii, res_iv;

  void do_real() {
    res_1_to_2  = 0.;
    res_1_to_2 += Phi(o,k,l_,S_,_448_);

    res_1_to_3  = 0.;
    res_1_to_3 += Rate_1_to_3(o,k,l_,Q_,S_,_B1_i );
    res_1_to_3 += Rate_1_to_3(o,k,S_,Q_,l_,_B1_ii);

    res_2_to_2  = 0.;
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,l_,S_,_45_i  );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,l_,Q_,S_,_45_ii );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,S_,l_,_45_iii);
    res_2_to_2 += Rate_2_to_2_tChan(o,k,S_,Q_,l_,_45_iii);

    res_2_to_2 += Rate_2_to_2_sChan(o,k,S_,Q_,l_,_B9_i);
    res_2_to_2 += Rate_2_to_2_sChan(o,k,l_,Q_,S_,_B9_ii);//*/

    res_3_to_1  = 0.;
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,l_,S_,_B27_i );
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,S_,l_,_B27_ii);

    res_3_to_1 += Rate_3_to_1_tChan(o,k,l_,S_,Q_,_B19_i );
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,S_,l_,_B19_ii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,S_,l_,Q_,_B19_iii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,l_,S_,_B19_iii);//*/
  }
    //res_virt  = 0.;
    //res_virt += Bubble(o,k,l_,Q_,l_,S_,_virt_i);
    //res_virt += Bubble(o,k,Q_,S_,l_,S_,_virt_ii);
    //
    //res_1_to_3 = Rate_1_to_3(o,k,S_,Q_,l_,_B1_ii);

    //res_virt = (m_S)*( Bubble(o,k,Q_,S_,l_,St_plus ,_virt_iii)
    //              - Bubble(o,k,Q_,S_,l_,St_minus,_virt_iii) )/(.5*dm_S);

    //res_virt -= (SQR(m_S)-SQR(M))*Triangle(o,k,l_,Q_,S_,l_,S_,_virt_iv);
    double ccc[3];

    printf("\n(i) single pole at m_l\n\n");
    _i = 1; _ii = 0; _iii = 0; _iv = 0;
    do_real();
    //res_1_to_3 = Rate_1_to_3(o,k,l_,Q_,S_,_B1_i);
    printf(" 3to1 = %.6f \n", G12*res_3_to_1);
    printf(" 1to3 = %.6f \n", G12*res_1_to_3);
    printf(" 2to2 = %.6f \n", G12*res_2_to_2);
    res_virt = Bubble_T(o,k,l_,Q_,S_,l_,_virt_i);
    printf(" virtual_T = %.6f \n", G12*res_virt);
    ccc[0] = 1.; ccc[1] = 0.; ccc[2] = 0.;
    res_virt_0 = Bubble_0(o,k,l_,Q_,S_,l_,ccc);
    printf(" virtual_0 = %.6f \n", G12*res_virt_0);
    printf(" sum of virtuals = %.6f \n", G12*(res_virt+res_virt_0));
    res_i = G12*(res_1_to_3+res_2_to_2+res_3_to_1+res_virt+res_virt_0);
    printf(" sum of all = %.6f \n", res_i );

    printf("\n(ii) poles at m_l & m_S\n\n");
    _i = 0; _ii = 1; _iii = 0; _iv = 0;
    do_real();
    //res_1_to_3 = Rate_1_to_3(o,k,l_,Q_,S_,_B1_i);
    printf(" 3to1 = %.6f \n", G12*res_3_to_1);
    printf(" 1to3 = %.6f \n", G12*res_1_to_3);
    printf(" 2to2 = %.6f \n", G12*res_2_to_2);
    res_virt = (SQR(m_S)-SQR(M))*Triangle_T(o,k,l_,Q_,S_,l_,S_,_virt_iv);
    printf(" virtual_T = %.6f \n", G12*res_virt);
    ccc[0] = 2.; ccc[1] = 1.; ccc[2] = 0.;
    res_virt_0 = (SQR(m_S)-SQR(M))*Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc);
    printf(" virtual_0 = %.6f \n", G12*res_virt_0);
    printf(" sum of virtuals = %.6f \n", G12*(res_virt+res_virt_0));
    res_ii = G12*(res_1_to_3+res_2_to_2+res_3_to_1+res_virt+res_virt_0);
    printf(" sum of all = %.6f \n", res_ii );

    printf("\n(iii) single pole at m_S\n\n");
    _i = 0; _ii = 0; _iii = 1; _iv = 0;
    do_real();
    printf(" 3to1 = %.6f \n", G12*res_3_to_1);
    printf(" 1to3 = %.6f \n", G12*res_1_to_3);
    printf(" 2to2 = %.6f \n", G12*res_2_to_2);
    res_virt = Bubble_T(o,k,Q_,S_,l_,S_,_virt_ii);
    printf(" virtual_T = %.6f \n", G12*res_virt);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 2.;
    res_virt_0 = Bubble_0(o,k,Q_,S_,l_,S_,ccc);
    printf(" virtual_0 = %.6f \n", G12*res_virt_0);
    printf(" sum of virtuals = %.6f \n", G12*(res_virt+res_virt_0));
    res_iii = G12*(res_1_to_3+res_2_to_2+res_3_to_1+res_virt+res_virt_0);
    printf(" sum of all = %.6f \n", res_iii );


    printf("\n(iv) double pole at m_S\n\n");
    _i = 0; _ii = 0; _iii = 0; _iv = 1;
    do_real();
    //res_1_to_3 = Rate_1_to_3(o,k,S_,Q_,l_,_B1_ii);
    printf(" 3to1 = %.6f \n", G12*res_3_to_1);
    printf(" 1to3 = %.6f \n", G12*res_1_to_3);
    printf(" 2to2 = %.6f \n", G12*res_2_to_2);
    res_virt = (m_S)*( Bubble_T(o,k,Q_,S_,l_,St_plus ,_virt_iii)
                     - Bubble_T(o,k,Q_,S_,l_,St_minus,_virt_iii) )/(2.*dm_S);
    printf(" virtual_T = %.6f \n", G12*res_virt);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 1.;
    res_virt_0 = (m_S)*( Bubble_0(o,k,Q_,S_,l_,St_plus,ccc)
                       - Bubble_0(o,k,Q_,S_,l_,St_minus,ccc) )/(2.*dm_S);
    printf(" virtual_0 = %.6f \n", G12*res_virt_0);
    printf(" sum of virtuals = %.6f \n", G12*(res_virt+res_virt_0));
    res_iv = G12*(res_1_to_3+res_2_to_2+res_3_to_1+res_virt+res_virt_0);
    printf(" sum of all = %.6f \n", res_iv );

    printf("\n Sum of (i)-(iv): %.6f \n", res_i+res_ii+res_iii+res_iv );

}

void test_C0() {
  int nn = 100;
  E = 'U';
  double M, ma, mb, mc, md, me; // ml, mQ, mS, ml, mS
  double M_i[6] = {.1e+3, .1e-2, .1e-2, .1e-1, .1e-2, .1e-1};
  double M_f[6] = {3.,    .1,    .01, .1e+1, .1,    .1e+1};
  double M_step[6];

  for (int i=0;i<6;i++) { M_step[i] = (M_f[i]-M_i[i])/((double)nn-1); }
  M  = M_i[0];
  ma = M_i[1];
  mb = M_i[2];
  mc = M_i[3];
  md = M_i[4];
  me = M_i[5];
  double ccc[3] = {1.,0.,1.};
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_B = 0., mu_Y = 2.*1e-2;
  double G12 = 2.*(SQR(g1)+3.*SQR(g2));

  //double m_l = .1, m_Q = .01, m_S = 1.;
  double mu_l = mu_L - .5*mu_Y, mu_Q = 0., mu_S = .5*mu_Y;
  double l_[3] = {0., mu_l, -1.};
  double Q_[3] = {0., mu_Q, +1.};
  double S_[3] = {0., mu_S, +1.};
  double o, k =1.;
  double C0_new, C0_old;


  for (int i=0;i<nn;i++) {
    //l_[0] = ma;
    //Q_[0] = mb;
    //S_[0] = mc;
    //o = sqrt(M*M+k*k);
    //C0 = G12*(SQR(mc)-SQR(M))*Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/o;
    //C0 = Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/o;
    C0_old = L3(SQR(M),SQR(ma),SQR(mb),SQR(mc),SQR(md),SQR(me))*SQR(OOFP);
    C0_new = L3_new(SQR(M),SQR(ma),SQR(mb),SQR(mc),SQR(md),SQR(me))*SQR(OOFP);

    printf("%.8e   %.8e   %.8e   %.8e   %.8e   %.8e   %.8e   %.8e\n",
            md,me,M,ma,mb,mc,C0_new,C0_old);
    M  += M_step[0];
    ma += M_step[1];
    mb += M_step[2];
    mc += M_step[3];
    md += M_step[4];
    me += M_step[5];

  }
}

void M_scan(double k, double m_l, double m_Q, double m_S) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_Y = 2.*1e-2;
  double G12 = 2.*(SQR(g1)+3.*SQR(g2));

  //double m_l = .1, m_Q = .01, m_S = 1.;
  //double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};
  double mlT = (1.+SQR(mu_l/M_PI))*G12/(32.); mlT = sqrt(mlT);

  double dm_S = 1e-7;
  double St_plus[3]  = {m_S+dm_S, mu_S, +1.};
  double St_minus[3] = {m_S-dm_S, mu_S, +1.};

  int N_M;
  double o, M, M_min, M_max, step;
  double res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1, tot;
  double res_virt;

  char *prefix=(char*)"out/M_scan_";
  char  suffix[50];
  char  filename[90];

  //printf( "2to2, HTL = %.8f\n", Rate_2_to_2_HTL(2.01,2.0,l_,S_,mlT) );

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"E=%c_{ml=%.2f,mQ=%.3f,mS=%.2f,k=%.2f}.dat",E,m_l,m_Q,m_S,k);
  strcat(filename,suffix);
  out=fopen(filename,"w");

  // Here are some parameters that can be changed:
  N_M=300; 
  M_min=0.01;
  M_max=10.;
  // don't change anything after that.

  M = M_min;
  step=pow(M_max/M_min,1./(N_M-1));


  printf(" Settings: k=%g, with M_min=%g, M_max=%g\n",k,M_min,M_max); 
  double frac;
  double ccc[3];

  for (int i=0; i<N_M; i++) {
    frac = (double)i/(double)(N_M-1);
    o = sqrt( SQR(M) + SQR(k) );

    //res_1_to_2  = 0.;
    //res_1_to_2 += Phi(o,k,l_,S_,_448_);

    res_1_to_3  = 0.;
    res_1_to_3 += Rate_1_to_3(o,k,l_,Q_,S_,_B1_i );
    res_1_to_3 += Rate_1_to_3(o,k,S_,Q_,l_,_B1_ii);//*/

    res_2_to_2  = 0.;
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,l_,S_,_45_i  );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,l_,Q_,S_,_45_ii );
    res_2_to_2 += Rate_2_to_2_tChan(o,k,Q_,S_,l_,_45_iii);
    res_2_to_2 += Rate_2_to_2_tChan(o,k,S_,Q_,l_,_45_iii);

    res_2_to_2 += Rate_2_to_2_sChan(o,k,S_,Q_,l_,_B9_i);
    res_2_to_2 += Rate_2_to_2_sChan(o,k,l_,Q_,S_,_B9_ii);//*/

    res_3_to_1  = 0.;
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,l_,S_,_B27_i );
    res_3_to_1 += Rate_3_to_1_sChan(o,k,Q_,S_,l_,_B27_ii);

    res_3_to_1 += Rate_3_to_1_tChan(o,k,l_,S_,Q_,_B19_i );
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,S_,l_,_B19_ii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,S_,l_,Q_,_B19_iii);
    res_3_to_1 += Rate_3_to_1_tChan(o,k,Q_,l_,S_,_B19_iii);//*/

    res_virt  = 0.;
    res_virt += Bubble_T(o,k,l_,Q_,S_,l_,_virt_i);
    res_virt += Bubble_T(o,k,Q_,S_,l_,S_,_virt_ii);

    res_virt += (m_S)*( Bubble_T(o,k,Q_,S_,l_,St_plus ,_virt_iii)
                      - Bubble_T(o,k,Q_,S_,l_,St_minus,_virt_iii) )/(2.*dm_S);

    res_virt += (SQR(m_S)-SQR(M))*Triangle_T(o,k,l_,Q_,S_,l_,S_,_virt_iv);

    ccc[0] = 1.; ccc[1] = 0.; ccc[2] = 0.;
    res_virt += Bubble_0(o,k,l_,Q_,S_,l_,ccc);
    ccc[0] = 2.; ccc[1] = 1.; ccc[2] = 0.;
    res_virt += (SQR(m_S)-SQR(M))*Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 2.;
    res_virt += Bubble_0(o,k,Q_,S_,l_,S_,ccc);
    ccc[0] = -1.; ccc[1] = -1.; ccc[2] = 1.;
    res_virt += (m_S)*( Bubble_0(o,k,Q_,S_,l_,St_plus,ccc)
                      - Bubble_0(o,k,Q_,S_,l_,St_minus,ccc) )/(2.*dm_S);//*/

    //res_2_to_2 += Rate_2_to_2_HTL(o,k,l_,S_,mlT);
    res_1_to_2  = Rate_1_to_2_HTL(o,k,l_,S_,mlT);
    res_1_to_2 += Rate_2_to_2_HTL(o,k,l_,S_,mlT);

    //tot = res_1_to_3+res_2_to_2+res_3_to_1+res_virt + Rate_1_to_2_HTL(o,k,l_,S_,mlT);
    res_1_to_3 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_2_to_2 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_3_to_1 *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    res_virt   *= G12;//2.*(SQR(g1)+3.*SQR(g2));
    tot = res_1_to_3 + res_2_to_2 + res_3_to_1 + res_1_to_2;
    //tot *= G12;
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" o = %.5e , [%2.2f%]\n", o , 100.*frac); 
    //fprintf( out, "%.10e    %.10e    %.10e    %.10e    %.10e    %.10e    %.10e\n", 
    //               M, res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1, res_virt, tot );
    fprintf( out, "%.10e    %.10e    %.10e    %.10e    %.10e    %.10e\n", 
                   M, res_1_to_2, res_1_to_3, res_2_to_2, res_3_to_1, tot );
    M *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}
