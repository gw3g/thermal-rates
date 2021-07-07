#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include "cubature.h"
/*
 *    code for rate computations
 *    @author  GJ
 *    @version 0.4
 *
 */

#define SQR(x) ((x)*(x))
#define SGN(x) (double) ((x>0)-(x<0))
#define THE(x) (double) (x>0)

#define OOFP 0.079577471545947667884441881686l
#define ZE2  1.644934066848226436472415166646l

double long mubar = 6.283185307179586476925286766559; // (renorm. scale)/T

char E = 'U';

//FILE *in;
void check_ML(double,double,double,double,double);
void virtual_check(double,double,double);   // omega,k,ma input

// integrating functions:
#include "cal.h"
#include "Theta.h"
#include "real.h"
#include "virt.h"
#include "htl.h"

#include "scan.h"
#include "rho.h"      // master SPF's
#include "hadronic.h" // top contrib

/*--------------------------------------------------------------------*/

int main () {
  //rho_j_virt(3.,1.,.1,.01,.4);
  //test_C0();
  //double c_[4] = {.1,.3,.001,.003};
  //double complex test = calG(1.+I*1e-5, .2, 3, c_);
  //printf("test = %g + I %g \n\n", test);
  E='K';
  //hadron_scan(.2,.5); // to scan over m_q
  M_scan(8.,.1,.01,1.);
  //M_scan2(1.,.1,.01,0.4);
  //check_ML(3.,0.96875,.1,.01,1.);
  //check_ML(.3,0.4,.1,.01,1.);
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

  //printf("1 -> 3 = %g\n\n", ans1);
  //printf("2 -> 2 (s) = %g\n", ans2s);
  //printf("2 -> 2 (t) = %g\n\n", ans2t);
  //printf("3 -> 1 (s) = %g\n", ans3s);
  //printf("3 -> 1 (t) = %g\n", ans3t);//*/

  //virtual_check(2.,1.1,.1);
  //virtual_check(2.,1.1,.01);
  //virtual_check(2.,1.1,.001);
  //virtual_check(2.,1.1,.0);

  //rho_f_scan(.0001);
  //rho_f_scan(.1);
  //rho_f_scan(1);
  //rho_f_scan(10);

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


void check_ML(double M, double k, double m_l, double m_Q, double m_S) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_B = 0., mu_Y = 2.*1e-2;

  //double m_l = .1, m_Q = .01, m_S = 1.;
  //double mu_l = mu_L - .5*mu_Y, mu_Q = 0., mu_S = .5*mu_Y;
  double mu_l = 0., mu_Q = 0., mu_S = 0.;
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

