#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include "cubature.h"
/*
 *    code for rate computations
 *    @author  GJ
 *    @version 0.1
 *
 */


#define SQR(x) ((x)*(x))
#define SGN(x) (double) ((x>0)-(x<0))
#define THE(x) (double) ((x>0))
#define OOFP 0.079577471545947667884441881686
char E = 'U';

double complex lam(double complex x, double complex y, double complex z) {
  double complex res = SQR(x) + SQR(y) + SQR(z) - 2.*( x*y + x*z + y*z );
  return res;
}
double complex n(int s, double complex x) {
  double complex denom = cexp(x) - ((double complex) s);
  return ((double complex) s)/denom;
}

double complex calF(double complex z, double *ab) {
  double a = ab[0], b = ab[1];
  double complex zz = z ;//+ I*1e-1;
  //printf("|z-a| = %g, |b| = %g\n", cabs(z-a), cabs(b) );
  double complex denom = csqrt( (zz-a)*(zz-a) - b*b  );
  return SGN(creal(z-a))/denom;
}

double complex calG(double complex z, double ad) {
  double complex res = (.5/ad)*clog( cabs( (z+ad)/(z-ad) ) );
  //printf("z = %g + I %g, ad = %g\n", creal(z), cimag(z), ad);
  //printf("res = %g + I %g \n", creal(res), cimag(res));
  return res;
}

// integrating functions:
#include "real.h"
#include "virt.h"

//FILE *in;
FILE *out;
void rho_f_scan(double);   // k input
void virtual_check(double,double,double);   // omega,k,ma input
void Rate_fig2_scan(double); // M


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

/*--------------------------------------------------------------------*/




int main () {

  Rate_fig2_scan(3.);

  double aa[3] = {.1,0.,-1.};
  double bb[3] = {.01,0.,+1.};
  double cc[3] = {10.,0.,+1.};
  double dd[3] = {0.,0.,-1.};

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
    res = Rate_1_to_2(o,k,aa,bb,cc,dd);

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


void Rate_fig2_scan(double M) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_B = 0., mu_Y = 2.*1e-2;
  double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double l_[3] = {0.10, mu_l, -1.};
  double Q_[3] = {0.01, mu_Q, +1.};
  double S_[3] = {1.0, mu_S, +1.};

  int N_k;
  double o, k, k_min, k_max, step;
  double res_1_to_3, res_2_to_2, res_3_to_1;

  char *prefix=(char*)"out/rate_1to3_lQS";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{M=%.1f}.dat",M);
  strcat(filename,suffix);
  out=fopen(filename,"w");

  // Here are some parameters that can be changed:
  N_k=20; 
  k_min=0.1;
  k_max=15.;
  // don't change anything after that.

  k = k_min;
  step=(k_max-k_min)/((double) N_k-1);


  printf(" Settings: M=%g, with k_min=%g, k_max=%g\n",M,k_min,k_max); 
  double frac;

  for (int i=0; i<N_k; i++) {
    frac = (double)i/(double)(N_k-1);
    o = sqrt( SQR(M) + SQR(k) );

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

    res_1_to_3 *= 2.*(SQR(g1)+3.*SQR(g2));
    res_2_to_2 *= 2.*(SQR(g1)+3.*SQR(g2));
    res_3_to_1 *= 2.*(SQR(g1)+3.*SQR(g2));
    //printf("omega = %g , rho_f = %g\n", o, ans);

    printf(" k = %.5e , [%2.2f%]\n", k , 100.*frac); 
    fprintf( out, "%.8e    %.8e    %.8e    %.8e\n", k, res_1_to_3, res_2_to_2, res_3_to_1 );
    k += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

