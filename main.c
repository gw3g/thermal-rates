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
  double complex denom = csqrt( (z-a)*(z-a) - b*b  );
  return SGN(creal(z-a))/denom;
}

// integrating functions:
#include "real.h"
#include "virt.h"

//FILE *in;
FILE *out;
void rho_f_scan(double);   // k input

double rho_f(double o, double k, void *A_, void *B_, void *C_) {
  double K2 = SQR(o) - SQR(k);
  double res = 
    Rate_1_to_3(o,k,A_,B_,C_)
   +Rate_2_to_2_sChan(o,k,A_,B_,C_)
   +Rate_2_to_2_sChan(o,k,B_,A_,C_)
   +Rate_2_to_2_sChan(o,k,C_,A_,B_)
   +Rate_3_to_1_sChan(o,k,A_,B_,C_)
   +Rate_3_to_1_sChan(o,k,B_,C_,A_)
   +Rate_3_to_1_sChan(o,k,C_,A_,B_);
  return ( res - pow(OOFP,3.)*K2/8.)*64.*M_PI;
}

/*--------------------------------------------------------------------*/




int main () {
  double aa[3] = {0.0,0.,-1.};
  double bb[3] = {0.0,0.,-1.};
  double cc[3] = {0.0,0.,-1.};
  double dd[3] = {0.0,0.,-1.};
/*
  double ans1 = Rate_1_to_3(2.1,1.,aa,bb,cc);
  double ans2 = Rate_2_to_2(2.1,1.,aa,bb,cc);
  double ans3 = Rate_3_to_1(2.1,1.,aa,bb,cc);
  printf("1 -> 3 = %g\n", ans1);
  printf("2 -> 2 = %g\n", ans2);
  printf("3 -> 1 = %g\n", ans3);//*/
  //
  /*
  double o = .12, ans;
  while (o<80.) {
    ans = rho_f(o,.1,aa,bb,cc);
    printf("omega = %g , rho_f = %g\n", o, ans);
    o*=1.3;
  }//*/

  double ans1 = Rate_1_to_2(2.,.1,aa,bb,cc,dd);
    printf("omega = %g , rho_virt = %g\n", 2.1, ans1);

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

