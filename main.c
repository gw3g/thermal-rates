#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "cubature.h"

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

int main () {
  double aa[3] = {0.0,0.,+1.};
  double bb[3] = {0.0,0.,+1.};
  double cc[3] = {0.0,0.,+1.};

  //double ans1 = Rate_1_to_3(2.1,1.,aa,bb,cc);
  //double ans2 = Rate_2_to_2(2.1,1.,aa,bb,cc);
  //double ans3 = Rate_3_to_1(2.1,1.,aa,bb,cc);
  //printf("1 -> 3 = %g\n", ans1);
  //printf("2 -> 2 = %g\n", ans2);
  //printf("3 -> 1 = %g\n", ans3);
  //
  double o = .12, ans;
  while (o<80.) {
    ans = rho_f(o,.1,aa,bb,cc);
    printf("omega = %g , rho_f = %g\n", o, ans);
    o*=1.3;
  }


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
