/*
 *  definitions of auxiliary functions:
 *
 *    lambda(x,y,z)    -- Kallen 
 *
 *    n_sigma(e)       -- distribution fncs.
 *    N_{sA,sB}(eA,eB) -- combination 1+n+n
 *
 *    F_{q;k,p1}(z)
 *    A_{alpha,beta,gamma}
 *    G_{pA;pD,k}(z)
 *    H_{pA;pD,k}(z1,z2)
 */

/*--------------------------------------------------------------------*/

double complex lam(double complex x, double complex y, double complex z) {
  double complex res = SQR(x) + SQR(y) + SQR(z) - 2.*( x*y + x*z + y*z );
  return res;
}

/*--------------------------------------------------------------------*/

double complex n(int s, double complex x) {
  double complex e = cexp(-x);
  double complex denom ;

  if (s==1) { denom = - expm1(-creal(x));  }
  else      { denom = 1. - ((double) s)*e; }

  return e*((double) s)/denom;
}

double complex calN(int sA, int sB,
                    double complex xA, double complex xB) {

  // only evaluate n(x) for x>0, to avoid
  // significant loss in the subtraction

  double pA = SGN(creal(xA));
  double pB = SGN(creal(xB));

  if ((pA>0)&&(pB>0)) return + 1.+ n(sA,+xA) + n(sB,+xB);
  if ((pA<0)&&(pB>0)) return     - n(sA,-xA) + n(sB,+xB);
  if ((pA>0)&&(pB<0)) return     + n(sA,+xA) - n(sB,-xB);
  if ((pA<0)&&(pB<0)) return - 1.- n(sA,-xA) - n(sB,-xB);
}

/*--------------------------------------------------------------------*/

double complex calF(double complex z, double *ab) {
  double a = ab[0], b = ab[1];
  double complex zz = z;// + I*1e-8; // not necessary?
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

  res = calA(al,be,ga)*res + calG(z1,pa*pd,nn,b); // or nn-1?

  return res;
}

