
/*--------------------------------------------------------------------*/

double Rate_1_to_3(double o, double k, void *A_, void *B_, void *C_) {
  double res[2], err[2];

  if (o<k) {printf("Below LC not yet implemented!"); return 0.;}
  double M = sqrt( SQR(o)-SQR(k) );

  double m1 = ((double *)A_)[0], m2 = ((double *)B_)[0], m3 = ((double *)C_)[0];
  double u1 = ((double *)A_)[1], u2 = ((double *)B_)[1], u3 = ((double *)C_)[1];
  int    s1 = ((double *)A_)[2], s2 = ((double *)B_)[2], s3 = ((double *)C_)[2];
  if (m1+m2+m3>M) { return 0.; }
  //printf(" m1 =  %g , m2 = %g , m3 = %g\n", m1, m2, m3 );
  //printf(" u1 =  %g , u2 = %g , u3 = %g\n", u1, u2, u3 );
  //printf(" s1 =  %d , s2 = %d , s3 = %d\n", s1, s2, s3 );

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2];

    double complex s12 = (1-_X_)*SQR(m1+m2) + _X_*SQR(M-m3) ;

    double complex lam_12 = csqrt( lam(s12,SQR(m1),SQR(m2)) ),
                   lam_3  = csqrt( lam(s12,SQR(M),SQR(m3)) );

    double complex q0  = .5*( o*(s12+SQR(M)-SQR(m3) ) + k*_Y_*(lam_3) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s12 );

    double complex e2  = .5*( q0*(s12+SQR(m2)-SQR(m1)) + q*_Z_*(lam_12) )/s12,
                   p2  = csqrt( SQR(e2) - SQR(m2) );

    double complex thermal_weight =  ( n(s1*s2,q0-u1-u2) - n(s3,q0-o+u3)  )
                                    *( n(s2,e2-u2)       - n(s1,e2-q0+u1) );//*/

    double complex jacobian =  ( SQR(M-m3) - SQR(m1+m2) )   // from X = [0,1]
                              *( .5*k*lam_3/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_12/s12 )          // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian);

    val[0] = creal(_inner); val[1] = cimag(_inner);
    //printf(" res = %g + I %g\n", val[0], val[1] );

    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  size_t calls = 1e7;
  {
    double    tol=1e-4;
    int       MaxEvls = calls;

    hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
    //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  }
  return res[0];

}

/*--------------------------------------------------------------------*/

double Rate_2_to_2(double o, double k, void *A_, void *B_, void *C_) {
  double res[2], err[2];

  if (o<k) {printf("Below LC not yet implemented!"); return 0.;}
  double M = sqrt( SQR(o)-SQR(k) );

  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s = fmax( SQR(M+M1), SQR(m1+m2) )/_X_;

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(M1)) ),
                   lam_12  = csqrt( lam(s,SQR(m1),SQR(m2)) );

    double complex q0  = .5*( o*(s+SQR(M)-SQR(M1) ) + k*_Y_*(lam_1) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s );

    double complex e2  = .5*( q0*(s+SQR(m2)-SQR(m1)) + q*_Z_*(lam_12) )/s,
                   p2  = csqrt( SQR(e2) - SQR(m2) );

    double complex thermal_weight =  ( n(t1,q0-o-U1) - n(s1*s2,q0-u1-u2)  )
                                    *( n(s2,e2-u2)       - n(s1,e2-q0+u1) );//*/

    double complex jacobian =  fmax( SQR(M+M1), SQR(m1+m2) )/SQR(_X_)   // from X = [0,1]
                              *( .5*k*lam_1/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_12/s )          // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian);

    val[0] = creal(_inner); val[1] = cimag(_inner);
    //printf(" res = %g + I %g\n", val[0], val[1] );

    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  size_t calls = 1e7;
  {
    double    tol=1e-5;
    int       MaxEvls = calls;

    hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
    //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  }
  return res[0];

}

/*--------------------------------------------------------------------*/

double Rate_3_to_1(double o, double k, void *A_, void *B_, void *C_) {
  double res[2], err[2];

  if (o<k) {printf("Below LC not yet implemented!"); return 0.;}
  double M = sqrt( SQR(o)-SQR(k) );

  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], U2 = ((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  if (m1>M+M1+M2) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s   = (1-_X_)*SQR(M1+M2) + _X_*SQR(m1-M) ;

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(m1)) ),
                   lam_12  = csqrt( lam(s,SQR(M1),SQR(M2)) );

    double complex q0  = .5*( -o*(s+SQR(M)-SQR(m1) ) + k*_Y_*(lam_1) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s );

    double complex E2  = .5*( q0*(s+SQR(M2)-SQR(M1)) + q*_Z_*(lam_12) )/s,
                   k2  = csqrt( SQR(E2) - SQR(M2) );

    double complex thermal_weight =  ( n(t1*t2,q0-U1-U2) - n(s1,q0+o-u1)  )
                                    *( n(t2,E2-U2)       - n(t1,E2-q0+U1) );//*/

    double complex jacobian =  ( SQR(m1-M) - SQR(M1+M2) )   // from X = [0,1]
                              *( .5*k*lam_1/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_12/s )            // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian);

    val[0] = creal(_inner);
    val[1] = cimag(_inner);
    //printf(" res = %g + I %g\n", val[0], val[1] );

    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  size_t calls = 1e7;
  {
    double    tol=1e-5;
    int       MaxEvls = calls;

    hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
    //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  }
  return res[0];

}

double rho_f(double o, double k, void *A_, void *B_, void *C_) {
  double K2 = SQR(o) - SQR(k);
  double res = 
    Rate_1_to_3(o,k,A_,B_,C_)

   +Rate_2_to_2(o,k,A_,B_,C_)
   +Rate_2_to_2(o,k,B_,A_,C_)
   +Rate_2_to_2(o,k,C_,A_,B_)

   +Rate_3_to_1(o,k,A_,B_,C_)
   +Rate_3_to_1(o,k,B_,C_,A_)
   +Rate_3_to_1(o,k,C_,A_,B_);
  return ( res - 4.*pow(OOFP,3.)*K2/8.)*64.*M_PI;
}

/*--------------------------------------------------------------------*/


