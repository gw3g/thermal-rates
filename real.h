/*--------------------------------------------------------------------*/

double    tol=1e-6;
int       MaxEvls = 1e7;

/*--------------------------------------------------------------------*/
// 1 \to 3 reactions in `s12-channel'

double Rate_1_to_3( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_ ) { double res[2], err[2];

  // reproducing eqs. (B.5)-(B.8)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double m1 = ((double *)A_)[0], m2 = ((double *)B_)[0], m3 = ((double *)C_)[0];
  double u1 = ((double *)A_)[1], u2 = ((double *)B_)[1], u3 = ((double *)C_)[1];
  int    s1 = ((double *)A_)[2], s2 = ((double *)B_)[2], s3 = ((double *)C_)[2];

  if (m1+m2+m3>M) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

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

    val[0] = creal(_inner); val[1] = cimag(_inner);//printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

/*--------------------------------------------------------------------*/
// 2 \to 2 reactions in s,t-channels

double Rate_2_to_2_sChan( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_) { double res[2], err[2];

  // reproducing eqs. (B.14)-(B.17)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

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

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

double Rate_2_to_2_tChan( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_) { double res[2], err[2];

  // reproducing eqs. (4.15)-(4.18)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex _inner(double complex t) {

      double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(m2)) ),
                     lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

      double complex q0  = .5*( o*(t+SQR(M)-SQR(m2) ) + k*_Y_*(lam_2) )/SQR(M); // [q0^+,q0^0]
      double complex q = csqrt( SQR(q0) - t );

      double complex e1p = .5*( q0*(t+SQR(m1)-SQR(M1)) + q*(lam_11) )/t,
                     e1m = .5*( q0*(t+SQR(m1)-SQR(M1)) - q*(lam_11) )/t, e1, jac_e1;

      if (creal(t)<0.) { e1 = (e1m-1.) + 1./_Z_; jac_e1 = 1./SQR(_Z_); }
      else { e1 = (1.-_Z_)*e1m + _Z_*e1p; jac_e1 = e1p-e1m; }

      double complex p1  = csqrt( SQR(e1) - SQR(m1) );

      double complex thermal_weight =  ( n(t1*s1,q0-u1+U1) - n(s1,q0-o+u2)  )
                                      *( n(t1,e1-q0-U1)    - n(s1,e1-u1)    );//*/

      double complex jacobian =
                                 ( .5*k*lam_2/SQR(M) )        // ..   Y = [-1,1]
                                *jac_e1          // ..   Z = [-1,1]
                                /(2.*q)                      ;//

      double prefactor = .5*pow(OOFP,3.)/k;

      return (prefactor)*(thermal_weight)*(jacobian);
    }

    double complex t1 = 1.-1./(1.-_X_), j1 = 1./SQR(1.-_X_); // (-inf,0]
    double complex j2 = fmin( SQR(M-m2), SQR(m1-M1) ), t2 = j2*_X_;  // [0,min(...))

    double complex _outer = _inner(t1)*j1;
    if ((M-m2)*(m1-M1) > 0.) { _outer+=_inner(t2)*j2; }

    val[0] = creal(_outer); val[1] = cimag(_outer); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0.0, -1., 0.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

/*--------------------------------------------------------------------*/
// 3 \to 1 reactions in s,t-channels

double Rate_3_to_1_sChan( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_) { double res[2], err[2];

  // reproducing eqs. (B.31)-(B.34)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], U2 = ((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  if (m1<M+M1+M2) { return 0.; }

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

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  printf("S-channel!! res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0]; // return the real part
}


double Rate_3_to_1_tChan( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_) { double res[2], err[2];

  // reproducing eqs. (B.23)-(B.26)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], U2 = ((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  if (m1<M+M1+M2) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex t   = (1-_X_)*SQR(M+M2) + _X_*SQR(m1-M1) ;

    double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(M2)) ),
                   lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

    double complex q0  = .5*( -o*(t+SQR(M)-SQR(M2) ) + k*_Y_*(lam_2) )/SQR(M);
    double complex q = csqrt( SQR(q0) - t );

    double complex E1  = .5*( q0*(t+SQR(M1)-SQR(m1)) + q*_Z_*(lam_11) )/t,
                   k1  = csqrt( SQR(E1) - SQR(M1) );

    double complex thermal_weight =  ( n(t1*s1,q0+u1-U2) - n(t2,q0+o+U2)  )
                                    *( n(t1,E1-U1)       - n(s1,E1-q0-u1) );//*/

    double complex jacobian =  ( SQR(m1-M1) - SQR(M+M2) )   // from X = [0,1]
                              *( .5*k*lam_2/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_11/t )            // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  printf("T-channel!! res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0]; // return the real part
}

/*--------------------------------------------------------------------*/
