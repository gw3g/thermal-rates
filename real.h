/*--------------------------------------------------------------------*/

double    TolReal = 1e-8;
int       MaxEvls = 1e8;

/*
 *  generic code: integrate n -> (4-n) reactions
 *
 *  K = (\omega,k)
 *  A_ = (m_a,\mu_a,\sig_a)
 *  [void *] argument of func() is a list of masses.
 *
 *  NB: should prefer K2 to M2 = M^2 (since M2 = M_2) !!
 */

/*--------------------------------------------------------------------*/
// 1 \to 3 reactions in `s12-channel'

double Rate_1_to_3(o,k,A_,B_,C_,func)
  double o,k;
  void *A_, *B_, *C_;
  double (*func)(double,double, // o, k
                 double complex,double complex,double complex,
                 void *); // pointer to masses
{
  double res[2], err[2]; // reproducing eqs. (B.5)-(B.8)

  if (o<k) { printf("Below light-cone not implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double m1 = ((double *)A_)[0], m2 = ((double *)B_)[0], m3 = ((double *)C_)[0];
  double u1 = ((double *)A_)[1], u2 = ((double *)B_)[1], u3 = ((double *)C_)[1];
  int    s1 = ((double *)A_)[2], s2 = ((double *)B_)[2], s3 = ((double *)C_)[2];

  double M_[3] = {m1,m2,m3}; // mass list

  if (m1+m2+m3>M) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s12 = (1-_X_)*SQR(m1+m2) + _X_*SQR(M-m3) ;

    double complex lam_12 = csqrt( lam(s12,SQR(m1),SQR(m2)) ),
                   lam_3  = csqrt( lam(s12,SQR(M),SQR(m3)) );

    double complex q0  = .5*( o*(s12+SQR(M)-SQR(m3) ) + k*_Y_*(lam_3) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s12 );

    double complex e2  = .5*( q0*(s12+SQR(m2)-SQR(m1)) + q*_Z_*(lam_12) )/s12;
    double complex e1  = q0 - e2, e3 = o - q0;

    double complex thermal_weight =  calN(s1*s2,s3,q0-u1-u2,o-q0-u3)
                                         *calN(s2,s1,e2-u2,q0-e2-u1);

    double complex jacobian =  ( SQR(M-m3) - SQR(m1+m2) )   // from X = [0,1]
                              *( .5*k*lam_3/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_12/s12 )          // ..   Z = [-1,1]
                              /(2.*q)                      ;//


    double prefactor = .5*pow(OOFP,3.)/k;
    double rate = func(o, k, s12, q0, e2, M_);
    //printf("s12 = %g, q0 = %g, e2 = %g, func = %g\n", creal(s12), creal(q0), creal(e2), func);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner);//printf(" [%g,%g,%g] res = %g + I %g\n",
                                                   //       _X_,_Y_,_Z_,val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

/*--------------------------------------------------------------------*/
// 2 \to 2 reactions in s,t-channels

double Rate_2_to_2_sChan(o,k,A_,B_,C_,func)
  double o,k;
  void *A_, *B_, *C_;
  double (*func)(double,double, // o, k
                 double complex,double complex,double complex,
                 void *); // pointer to masses
{
  double res[2], err[2]; // reproducing eqs. (B.14)-(B.17)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 =-((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];
  double M_[3] = {M1,m1,m2};
  double s_min = fmax( SQR(M+M1), SQR(m1+m2) );
  //printf("2->2, s channel:  M = %g, M1 = %g, m1 = %g, m2 = %g \n",M,M1,m1,m2); 
  //printf(" s_min = %g\n",s_min); 

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s = s_min - SQR(M)*(1.-1./_X_);

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(M1)) ),
                   lam_12  = csqrt( lam(s,SQR(m1),SQR(m2)) );

    double alpha = 1., _Y2_ = tanh(alpha*_Y_)/tanh(alpha),
                       _Z2_ = tanh(alpha*_Z_)/tanh(alpha);

    double complex q0  = .5*( o*(s+SQR(M)-SQR(M1) ) + k*_Y2_*(lam_1) )/(SQR(M));
    double complex q = csqrt( SQR(q0) - s );
    //printf(" q0 = %g, (s+SQR(M)-SQR(M1) = %g, _Y_ = %g \n",creal(q0), creal(s+SQR(M)-SQR(M1)),_Y_ );

    double complex e2  = .5*( q0*(s+SQR(m2)-SQR(m1)) + q*_Z2_*(lam_12) )/s;

    double complex thermal_weight = calN(t1,s1*s2,q0-o-U1,u1+u2-q0)*calN(s2,s1,e2-u2,q0-e2-u1);

    double complex jacobian =  SQR(M/_X_) // from X = [0,1]
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Y_)))
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Z_)))
                              *( .5*lam_1/SQR(M) )                  // ..   Y = [-1,1]
                              *( .5*lam_12/s )                      // ..   Z = [-1,1]
                              /(2.)                      ;          //

    double prefactor = .5*pow(OOFP,3.);
    double rate = func(o, k, s-I*1e-9, q0, e2, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner);
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  return res[0]; // return the real part
}

double Rate_2_to_2_tChan(o,k,A_,B_,C_,func)
  double o,k;
  void *A_, *B_, *C_;
  double (*func)(double,double, // o, k
                 double complex,double complex,double complex,
                 void *); // pointer to masses
{
  double res[2], err[2]; // reproducing eqs. (4.15)-(4.18)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 =-((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];
  double M_[3] = {M1,m1,m2};
  //printf("2->2, t channel:  M = %g, M1 = %g, m1 = %g, m2 = %g \n",M,M1,m1,m2); 
  double t_max = fmin( SQR(M-m2), SQR(m1-M1) );
  /*if ((M-m2)*(m1-M1) > 0.) { 
    printf(" YES: M = %g, M1 = %g, m1 = %g, m2 = %g \n",M,M1,m1,m2); 
    printf(" t_max = %g\n",t_max); 
  }//*/

  int integrand1(unsigned dim,  const double *x, void *data_,
                 unsigned fdim, double *val) { // t<0

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex t = SQR(M)*(1.-1./(_X_));

    double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(m2)) ),
                   lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

    double alpha = 1., _Y2_ = tanh(alpha*_Y_)/tanh(alpha);

    double complex q0  = .5*( o*(t+SQR(M)-SQR(m2) ) + k*_Y2_*(lam_2) )/SQR(M); // [q0^+,q0^0]
    double complex q = csqrt( SQR(q0) - t );

    double complex e1m = .5*( q0*(t+SQR(m1)-SQR(M1)) - q*(lam_11) )/t,
                   e1 = (e1m)-fabs(e1m)*(1.-1./(_Z_));

    double complex thermal_weight;
    if (fabs(creal(q0-u1+U1))>1e-2) {
                   thermal_weight =  - calN(t1*s1,s2,q0-u1+U1,o-q0-u2)*calN(t1,s1,q0-e1+U1,e1-u1);
    } else {
                   thermal_weight =  ( 1. + (1.+n(s2,o-q0-u2))/n(t1*s1,q0-u1+U1) )
                                    *( 1. + n(t1,e1-q0-U1) )*n(s1,e1-u1); }
                   //printf(" q0-u1+U1 = %g , thermal_weight = %g \n", creal(q0-u1+U1), creal(thermal_weight) ); }

    double complex jacobian = ( SQR(M/_X_) )                //      X = [0,1]
                              *( .5*lam_2/SQR(M) )          // ..   Y = [-1,1]
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Y_)))
                              *( fabs(e1m)/SQR(_Z_) )         // ..   Z = [0,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.);
      //if ((M-m2)*(m1-M1) > 0.) { printf(" YES \n" ); }
    double rate = func(o, k, t, q0, e1, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  int integrand2(unsigned dim,  const double *x, void *data_,
                 unsigned fdim, double *val) { // t>0

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex t = t_max*_X_;

    double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(m2)) ),
                   lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

    double alpha = 1., _Y2_ = tanh(alpha*_Y_)/tanh(alpha),
                       _Z2_ = tanh(alpha*(2.*_Z_-1.))/tanh(alpha);

    double complex q0 = .5*( o*(t+SQR(M)-SQR(m2) ) + k*_Y2_*(lam_2) )/SQR(M); // [q0^+,q0^0]
    double complex q  = csqrt( SQR(q0) - t );

    double complex e1 = .5*( q0*(t+SQR(m1)-SQR(M1)) + q*_Z2_*(lam_11) )/t;

    double complex thermal_weight;

    if (fabs(q0-u1+U1)>1e-2) {
                   thermal_weight = - calN(t1*s1,s2,q0-u1+U1,o-q0-u2)*calN(t1,s1,q0-e1+U1,e1-u1);
    } else {
                   thermal_weight =  ( 1. + (1.+n(s2,o-q0-u2))/n(t1*s1,q0-u1+U1) )
                                    *( 1. + n(t1,e1-q0-U1) )*n(s1,e1-u1); 
    }

    double complex jacobian = t_max                      //      X = [0,1]
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Y_)))
                              *2.*alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*(2.*_Z_-1.))))
                              *( .5*lam_2/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*lam_11/t )                 // ..   Z = [0,1]
                              /(2.)                      ;//

    double prefactor = .5*pow(OOFP,3.);
    double rate = func(o, k, t, q0, e1, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); 
    return 0;
  }

  double xl[3] = { 0., -1., 0.};
  double xu[3] = { 1., +1., +1.};

  double temp = 0.;
  hcubature(2, integrand1, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  temp += res[0];
  if ((M-m2)*(m1-M1) > 0.) {
    hcubature(2, integrand2, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
    temp += res[0];
  }//*/
  return (temp);
}

/*--------------------------------------------------------------------*/
// 3 \to 1 reactions in s,t-channels

double Rate_3_to_1_sChan(o,k,A_,B_,C_,func)
  double o,k;         // K = (\omega,k)
  void *A_, *B_, *C_; // (mass, chem. potential, boson/fermion)
  double (*func)(double,double, // o, k
                 double complex,double complex,double complex,
                 void *); // pointer to masses
{
  double res[2], err[2]; // reproducing eqs. (B.31)-(B.34)

  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 =-((double *)A_)[1], U2 =-((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  double M_[3] = {M1,M2,m1};

  if (m1<M+M1+M2) return 0.;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s   = (1-_X_)*SQR(M1+M2) + _X_*SQR(m1-M) ;

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(m1)) ),
                   lam_12  = csqrt( lam(s,SQR(M1),SQR(M2)) );

    double alpha = 1., _Y2_ = tanh(alpha*_Y_)/tanh(alpha),
                       _Z2_ = tanh(alpha*_Z_)/tanh(alpha);

    double complex q0  = .5*( -o*(s+SQR(M)-SQR(m1) ) + k*_Y2_*(lam_1) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s );

    double complex E2  = .5*( q0*(s+SQR(M2)-SQR(M1)) + q*_Z2_*(lam_12) )/s;

    double complex thermal_weight =  calN(t1*t2,s1,q0-U1-U2,u1-o-q0)*calN(t2,t1,E2-U2,q0-E2-U1);

    double complex jacobian =  ( SQR(m1-M) - SQR(M1+M2) )   // from X = [0,1]
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Y_)))
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Z_)))
                              *( .5*lam_1/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*lam_12/s )            // ..   Z = [-1,1]
                              /(2.)                      ;//

    double prefactor = .5*pow(OOFP,3.);
    double rate = func(o, k, s, q0, E2, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  return res[0]; // return the real part
}


double Rate_3_to_1_tChan(o,k,A_,B_,C_,func)
  double o,k;         // K = (\omega,k)
  void *A_, *B_, *C_; // (mass, chem. potential, boson/fermion)
  double (*func)(double,double, // o, k
                 double complex,double complex,double complex,
                 void *); // pointer to masses
{
  double res[2], err[2]; // reproducing eqs. (B.23)-(B.26)

  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 =-((double *)A_)[1], U2 =-((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  double M_[3] = {M1,M2,m1};

  if (m1<M+M1+M2) return 0.;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex t   = (1-_X_)*SQR(M+M2) + _X_*SQR(m1-M1) ;

    double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(M2)) ),
                   lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

    double alpha = 1., _Y2_ = tanh(alpha*_Y_)/tanh(alpha),
                       _Z2_ = tanh(alpha*_Z_)/tanh(alpha);

    double complex q0  = .5*( -o*(t+SQR(M)-SQR(M2) ) + k*_Y2_*(lam_2) )/SQR(M);
    double complex q = csqrt( SQR(q0) - t );

    double complex E1  = .5*( q0*(t+SQR(M1)-SQR(m1)) + q*_Z2_*(lam_11) )/t;

    double complex thermal_weight = calN(t1*s1,t2,q0+u1-U1,-q0-o-U2)*calN(t1,s1,E1-U1,q0+u1-E1);

    double complex jacobian =  ( SQR(m1-M1) - SQR(M+M2) ) // from X = [0,1]
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Y_)))
                              *alpha*cosh(alpha)/(sinh(alpha)*SQR(cosh(alpha*_Z_)))
                              *( .5*lam_2/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*lam_11/t )            // ..   Z = [-1,1]
                              /(2.)                      ;//

    double prefactor = .5*pow(OOFP,3.);
    double rate = func(o, k, t, q0, E1, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); 
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  return res[0]; // return the real part
}

/*--------------------------------------------------------------------*/
