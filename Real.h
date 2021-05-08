/*--------------------------------------------------------------------*/

double    TolReal=1e-8;
int       MaxEvls = 1e8;

#include "Theta.h"
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
  double (*func)(double,double,double complex,double complex,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (B.5)-(B.8)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

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

    double complex thermal_weight =  ( n(s1*s2,q0-u1-u2) - n(s3,q0-o+u3)  )
                                    *( n(s2,e2-u2)       - n(s1,e2-q0+u1) );//*/

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
  double (*func)(double,double,double complex,double complex,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (B.14)-(B.17)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];
  double M_[3] = {M1,m1,m2};

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s = fmax( SQR(M+M1), SQR(m1+m2) )/_X_;

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(M1)) ),
                   lam_12  = csqrt( lam(s,SQR(m1),SQR(m2)) );

    double complex q0  = .5*( o*(s+SQR(M)-SQR(M1) ) + k*_Y_*(lam_1) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s );

    double complex e2  = .5*( q0*(s+SQR(m2)-SQR(m1)) + q*_Z_*(lam_12) )/s;

    double complex thermal_weight =  ( n(t1,q0-o-U1) - n(s1*s2,q0-u1-u2) )
                                    *( n(s2,e2-u2)   - n(s1,e2-q0+u1)    );//*/

    double complex jacobian =  fmax( SQR(M+M1), SQR(m1+m2) )/SQR(_X_) // from X = [0,1]
                              *( .5*k*lam_1/SQR(M) )                  // ..   Y = [-1,1]
                              *( .5*q*lam_12/s )                      // ..   Z = [-1,1]
                              /(2.*q)                      ;          //

    double prefactor = .5*pow(OOFP,3.)/k;
    double rate = func(o, k, s, q0, e2, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

double Rate_2_to_2_tChan(o,k,A_,B_,C_,func)
  double o,k;
  void *A_, *B_, *C_;
  double (*func)(double,double,double complex,double complex,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (4.15)-(4.18)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], m1 = ((double *)B_)[0], m2 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], u1 = ((double *)B_)[1], u2 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], s1 = ((double *)B_)[2], s2 = ((double *)C_)[2];
  double M_[3] = {M1,m1,m2};

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

      if (creal(t)<0.) { e1 = (e1m)*(1.)/_Z_; jac_e1 = e1m/SQR(_Z_); }
      //double LAM = fmax(20.,q0+20.);
      //if (creal(t)<0.)      { e1 = (1.-_Z_)*e1m + _Z_*LAM; jac_e1 = LAM-e1m; }
      else if (creal(t)>0.) { e1 = (1.-_Z_)*e1m + _Z_*e1p; jac_e1 = e1p-e1m; }

      double complex thermal_weight =  ( n(t1*s1,q0-u1+U1) - n(s2,q0-o+u2)  )
                                      *( n(t1,e1-q0-U1)    - n(s1,e1-u1)    );//*/

      double complex jacobian = ( .5*k*lam_2/SQR(M) )         // ..   Y = [-1,1]
                                *jac_e1                       // ..   Z = [-1,1]
                                /(2.*q)                      ;//

      double prefactor = .5*pow(OOFP,3.)/k;
      //if ((M-m2)*(m1-M1) > 0.) { printf(" YES \n" ); }
      double rate = func(o, k, t, q0, e1, M_);

      return (prefactor)*(thermal_weight)*(jacobian)*(rate);
    }

    double complex t_1 = 1.-1./(_X_), j_1 = 1./SQR(_X_); // (-inf,0]
    double complex j_2 = fmin( SQR(M-m2), SQR(m1-M1) ), t_2 = j_2*_X_;  // [0,min(...))
    //printf(" j_2 = %.4f, t_2 = %.4f\n",j_2, t_2);

    double complex _outer = _inner(t_1)*j_1; //printf("(t<0) t_1 = %g\n", t_1);
    //if ((M-m2)*(m1-M1) > 0.) { printf("% M-m2 = %.4f and m1-M1 = %.4f\n------\n", M-m2,m1-M1); _outer+=_inner(t_2)*j_2; }//printf("(t>0) t_2 = %g\n", t_2);}
    if ((M-m2)*(m1-M1) > 0.) { _outer+=_inner(t_2)*j_2; }//printf("(t>0) t_2 = %g\n", t_2);}

    val[0] = creal(_outer); val[1] = cimag(_outer); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., 0.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

/*--------------------------------------------------------------------*/
// 3 \to 1 reactions in s,t-channels

double Rate_3_to_1_sChan(o,k,A_,B_,C_,func)
  double o,k;
  void *A_, *B_, *C_;
  double (*func)(double,double,double complex,double complex,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (B.31)-(B.34)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], U2 = ((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  double M_[3] = {M1,M2,m1};

  if (m1<M+M1+M2) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex s   = (1-_X_)*SQR(M1+M2) + _X_*SQR(m1-M) ;

    double complex lam_1   = csqrt( lam(s,SQR(M ),SQR(m1)) ),
                   lam_12  = csqrt( lam(s,SQR(M1),SQR(M2)) );

    double complex q0  = .5*( -o*(s+SQR(M)-SQR(m1) ) + k*_Y_*(lam_1) )/SQR(M);
    double complex q = csqrt( SQR(q0) - s );

    double complex E2  = .5*( q0*(s+SQR(M2)-SQR(M1)) + q*_Z_*(lam_12) )/s;

    double complex thermal_weight =  ( n(t1*t2,q0-U1-U2) - n(s1,q0+o-u1)  )
                                    *( n(t2,E2-U2)       - n(t1,E2-q0+U1) );//*/

    double complex jacobian =  ( SQR(m1-M) - SQR(M1+M2) )   // from X = [0,1]
                              *( .5*k*lam_1/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_12/s )            // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;
    double rate = func(o, k, s, q0, E2, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  //printf("S-channel!! res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0]; // return the real part
}


double Rate_3_to_1_tChan(o,k,A_,B_,C_,func)
  double o,k;         // K = (\omega,k)
  void *A_, *B_, *C_; // (mass, chem. potential, boson/fermion)
  double (*func)(double,double,double complex,double complex,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (B.23)-(B.26)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double M1 = ((double *)A_)[0], M2 = ((double *)B_)[0], m1 = ((double *)C_)[0];
  double U1 = ((double *)A_)[1], U2 = ((double *)B_)[1], u1 = ((double *)C_)[1];
  int    t1 = ((double *)A_)[2], t2 = ((double *)B_)[2], s1 = ((double *)C_)[2];
  double M_[3] = {M1,M2,m1};

  if (m1<M+M1+M2) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1], _Z_ = x[2]; // integration variables

    double complex t   = (1-_X_)*SQR(M+M2) + _X_*SQR(m1-M1) ;

    double complex lam_2   = csqrt( lam(t,SQR(M ),SQR(M2)) ),
                   lam_11  = csqrt( lam(t,SQR(m1),SQR(M1)) );

    double complex q0  = .5*( -o*(t+SQR(M)-SQR(M2) ) + k*_Y_*(lam_2) )/SQR(M);
    double complex q = csqrt( SQR(q0) - t );

    double complex E1  = .5*( q0*(t+SQR(M1)-SQR(m1)) + q*_Z_*(lam_11) )/t;

    double complex thermal_weight =  ( n(t1*s1,q0+u1-U1) - n(t2,q0+o+U2)  )
                                    *( n(t1,E1-U1)       - n(s1,E1-q0-u1) );//*/

    double complex jacobian =  ( SQR(m1-M1) - SQR(M+M2) )   // from X = [0,1]
                              *( .5*k*lam_2/SQR(M) )        // ..   Y = [-1,1]
                              *( .5*q*lam_11/t )            // ..   Z = [-1,1]
                              /(2.*q)                      ;//

    double prefactor = .5*pow(OOFP,3.)/k;
    double rate = func(o, k, t, q0, E1, M_);

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner); //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0., -1., -1.};
  double xu[3] = { 1., +1., +1.};

  hcubature(2, integrand, NULL, 3, xl, xu, MaxEvls, 0, TolReal, ERROR_INDIVIDUAL, res, err);
  //printf("T-channel!! res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0]; // return the real part
}

/*--------------------------------------------------------------------*/
