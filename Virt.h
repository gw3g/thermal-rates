/*
 *  generic code: integrate 1 -> 2 reactions
 *
 *  K = (\omega,k)
 *  A_ = (m_a,\mu_a,\sig_a)
 *  [void *] argument of func() is a list of masses.
 *
 */

/*--------------------------------------------------------------------*/
// 1 \to 2 BORN reactions 

double Phi(o,k,A_,B_,func)
  double o,k;
  void *A_, *B_;
  double (*func)(double,double,double complex,void *);
{
  double res[2], err[2]; // reproducing eqs. (B.5)-(B.8)

  double K2 = SQR(o)-SQR(k) ;
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2];
  double M_[2] = {mA,mB};

  double complex lam_AB = csqrt( lam(K2,SQR(mA),SQR(mB)) );
  if (fabs(cimag(lam_AB))>1e-7) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];
    double complex eA = .5*( o*(K2+SQR(mA)-SQR(mB)) + k*_X_*lam_AB )/K2, eB = o-eA;
    double complex thermal_weight =  ( 1. + n(sA,eA-uA) + n(sB,eB-uB) );
    double complex jacobian = .5*k*lam_AB/K2;
    double complex rate = func(o,k,eA,M_);
    double prefactor = OOFP/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);

    val[0] = creal(_inner); val[1] = cimag(_inner);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[1] = { -1. };
  double xu[1] = { +1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  return res[0];
}

double Bubble(o,k,A_,B_,C_,D_,func)
  double o,k;
  void *A_, *B_, *C_, *D_;
  double (*func)(double,double,double *,double *,double *,double *,double complex *);
{
  double res[2], err[2];

  double M = sqrt( SQR(o)-SQR(k) );
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0], mC = ((double *)C_)[0], mD = ((double *)D_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1], uC = ((double *)C_)[1], uD = ((double *)D_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2], sC = ((double *)C_)[2], sD = ((double *)D_)[2];

  double complex lam_DC = csqrt( lam(SQR(M),SQR(mD),SQR(mC)) );
  if (fabs(cimag(lam_DC))>1e-7) { return 0.; }

  double complex eC_p  = .5*( o*(SQR(M)+SQR(mC)-SQR(mD)) + k*(lam_DC) )/SQR(M) ;
  double complex eC_m  = .5*( o*(SQR(M)+SQR(mC)-SQR(mD)) - k*(lam_DC) )/SQR(M) ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1];// _Z_ = x[2]; // integration variables

    double complex eC = (1-_X_)*eC_m + _X_*eC_p , eD = o - eC,
                   pC = csqrt( SQR(eC) - SQR(mC) ),
                   pD = csqrt( SQR(eD) - SQR(mD) );

    double complex eA = mA + (1.-_Y_)/_Y_, pA = csqrt( SQR(eA) - SQR(mA) );
    double complex eB = mB + (1.-_Y_)/_Y_, pB = csqrt( SQR(eB) - SQR(mB) );

    double PA[3] = {creal(eA),creal(pA),mA};
    double PB[3] = {creal(eB),creal(pB),mB};
    double PC[3] = {creal(eC),creal(pC),mC};
    double PD[3] = {creal(eD),creal(pD),mD};

    double complex _inner[4];

    func(o,k,PA,PB,PC,PD,_inner);

    double complex tempA = n(sA,eA-uA)*_inner[0] - n(sA,eA+uA)*_inner[1], // Eq.(4.25)
                   tempB = n(sB,eB-uB)*_inner[2] - n(sB,eB+uB)*_inner[3];

    double complex thermal_weight = 1. + n(sC,eC-uC) + n(sD,eD-uD) ;

    double complex jacobian =  ( eC_p - eC_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              *SGN(creal(eC*eD))            ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0.,  0.};
  double xu[3] = { 1., +1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


