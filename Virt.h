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

  double xl[2] = { 0.,  0.};
  double xu[2] = { 1., +1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


double Triangle(o,k,A_,B_,C_,D_,E_,func)
  double o,k;
  void *A_, *B_, *C_, *D_, *E_;
  double (*func)(double,double,double *,double *,double *,double *, double *,double complex *);
{
  double res[2], err[2];

  double M = sqrt( SQR(o)-SQR(k) );
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0], mC = ((double *)C_)[0], mD = ((double *)D_)[0], mE = ((double *)E_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1], uC = ((double *)C_)[1], uD = ((double *)D_)[1], uE = ((double *)E_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2], sC = ((double *)C_)[2], sD = ((double *)D_)[2], sE = ((double *)E_)[2];

  double complex lam_DE = csqrt( lam(SQR(M),SQR(mD),SQR(mE)) );
  if (fabs(cimag(lam_DE))>1e-7) { return 0.; }

  double complex eD_p  = .5*( o*(SQR(M)+SQR(mD)-SQR(mE)) + k*(lam_DE) )/SQR(M) ;
  double complex eD_m  = .5*( o*(SQR(M)+SQR(mD)-SQR(mE)) - k*(lam_DE) )/SQR(M) ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1];// _Z_ = x[2]; // integration variables

    double complex eD = (1-_X_)*eD_m + _X_*eD_p , eE = o - eD,
                   pD = csqrt( SQR(eD) - SQR(mD) ),
                   pE = csqrt( SQR(eE) - SQR(mE) );

    double complex eA = mA + (1.-_Y_)/_Y_, pA = csqrt( SQR(eA) - SQR(mA) );
    double complex eB = mB + (1.-_Y_)/_Y_, pB = csqrt( SQR(eB) - SQR(mB) );
    double complex eC = mC + (1.-_Y_)/_Y_, pC = csqrt( SQR(eC) - SQR(mC) );

    double PA[3] = {creal(eA),creal(pA),mA};
    double PB[3] = {creal(eB),creal(pB),mB};
    double PC[3] = {creal(eC),creal(pC),mC};
    double PD[3] = {creal(eD),creal(pD),mD};
    double PE[3] = {creal(eE),creal(pE),mE};

    double complex _inner[8];

    func(o,k,PA,PB,PC,PD,PE,_inner);

    double complex tempA = + n(sA,eA-uA)*_inner[0] + n(sA,eA+uA)*_inner[1], // Eq.(4.24)
                   tempB = - n(sB,eB-uB)*( _inner[2] + _inner[4] )
                           - n(sB,eB+uB)*( _inner[3] + _inner[5] ),
                   tempC = + n(sC,eC-uC)*_inner[6] + n(sC,eC+uC)*_inner[7];

    double complex thermal_weight = 1. + n(sD,eD-uD) + n(sE,eE-uE) ;

    double complex jacobian =  ( eD_p - eD_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              *SGN(creal(eD*eE))            ;//

    double prefactor = .25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB+tempC);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0.,  0.};
  double xu[2] = { 1., +1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


/*--------------------------------------------------------------------*/
// Vacuum functions

double L1(double M2) {
  return -M2*( log(fabs( SQR(mubar)/M2)) + 1. );
}

double L2(double M2a, double M2b, double M2d) {
  double res[1], err[1];

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {
    double _X_ = x[0];
    val[0] = log(fabs( SQR(mubar)/(M2a*_X_+M2b*(1.-_X_)-M2d*_X_*(1.-_X_)) ));
    return 0;
  }

  double xl[1] = { 0. };
  double xu[1] = { 1. };

  hcubature(1, integrand, NULL, 1, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  return -res[0];
}

double L3(double K2,
          double M2a, double M2b, double M2c, double M2d, double M2e) {
  double res[1], err[1];

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {
    double _Y_ = x[0];

    double al  = SQR(_Y_)*M2d,
           be  = _Y_*(M2a-M2b-M2d*_Y_+(M2e-K2)*(1.-_Y_)),
           ga  = M2b*_Y_+M2c*(1.-_Y_)-M2e*_Y_*(1.-_Y_);

    val[0] = _Y_*calA(al,be,ga);
    return 0;
  }

  double xl[1] = { 0. };
  double xu[1] = { 1. };

  hcubature(1, integrand, NULL, 1, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  return res[0];
}


double Bubble_0(o,k,A_,B_,C_,D_,c)
  double o,k;
  void *A_, *B_, *C_, *D_;
  double *c;
{
  double res[2], err[2];

  double M = sqrt( SQR(o)-SQR(k) );
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0], mC = ((double *)C_)[0], mD = ((double *)D_)[0];
  double uC = ((double *)C_)[1], uD = ((double *)D_)[1];
  int    sC = ((double *)C_)[2], sD = ((double *)D_)[2];

  double L2ab = L2(SQR(mA),SQR(mB),SQR(mD)),
         L1a  = L1(SQR(mA)),
         L1b  = L1(SQR(mB));

  double ca = ((double *)c)[0];
  double cb = ((double *)c)[1];
  double cK = ((double *)c)[2];

  double complex lam_DC = csqrt( lam(SQR(M),SQR(mD),SQR(mC)) );
  if (fabs(cimag(lam_DC))>1e-7) { return 0.; }

  double complex eC_p  = .5*( o*(SQR(M)+SQR(mC)-SQR(mD)) + k*(lam_DC) )/SQR(M) ;
  double complex eC_m  = .5*( o*(SQR(M)+SQR(mC)-SQR(mD)) - k*(lam_DC) )/SQR(M) ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];// integration variables

    double complex eC = (1-_X_)*eC_m + _X_*eC_p , eD = o - eC,
                   pC = csqrt( SQR(eC) - SQR(mC) ),
                   pD = csqrt( SQR(eD) - SQR(mD) );

    double complex temp = cb*L2ab
                        + .5*( (ca-cb)/SQR(mD) )*( (SQR(mD)+SQR(mA)-SQR(mB))*L2ab
                                                 + L1b - L1a );

    double E_PD, E_K;
    if (E=='K') {
      E_PD = .5*( SQR(M) + SQR(mD) - SQR(mC) );
      E_K = SQR(M);
    }
    if (E=='U') {
      E_PD = eD;
      E_K  = o;
    }
    temp *= E_PD; temp += cK*E_K*L2ab;

    double complex thermal_weight = 1. + n(sC,eC-uC) + n(sD,eD-uD) ;

    double complex jacobian =  ( eC_p - eC_m )              // from X = [0,1]
                              *SGN(creal(eC*eD))            ;//

    double prefactor = -.25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(temp);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }//*/

  double xl[1] = { 0. };
  double xu[1] = { 1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


double Triangle_0(o,k,A_,B_,C_,D_,E_,c)
  double o,k;
  void *A_, *B_, *C_, *D_, *E_;
  double *c;
{
  double res[2], err[2];

  double M = sqrt( SQR(o)-SQR(k) );
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0], mC = ((double *)C_)[0], mD = ((double *)D_)[0], mE = ((double *)E_)[0];
  double uD = ((double *)D_)[1], uE = ((double *)E_)[1];
  int    sD = ((double *)D_)[2], sE = ((double *)E_)[2];

  double L2bc    = L2(SQR(mB),SQR(mC),SQR(mE)),
         L2ac    = L2(SQR(mA),SQR(mC),SQR(M)),
         L2ab    = L2(SQR(mA),SQR(mB),SQR(mD)),
         L3abc   = L3(SQR(M),SQR(mA),SQR(mB),SQR(mC),SQR(mD),SQR(mE));

  double Ib = (SQR(mD)+SQR(mA)-SQR(mB))*L3abc + L2bc - L2ac,
         Ic = (SQR(M) +SQR(mA)-SQR(mC))*L3abc + L2bc - L2ab;

  double al_b = - 2.*SQR(M)*Ib + (SQR(mD)-SQR(mE)+SQR(M))*Ic,
         al_c = + (SQR(mD)-SQR(mE)+SQR(M))*Ib - 2.*SQR(mD)*Ic;

  double ca = ((double *)c)[0];
  double cb = ((double *)c)[1];
  double cc = ((double *)c)[2];

  double complex lam_DE = csqrt( lam(SQR(M),SQR(mD),SQR(mE)) );
  if (fabs(cimag(lam_DE))>1e-7) { return 0.; }

  al_b *= 1./lam_DE;
  al_c *= 1./lam_DE;

  double complex eD_p  = .5*( o*(SQR(M)+SQR(mD)-SQR(mE)) + k*(lam_DE) )/SQR(M) ;
  double complex eD_m  = .5*( o*(SQR(M)+SQR(mD)-SQR(mE)) - k*(lam_DE) )/SQR(M) ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];// _Z_ = x[2]; // integration variables

    double complex eD = (1-_X_)*eD_m + _X_*eD_p , eE = o - eD,
                   pD = csqrt( SQR(eD) - SQR(mD) ),
                   pE = csqrt( SQR(eE) - SQR(mE) );

    double E_PD, E_K;
    if (E=='K') {
      E_PD = .5*( SQR(M) + SQR(mD) - SQR(mE) );
      E_K = SQR(M);
    }
    if (E=='U') {
      E_PD = eD;
      E_K  = o;
    }
    double complex temp = al_b*E_PD+al_c*E_K;
    temp *= (ca-cb-cc);
    temp += (cb*E_PD+cc*E_K)*L3abc;


    double complex thermal_weight = 1. + n(sD,eD-uD) + n(sE,eE-uE) ;

    double complex jacobian =  ( eD_p - eD_m )              // from X = [0,1]
                              *SGN(creal(eD*eE))            ;//

    double prefactor = .25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(temp);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0. };
  double xu[2] = { 1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}
