/*--------------------------------------------------------------------*/
// 1 \to 2 reactions 

double Rate_1_to_2( double o, double k, // K = (\omega,k)
                    void *A_, void *B_, void *C_, void *D_) 
                  { double res[2], err[2];

  // reproducing eqs. (4.20)-(4.28)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

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

    double complex cos_kpD = ( o*eD + .5*(SQR(mC)-SQR(mD)-SQR(M)) );

    double complex eA = mA + (1.-_Y_)/_Y_, pA = csqrt( SQR(eA) - SQR(mA) );
    double complex eB = mB + (1.-_Y_)/_Y_, pB = csqrt( SQR(eB) - SQR(mB) );
    //double _ab[2];

    //_ab[0] = creal(pA*pD*cos_kpD*_Z_);
    //_ab[1] = creal(pA*pD*csqrt(1.-SQR(cos_kpD))*csqrt(1.-SQR(_Z_)));
    double complex tempA = 
    pA*(  - (.0+n(sA,eA-uA))*calG(eA*eD+.5*(SQR(mB)-SQR(mA)-SQR(mD)), pA*pD)
          - (.0+n(sA,eA+uA))*calG(eA*eD-.5*(SQR(mB)-SQR(mA)-SQR(mD)), pA*pD)  
          );

    //_ab[0] = creal(pB*pD*cos_kpD*_Z_);
    //_ab[1] = creal(pB*pD*csqrt(1.-SQR(cos_kpD))*csqrt(1.-SQR(_Z_)));
    double complex tempB = 
    pB*(  - (.0+n(sB,eB-uB))*calG(eB*eD+.5*(SQR(mA)-SQR(mB)-SQR(mD)), pB*pD)
          - (.0+n(sB,eB+uB))*calG(eB*eD-.5*(SQR(mA)-SQR(mB)-SQR(mD)), pB*pD)  );

    double complex thermal_weight = 1. + n(sC,eC-uC) + n(sD,eD-uD) ;

    double complex jacobian =  ( eC_p - eC_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              //*( .5 )                       // ..   Z = [-1,1]
                              *SGN(creal(eC*eD))            ;//

    double prefactor = .5*pow(OOFP,3.)/k;
    printf(" res = %g + I %g\n", creal(tempB), cimag(tempB) );

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB);

    val[0] = creal(_inner); val[1] = cimag(_inner);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[3] = { 0.,  0.};
  double xu[3] = { 1., +1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, tol, tol, ERROR_INDIVIDUAL, res, err);
  printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

