/*
 *  generic code: integrate 1 -> 2 reactions
 *
 *  K = (\omega,k)
 *  A_ = (m_a,\mu_a,\sig_a)
 *  [void *] argument of func() is a list of masses.
 *
 */
double    TolVirt=1e-9;

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
  if (fabs(cimag(lam_AB))>0.) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];
    double complex eA = .5*( o*(K2+SQR(mA)-SQR(mB)) + k*_X_*lam_AB )/K2, eB = o-eA;
    //double complex thermal_weight =  ( 1. + n(sA,eA-uA) + n(sB,eB-uB) );
    double complex thermal_weight =  calN(sA,sB,eA-uA,eB-uB);
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

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, TolVirt, 0, ERROR_INDIVIDUAL, res, err);
  return res[0];
}

/*--------------------------------------------------------------------*/
// B-operator: thermal part

double Bubble_T(o,k,A_,B_,C_,D_,func)
  double o,k;
  void *A_, *B_, *C_, *D_;
  double (*func)(double,double,double *,double *,double *,double *,double complex *);
{
  double res[2], err[2]; // eq. (2.29)

  double K2 = SQR(o)-SQR(k) ;
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0],
         mC = ((double *)C_)[0], mD = ((double *)D_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1],
         uC = ((double *)C_)[1], uD = ((double *)D_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2],
         sC = ((double *)C_)[2], sD = ((double *)D_)[2];

  double complex lam_DC = csqrt( lam(K2,SQR(mD),SQR(mC)) );

  if (fabs(cimag(lam_DC))>0.) { return 0.; } // catch boundary cases

  double complex eC_p  = .5*( o*(K2+SQR(mC)-SQR(mD)) + k*(lam_DC) )/K2 ,
                 eC_m  = .5*( o*(K2+SQR(mC)-SQR(mD)) - k*(lam_DC) )/K2 ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1]; // integration variables

    double complex eC = (1-_X_)*eC_m + _X_*eC_p , eD = o - eC,
                   pC = csqrt( SQR(eC) - SQR(mC) ),
                   pD = csqrt( SQR(eD) - SQR(mD) );

    double complex eA = mA + (1.-_Y_)/_Y_, pA = csqrt( SQR(eA) - SQR(mA) ),
                   eB = mB + (1.-_Y_)/_Y_, pB = csqrt( SQR(eB) - SQR(mB) );

    double PA[3] = {creal(eA),creal(pA),mA};
    double PB[3] = {creal(eB),creal(pB),mB};
    double PC[3] = {creal(eC),creal(pC),mC};
    double PD[3] = {creal(eD),creal(pD),mD};

    double complex _inner[4];

    func(o,k,PA,PB,PC,PD,_inner);

    double complex tempA = n(sA,eA-uA)*_inner[0] - n(sA,eA+uA)*_inner[1], // Eq.(4.25)
                   tempB = n(sB,eB-uB)*_inner[2] - n(sB,eB+uB)*_inner[3];

    //double complex thermal_weight = 1. + n(sC,eC-uC) + n(sD,eD-uD) ;
    double complex thermal_weight = calN(sC,sD,eC-uC,eD-uD);

    double complex jacobian =  ( eC_p - eC_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              *SGN(creal(eC*eD))           ;//

    double prefactor = .5*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0., 0.};
  double xu[2] = { 1., 1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, 0, TolVirt, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

/*--------------------------------------------------------------------*/
// C-operator: thermal part

double Triangle_T(o,k,A_,B_,C_,D_,E_,func)
  double o,k;
  void *A_, *B_, *C_, *D_, *E_;
  double (*func)(double,double,double *,double *,double *,double *, double *,double complex *);
{
  double res[2], err[2]; // eq.(2.36)

  double K2 = SQR(o)-SQR(k);
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0],
         mC = ((double *)C_)[0], mD = ((double *)D_)[0], mE = ((double *)E_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1],
         uC = ((double *)C_)[1], uD = ((double *)D_)[1], uE = ((double *)E_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2],
         sC = ((double *)C_)[2], sD = ((double *)D_)[2], sE = ((double *)E_)[2];

  double complex lam_DE = csqrt( lam(K2,SQR(mD),SQR(mE)) );

  if (fabs(cimag(lam_DE))>0.) { return 0.; } // catch boundary cases

  double complex eD_p  = .5*( o*(K2+SQR(mD)-SQR(mE)) + k*(lam_DE) )/K2 ,
                 eD_m  = .5*( o*(K2+SQR(mD)-SQR(mE)) - k*(lam_DE) )/K2 ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1]; // integration variables

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

    //double complex thermal_weight = 1. + n(sD,eD-uD) + n(sE,eE-uE) ;
    double complex thermal_weight = calN(sD,sE,eD-uD,eE-uE);

    double complex jacobian =  ( eD_p - eD_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              *SGN(creal(eD*eE))           ;//

    double prefactor = -.25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB+tempC);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0., 0.};
  double xu[2] = { 1., 1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, 0, TolVirt, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}

// TEST!
double Triangle_0_new(o,k,A_,B_,C_,D_,E_,func)
  double o,k;
  void *A_, *B_, *C_, *D_, *E_;
  double (*func)(double,double,double *,double *,double *,double *, double *,double complex *);
{
  double res[2], err[2]; // eq.(2.36)

  double K2 = SQR(o)-SQR(k);
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0],
         mC = ((double *)C_)[0], mD = ((double *)D_)[0], mE = ((double *)E_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1],
         uC = ((double *)C_)[1], uD = ((double *)D_)[1], uE = ((double *)E_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2],
         sC = ((double *)C_)[2], sD = ((double *)D_)[2], sE = ((double *)E_)[2];

  double complex lam_DE = csqrt( lam(K2,SQR(mD),SQR(mE)) );

  if (fabs(cimag(lam_DE))>0.) { return 0.; } // catch boundary cases

  double complex eD_p  = .5*( o*(K2+SQR(mD)-SQR(mE)) + k*(lam_DE) )/K2 ,
                 eD_m  = .5*( o*(K2+SQR(mD)-SQR(mE)) - k*(lam_DE) )/K2 ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1]; // integration variables

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

    double complex tempA = + _inner[0] + _inner[1], // Eq.(4.24)
                   tempB = - ( _inner[2] + _inner[4] )
                           - ( _inner[3] + _inner[5] ),
                   tempC = + _inner[6] + _inner[7];

    //double complex thermal_weight = 1. + n(sD,eD-uD) + n(sE,eE-uE) ;
    double complex thermal_weight = calN(sD,sE,eD-uD,eE-uE);

    double complex jacobian =  ( eD_p - eD_m )              // from X = [0,1]
                              *( 1./SQR(_Y_) )              // ..   Y = [0,1]
                              *SGN(creal(eD*eE))           ;//

    double prefactor = -.25*pow(OOFP,3.)/k;

    double complex _outer = .5*(prefactor)*(thermal_weight)*(jacobian)*(tempA+tempB+tempC);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0., 0.};
  double xu[2] = { 1., 1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, TolVirt, 0, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


/*--------------------------------------------------------------------*/
// Vacuum functions

/*double Li2(double x) { // Real part only, for x>1
  double omx = 1.-x;
  if (x<-1.)         return +Li2(1./omx)-log(omx)*log(-x)+.5*SQR(log(omx))-ZE2;
  if (-1.<x && x<0.) return -Li2(-x/omx)-.5*SQR(log(omx));
  if ( .5<x && x<1.) return -Li2(omx)-log(omx)*log(x)+ZE2;
  if ( 1.<x && x<2.) return +Li2(-omx/x)-log(-omx)*log(x)+.5*SQR(log(x))+ZE2;
  if (x>2.)          return -Li2(1./x)-.5*SQR(log(x))+.5*ZE2;

  double tmp = 0.0, i = 1.0;
  for (; i < 50.; i += 1.0) { tmp += pow(x,i)/SQR(i); }
  return tmp;
}//*/

double GRID[16][2] = { // Gaussian-quadrature: 16-pt [-1,1]
                       {0.1894506104550685,-0.0950125098376374},
                       {0.1894506104550685,+0.0950125098376374},
                       {0.1826034150449236,-0.2816035507792589},
                       {0.1826034150449236,+0.2816035507792589},
                       {0.1691565193950025,-0.4580167776572274},
                       {0.1691565193950025,+0.4580167776572274},
                       {0.1495959888165767,-0.6178762444026438},
                       {0.1495959888165767,+0.6178762444026438},
                       {0.1246289712555339,-0.7554044083550030},
                       {0.1246289712555339,+0.7554044083550030},
                       {0.0951585116824928,-0.8656312023878318},
                       {0.0951585116824928,+0.8656312023878318},
                       {0.0622535239386479,-0.9445750230732326},
                       {0.0622535239386479,+0.9445750230732326},
                       {0.0271524594117541,-0.9894009349916499},
                       {0.0271524594117541,+0.9894009349916499}
                     };


double complex Li2(double complex z) { // Osacar, Palacian & Palacios (1995)
  double r = cabs(z);
  double complex omz = 1.-z, res;
  //printf( " z = %g + I %g \n", z);
  if (cabs(omz)<1e-8) return ZE2;
  if        (r<.5) { // circle
    res = 0.;
    for (int i=1;i<30;i++) {
      res += cpow(z,i)/SQR(i*(i+1.)*(i+2.));
    }
    res *= 4.*SQR(z);
    res += 4.*z + 5.75*SQR(z) + 3.*(1.-SQR(z))*clog(omz);
    return res/(SQR(z)+4.*z+1.);
  } else if (r<1.) { // crown
    if (cabs(omz)<.5) {
      return - Li2(omz) - clog(z)*clog(omz) + ZE2;
    } else {
      res = 0.;
      double ai, ti;
      for (int i=0;i<16;i++) {
        ai = .5*GRID[i][0]; ti = .5*(GRID[i][1]+1.);
        res += ai*clog(1.-z*ti)/ti;
      }
      return res;
    }
  } else {
    return - Li2(1./z) - .5*SQR(clog(-z)) - ZE2;
  }
}

double complex R(double complex y0, double complex y1) {
  double complex omy0 = 1.-y0, omy1 = 1.-y1, y1my0=y1-y0;
  //printf( " y0 = %g + I %g ... 1-y0 = %g + I %g \n", y0, omy0);
  double complex tmp = Li2(-omy1/y1my0) - Li2(y1/y1my0);
  if ( (cabs(omy0)>1e-12) && (cabs(omy1)>1e-12) ) {tmp += +clog(omy0/y1my0)*clog(-omy1/y1my0);}
  if ( (cabs(  y0)>1e-12) && (cabs(y1  )>1e-12) ) {tmp += -clog(-y0/y1my0)*clog(y1/y1my0);}
  //printf( " result = %g + I %g \n", tmp);
  return (tmp);
}

double complex S3(double complex y0,
                  double complex a, double complex b, double complex c) {
  double complex D = SQR(b)-4.*a*c;
  double complex y1 = ( .5*(-b-csqrt(D))/a ),
                 y2 = ( .5*(-b+csqrt(D))/a );

  //printf( " y1 = %g + I %g ... y2 = %g + I %g \n", y1, y2);
  return ( R(y0,y1) + R(y0,y2) );
}


double L1(double M2) {
  return -M2*( log(fabs(SQR(mubar)/M2)) + 1. );
}

double L2_old(double M2a, double M2b, double M2d) {
  double res[1], err[1];

  //if (sqrt(M2d)>sqrt(M2a)+sqrt(M2b)) {
  double complex lam_ABD = csqrt( lam(M2a,M2b,M2d) );

    double x_plus  = fmin( fmax( .5*(M2b+M2d-M2a+creal(lam_ABD))/M2d, 0. ), 1.),
           x_minus = fmin( fmax( .5*(M2b+M2d-M2a-creal(lam_ABD))/M2d, 0. ), 1.);

    int integrand(unsigned dim,  const double *x, void *data_,
                  unsigned fdim, double *val) {
      double _X_ = x[0]*x_minus,
             _Y_ = (1.-x[0])*x_minus + x[0]*x_plus,
             _Z_ = (1.-x[0])*x_plus  + x[0];
      double log1 = log(fabs( SQR(mubar)/(M2a*_X_+M2b*(1.-_X_)-M2d*_X_*(1.-_X_)) )),
             log2 = log(fabs( SQR(mubar)/(M2a*_Y_+M2b*(1.-_Y_)-M2d*_Y_*(1.-_Y_)) )),
             log3 = log(fabs( SQR(mubar)/(M2a*_Z_+M2b*(1.-_Z_)-M2d*_Z_*(1.-_Z_)) ));

      val[0] = x_minus*log1 + (x_plus-x_minus)*log2 + (1.-x_plus)*log3;
      return 0;
    }//*/
  /*} else {
    int integrand(unsigned dim,  const double *x, void *data_,
                  unsigned fdim, double *val) {
      double _X_ = x[0];
      val[0] = log(fabs( SQR(mubar)/(M2a*_X_+M2b*(1.-_X_)-M2d*_X_*(1.-_X_)) ));
      return 0;
    }//*/
  //}

  double xl[1] = { 0. };
  double xu[1] = { 1. };

  hcubature(1, integrand, NULL, 1, xl, xu, MaxEvls, 1e-10, 0, ERROR_INDIVIDUAL, res, err);
  return -res[0];
  /*  double complex x1 = .5*(M2b+M2d-M2a+creal(lam_ABD))/M2d,
                   x2 = .5*(M2b+M2d-M2a-creal(lam_ABD))/M2d;
  return - log(fabs(SQR(mubar)/M2d)) 
    + creal( clog(1.-x1)-x1*clog((x1-1.)/x1)-1.)
    + creal( clog(1.-x2)-x2*clog((x2-1.)/x2)-1.);//*/
}

double L2(double M2a, double M2b, double M2d) {
  double res[1], err[1];

  double complex lam_ABD = csqrt( lam(M2a,M2b,M2d) );
  double complex x1 = .5*(M2b+M2d-M2a+lam_ABD)/M2d,
                 x2 = .5*(M2b+M2d-M2a-lam_ABD)/M2d;

  return - log(fabs(SQR(mubar)/M2d))
         + creal( clog(1.-x1)-x1*clog((x1-1.)/x1)-1.)
         + creal( clog(1.-x2)-x2*clog((x2-1.)/x2)-1.);//*/
}


double L3_new(double K2,
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

  hcubature(1, integrand, NULL, 1, xl, xu, MaxEvls, 1e-10, 0, ERROR_INDIVIDUAL, res, err);
  return res[0];
}//*/


double L3(double K2,
          double M2a, double M2b, double M2c, double M2d, double M2e) {
  double res[1], err[1];
  double a = M2e, b = M2d, c = K2-M2d-M2e, d = M2b-M2c-M2e, e = M2a-M2b+M2e-K2, f = M2c;
  double complex lam_DE = csqrt( lam(M2d,M2e,K2) );
  //printf(" lam = %g + I %g \n", lam_DE );

  double complex al,
                 al1 = .5*( M2d+M2e-K2 + lam_DE )/M2d ,
                 al2 = .5*( M2d+M2e-K2 - lam_DE )/M2d ;
  //printf(" al1 = %g + I %g \n", al );
  //printf(" al2 = %g + I %g \n", al2 );
  if ( cabs(al1) > cabs(al2) ) { al = al1; }
  else { al = al2; }

  double complex y0 = - (d+e*al)/(c+2.*b*al),
                 y1 = y0 + al,
                 y2 = y0/(1.-al),
                 y3 = -y0/al;

  //printf(" y0 = %g + I %g \n", y0 );
  //printf(" y1 = %g + I %g \n", y1 );
  //printf(" y2 = %g + I %g \n", y2 );
  //printf(" y3 = %g + I %g \n", y3 );

  double complex tmp = + S3(y1,b,c+e,a+d+f) - S3(y2,a+b+c,d+e,f) + S3(y3,a,d,f);
  //printf(" res = %g \n", tmp/(c+2.*b*al) );
  return creal( tmp/(c+2.*b*al) );
}//*/



double Bubble_0(o,k,A_,B_,C_,D_,c)
  double o,k;
  void *A_, *B_, *C_, *D_;
  double *c; // Phi(Pa,Pb) = ca*Pa + cb*Pc + cK*K
{
  double res[2], err[2];

  double K2 = SQR(o)-SQR(k);
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0],
         mC = ((double *)C_)[0], mD = ((double *)D_)[0];
  double uC = ((double *)C_)[1], uD = ((double *)D_)[1];
  int    sC = ((double *)C_)[2], sD = ((double *)D_)[2];

  double L2ab = L2(SQR(mA),SQR(mB),SQR(mD)),
         L1a  = L1(SQR(mA)),
         L1b  = L1(SQR(mB));

  double ca = ((double *)c)[0];
  double cb = ((double *)c)[1];
  double cK = ((double *)c)[2];

  double complex lam_DC = csqrt( lam(K2,SQR(mD),SQR(mC)) );

  if (fabs(cimag(lam_DC))>0.) { return 0.; }

  double complex eC_p  = .5*( o*(K2+SQR(mC)-SQR(mD)) + k*(lam_DC) )/K2 ,
                 eC_m  = .5*( o*(K2+SQR(mC)-SQR(mD)) - k*(lam_DC) )/K2 ;

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
      E_PD = .5*( K2 + SQR(mD) - SQR(mC) );
      E_K  = K2;
    }
    if (E=='U') {
      E_PD = eD;
      E_K  = o;
    }
    temp *= E_PD; temp += cK*E_K*L2ab;

    //double complex thermal_weight = 1. + n(sC,eC-uC) + n(sD,eD-uD) ;
    double complex thermal_weight = calN(sC,sD,eC-uC,eD-uD);

    double complex jacobian =  ( eC_p - eC_m )              // from X = [0,1]
                              *SGN(creal(eC*eD))           ;//

    double prefactor = -.25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(temp);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }//*/

  double xl[1] = { 0. };
  double xu[1] = { 1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, TolVirt, 0, ERROR_INDIVIDUAL, res, err);
  //printf(" res = %g + I %g    ... err = %g + I %g \n", res[0], res[1], err[0], err[1] );
  return res[0];
}


double Triangle_0(o,k,A_,B_,C_,D_,E_,c)
  double o,k;
  void *A_, *B_, *C_, *D_, *E_;
  double *c; // Phi(Pa,Pb,Pc) = ca*Pa + cb*Pb + cc*Pc
{
  double res[2], err[2];

  double K2 = SQR(o)-SQR(k);
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0],
         mC = ((double *)C_)[0], mD = ((double *)D_)[0], mE = ((double *)E_)[0];
  double uD = ((double *)D_)[1], uE = ((double *)E_)[1];
  int    sD = ((double *)D_)[2], sE = ((double *)E_)[2];

  double L2bc    = L2(SQR(mB),SQR(mC),SQR(mE)),
         L2ac    = L2(SQR(mA),SQR(mC),K2),
         L2ab    = L2(SQR(mA),SQR(mB),SQR(mD)),
         L3abc   = L3(K2,SQR(mA),SQR(mB),SQR(mC),SQR(mD),SQR(mE));

  double Ib = (SQR(mD)+SQR(mA)-SQR(mB))*L3abc + L2bc - L2ac,
         Ic = (K2 +SQR(mA)-SQR(mC))*L3abc + L2bc - L2ab;

  double al_b = - 2.*K2*Ib + (SQR(mD)-SQR(mE)+K2)*Ic,
         al_c = + (SQR(mD)-SQR(mE)+K2)*Ib - 2.*SQR(mD)*Ic;

  double ca = ((double *)c)[0];
  double cb = ((double *)c)[1];
  double cc = ((double *)c)[2];

  double complex lam_DE = csqrt( lam(K2,SQR(mD),SQR(mE)) );

  if (fabs(cimag(lam_DE))>0.) { return 0.; }

  al_b *= 1./SQR(lam_DE);
  al_c *= 1./SQR(lam_DE);

  double complex eD_p  = .5*( o*(K2+SQR(mD)-SQR(mE)) + k*(lam_DE) )/K2 ;
  double complex eD_m  = .5*( o*(K2+SQR(mD)-SQR(mE)) - k*(lam_DE) )/K2 ;

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];// _Z_ = x[2]; // integration variables

    double complex eD = (1-_X_)*eD_m + _X_*eD_p , eE = o - eD,
                   pD = csqrt( SQR(eD) - SQR(mD) ),
                   pE = csqrt( SQR(eE) - SQR(mE) );

    double E_PD, E_K;
    if (E=='K') {
      E_PD = .5*( K2 + SQR(mD) - SQR(mE) );
      E_K = K2;
    }
    if (E=='U') {
      E_PD = eD;
      E_K  = o;
    }
    double complex temp = al_b*E_PD+al_c*E_K;
    temp *= (ca-cb-cc);
    temp += (cb*E_PD+cc*E_K)*L3abc;

    //double complex thermal_weight = 1. + n(sD,eD-uD) + n(sE,eE-uE) ;
    double complex thermal_weight = calN(sD,sE,eD-uD,eE-uE);

    double complex jacobian =  ( eD_p - eD_m )              // from X = [0,1]
                              *SGN(creal(eD*eE))           ;//

    double prefactor = -.25*pow(OOFP,3.)/k;

    double complex _outer = (prefactor)*(thermal_weight)*(jacobian)*(temp);

    val[0] = creal(_outer); val[1] = cimag(_outer);
    //printf(" res = %g + I %g\n", val[0], val[1] );
    return 0;
  }

  double xl[2] = { 0. };
  double xu[2] = { 1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, TolVirt, 0, ERROR_INDIVIDUAL, res, err);
  return res[0];
}
