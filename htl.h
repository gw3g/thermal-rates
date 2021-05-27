/*
 *  code for HTL corrections
 *
 */
double    TolHTL=1e-5;

double complex Q0(double complex x) { return .5*clog( (1.+x)/(1.-x) ); }

double rho_0(double complex q0, double complex q, double m) {
  double complex LL = Q0(q/q0)/q, res;
  res  = 1.-.5*SQR(m)*LL/q0;
  res *= 1./( SQR(q0)-SQR(q)-SQR(m)
            + .25*SQR(SQR(m))*(SQR(q*LL)-SQR(1.-q0*LL))/SQR(q) );
  return cimag(res);
}

double rho_s(double complex q0, double complex q, double m) {
  double complex LL = Q0(q/q0)/q, res;
  res  = 1.+.5*SQR(m)*(1.-q0*LL)/SQR(q);
  res *= 1./( SQR(q0)-SQR(q)-SQR(m)
            + .25*SQR(SQR(m))*(SQR(q*LL)-SQR(1.-q0*LL))/SQR(q) );
  return cimag(res);
}


/*--------------------------------------------------------------------*/
// 1+n \to 2+n

double Rate_1_to_2_HTL(o,k,A_,B_,mlT)
  double o,k,mlT;
  void *A_, *B_;
{
  double res[2], err[2]; // reproducing eqs. (B.5)-(B.8)

  double K2 = SQR(o)-SQR(k) ;
  double mA = ((double *)A_)[0], mB = ((double *)B_)[0];
  double uA = ((double *)A_)[1], uB = ((double *)B_)[1];
  int    sA = ((double *)A_)[2], sB = ((double *)B_)[2];

  double complex lam_AB = csqrt( lam(K2,SQR(mA),SQR(mB)) );
  if (fabs(cimag(lam_AB))>0.) { return 0.; }

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0];
    double complex eA = .5*( o*(K2+SQR(mA)-SQR(mB)) + k*_X_*lam_AB )/K2, eB = o-eA;
    double complex pA = csqrt(SQR(eA)-SQR(mA)),
                   pB = csqrt(SQR(eB)-SQR(mB));
    //double complex thermal_weight =  ( 1. + n(sA,eA-uA) + n(sB,eB-uB) );
    double complex thermal_weight =  calN(sA,sB,eA-uA,eB-uB);
    double complex jacobian = .5*k*lam_AB/K2;
    double complex rate;
    double complex kpA = .5*(SQR(k)+SQR(pA)-SQR(pB));
    if (E=='K') {
      rate = (o-kpA*eA/SQR(pA))*log(fabs( (eA-pA)/(eA+pA) ))/pA - 2.*kpA/SQR(pA) ;
    }
    if (E=='U') {
      rate = log(fabs( (eA-pA)/(eA+pA) ))/pA ;
    }
    rate *= SQR(mlT);
    double prefactor = .25*OOFP/k;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate)*SGN(creal(eA*eB));
    //printf(" inner = %g \n", creal(_inner));

    val[0] = creal(_inner); val[1] = cimag(_inner);
    return 0;
  }

  double xl[1] = { -1. };
  double xu[1] = { +1. };

  hcubature(2, integrand, NULL, 1, xl, xu, MaxEvls, 0, TolHTL, ERROR_INDIVIDUAL, res, err);
  return -res[0];
}

/*--------------------------------------------------------------------*/
// 2 \to 2 w/ soft t-channel

double Rate_2_to_2_HTL(o,k,A_,B_,mlT) 
  double o, k, mlT;
  void *A_, *B_;// *C_;
{
  double res[2], err[2]; // reproducing eqs. (5.11)-(5.12)

  // TODO (!)
  if (o<k) { printf("Below LC not yet implemented!"); return 0.; }

  double M = sqrt( SQR(o)-SQR(k) );
  double ml = ((double *)A_)[0], mS = ((double *)B_)[0];
  double ul = ((double *)A_)[1], uS = ((double *)B_)[1];
  int    sl = ((double *)A_)[2], sS = ((double *)B_)[2];

  int integrand(unsigned dim,  const double *x, void *data_,
                unsigned fdim, double *val) {

    double _X_ = x[0], _Y_ = x[1]; // integration variables

    double complex t = SQR(M)*((_X_-1.)/(_X_)); // (-inf,0]
    //double complex t = -SQR(M)*1e3*_X_; // (-inf,0]

    double complex lam_S   = csqrt( lam(t,SQR(M ),SQR(mS)) );

    double complex q0  = .5*( o*(t+SQR(M)-SQR(mS) ) + k*_Y_*(lam_S) )/SQR(M); // [q0^+,q0^0]
    double complex q = csqrt( SQR(q0) - t );
    double complex qk = q0*o + .5*( SQR(mS) - SQR(M) - t);


    //double complex thermal_weight =  1. + n(sl,q0-ul) + n(sS,o-q0-uS)  ;
    double complex thermal_weight =  calN(sl,sS,q0-ul,o-q0-uS);
    //double complex LAM = ( -n(sl,ul)+n(sS,o-uS) )/thermal_weight;
    //double complex thermal_weight = ( -n(sl,ul)+n(sS,o-uS) );

    double complex jacobian = ( .5*k*lam_S/SQR(M) )        // ..   Y = [-1,1]
                              *SQR(M/_X_) ;               //      X = [0,1]
                              //*SQR(M)*1e3 ;               //      X = [0,1]

    double prefactor = 4.*pow(OOFP,2.)/k;
    double R0 = rho_0(q0+I*1e-9, q, mlT),
           RS = rho_s(q0+I*1e-9, q, mlT);

    double complex rate;
    if (E=='K') {
      rate  = + o*q0*( R0 - .25*M_PI*SQR(mlT)/(q0*q*(t-SQR(ml))) );
      rate += - qk*(RS - .25*M_PI*SQR(mlT)*q0/(SQR(q)*q*(t-SQR(ml))) );
    }
    if (E=='U') {
      rate  = + q0*( R0 - .25*M_PI*SQR(mlT)/(q0*q*(t-SQR(ml))) );
    }
    //rate *= LAM;

    double complex _inner = (prefactor)*(thermal_weight)*(jacobian)*(rate);
    //printf(" q0 = %g, t = %g, rate = %g \n", creal(q0), creal(t), creal(rate));

    val[0] = creal(_inner); val[1] = cimag(_inner);
    return 0;
  }

  double xl[3] = { 0., -1.};
  double xu[3] = { 1., +1.};

  hcubature(2, integrand, NULL, 2, xl, xu, MaxEvls, 0, TolHTL, ERROR_INDIVIDUAL, res, err);
  return -res[0];
}

/*--------------------------------------------------------------------*/

void Rate_fig3_scan(double M, double m_l, double m_Q, double m_S) {
  double g1 = 1./3., g2 = 2./3., mu_L = 1e-3, mu_Y = 2.*1e-2;

  //double m_l = .1, m_Q = .01, m_S = 1.;
  double mu_l = mu_L - .5*mu_Y, mu_Q = 0, mu_S = .5*mu_Y;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};

  double mlT = (1.+SQR(mu_l/M_PI))*(SQR(g1)+3.*SQR(g2))/(16.); mlT = sqrt(mlT);

  int N_k;
  double o, k, k_min, k_max, step;
  double res_2_to_2;
  double res_1_to_2;

  char *prefix=(char*)"out/htl_";
  char  suffix[50];
  char  filename[90];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"E=%c_{ml=%.2f,mQ=%.3f,mS=%.2f,M=%.2f}_LamS.dat",E,m_l,m_Q,m_S,M);
  strcat(filename,suffix);
  out=fopen(filename,"w");

  // Here are some parameters that can be changed:
  N_k=100; 
  k_min=0.01;
  k_max=15.;
  // don't change anything after that.

  k = k_min;
  step=(k_max-k_min)/((double) N_k-1);


  printf(" Settings: M=%g, with k_min=%g, k_max=%g\n",M,k_min,k_max); 
  double frac;

  for (int i=0; i<N_k; i++) {
    frac = (double)i/(double)(N_k-1);
    o = sqrt( SQR(M) + SQR(k) );

    res_2_to_2 = Rate_2_to_2_HTL(o,k,l_,S_,mlT);
    res_1_to_2 = Rate_1_to_2_HTL(o,k,l_,S_,mlT);

    printf(" k = %.5e , [%2.2f%]\n", k , 100.*frac); 
    fprintf( out, "%.10e    %.10e    %.10e\n", 
                   k, res_1_to_2, res_2_to_2);
    k += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

