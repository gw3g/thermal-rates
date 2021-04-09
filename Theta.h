

double _id(double o, double k, // K = (\omega,k)
           double complex s, double complex q0, double complex e,
           void *M_)
{
  return 1.;
}

/*--------------------------------------------------------------------*/
// eq. (B.1)

double _B1_i(double o, double k,
              double complex s12, double complex q0, double complex e2,
              void *M_)
{
  double m1=((double*)M_)[0], m2=((double*)M_)[1], m3=((double*)M_)[2];
  //     % = m_\ell           % = m_Q              % = m_S

  double complex  e1  = q0 - e2, e3 = o - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s12 ), res;

  double ms = m3, ml = m1; // TODO ?

  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );
  double complex p3 = csqrt( SQR(e3) - SQR(m3) );

  double complex cos_qp2 = ( q0*e2 + .5*(SQR(m1)- SQR(m2) - s12 ) )/(q*p2), 
                 sin_qp2 = csqrt( 1. - SQR(cos_qp2) );
  double complex cos_qk  = ( q0*o + .5*( SQR(m3) - K2      - s12 ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*p2*cos_qk*cos_qp2), b = creal(k*p2*sin_qk*sin_qp2);
  double ab[2] = {a,b};
  double complex F_ = calF(o*e2-.5*(SQR(ms)+s12-SQR(m1)-SQR(m3)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = -.5*(  K2 - SQR(m3) + s12 + 2.*(a-o*e2) + (K2-SQR(m3))*(1.
              - (K2-.5*SQR(m3)+.5*s12-.5*SQR(ms)+.5*SQR(m1))*creal(F_) ) 
              )/(s12-SQR(ml));
  } else
  if (E=='U') {
    res = -.5*( 2.*e1 - (K2-SQR(m3))*(e2+2.*e1)*creal(F_)
              )/(s12-SQR(ml));
  }
  return creal(res);
}

double _B1_ii(double o, double k,
              double complex s12, double complex q0, double complex e2,
              void *M_)
{
  double m1=((double*)M_)[0], m2=((double*)M_)[1], m3=((double*)M_)[2];
  //     % = m_S              % = m_Q              % = m_\ell

  double complex  e1  = q0 - e2, e3 = o - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s12 ), res;

  double ms = m1, ml = m3; // TODO ?

  // Here are the explicit integrands:
  if (E=='K') {
    res = -.5*(  3.*K2 + SQR(m3) - s12 - 2.*SQR(ms)
              + 2.*SQR(ms)*( K2 + SQR(m3) - SQR(ms) )/(s12-SQR(ms)) )/(s12-SQR(ms));
  } else
  if (E=='U') {
    res = ( - (e3+o) - 2.*SQR(ms)*e3/(s12-SQR(ms)) )/(s12-SQR(ms));
  }
  return creal(res);
}

/*--------------------------------------------------------------------*/
// eq. (4.5) 2 -> 2 reactions (t-channel)

double _45_i(double o, double k,
             double complex t, double complex q0, double complex e1,
             void *M_)
{
  double M1=((double*)M_)[0], m1=((double*)M_)[1], m2=((double*)M_)[2];
  //     % = m_Q              % = m_l              % = m_S

  double complex  E1  = e1 - q0, e2 = o - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = m2, ml = m1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );

  double complex cos_qp1 = ( q0*e1 + .5*(SQR(M1)- SQR(m1) - t ) )/(q*p1), 
                 sin_qp1 = csqrt( 1. - SQR(cos_qp1) );
  double complex cos_qk  = ( q0*o + .5*( SQR(m2) - K2 - t ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*p1*cos_qk*cos_qp1), b = creal(k*p1*sin_qk*sin_qp1);
  double ab[2] = {a,b};
  double complex F_ = calF(o*e1 + .5*(SQR(ms)-K2-SQR(m1)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( 2.*(a-o*e1) + ( SQR(ms) - K2 )
              - .5*( SQR(ms) - K2 )*( SQR(m2) - 2.*K2 - SQR(m1) - t + SQR(ms) )*creal(F_) 
              )/(t-SQR(ml));
  } else
  if (E=='U') {
    res = ( - e1 - .5*( SQR(ms)-K2 )*( E1-2.*e1 )*creal(F_) )/(t-SQR(ml));
  }
  return creal(res);
}

double _45_ii(double o, double k,
              double complex t, double complex q0, double complex e1,
              void *M_)
{
  double M1=((double*)M_)[0], m1=((double*)M_)[1], m2=((double*)M_)[2];
  //     % = m_l              % = m_Q              % = m_S

  double complex  E1  = e1 - q0, e2 = o - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = m2, ml = M1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );

  double complex cos_qp1 = ( q0*e1 + .5*(SQR(M1)- SQR(m1) - t ) )/(q*p1), 
                 sin_qp1 = csqrt( 1. - SQR(cos_qp1) );
  double complex cos_qk  = ( q0*o + .5*( SQR(m2) - K2 - t ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*p1*cos_qk*cos_qp1), b = creal(k*p1*sin_qk*sin_qp1);
  double ab[2] = {a,b};
  double complex F_ = calF(o*e1 + .5*(SQR(m2)+SQR(M1)-t-SQR(ms)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( - K2 + SQR(m2) -t - 2.*(a-o*e1) - ( K2 - SQR(ms) )
              + .5*( K2 - SQR(ms) )*( 2.*K2 - SQR(m2) + t + SQR(M1) - SQR(ms) )*creal(F_) 
              )/(t-SQR(ml));
  } else
  if (E=='U') {
    res = ( + E1 + .5*( K2-SQR(ms) )*( e1-2.*E1 )*creal(F_) )/(t-SQR(ml));
  }
  return creal(res);
}

double _45_iii(double o, double k,
               double complex t, double complex q0, double complex e1,
               void *M_)
{
  double M1=((double*)M_)[0], m1=((double*)M_)[1], m2=((double*)M_)[2];
  //     % = m_Q/S            % = m_S/Q            % = m_l

  double complex  E1  = e1 - q0, e2 = o - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = fmax(M1,m1), ml = m2; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( - 2.*K2 - SQR(m2) + t + 2.*SQR(ms)
                - 2.*SQR(ms)*( K2 + SQR(m2) - SQR(ms) )/(t-SQR(ms)) )/(t-SQR(ms));
  } else
  if (E=='U') {
    res = ( - e2 - o - 2.*SQR(ms)*e2/(t-SQR(ms)) )/(t-SQR(ms));
  }
  return creal(res);
}

// eq. (B.9) 2 -> 2 reactions (s-channel)

double _B9_i(double o, double k,
             double complex s, double complex q0, double complex e2,
             void *M_)
{
  double M1=((double*)M_)[0], m1=((double*)M_)[1], m2=((double*)M_)[2];
  //     % = m_S              % = m_Q              % = m_l

  double complex  e1  = q0 - e2, E1 = q0 - o, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s ), res;

  double ms = M1, ml = m2; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );

  double complex cos_qp2 = ( q0*e2 + .5*(SQR(m1)- SQR(m2) - s ) )/(q*p2), 
                 sin_qp2 = csqrt( 1. - SQR(cos_qp2) );
  double complex cos_qk  = ( q0*o + .5*( SQR(M1) - K2 - s ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*p2*cos_qk*cos_qp2), b = creal(k*p2*sin_qk*sin_qp2);
  double ab[2] = {a,b};
  double complex F_ = calF(o*e2+.5*(SQR(ms)-K2-SQR(m2)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( + 2.*(a-o*e2) - ( K2 - SQR(ms) )
                - .5*( K2 - SQR(ms) )*( 2.*K2 - SQR(M1) + s + SQR(m2) - SQR(ms) )*creal(F_) 
              )/(s-SQR(ml));
  } else
  if (E=='U') {
    res = ( - e2 - .5*( K2-SQR(ms) )*( e1+2.*e2 )*creal(F_) )/(s-SQR(ml));
  }
  return creal(res);
}

double _B9_ii(double o, double k,
              double complex s, double complex q0, double complex e2,
              void *M_)
{
  double M1=((double*)M_)[0], m1=((double*)M_)[1], m2=((double*)M_)[2];
  //     % = m_l              % = m_Q              % = m_S

  double complex  e1  = q0 - e2, E1 = q0 - o, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s ), res;

  double ms = m2, ml = M1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );
  double complex p2 = csqrt( SQR(e2) - SQR(m2) );

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( + s - 2.*K2 - SQR(M1)
                + 2.*SQR(ms)*(s-K2-SQR(M1))/(s-SQR(ms))
              )/(s-SQR(ms));
  } else
  if (E=='U') {
    res = ( E1 - o + 2.*SQR(ms)*E1/(s-SQR(ms)) )/(s-SQR(ms));
  }
  return creal(res);
}


/*--------------------------------------------------------------------*/
// eq. (B.19)

double _B19_i(double o, double k,
              double complex t, double complex q0, double complex E1,
              void *M_)
{
  double M1=((double*)M_)[0], M2=((double*)M_)[1], m1=((double*)M_)[2];
  //     % = m_l              % = m_S              % = m_Q

  double complex  E2  = - q0 - o, e1 = E1 - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = M2, ml = M1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex k2 = csqrt( SQR(E2) - SQR(M2) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );

  double complex cos_qk1 = ( q0*E1 + .5*(SQR(m1)- SQR(M1) - t ) )/(q*k1), 
                 sin_qk1 = csqrt( 1. - SQR(cos_qk1) );
  double complex cos_qk  = ( q0*o + .5*( t + K2 - SQR(M2) ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*k1*cos_qk*cos_qk1), b = creal(k*k1*sin_qk*sin_qk1);
  double ab[2] = {a,b};
  double complex F_ = calF(o*E1+.5*(K2+SQR(M1)-SQR(ms)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( - 2.*(a-o*E1) + ( K2-SQR(ms) )
              + .5*( K2-SQR(ms) )*( t - SQR(M2) - SQR(M1) + SQR(ms) )*creal(F_) 
              )/(t-SQR(ml));
  } else
  if (E=='U') {
    res = ( E1 + .5*( K2-SQR(ms) )*( e1-2.*E1 )*creal(F_) )/(t-SQR(ml));
  }
  return creal(res);
}

double _B19_ii(double o, double k,
               double complex t, double complex q0, double complex E1,
               void *M_)
{
  double M1=((double*)M_)[0], M2=((double*)M_)[1], m1=((double*)M_)[2];
  //     % = m_Q              % = m_S              % = m_l

  double complex  E2  = - q0 - o, e1 = E1 - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = M2, ml = m1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex k2 = csqrt( SQR(E2) - SQR(M2) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );

  double complex cos_qk1 = ( q0*E1 + .5*(SQR(m1)- SQR(M1) - t ) )/(q*k1), 
                 sin_qk1 = csqrt( 1. - SQR(cos_qk1) );
  double complex cos_qk  = ( q0*o + .5*( t + K2 - SQR(M2) ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*k1*cos_qk*cos_qk1), b = creal(k*k1*sin_qk*sin_qk1);
  double ab[2] = {a,b};
  double complex F_ = calF(o*E1+.5*(t+SQR(ms)-SQR(m1)-SQR(M2)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( SQR(M2) - K2 - t + 2.*(a-o*E1) + ( SQR(ms) - K2 )
              - .5*( SQR(ms) - K2 )*( - t - 2.*K2 + SQR(M2) - SQR(m1) + SQR(ms) )*creal(F_) 
              )/(t-SQR(ml));
  } else
  if (E=='U') {
    res = ( - e1 - .5*( SQR(ms) - K2 )*( E1-2.*e1 )*creal(F_) )/(t-SQR(ml));
  }
  return creal(res);
}

double _B19_iii(double o, double k,
                double complex t, double complex q0, double complex E1,
                void *M_)
{
  double M1=((double*)M_)[0], M2=((double*)M_)[1], m1=((double*)M_)[2];
  //     % = m_S/Q            % = m_l              % = m_Q/S

  double complex  E2  = - q0 - o, e1 = E1 - q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - t ), res;

  double ms = fmax(M1,m1), ml = M2; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex k2 = csqrt( SQR(E2) - SQR(M2) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );

  // Here are the explicit integrands:
  if (E=='K') {
    res = +.5*( t -2.*SQR(M2) - 2.*K2 + 2.*SQR(ms) 
              - 2.*SQR(ms)*(K2+SQR(M2)-SQR(ms))/(t-SQR(ms))
              )/(t-SQR(ms));
  } else
  if (E=='U') {
    res = ( (E2-o) + 2.*SQR(ms)*E2/(t-SQR(ms)) )/(t-SQR(ms));
  }
  return creal(res);
}



/*--------------------------------------------------------------------*/
// eq. (B.27)

double _B27_i(double o, double k,
              double complex s, double complex q0, double complex E2,
              void *M_)
{
  double M1=((double*)M_)[0], M2=((double*)M_)[1], m1=((double*)M_)[2];
  //     % = m_Q              % = m_l              % = m_S

  double complex  E1  = q0 - E2, e1 = o + q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s ), res;

  double ms = m1, ml = M2; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex k2 = csqrt( SQR(E2) - SQR(M2) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );

  double complex cos_qk2 = ( q0*E2 + .5*(SQR(M1)- SQR(M2) - s ) )/(q*k2), 
                 sin_qk2 = csqrt( 1. - SQR(cos_qk2) );
  double complex cos_qk  = ( q0*o + .5*( s + K2 - SQR(m1) ) )/(q*k),
                 sin_qk  = csqrt( 1. - SQR(cos_qk ) );

  double a = creal(k*k2*cos_qk*cos_qk2), b = creal(k*k2*sin_qk*sin_qk2);
  double ab[2] = {a,b};
  double complex F_ = calF(o*E2+.5*(K2+SQR(M2)-SQR(ms)), ab) ;

  // Here are the explicit integrands:
  if (E=='K') {
    res = -.5*( 2.*(a-o*E2) - (SQR(ms)-K2)
              + .5*( SQR(ms)-K2 )*( 2.*K2 + SQR(M2) - SQR(m1) - SQR(ms) + s )*creal(F_) 
              )/(s-SQR(ml));
  } else
  if (E=='U') {
    res = ( E2 + .5*( SQR(ms) - K2 )*( E1 + 2.*E2 )*creal(F_) )/(s-SQR(ml));
  }
  return creal(res);
}

double _B27_ii(double o, double k,
               double complex s, double complex q0, double complex E2,
               void *M_)
{
  double M1=((double*)M_)[0], M2=((double*)M_)[1], m1=((double*)M_)[2];
  //     % = m_Q              % = m_S              % = m_l

  double complex  E1  = q0 - E2, e1 = o + q0, // resolve energies
                  K2 = SQR(o) - SQR(k), q = csqrt( SQR(q0) - s ), res;

  double ms = M2, ml = m1; // TODO ?

  double complex k1 = csqrt( SQR(E1) - SQR(M1) );
  double complex k2 = csqrt( SQR(E2) - SQR(M2) );
  double complex p1 = csqrt( SQR(e1) - SQR(m1) );

  // Here are the explicit integrands:
  if (E=='K') {
    res = -.5*( 
        3.*K2 + SQR(m1) - s - 2.*SQR(ms) + 2.*SQR(ms)*(K2+SQR(m1)-SQR(ms))/(s-SQR(ms))
              )/(s-SQR(ms));
  } else
  if (E=='U') {
    res = ( - (e1+o) - 2.*SQR(ms)*e1/(s-SQR(ms)) )/(s-SQR(ms));
  }
  return creal(res);
}

/*--------------------------------------------------------------------*/
// eq. (4.48)

double _448_(double o, double k, double complex eA, void *M_)
{
  double mA=((double*)M_)[0], mB=((double*)M_)[1];
  //     % = m_l              % = m_S

  double complex  eB  = o - eA,
                  K2 = SQR(o) - SQR(k), res;

  // Here are the explicit integrands:
  if (E=='K') {
    res = - ( SQR(mB) - SQR(mA) - K2);
  } else
  if (E=='U') {
    res = 2.*eA;
  }
  return .5*creal(res)*SGN(creal(eA*eB));
}

/*--------------------------------------------------------------------*/
// Definitions of angular ave. for virtual integrands:
//

void _virt_i(double o, double k,
             double *PA, double *PB, double *PC, double *PD,
             double complex *_integrands) {

  double eA = PA[0], pA = PA[1], mA = PA[2],
         eB = PB[0], pB = PB[1], mB = PB[2],
         eC = PC[0], pC = PC[1], mC = PC[2],
         eD = PD[0], pD = PD[1], mD = PD[2];

  double K2 = SQR(o)-SQR(k);
  double k_pD = ( o*eD + .5*(SQR(mC)-SQR(mD)-K2) );

  int nn; double coeffs[4][10];
  if (E=='K') {
    nn = 1;
    coeffs[0][0] = o*eA ; coeffs[0][1] = -k_pD/SQR(pD);
    coeffs[1][0] = -o*eA; coeffs[1][1] = +k_pD/SQR(pD);
    coeffs[2][0] = o*eD-k_pD-o*eB; coeffs[2][1] = +k_pD/SQR(pD);
    coeffs[3][0] = o*eD-k_pD+o*eB; coeffs[3][1] = -k_pD/SQR(pD);
  }
  if (E=='U') {
    nn = 0;
    coeffs[0][0] = eA;
    coeffs[1][0] = -eA;
    coeffs[2][0] = eD-eB;
    coeffs[3][0] = eD+eB;
  }
  double p_AD = pA*pD, p_BD = pB*pD;
  double e_AD = eA*eD, e_BD = eB*eD;
  double m_BAD = .5*( SQR(mB) - SQR(mA) - SQR(mD) );
  double m_ABD = .5*( SQR(mA) - SQR(mB) - SQR(mD) );

  _integrands[0] = - pA*calG(e_AD+m_BAD, p_AD, nn, coeffs[0]);
  _integrands[1] = - pA*calG(e_AD-m_BAD, p_AD, nn, coeffs[1]);
  _integrands[2] = - pB*calG(e_BD+m_ABD, p_BD, nn, coeffs[2]);
  _integrands[3] = - pB*calG(e_BD-m_ABD, p_BD, nn, coeffs[3]);

  //printf("pA*calG = %g\n", creal(_integrands[0]) );
}

void _virt_ii(double o, double k,
             double *PA, double *PB, double *PC, double *PD,
             double complex *_integrands) {

  double eA = PA[0], pA = PA[1], mA = PA[2],
         eB = PB[0], pB = PB[1], mB = PB[2],
         eC = PC[0], pC = PC[1], mC = PC[2],
         eD = PD[0], pD = PD[1], mD = PD[2];

  double K2 = SQR(o)-SQR(k);
  double k_pD = ( o*eD + .5*(SQR(mC)-SQR(mD)-K2) );

  int nn; double coeffs[4][10];
  if (E=='K') {
    nn = 0;
    coeffs[0][0] = 2.*K2-o*eD+k_pD;
    coeffs[1][0] = 2.*K2-o*eD+k_pD;
    coeffs[2][0] = 2.*K2-o*eD+k_pD;
    coeffs[3][0] = 2.*K2-o*eD+k_pD;
  }
  if (E=='U') {
    nn = 0;
    coeffs[0][0] = eC+o;
    coeffs[1][0] = eC+o;
    coeffs[2][0] = eC+o;
    coeffs[3][0] = eC+o;
  }
  double p_AD = pA*pD, p_BD = pB*pD;
  double e_AD = eA*eD, e_BD = eB*eD;
  double m_BAD = .5*( SQR(mB) - SQR(mA) - SQR(mD) );
  double m_ABD = .5*( SQR(mA) - SQR(mB) - SQR(mD) );

  _integrands[0] = - pA*calG(e_AD+m_BAD, p_AD, nn, coeffs[0]);
  _integrands[1] = - pA*calG(e_AD-m_BAD, p_AD, nn, coeffs[1]);
  _integrands[2] = - pB*calG(e_BD+m_ABD, p_BD, nn, coeffs[2]);
  _integrands[3] = - pB*calG(e_BD-m_ABD, p_BD, nn, coeffs[3]);

  //printf("pA*calG = %g\n", creal(_integrands[0]) );
}

void _virt_iii(double o, double k,
             double *PA, double *PB, double *PC, double *PD,
             double complex *_integrands) {

  double eA = PA[0], pA = PA[1], mA = PA[2],
         eB = PB[0], pB = PB[1], mB = PB[2],
         eC = PC[0], pC = PC[1], mC = PC[2],
         eD = PD[0], pD = PD[1], mD = PD[2];

  double K2 = SQR(o)-SQR(k);
  double k_pD = ( o*eD + .5*(SQR(mC)-SQR(mD)-K2) );

  int nn; double coeffs[4][10];
  if (E=='K') {
    nn = 0;
    coeffs[0][0] = K2-o*eD+k_pD;
    coeffs[1][0] = K2-o*eD+k_pD;
    coeffs[2][0] = K2-o*eD+k_pD;
    coeffs[3][0] = K2-o*eD+k_pD;
  }
  if (E=='U') {
    nn = 0;
    coeffs[0][0] = eC;
    coeffs[1][0] = eC;
    coeffs[2][0] = eC;
    coeffs[3][0] = eC;
  }
  double p_AD = pA*pD, p_BD = pB*pD;
  double e_AD = eA*eD, e_BD = eB*eD;
  double m_BAD = .5*( SQR(mB) - SQR(mA) - SQR(mD) );
  double m_ABD = .5*( SQR(mA) - SQR(mB) - SQR(mD) );

  _integrands[0] = - pA*2.*calG(e_AD+m_BAD, p_AD, nn, coeffs[0]);
  _integrands[1] = - pA*2.*calG(e_AD-m_BAD, p_AD, nn, coeffs[1]);
  _integrands[2] = - pB*2.*calG(e_BD+m_ABD, p_BD, nn, coeffs[2]);
  _integrands[3] = - pB*2.*calG(e_BD-m_ABD, p_BD, nn, coeffs[3]);

  //printf("pA*calG = %g\n", creal(_integrands[0]) );
}

void _virt_iv(double o, double k,
              double *PA, double *PB, double *PC, double *PD, double *PE,
              double complex *_integrands) {

  double eA = PA[0], pA = PA[1], mA = PA[2],
         eB = PB[0], pB = PB[1], mB = PB[2],
         eC = PC[0], pC = PC[1], mC = PC[2],
         eD = PD[0], pD = PD[1], mD = PD[2],
         eE = PE[0], pE = PE[1], mE = PE[2];

  double K2 = SQR(o)-SQR(k);
  double K_PD = ( -.5*(SQR(mE)-SQR(mD)-K2) );
  double K_PE = ( -.5*(SQR(mD)-SQR(mE)-K2) );

  int nn; double coeffs[6][10];
  if (E=='K') {
    nn = 1;
    coeffs[0][0] = +K_PD-3.*o*eA,       coeffs[0][1] = +3.; // for eA-integral
    coeffs[1][0] = +K_PD+3.*o*eA,       coeffs[1][1] = -3.;
    coeffs[2][0] = -2.*K_PD+3.*o*eB,    coeffs[2][1] = -3.; // for eB-integral
    coeffs[3][0] = -2.*K_PD-3.*o*eB,    coeffs[3][1] = +3.;
    coeffs[4][0] = -K_PE-2.*K2+3.*o*eC, coeffs[4][1] = -3.; // for eC-integral
    coeffs[5][0] = -K_PE-2.*K2-3.*o*eC, coeffs[4][1] = +3.;
  }
  if (E=='U') {
    nn = 0;
    coeffs[0][0] = +eD-3.*eA; // for eA-integral
    coeffs[1][0] = +eD+3.*eA;
    coeffs[2][0] = +3.*eB-2.*eD; // for eB-integral
    coeffs[3][0] = -3.*eB-2.*eD;
    coeffs[4][0] = +3.*eC-eE-2.*o; // for eC-integral
    coeffs[5][0] = -3.*eC-eE-2.*o;
  }
  double p_AD = pA*pD, p_BD = pB*pD, p_BE = pB*pE, p_CE = pC*pE;
  double e_AD = eA*eD, e_BD = eB*eD, e_BE = eB*eE, e_CE = eC*eE;

  double m_BAD = .5*( SQR(mB) - SQR(mA) - SQR(mD) );
  double m_ABD = .5*( SQR(mA) - SQR(mB) - SQR(mD) );

  double m_BCE = .5*( SQR(mB) - SQR(mC) - SQR(mE) );
  double m_CBE = .5*( SQR(mC) - SQR(mB) - SQR(mE) );

  double m_CAM = .5*( SQR(mC) - SQR(mA) - K2 );
  double m_ACM = .5*( SQR(mA) - SQR(mC) - K2 );

  double m_AECD = .5*( SQR(mA) + SQR(mE) - SQR(mC) - SQR(mD) );

  _integrands[0] = - pA*calH(e_AD+m_BAD, o*eA+m_CAM,  pA, pD, pE, k, nn, coeffs[0]);
  _integrands[1] = - pA*calH(e_AD-m_BAD, o*eA-m_CAM,  pA, pD, pE, k, nn, coeffs[1]);

  _integrands[2] = - pB*calH(e_BD+m_ABD, o*eB+m_AECD, pB, pD, pE, k, nn, coeffs[2]);
  _integrands[3] = - pB*calH(e_BD-m_ABD, o*eB-m_AECD, pB, pD, pE, k, nn, coeffs[3]);

  _integrands[4] = - pB*calH(e_BE-m_CBE, o*eB+m_AECD, pB, pE, pD, k, nn, coeffs[2]);
  _integrands[5] = - pB*calH(e_BE+m_CBE, o*eB-m_AECD, pB, pE, pD, k, nn, coeffs[3]);

  _integrands[6] = - pC*calH(e_CE+m_BCE, o*eC+m_ACM,  pC, pE, pD, k, nn, coeffs[4]);
  _integrands[7] = - pC*calH(e_CE-m_BCE, o*eC-m_ACM,  pC, pE, pD, k, nn, coeffs[5]);

  //printf("pA*calG = %g\n", creal(_integrands[0]) );
}


