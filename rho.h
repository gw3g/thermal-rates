/*
 *  master spectral functions, as a cross check
 *
 */

#include <gsl/gsl_sf_bessel.h>
double k2av(double M) {
  // worth noting: k0=3*k for k0~26.
  return 3.*M*gsl_sf_bessel_Kn(3,M)/gsl_sf_bessel_Kn(2,M);
}

/*--------------------------------------------------------------------*/

double rho_f(double o, double k, void *A_, void *B_, void *C_) {
  double K2 = SQR(o) - SQR(k);
  double res = 
    Rate_1_to_3(o,k,A_,B_,C_,_id)
   +Rate_2_to_2_sChan(o,k,A_,B_,C_,_id)
   +Rate_2_to_2_sChan(o,k,B_,A_,C_,_id)
   +Rate_2_to_2_sChan(o,k,C_,A_,B_,_id)
   +Rate_3_to_1_sChan(o,k,A_,B_,C_,_id)
   +Rate_3_to_1_sChan(o,k,B_,C_,A_,_id)
   +Rate_3_to_1_sChan(o,k,C_,A_,B_,_id);
  return ( res - pow(OOFP,3.)*K2/8.)*64.*M_PI;
}//*/

void rho_f_scan(double k) {
  double aa[3] = {0.0,0.,+1.},
         bb[3] = {0.0,0.,+1.},
         cc[3] = {0.0,0.,+1.};

  int N_o;
  double o, o_min, o_max, step;
  double res;

  char *prefix=(char*)"out/rho_f_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{k=%.1f}.dat",k);
  strcat(filename,suffix);
  out=fopen(filename,"w");

  N_o=200; 
  o_min=fmax(k*1.01,.01);
  o_max=100.;

  o = o_min;
  step=pow(o_max/o_min,1./(double)(N_o-1));

  printf(" Settings: k=%g, with o_min=%g, o_max=%g\n",k,o_min,o_max); 
  double frac;

  for (int i=0; i<N_o; i++) {
    frac = (double)i/(double)(N_o-1);

    res = rho_f(o,k,aa,bb,cc);

    printf(" omega = %.5e , [%2.2f%]\n", o , 100.*frac); 
    fprintf( out, "%.8e    %.8e\n", o, res );
    o *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}//*/


/*--------------------------------------------------------------------*/

double rho_j_real(double M, double k, double m_l, double m_Q, double m_S) {

  double o = sqrt( SQR(M) + SQR(k) );

  double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.};
  double Q_[3] = {m_Q, mu_Q, +1.};
  double S_[3] = {m_S, mu_S, +1.};

  double res =
    Rate_1_to_3(o,k,l_,Q_,S_,_j_1_3_s)
   +Rate_2_to_2_tChan(o,k,Q_,l_,S_,_j_2_2_t1)
   +Rate_2_to_2_tChan(o,k,l_,Q_,S_,_j_2_2_t2)
   +Rate_2_to_2_sChan(o,k,S_,Q_,l_,_j_2_2_s)
   +Rate_3_to_1_tChan(o,k,l_,S_,Q_,_j_3_1_t1)
   +Rate_3_to_1_tChan(o,k,Q_,S_,l_,_j_3_1_t2)
   +Rate_3_to_1_sChan(o,k,Q_,l_,S_,_j_3_1_s);
    ;
  return  res*32.*M_PI ;
}//*/

double rho_j_virt(double M, double k, double m_l, double m_Q, double m_S) {

  double o = sqrt( SQR(M) + SQR(k) );

  double mu_l = 0., mu_Q = 0., mu_S = 0.;
  double l_[3] = {m_l, mu_l, -1.},
         Q_[3] = {m_Q, mu_Q, +1.},
         S_[3] = {m_S, mu_S, +1.};

  double ccc[3];
  ccc[0] = 1.; ccc[1] = 0.; ccc[2] = 1.;
  E = 'K';

  // CHECK (1/2+n) -> 1/2
  //double temp1 = Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/SQR(M); //*/
  //double temp2 = Triangle_0_new(o,k,l_,Q_,S_,l_,S_,_j_virt); //*/
  //printf(" temp1 = %.8f \n", temp1);
  //printf(" temp2 = %.8f \n", temp2);
  //return 0.;

  double res =
     Triangle_T(o,k,l_,Q_,S_,l_,S_,_j_virt)
   + Triangle_0(o,k,l_,Q_,S_,l_,S_,ccc)/SQR(M)
    ;
  return  res*32.*M_PI ;
}//*/

void rho_j_scan1() {
  int N_M;
  double M, M_min, M_max, step, o, k;
  double res_r, res_v;

  char *prefix=(char*)"out/rho_j_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"1.dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");

  N_M=30; 
  M_min=.1;
  M_max=100.;

  M = M_min;
  step=pow(M_max/M_min,1./(double)(N_M-1));

  printf(" Settings: M_min=%g, M_max=%g\n",M_min,M_max); 
  double frac;

  for (int i=0; i<N_M; i++) {
    frac = (double)i/(double)(N_M-1);
    k = sqrt(k2av(M));

    res_r = rho_j_real(M,k,.001,.001,.01);
    res_v = rho_j_virt(M,k,.001,.001,.01);

    printf(" M = %.5e , [%2.2f%]\n", M , 100.*frac); 
    fprintf( out, "%.8e    %.8e    %.8e\n", M, res_r, res_v );
    M *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}//*/

void rho_j_scan2() {
  int N_M;
  double mQ, M_min, M_max, step, o, k;
  double res_r, res_v;

  char *prefix=(char*)"out/rho_j_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"4.dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");

  N_M=30; 
  M_min=.0001;
  M_max=.01;

  mQ = M_min;
  step=pow(M_max/M_min,1./(double)(N_M-1));
  double M=3.;
  k = sqrt(k2av(M));
  o = sqrt( SQR(M)+SQR(k) );

  printf(" Settings: M_min=%g, M_max=%g\n",M_min,M_max); 
  double frac;

  for (int i=0; i<N_M; i++) {
    frac = (double)i/(double)(N_M-1);

    res_r = rho_j_real(M,k,.001,mQ,.01); // l,Q,S
    res_v = rho_j_virt(M,k,.001,mQ,.01);

    printf(" M = %.5e , [%2.2f%]\n", M , 100.*frac); 
    fprintf( out, "%.8e    %.8e    %.8e\n", mQ, res_r, res_v );
    mQ *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}//*/

