#include "hyperbolic_geodesics.h"

void companion_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *p){
  size_t i,j;
  double complex aux;

  for(i=0; i<nb_vectors; ++i){
    aux = v_all[i];
    for(j=0; j<nb_coordinates-1; ++j)
      v_all[i + nb_vectors*j] = v_all[i+nb_vectors*(j+1)];
    v_all[i+nb_vectors*(nb_coordinates-1)] = 0+0*I;

    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors*j] -= conj(p[nb_coordinates-j-1])*aux;
  }
}

void companion(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *p){
  size_t i,j;
  double complex aux;

  for(i=0; i<nb_vectors; ++i){
    aux = v_all[i+nb_vectors*(nb_coordinates-1)];
    for(j=nb_coordinates-1; j>0; --j)
      v_all[i+nb_vectors*j] = v_all[i+nb_vectors*(j-1)];
    v_all[i] = 0+0*I;

    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors*j] -= p[j]*aux;
  }
}

void monodromy(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all,
	       double complex *p_zero, double complex *p_infinity){
  if (n == 0)
    // M_inf
    companion(nb_vectors, nb_coordinates, v_all, p_infinity);
  if (n == 1)
    // M_0
    companion_inverse(nb_vectors, nb_coordinates, v_all, p_zero);
  if (n == 2){
    // M_1 = M_0^-1 * M_inf^-1
    companion_inverse(nb_vectors, nb_coordinates, v_all, p_infinity);
    companion(nb_vectors, nb_coordinates, v_all, p_zero);
  }
}


void monodromy_inverse(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all,
		       double complex *p_zero, double complex *p_infinity){
  if (n == 0)
    // M_inf^-1
    companion_inverse(nb_vectors, nb_coordinates, v_all, p_infinity);
  if (n == 1)
    // M_0^-1
    companion(nb_vectors, nb_coordinates, v_all, p_zero);
  if (n == 2){
    // M_1^-1 = M_inf*M_0
    companion_inverse(nb_vectors, nb_coordinates, v_all, p_zero);
    companion(nb_vectors, nb_coordinates, v_all, p_infinity);
  }
}

void monodromy_T(size_t nb_vectors, double complex* v_all){
  size_t nb_coordinates = 4;
  size_t i,j;
  double complex res[4] = {0,0,0,0};

  for(i=0; i<nb_vectors; ++i) {
    res[0] = v_all[i];
    res[1] = v_all[i] + v_all[i + nb_vectors * 1];
    res[2] = v_all[i]/2 + v_all[i + nb_vectors * 1] + v_all[i + nb_vectors * 2];
    res[3] = v_all[i]/6 + v_all[i + nb_vectors * 1]/2 + v_all[i + nb_vectors * 2] + v_all[i + nb_vectors * 3];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] = res[j];
  }
}

void monodromy_T_inverse(size_t nb_vectors, double complex* v_all){
  size_t nb_coordinates = 4;
  size_t i,j;
  double complex res[4] = {0,0,0,0};

  for(i=0; i<nb_vectors; ++i) {
    res[0] = v_all[i];
    res[1] = -v_all[i] + v_all[i + nb_vectors * 1];
    res[2] = v_all[i]/2 - v_all[i + nb_vectors * 1] + v_all[i + nb_vectors * 2];
    res[3] = -v_all[i]/6 + v_all[i + nb_vectors * 1]/2 - v_all[i + nb_vectors * 2] + v_all[i + nb_vectors * 3];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] = res[j];
  }
}

void monodromy_S(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  size_t i;

  for(i=0; i<nb_vectors; ++i)
    v_all[i] -= C*v_all[i + nb_vectors * 1]/12 + d*v_all[i + nb_vectors * 3];
}

void monodromy_S_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  size_t i;

  for(i=0; i<nb_vectors; ++i)
    v_all[i] += C*v_all[i + nb_vectors * 1]/12 + d*v_all[i + nb_vectors * 3];
}

void monodromy_ST(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_T(nb_vectors, v_all);
  monodromy_S(nb_vectors, v_all, C, d);
}

void monodromy_ST_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_S_inverse(nb_vectors, v_all, C, d);
  monodromy_T_inverse(nb_vectors, v_all);
}

void monodromy_TS(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_S(nb_vectors, v_all, C, d);
  monodromy_T(nb_vectors, v_all);
}

void monodromy_TS_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_T_inverse(nb_vectors, v_all);
  monodromy_S_inverse(nb_vectors, v_all, C, d);
}


void monodromy_CY_K(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id
  if (n == 0)
    monodromy_T(nb_vectors, v_all);                         // M0
  if (n == 1)
    monodromy_S(nb_vectors, v_all, C, d);                   // M1
  if (n == 2)
    monodromy_TS_inverse(nb_vectors, v_all, C, d);          // Minf
}

void monodromy_CY_inverse_K(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  if (n == 0)
    monodromy_T_inverse(nb_vectors, v_all);
  if (n == 1)
    monodromy_S_inverse(nb_vectors, v_all, C, d);
  if (n == 2)
    monodromy_TS(nb_vectors, v_all, C, d);
}


void monodromy_CY_zero_inv(size_t nb_vectors, double complex* v_all){
  size_t i,j;
  double complex res[4] = {0,0,0,0};

  for(i=0; i<nb_vectors; ++i){
    res[0] = 0;
    for(j=0; j<3; ++j)
      res[j+1] = v_all[i+nb_vectors*j];

    res[0] -= 1*v_all[i+nb_vectors*3];
    res[1] += 4*v_all[i+nb_vectors*3];
    res[2] -= 6*v_all[i+nb_vectors*3];
    res[3] += 4*v_all[i+nb_vectors*3];

    for(j=0; j<4; ++j)
      v_all[i + nb_vectors*j] = res[j];
  }
}

void monodromy_CY_zero(size_t nb_vectors, double complex* v_all){
  size_t i,j;
  double complex res[4] = {0,0,0,0};

  for(i=0; i<nb_vectors; ++i){
    for(j=0; j<3; ++j)
      res[j] = v_all[i+nb_vectors*(j+1)];
    res[3] = 0;

    res[0] += 4*v_all[i];
    res[1] -= 6*v_all[i];
    res[2] += 4*v_all[i];
    res[3] -= 1*v_all[i];

    for(j=0; j<4; ++j)
      v_all[i + nb_vectors*j] = res[j];
  }
}

void monodromy_CY_inf_inv(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  size_t i,j;
  double complex res[4] = {0,0,0,0};

  for(i=0; i<nb_vectors; ++i){
    for(j=0; j<3; ++j)
      res[j] = v_all[i+nb_vectors*(j+1)];
    res[3] = 0;

    res[0] += a*v_all[i];
    res[1] += b*v_all[i];
    res[2] += a*v_all[i];
    res[3] -= 1*v_all[i];

    for(j=0; j<4; ++j)
      v_all[i + nb_vectors*j] = res[j];
  }
}

void monodromy_CY_inf(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  size_t i,j;
  double complex res[4] = {0,0,0,0};

   for(i=0; i<nb_vectors; ++i){
    res[0] = 0;
    for(j=0; j<3; ++j)
      res[j+1] = v_all[i+nb_vectors*j];

    res[0] -= v_all[i+nb_vectors*3];
    res[1] += a*v_all[i+nb_vectors*3];
    res[2] += b*v_all[i+nb_vectors*3];
    res[3] += a*v_all[i+nb_vectors*3];

    for(j=0; j<4; ++j)
      v_all[i + nb_vectors*j] = res[j];
   }
}

void monodromy_CY_one(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  monodromy_CY_inf_inv(nb_vectors, v_all, a, b);
  monodromy_CY_zero_inv(nb_vectors, v_all);
}

void monodromy_CY_one_inv(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  monodromy_CY_zero(nb_vectors, v_all);
  monodromy_CY_inf(nb_vectors, v_all, a, b);
}


void monodromy_CY(size_t n, size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id
  if (n == 0)
    monodromy_CY_zero(nb_vectors, v_all);             // M0
  if (n == 1)
    monodromy_CY_one(nb_vectors, v_all, a, b);             // M1
  if (n == 2)
    monodromy_CY_inf(nb_vectors,v_all, a, b);             // Minf
}

void monodromy_CY_inverse(size_t n, size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id
  if (n == 0)
    monodromy_CY_zero_inv(nb_vectors, v_all);               // M0
  if (n == 1)
    monodromy_CY_one_inv(nb_vectors, v_all, a, b);          // M1
  if (n == 2)
    monodromy_CY_inf_inv(nb_vectors, v_all, a, b);          // Minf
}
