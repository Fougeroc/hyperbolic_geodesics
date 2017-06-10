#include "hyperbolic_geodesics.h"

void monodromy_zero(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha){
  size_t i,j;
  for(i=0; i<nb_vectors; ++i)
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] *= e_alpha[j];
}

void monodromy_zero_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha){
  size_t i,j;
  for(i=0; i<nb_vectors; ++i)
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] /= e_alpha[j];
}


void monodromy_one(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // We assume here that we took a basis such that Minf is diagonal and M1 = Id + M0^-1*(1,...,1).transpose()*w^t
  size_t i,j;
  double complex scalar;

  for(i=0; i<nb_vectors; ++i){
    scalar = 0;
    for(j=0; j<nb_coordinates; ++j)
      scalar += v_all[i + nb_vectors * j]*w[j];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] += scalar/e_alpha[j];
  }
}

void monodromy_one_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // We have M1.v = v' = v + <w,v>*(e^-a1,..,e^-an) so <w,v'> = <w,v> + <w,v>*<w,1/alpha> and <w,v> = <w,v'>/(1+<w,1/alpha>)
  // Hence v = v' - <w,v'>/(1+<w,1/alpha>)*(e^-a1,..,e^-an)
  size_t i,j;
  double complex scalar, sum_w, coeff;

  sum_w = 0;
  for(j=0; j<nb_coordinates; ++j)
    sum_w += w[j]/e_alpha[j];
  coeff = 1+sum_w;

  for(i=0; i<nb_vectors; ++i){
    scalar = 0;
    for(j=0; j<nb_coordinates; ++j)
      scalar += v_all[i + nb_vectors * j]*w[j];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] -= scalar/(coeff*e_alpha[j]);
  }
}


void monodromy_infinity(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // MInfty = (M0*M1)^-1 = M1^-1*M0^-1
  monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
  monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
}

void monodromy_infinity_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // MInfty^-1 = M0*M1
  monodromy_one(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
}


void monodromy(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one(nb_vectors, nb_coordinates, v_all, e_alpha, w);
}

void monodromy_inverse(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
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
    monodromy_T(nb_vectors, v_all);                         // MInfty
  if (n == 1)
    monodromy_S(nb_vectors, v_all, C, d);                   // M0
  if (n == 2)
    monodromy_TS_inverse(nb_vectors, v_all, C, d);          // M1
}

void monodromy_CY_inverse_K(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  if (n == 0)
    monodromy_T_inverse(nb_vectors, v_all);
  if (n == 1)
    monodromy_S_inverse(nb_vectors, v_all, C, d);
  if (n == 2)
    monodromy_TS(nb_vectors, v_all, C, d);
}


void monodromy_CY_zero(size_t nb_vectors, double complex* v_all){
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

void monodromy_CY_zero_inv(size_t nb_vectors, double complex* v_all){
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

void monodromy_CY_inf(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
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

void monodromy_CY_inf_inv(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
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
  monodromy_CY_zero_inv(nb_vectors, v_all);
  monodromy_CY_inf_inv(nb_vectors, v_all,a,b);
}

void monodromy_CY_one_inv(size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  monodromy_CY_inf(nb_vectors, v_all,a,b);
  monodromy_CY_zero(nb_vectors, v_all);
}


void monodromy_CY(size_t n, size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id 
  if (n == 0)
    monodromy_CY_inf(nb_vectors, v_all, a, b);             // MInfty
  if (n == 1)
    monodromy_CY_zero(nb_vectors, v_all);                 // M0
  if (n == 2)
    monodromy_CY_one(nb_vectors,v_all, a, b);             // M1
}

void monodromy_CY_inverse(size_t n, size_t nb_vectors, double complex* v_all, double complex a, double complex b){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id 
  if (n == 0)
    monodromy_CY_inf_inv(nb_vectors, v_all, a, b);            // MInfty
  if (n == 1)
    monodromy_CY_zero_inv(nb_vectors, v_all);                 // M0
  if (n == 2)
    monodromy_CY_one_inv(nb_vectors, v_all, a, b);            // M1
}
