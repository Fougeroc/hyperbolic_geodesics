#include "hyperbolic_geodesics.h"

inline void monodromy_zero(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha){
  size_t i,j;
  for(i=0; i<nb_vectors; ++i)
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] *= e_alpha[j];
}

inline void monodromy_zero_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha){
  size_t i,j;
  for(i=0; i<nb_vectors; ++i)
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] /= e_alpha[j];
}


inline void monodromy_one(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
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

inline void monodromy_one_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
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


inline void monodromy_infinity(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // MInfty = (M0*M1)^-1 = M1^-1*M0^-1
  monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
  monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
}

inline void monodromy_infinity_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  // MInfty^-1 = M0*M1
  monodromy_one(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
}


inline void monodromy(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one(nb_vectors, nb_coordinates, v_all, e_alpha, w);
}

inline void monodromy_inverse(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
}


inline void monodromy_T(size_t nb_vectors, double complex* v_all){
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

inline void monodromy_T_inverse(size_t nb_vectors, double complex* v_all){
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

inline void monodromy_S(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  size_t i;

  for(i=0; i<nb_vectors; ++i)
    v_all[i] -= C*v_all[i + nb_vectors * 1]/12 + d*v_all[i + nb_vectors * 3];
}

inline void monodromy_S_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  size_t i;

  for(i=0; i<nb_vectors; ++i)
    v_all[i] += C*v_all[i + nb_vectors * 1]/12 + d*v_all[i + nb_vectors * 3];
}

inline void monodromy_ST(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_T(nb_vectors, v_all);
  monodromy_S(nb_vectors, v_all, C, d);
}

inline void monodromy_ST_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_S_inverse(nb_vectors, v_all, C, d);
  monodromy_T_inverse(nb_vectors, v_all);
}

inline void monodromy_TS(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_S(nb_vectors, v_all, C, d);
  monodromy_T(nb_vectors, v_all);
}

inline void monodromy_TS_inverse(size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  monodromy_T_inverse(nb_vectors, v_all);
  monodromy_S_inverse(nb_vectors, v_all, C, d);
}


inline void monodromy_CY(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  //Minf*M0*M1 = Id, M0*M1*Minf = Id, A*B*C = Id 
  if (n == 0)
    monodromy_T(nb_vectors, v_all);                         // MInfty
  if (n == 1)
    monodromy_S(nb_vectors, v_all, C, d);                   // M0
  if (n == 2)
    monodromy_TS_inverse(nb_vectors, v_all, C, d);          // M1
}

inline void monodromy_CY_inverse(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d){
  if (n == 0)
    monodromy_T_inverse(nb_vectors, v_all);
  if (n == 1)
    monodromy_S_inverse(nb_vectors, v_all, C, d);
  if (n == 2)
    monodromy_TS(nb_vectors, v_all, C, d);
}
