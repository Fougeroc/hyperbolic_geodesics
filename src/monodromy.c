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


inline void monodromy_one(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *w){
  size_t i,j;
  double complex sum;

  for(i=0; i<nb_vectors; ++i){
    sum = 0;
    for(j=0; j<nb_coordinates; ++j)
      sum += v_all[i + nb_vectors * j];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] += w[j]*sum;
  }
}

inline void monodromy_one_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *w){
  size_t i,j;
  double complex sum, sum_w, coeff;

  sum_w = 0;
  for(j=0; j<nb_coordinates; ++j)
    sum_w += w[j];
  coeff = 1 / (1 + sum_w);

  for(i=0; i<nb_vectors; ++i){
    sum = 0;
    for(j=0; j<nb_coordinates; ++j)
      sum += v_all[i + nb_vectors * j];
    for(j=0; j<nb_coordinates; ++j)
      v_all[i + nb_vectors * j] -= coeff*w[j]*sum;
  }
}


inline void monodromy_infinity(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, w);
  monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
}

inline void monodromy_infinity_inverse(size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
  monodromy_one(nb_vectors, nb_coordinates, v_all, w);
}


inline void monodromy(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one(nb_vectors, nb_coordinates, v_all, w);
}


inline void monodromy_inverse(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex* v_all, double complex *e_alpha, double complex *w){
  if (n == 0)
    monodromy_infinity_inverse(nb_vectors, nb_coordinates, v_all, e_alpha, w);
  if (n == 1)
    monodromy_zero_inverse(nb_vectors, nb_coordinates, v_all, e_alpha);
  if (n == 2)
    monodromy_one_inverse(nb_vectors, nb_coordinates, v_all, w);
}
