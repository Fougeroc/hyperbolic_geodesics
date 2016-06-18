#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

double log(double x);


/****************************/
/* custom random functions */
/***************************/

inline double drand(void);
inline long double ldrand(void);
inline double gauss_rand(void);
inline double gauss_map(double x);

/**********************************/
/* data and functions for vectors */
/**********************************/

void set_vectors_caracteristics(size_t nb_vectors, size_t nb_coordinates);
int init_simulation(size_t nb_v, size_t nb_c);
void free_simulation(void);

void print_vectors(double complex *v_all);
void set_random_vectors(double complex *v_all);

void orthogonalize_GS(double complex *v_all, double *theta);
void ortho_plus_check(double complex *v_all, double complex *v_buffer, double *theta);
void check_orthogonality(double complex *v_all);
inline double first_norm(double complex *v_all);
inline double max_norm(double complex *v_all);
inline double min_norm(double complex *v_all);
inline int test_norm(double complex *v_all);


/**********************************************/
/* monodromy functions and lyapunov exponents */
/**********************************************/

inline void monodromy(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex *v_all, double complex *e_alpha, double complex *w);
inline void monodromy_inverse(size_t n, size_t nb_vectors, size_t nb_coordinates, double complex *v_all, double complex *e_alpha, double complex *w);

void lyap_exp(size_t nb_iteration, size_t nb_vectors, size_t nb_coordinates, double complex *e_alpha, double complex *w, double *theta);

inline void monodromy_CY(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d);
inline void monodromy_CY_inverse(size_t n, size_t nb_vectors, double complex* v_all, double complex C, double complex d);

void lyap_exp_CY(size_t nb_iteration, size_t nb_vectors, double complex C, double complex d, double *theta);
