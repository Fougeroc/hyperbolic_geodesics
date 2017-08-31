/* Some highly non optimized linear algebra routines */
/* 1) random orthogonal vectors                      */
/* 2) orthogonalization                              */

#include "hyperbolic_geodesics.h"
#include <time.h>

/* 0 <= i < nb_vectors   */
/* 0 <= j < nb_coordinates       */
/*  index at pos v_i[j] = v[i + nb_vectors * j]  */

/* global variable used in GS orthogonalisation                                 */
/* they must be dynamic so we keep them here in order to avoid many malloc/free */
double complex * scal     = NULL;
double complex * scal_new = NULL;
size_t scal_size  = 0;
int nb_vectors = 0, nb_coordinates = 0;

void print_vectors(double complex * v_all)
{
	size_t i,j;

	for(i=0; i<nb_vectors; ++i)
	{
		printf("v[%zu]=",i);
		for(j=0; j < nb_coordinates; ++j)
		{
		  printf(" %f + %f*i", creal(v_all[i + nb_vectors * j]), cimag(v_all[i + nb_vectors * j]));
		  printf(",");
		}
		printf("\n");
	}
}

void set_random_vectors(double complex * v_all)
/* set random orthogonal frame */
/* warning: before calling this function call init_GS(dim) */
{
	size_t i,j;

	for(i=0; i<nb_vectors; ++i)
	{
		for(j=0; j<nb_coordinates; j++)
		  v_all[i + nb_vectors * j] = (drand() - .5) + (drand() - .5)*I;
	}
	orthogonalize_GS(v_all, NULL);
}

void set_identity_vectors(double complex * v_all)
/* set random orthogonal frame */
/* warning: before calling this function call init_GS(dim) */
{
	size_t i,j;

	for(i=0; i<nb_vectors; ++i)
	{
	  for(j=0; j<nb_coordinates; j++){
	    if (i == j) {
	      v_all[i + nb_vectors * j] = 1;
	    }
	    else {
	      v_all[i + nb_vectors * j] = 0;
	    };
	  }
	}
	//	orthogonalize_GS(v_all, NULL);
}


int init_simulation(size_t nb_v, size_t nb_c)
/* allocate scal and scal_new in order that it contains nb_vectors elements */
{
  nb_vectors = nb_v;
  nb_coordinates = nb_c;

  if(scal == NULL || scal_new == NULL) {
    scal = (double complex *) malloc(nb_vectors * sizeof(double complex));
    scal_new = (double complex *) malloc(nb_vectors * sizeof(double complex));
  }

  srand(time(NULL));
  return 0;
}

void free_simulation(void)
{
	scal_size = 0;
	if(scal != NULL)
	{
		free(scal);
		free(scal_new);
		scal = NULL;
		scal_new = NULL;
	}
}

void orthogonalize_GS(double complex *v_all, double *theta)
/* orthogonalization using Gram-Schmidt process                                                */
/* warning: this method should be called AFTER being sure that init_GS has been called         */
/* it theta is not NULL it is updated by 2*log of the diagonal from position 1 to nb_vectors   */
/* warning: it is 2*log(entry) and not log(entry)                                              */
/* element at position (i,j,k): (qcc->labels)[j].v[i + nb_vectors*k];                          */
{
	size_t i,ii,j;

	double norm,sqnorm;
	double complex c;
	double complex * tmp = NULL;

	/* some check that we will remove */
	if(nb_vectors < 1)
	{
		fprintf(stderr, "calling GS with nb_vectors < 1 is not possible.\n");
		exit(EXIT_FAILURE);
	}

	if(scal == NULL || scal_new == NULL)
	{
		fprintf(stderr, "you must call init_GS before calling orthogonalize_GS.\n");
		exit(EXIT_FAILURE);
	}

	if (nb_vectors == 1) {
	  sqnorm = 0.;
	  for(j=0; j<nb_coordinates; ++j)
	    sqnorm  += creal(v_all[j] * conj(v_all[j]));
	  if (theta != NULL) theta[0] += log(sqnorm);
	  norm = sqrt(sqnorm);
	  for(j=0; j<nb_coordinates; ++j)
	    v_all[j] /= norm;
	}

	if (nb_vectors > 1) {
	  /* put v in the holonomy free subspace */
	  /* compute <v_0,v_1> in scal_new[0] */
	  /* compute <v_0,v_0> in sqnorm */
	  scal_new[0] = 0.;
	  sqnorm = 0.;
	  for(j=0; j<nb_coordinates; ++j)
	    {
	      scal_new[0] += v_all[0 + nb_vectors * j] * conj(v_all[1 + nb_vectors * j]);
	      sqnorm  += creal(v_all[0 + nb_vectors * j] * conj(v_all[0 + nb_vectors * j]));
	    }

	  /* vector by vector orthonormalisation */
	  for(i=1; i < nb_vectors-1; ++i)
	    {
	      /* sqnorm contains <v_(i-1),v_(i-1)> */
	      /* v_0, v_1, ..., v_(i-2) are normalized */
	      /* v_0, v_1, ..., v_(i-1) are orthogonal */
	      /* scal_new contains the i scalar products <v_0,v_i>, <v_1,v_i>, ..., <v_(i-1),v_i> */
	      tmp = scal; scal = scal_new; scal_new = tmp;
	      for(ii=0; ii<=i; ++ii)
		scal_new[ii] = 0.;

	      c = scal[i-1]/sqnorm;
	      if(theta != NULL) theta[i-1] += log(sqnorm);
	      norm = sqrt(sqnorm);
	      sqnorm = 0.;

	      /* c = <v_i,v_(i-1)> / <v_(i-1,v_(i-1)> */
	      for(j=0; j < nb_coordinates; ++j)
		{
		  /* subtract the part of the span of v_0, v_1, ..., v_(i-2) */
		  for(ii=0; ii<i-1; ++ii)
		    v_all[i + nb_vectors * j] -= conj(scal[ii]) * v_all[ii + nb_vectors * j];

		  /* subtract the part of the span of v_(i-1) */
		  v_all[i + nb_vectors * j] -= conj(c) * v_all[ii + nb_vectors * j];

		  /* normalize v_(i-1) */
		  v_all[(i-1) + nb_vectors * j] /= norm;

		  /* compute scalar products and norms for next loop */
		  /* sqnorm = <v_i, v_i> */
		  /* scal_new[ii] = <v_(i+1), v_ii>*/
		  sqnorm += creal(v_all[i + nb_vectors * j] * conj(v_all[i + nb_vectors * j]));
		  for(ii=0; ii<=i; ++ii)
		    scal_new[ii] += v_all[ii + nb_vectors * j] * conj(v_all[i+1 + nb_vectors * j]);
		}
	    }

	  /* here i = NB_VECTORS-1 */
	  /* renormalize v_(i-1)   */
	  /* put v_i in the orthogonal of the span of v_0, v_1, ..., v_(i-1) */
	  c = scal_new[i-1] / sqnorm;
	  if(theta != NULL) theta[i-1] += log(sqnorm);
	  norm = sqrt(sqnorm);
	  sqnorm = .0;

	  for(j=0; j< nb_coordinates; ++j)
	    {
	      for(ii=0; ii<i-1; ++ii)
		v_all[i + nb_vectors * j] -= conj(scal_new[ii]) * v_all[ii + nb_vectors * j];

	      v_all[i + nb_vectors * j] -= conj(c) * v_all[i-1 + nb_vectors * j];

	      v_all[i-1 + nb_vectors * j] /= norm;

	      sqnorm += creal(v_all[i + nb_vectors * j] * conj(v_all[i + nb_vectors * j]));
	    }

	  /* we renormalize v_i*/
	  if(theta != NULL) theta[i] += log(sqnorm);
	  norm = sqrt(sqnorm);
	  for(j=0; j < nb_coordinates; ++j)
	    v_all[i + nb_vectors * j] /= norm;
	}
}



void check_orthogonality(double complex *v_all)
{
	double complex s;
	size_t i1,i2,k;

	for(i1=0; i1<nb_vectors; ++i1)
	  for(i2=i1; i2<nb_vectors; ++i2)
	    {
	      s = 0;
	      for(k=0; k < nb_coordinates; ++k)
		s += v_all[i1 + nb_vectors*k] * conj(v_all[i2 + nb_vectors*k]);

	      if (i1==i2 && creal(s - 1) > 0.01){
		fprintf(stderr, "Wrong normalisation\nNorm : %lf\nDiff : %lf > %lf\n", creal(s), creal(s-1), 1./0x4000);
		print_vectors(v_all);
		exit(EXIT_FAILURE);
	      }

	      if (i1!=i2 && creal(s) > 0.01) {
		fprintf(stderr, "Wrong orthogonalisation\nScalar product <v[i%zu], v[i%zu]> : %lf\nnb_vector : %i nb_coordinates : %i\n\n",
			i1, i2, creal(s), nb_vectors, nb_coordinates);
		print_vectors(v_all);
		exit(EXIT_FAILURE);
	      }
	    }
}



inline double first_norm(double complex *v_all)
{
	double sqnorm = 0.;
	size_t j;

	for(j=0; j < nb_coordinates; ++j)
	  sqnorm += creal(v_all[0 + nb_vectors*j] * conj(v_all[0 + nb_vectors*j]));

	return sqnorm;
}



inline double max_norm(double complex *v_all)
{
	double sqnorm = 0.;
	double max = 0;
	size_t i,j;

	for (i=0; i<nb_vectors; ++i){
	  for(j=0; j < nb_coordinates; ++j)
	    sqnorm += creal(v_all[i + nb_vectors*j] * conj(v_all[i + nb_vectors*j]));
	  if (sqnorm > max)
	    max = sqnorm;
	}

	return max;
}


inline double min_norm(double complex *v_all)
{
	double sqnorm = 0.;
	double min = 1;
	size_t i,j;

	for (i=0; i<nb_vectors; ++i){
	  for(j=0; j < nb_coordinates; ++j)
	    sqnorm += creal(v_all[i + nb_vectors*j] * conj(v_all[i + nb_vectors*j]));
	  if (sqnorm < min)
	    min = sqnorm;
	}

	return min;
}


inline int test_norm(double complex *v_all)
{
	double sqnorm = 0.;
	double max = 0, min = 1;
	size_t i,j;

	for (i=0; i<nb_vectors; ++i){
	  for(j=0; j < nb_coordinates; ++j)
	    sqnorm += creal(v_all[i + nb_vectors*j] * conj(v_all[i + nb_vectors*j]));
	  if (sqnorm > max)
	    max = sqnorm;
	  if (sqnorm < min)
	    min = sqnorm;
	}

	return ((max > 1024.) || (min < .0001));
}



void ortho_plus_check(double complex *v_all, double complex *v_buffer, double* theta)
{
	double complex s;
	size_t i1,i2,k;

	for (i1=0; i1<nb_vectors; ++i1)
	  for (i2=0; i2<nb_coordinates; ++i2)
	    v_buffer[i1 + nb_vectors * i2] = v_all[i1 + nb_vectors * i2];

	orthogonalize_GS(v_all, theta);

	for(i1=0; i1<nb_vectors; ++i1)
	  for(i2=i1; i2<nb_vectors; ++i2)
	    {
	      s = 0;
	      for(k=0; k < nb_coordinates; ++k)
		s += v_all[i1 + nb_vectors*k] * conj(v_all[i2 + nb_vectors*k]);

	      if (i1==i2 && creal(s - 1) > 0.01){
		fprintf(stderr, "Wrong normalisation\nNorm : %lf\nDiff : %lf > %lf\n", creal(s), creal(s-1), 1./0x4000);
		print_vectors(v_buffer);
		print_vectors(v_all);
		exit(EXIT_FAILURE);
	      }

	      if (i1!=i2 && creal(s) > 0.01) {
		fprintf(stderr, "Wrong orthogonalisation\nScalar product <v[i%zu], v[i%zu]> : %lf\nnb_vector : %i nb_coordinates : %i\n\n",
			i1, i2, creal(s), nb_vectors, nb_coordinates);
		print_vectors(v_buffer);
		print_vectors(v_all);
		exit(EXIT_FAILURE);
	      }
	    }
}
