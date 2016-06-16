#include "hyperbolic_geodesics.h"

void lyap_exp(size_t nb_iteration, size_t nb_vectors, size_t nb_coordinates, double complex *e_alpha, double complex *w, double *theta)
/*  - nb_iteration : number of iteration of the continued fraction algorithm to simulate geodesics
    - nb_vectors : number of vector we are applying the monodromy matrices to. (also the number of lyapunov exponents we will get)
    - nb_coordinates : size of the monodromy matrices and vectors
    - e_alpha : the eigenvalues of the monodromy around 0
    - w  : the vector such that monodromy around 1 is Id + w * [1]^t 
    - theta : a pointer to which we store the computed lyapunov exponents */
{
  long double x = gauss_rand(), y;
  unsigned long aux;
  int penalty = 0;
  size_t i, it;
  int current_letter = 0;

  double complex *v_all = malloc(nb_vectors * nb_coordinates * sizeof(double complex));

  set_random_vectors(v_all);
  orthogonalize_GS(v_all, theta);

  for(i=0; i < nb_vectors; ++i) theta[i] = 0.;

  for (it=0; it<nb_iteration; ++it)
    {
      y = 1./x + ((long double) (.5-drand())) / 0x40000000000;
      aux = floor(y);
      x = y - aux;

      // if aux is too big, we make it smaller for the code to end.
      if (aux > 100000000) aux = 100000000;


      if (current_letter < 0) current_letter += 3;

      aux -= penalty;
      penalty = 0;

      if (aux % 2 == 0) {
	if (it % 2 == 0) {
	  for (i=0; i<aux/2; ++i) {
	    monodromy(current_letter, nb_vectors, nb_coordinates, v_all, e_alpha, w);	      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }

	  
	  current_letter = (current_letter + 1) % 3;
	}

	else {
	  for (i=0; i<aux/2; ++i) {
	    monodromy_inverse(current_letter, nb_vectors, nb_coordinates, v_all, e_alpha, w);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }
	  
	  current_letter = (current_letter - 1) % 3;
	}
      }

      else {
	if (it % 2 == 0) {
	  for (i=0; i<(aux+1)/2; ++i) {
	    monodromy(current_letter, nb_vectors, nb_coordinates, v_all, e_alpha, w);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }

	  current_letter = (current_letter - 1) % 3;
	}

	else {
	  for (i=0; i<(aux+1)/2; ++i) {
	    monodromy_inverse(current_letter, nb_vectors, nb_coordinates, v_all, e_alpha, w);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }

	  current_letter = (current_letter + 1) % 3;
	}

	  penalty += 1;
      }
    }
    
  orthogonalize_GS(v_all, theta);

  for(i=0; i<nb_vectors; i++) theta[i] /=  2.37313822083125 * nb_iteration;

  free(v_all);
}


void lyap_exp_CY(size_t nb_iteration, size_t nb_vectors, double complex C, double complex d, double *theta)
/*  - nb_iteration : number of iteration of the continued fraction algorithm to simulate geodesics
    - nb_vectors : number of vector we are applying the monodromy matrices to. (also the number of lyapunov exponents we will get)
    - C : parameters for CY familly
    - d : parameters for CY familly
    - theta : a pointer to which we store the computed lyapunov exponents */
{
  long double x = gauss_rand(), y;
  unsigned long aux;
  int penalty = 0;
  size_t i, it;
  int current_letter = 0;

  size_t nb_coordinates = 4;
  double complex *v_all = malloc(nb_vectors * nb_coordinates * sizeof(double complex));

  set_random_vectors(v_all);
  orthogonalize_GS(v_all, theta);

  for(i=0; i < nb_vectors; ++i) theta[i] = 0.;

  for (it=0; it<nb_iteration; ++it)
    {
      y = 1./x + ((long double) (.5-drand())) / 0x40000000000;
      aux = floor(y);
      x = y - aux;

      // if aux is too big, we make it smaller for the code to end.
      if (aux > 100000000) aux = 100000000;

      if (current_letter < 0) current_letter += 3;

      aux -= penalty;
      penalty = 0;

      if (aux % 2 == 0) {
	if (it % 2 == 0) {
	  for (i=0; i<aux/2; ++i) {
	    monodromy_CY(current_letter, nb_vectors, v_all, C, d);	      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }
	  
	  current_letter = (current_letter + 1) % 3;
	}

	else {
	  for (i=0; i<aux/2; ++i) {
	    monodromy_CY_inverse(current_letter, nb_vectors, v_all, C, d);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }
	  
	  current_letter = (current_letter - 1) % 3;
	}
      }

      else {
	if (it % 2 == 0) {
	  for (i=0; i<(aux+1)/2; ++i) {
	    monodromy_CY(current_letter, nb_vectors, v_all, C, d);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }

	  current_letter = (current_letter - 1) % 3;
	}

	else {
	  for (i=0; i<(aux+1)/2; ++i) {
	    monodromy_CY_inverse(current_letter, nb_vectors, v_all, C, d);		      
	    if (i % 100 == 1 && test_norm(v_all))
	      orthogonalize_GS(v_all, theta);
	  }

	  current_letter = (current_letter + 1) % 3;
	}

	  penalty += 1;
      }
    }
    
  orthogonalize_GS(v_all, theta);

  for(i=0; i<nb_vectors; i++) theta[i] /=  2.37313822083125 * nb_iteration;

  free(v_all);
}
