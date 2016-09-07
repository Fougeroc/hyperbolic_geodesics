#include "hyperbolic_geodesics.h"

/* random (?) uniform double in [0,1] */
double drand(void)
{
	return (((double) random() + (double) random() / ( (double) RAND_MAX) ) / ((double) RAND_MAX+1) );
}

/* random (?) uniform long double in [0,1] */
long double ldrand(void)
{
	return ((long double) random() + (((long double) random() + (long double) random() / ((long double) RAND_MAX)) / ((long double) RAND_MAX))) / ((long double) RAND_MAX+1);
}


double gauss_rand(void)
{
  double y = drand();
  return( exp( y * log(2.) ) - 1. );
}

double gauss_map(double x)
{
  double y = 1/x;
  return( y - floor(y) );
}
