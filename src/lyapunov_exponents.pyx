r"""
Python bindings for various computation of Lyapunov exponents.
"""

from libc.stdlib cimport malloc,free

cdef extern from "hyperbolic_geodesics.h":
    void lyap_exp(size_t nb_iteration, size_t nb_vectors, size_t nb_coordinates, double complex *e_alpha, double complex *w, double *theta)
    int init_simulation(size_t nb_v, size_t nb_c)
    void free_simulation()


def lyapunov_exponents(e_alpha_py, w_py, nb_coordinates, nb_vectors,  nb_experiments, nb_iterations):
    r"""
    Compute the Lyapunov exponents of the geodesic flow in the hypergeometric function
    space.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``a_py`` -- roots of thenumerator polynom

    - ``w_py`` -- vector associated to monodromy of one

    - ``nb_vectors`` -- the number of vectors to use

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the induction to perform

    """
    
    cdef double complex *e_alpha, *w
    cdef double *theta

    # convert the data of into C values
    e_alpha = <double complex *> malloc(nb_coordinates * sizeof(double complex))
    w = <double complex *> malloc(nb_coordinates * sizeof(double complex))
   
    for j from 0 <= j < nb_coordinates:
        e_alpha[j] = e_alpha_py[j]
        w[j] = w_py[j]

    theta = <double *> malloc(nb_vectors * sizeof(double))

    res = [[] for _ in xrange(nb_vectors)]    

    init_simulation(nb_vectors, nb_coordinates); 
    for i in xrange(nb_experiments):
        lyap_exp(nb_iterations, nb_vectors, nb_coordinates, e_alpha, w, theta)
        for j in xrange(nb_vectors):
            res[j].append(theta[j])

    free_simulation()

    free(e_alpha)
    free(w)
    free(theta)

    return res
