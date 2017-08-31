r"""
Python bindings for various computation of Lyapunov exponents.
"""

from libc.stdlib cimport malloc,free

cdef extern from "hyperbolic_geodesics.h":
    void lyap_exp(size_t nb_iteration, size_t nb_vectors, size_t nb_coordinates, double complex *p_zero, double complex *p_infinity, double *theta)
    void lyap_exp_CY(size_t nb_iteration, size_t nb_vectors, double complex C, double complex d, double *theta)
    int init_simulation(size_t nb_v, size_t nb_c)
    void free_simulation()


def lyapunov_exponents(p_zero_py, p_infinity_py, nb_coordinates, nb_vectors,  nb_experiments, nb_iterations, CY=[]):
    r"""
    Compute the Lyapunov exponents of the geodesic flow in the hypergeometric function
    space.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``p_zero_py`` -- coefficient of the characteristic polynomial around zero

    - ``p_infinity_py`` -- coefficient of the characteristic polynomial around inifinity

    - ``nb_vectors`` -- the number of vectors to use

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the induction to perform

    """

    cdef double complex *p_zero, *p_infinity
    cdef double *theta

    # convert the data of into C values
    p_zero = <double complex *> malloc(nb_coordinates * sizeof(double complex))
    p_infinity = <double complex *> malloc(nb_coordinates * sizeof(double complex))

    for j from 0 <= j < nb_coordinates:
        p_zero[j] = <double complex> p_zero_py[j]
        p_infinity[j] = <double complex> p_infinity_py[j]

    theta = <double *> malloc(nb_vectors * sizeof(double))

    res = [[] for _ in xrange(nb_vectors)]

    init_simulation(nb_vectors, nb_coordinates);
    for i in xrange(nb_experiments):
        if CY:
            [a, b] = CY
            lyap_exp_CY(nb_iterations, nb_vectors, a, b, theta)
        else:
            lyap_exp(nb_iterations, nb_vectors, nb_coordinates, p_zero, p_infinity, theta)
        for j in xrange(nb_vectors):
            res[j].append(theta[j])

    free_simulation()

    free(p_zero)
    free(p_infinity)
    free(theta)

    return res
