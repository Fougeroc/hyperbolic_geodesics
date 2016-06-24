from pylab import *
from matplotlib import *

def lin_space(a, b, s):
    delta = (b-a)/(s-1.)
    return [a + delta*k for k in range(s)]

E = lambda z: exp(2*1j*pi*z)

class Experiment(object):
    r"""
    Class of hypergeometric flat bundle on the sphere minus three points.

    Has three attributes :

    - ``_dimension`` -- dimension of the flat bundle

    - ``_alpha`` -- roots of the denominator polynom - also eigenvalues of monodromy around 0

    - ``_beta`` -- roots of the numerator polynom - also minus eigenvalues of monodromy around infinity
    """
    _dimension = None
    _beta = None
    _alpha = None
    _zero_diagonalizable = True

    def __repr__(self):
        return str(self._alpha) + ", " + str(self._beta)

    def __copy__(self):
        r"""
        Returns a copy of self.

        EXAMPLES::

            sage: p = Experiment([0, 0.5], [0.2, 0.3])
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
        """
        from copy import copy
        q = self.__class__()

        q.__hash__ = self.__hash__
        q._beta = copy(self._beta)
        q._alpha = copy(self._alpha)
        q._zero_diagonalizable = copy(self._zero_diagonalizable)
        q._dimension = self._dimension

        return q

    def __len__(self) :
        return self._dimension

    def __init__(self, alpha, beta, zero_diagonalizable=True, x_plot=None, y_plot=None):
        if len(alpha) <> len(beta):
            raise ValueError("The two parameter lists must be of the same length")

        n = len(alpha)
        if not(set([E(beta[j]) for j in range(n)]).isdisjoint(set([E(alpha[j]) for j in range(n)]))):
            print "Warning the sets of eigenvalues are not disjoint, the flat bundle won't be minimal, adding a little bit of noise"
            print alpha, beta
            alpha = map(lambda x: x + random()/2**15, alpha)
            beta = map(lambda x: x + random()/2**15, beta)

        self._dimension = len(alpha)
        self._alpha = [alpha[k] if 0<=alpha[k]<=1 else alpha[k] - floor(alpha[k]) for k in xrange(self._dimension)]
        self._beta = [beta[k] if 0<=beta[k]<=1 else beta[k] - floor(beta[k]) for k in xrange(self._dimension)]
        self._alpha.sort()
        self._beta.sort()
        self._zero_diagonalizable = zero_diagonalizable
        self.x_plot = x_plot
        self.y_plot = y_plot

    def hypergeometric_lyap_exp(self, nb_vectors=None, nb_experiments=10 ,
                                nb_iterations=10**4, verbose=False, output_file=None, return_error=False):
        r"""
        Compute the Lyapunov exponents of the geodesic flow in the hypergeometric function
        space.
        
        INPUT:

        - ``nb_vectors`` -- the number of vectors to use
        
        - ``nb_experiments`` -- number of experimets
        
        - ``nb_iterations`` -- the number of iterations of the induction to perform
        
        - ``output_file`` -- place where we print the results
        
        - ``verbose`` -- do we print the result with an extensive number of information or not
    
        """
        import time
        import lyapunov_exponents    # the cython bindings
        from math import sqrt
        
        if output_file is None:
            from sys import stdout
            output_file = stdout
            closed = True
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")
            closed = False

        if nb_vectors <> None and nb_vectors <= 0:
            raise ValueError("the number of vectors must be positive")
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        #recall that the lyapunov exponents are symmetric
        if nb_vectors == None:
            nb_vectors = self._dimension//2

        if max(self.hodge_numbers()) == self._dimension:
            if verbose:
                output_file.write("by signature of the invariant hermitian form, all lyapunov exponents are 0\n")
            return [0]*nb_vectors, [0]*nb_vectors, 0, 0

        e_alpha, w = self.compute_monodromy()

        t0 = time.time()
        res = lyapunov_exponents.lyapunov_exponents(e_alpha, w, self._dimension, nb_vectors,  nb_experiments, nb_iterations)
        t1 = time.time()

        res_final = []
        std_final = []
        s_m, s_d = 0, 0

        if verbose:
            from math import floor, log
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2**%d)\n"%(nb_iterations, floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
        for i in xrange(nb_vectors):
            m,d = mean_and_std_dev(res[i])
            s_m += m
            s_d += d**2
            if verbose:
                output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    i,m,d, 2.576*d/sqrt(nb_experiments)))
            res_final.append(m)
            std_final.append(2.576*d/sqrt(nb_experiments))

        s_d = sqrt(s_d)
        s_d_final = 2.576*s_d/sqrt(nb_experiments)
        if verbose:
            output_file.write("sum_theta        : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n\n"%(
                s_m,s_d, 2.576*s_d/sqrt(nb_experiments)))
        
        if not closed :
            output_file.close()
            print "file closed"

        if return_error:
            return (res_final, std_final, s_m, s_d_final)
        else:
            return res_final

    def compute_monodromy(self):
        r"""
        Return vectors e_alpha and v, associated to monodromy (When winding around
        singular points counterclockwise).

        We conjugate all the matrices such that the monodromy matrix around zero
        is diagonal, and the one around 1 is identity + column_matrix(1, .., 1)*w.transpose().

        Then calculus will be easy since we just make scalar product with v and add it in each
        coordinate.
        Those vector will be used afterward in the C algorithm.

        We have a formula for these vectors (see [KF]).

        WARNING ::

            We expect alpha and beta to be disjoint and to have disjoint values.
        """
        from cmath import exp, pi
        import numpy as np

        n = self._dimension
        alpha, beta = self._alpha, self._beta
        
        if not(set([E(beta[j]) for j in range(n)]).isdisjoint(set([E(alpha[j]) for j in range(n)]))):
            print "Warning the sets of eigenvalues are not disjoint, the flat bundle is not minimal"
            print alpha, beta
            raise ValueError

        if self._zero_diagonalizable:            
            N = np.mat([[1/(E(beta[j]) - E(alpha[i])) for j in range(n)] for i in range(n)]) # list of lines : M_(i,j) = M[i][j]
            w = [(np.mat([1]*n)*N.I)[0,k] for k in range(n)]
            
            return [E(a) for a in alpha], w

        else:
            N = np.mat([[E(-beta[i])/(E(beta[j]) - E(alpha[i])) for j in range(n)] for i in range(n)]) # list of lines : M_(i,j) = M[i][j]
            w = [(np.mat([1]*n)*N.I)[0,k] for k in range(n)]
            
            return [E(-b) for b in beta], w

    def monodromy_matrices(self):
        r"""
        Return monodromy matrices of the hypergeometric equation.

        OUTPUT::
            -- M0 - Monodromy around 0
            -- M1 - Monodromy around 1
            -- MInfty - Monodromy around infinity
        """
        import sage.all 
        from sage.matrix.special import identity_matrix, diagonal_matrix, column_matrix
        from sage.rings.complex_double import CDF

        if self._zero_diagonalizable:
            e_alpha, w = self.compute_monodromy()
            M0 = diagonal_matrix(CDF,e_alpha)
            M1 = identity_matrix(CDF,self._dimension) + M0.inverse()*column_matrix(CDF,[1]*self._dimension)*column_matrix(CDF,w).transpose()
            MInfty = (M0*M1).inverse()

        else:
            e_beta, w = self.compute_monodromy()
            MInfty = diagonal_matrix(CDF, e_beta)
            M1 = identity_matrix(CDF,self._dimension) + column_matrix(CDF,[1]*self._dimension)*column_matrix(CDF,w).transpose()
            M0 = (M1*MInfty).inverse()

        return (M0, M1, MInfty)

    def herm_signature(self):
        from sage.all import exp, vector, matrix
        from sage.misc.functional import numerical_approx

        n = self._dimension

        alpha_exp = [numerical_approx(exp(2*1j*pi*self._alpha[i]),30) for i in xrange(self._dimension)]
        beta_exp = [numerical_approx(exp(2*1j*pi*self._beta[i]),30) for i in xrange(self._dimension)]
        v = vector([1]*n)*matrix(n,n, lambda i, j: 1/(alpha_exp[j] + beta_exp[i])).inverse()
        d_val = [1/(-beta_exp[i]/v[i]).real_part() for i in xrange(self._dimension)]
        d_val.sort()
    
        k = 0
        while (k < n) and (d_val[k] < 0):
            k += 1
        
        return(k, n-k)


    def hodge_numbers(self):
        from copy import copy
        from math import floor
        #we need to be sure that every value is positive
        alpha, beta = copy(self._alpha), copy(self._beta)
        for x in alpha: 
            if x<0: raise NameError('one value in hodge test is not positive')
        for x in beta: 
            if x<0: raise NameError('one value in hodge test is not positive')
        alpha.append(float('Inf')), beta.append(float('Inf'))
        alpha.sort(), beta.sort()

        haar = [0]*2*self._dimension
        i, j = 0, 0
        k, switch = (self._dimension, 0) if alpha[0]>beta[0] else (self._dimension-1, 1)

        while (i < self._dimension) or (j < self._dimension):
        # We are aranging alpha and beta sorted lists. We know beta[0] == 0 is the first to be set.
        # Then suppose we have set beta[i] if switch == 0 and alpha[j] if switch == 0,
        # We add some conditions on the limit.

            haar[k] += 1

            if not switch :
                if beta[i+1] > alpha[j]:
                    switch = 1
                    i += 1
                    # place alpha[j]
                else:
                    i += 1
                    k += 1
                    # place beta[i]
            else:
                if alpha[j+1] > beta[i]:
                    switch = 0
                    j += 1
                    # place beta[i]
                else:
                    j += 1
                    k -= 1
                    # place alpha[j]
    
        s = 1
        res = 0
        for x in haar:
            res += s * x/2
            s = -s

        return [haar[k]/2 for k in xrange(2*self._dimension)]




def mean_and_std_dev(l):
    r"""
    Return the mean and standard deviation of the floatting point numbers in
    the list l.

    The implementation is very naive and should not be used for large list
    (>1000) of numbers.

    .. NOTE::

    mean and std are implemented in Sage but are quite buggy!
    """
    from math import sqrt
    m = sum(l) / len(l)
    if len(l) == 1:
        d = 0
    else:
        d = sum((x-m)**2 for x in l) / (len(l)-1)
    return m,sqrt(d)




class Problem(object):
    r"""
    Has attributes :
    
    - ``dim``


    """
    def __repr__(self):
        return "Problem in dimension " + str(self.dim)



class TorusPlanarSection(Problem):
    r"""
    - ``_zones`` -- list of (test for a zone, name for it)

    - ``dim`` -- dimention of the studied flat bundle
 
    - ``section_f`` -- function that take two real parameter and map to the section of the torus

    - ``_xmin`` -- Boundary of the studied domain 

    - ``_xmax`` -- Boudary of the studied domain

    - ``_ymin`` -- Boundary of the studied domain 

    - ``_ymax`` -- Boudary of the studied domain
    """
    def __copy__(self):
        r"""
        Returns a copy of self.

        EXAMPLES::

            sage: p = Experiment([0, 0.5], [0.2, 0.3])
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
        """
        from copy import copy
        q = self.__class__()
        q.__hash__ = self.__hash__
        q.dim = self.dim
        q.section_f = self.section_f
        return q

    def __init__(self, dim, section_f, section_name, xmin=None, xmax=None, ymin=None, ymax=None):
        self.dim = dim
        self.section_f = section_f
        self.section_name = section_name
        self._xmin, self._ymin, self._xmax, self._ymax = xmin, ymin, xmax, ymax

    def zone(self, zone_test=(lambda x,y: True), zone_name='all', xmin=None, xmax=None, ymin=None, ymax=None):
        if not type(zone_name)==str:
            raise NameError('the name of the zone must be a str')
        if xmin and xmax and ymin and ymax:
            return TorSecZone(self.dim, self.section_f, self.section_name, zone_test, zone_name, xmin, xmax, ymin, ymax)
        else:
            return TorSecZone(self.dim, self.section_f, self.section_name, zone_test, zone_name, self._xmin, self._xmax, self._ymin, self._ymax)

    def zone_list(self, zones, xmin=None, xmax=None, ymin=None, ymax=None):
        return [self.zone(z[0], z[1], xmin, xmax, ymin, ymax) for z in zones]

    def f_name(self, nsteps_x, nsteps_y, nb_iterations, nb_experiments):
        res = 'res/k' + str(self.dim) + '_' + self.section_name
        res += '_(%0.1e<%0.1e:%s),(%0.1e<%0.1e:%s)'%(float(self._xmin),float(self._xmax),nsteps_x,float(self._ymin),float(self._ymax),nsteps_y) 
        res += '_it%0.1e_exp%s'%(nb_iterations, nb_experiments)
        return(res)

    def compute_discretized(self, nsteps_x, nsteps_y, nb_iterations=10**6, nb_experiments=10):
        r"""
        Compute discretized values of the given torus section.

        INPUT::
            - nsteps_x -- 
            - nsteps_y --
            - nb_iterations --
            - nb_experiments --

        OUTPUT::
            - None
        """
        import os.path, csv
        from csv import DictReader, DictWriter

        x_list = lin_space(self._xmin, self._xmax, nsteps_x)
        y_list = lin_space(self._ymin, self._ymax, nsteps_y)

        tab = [[self.section_f(x,y) for y in y_list] for x in x_list]
        fieldnames = ['Exp', 'i', 'j', 'x_plot', 'y_plot', 'res', 'std', 'sum', 'sum_std']

        f_name = self.f_name(nsteps_x, nsteps_y, nb_iterations, nb_experiments)

        if os.path.isfile(f_name) :
            with open(f_name, 'r') as csvfile :
                reader= csv.DictReader(csvfile, delimiter=',')
                test = 1
                for row in reader :
                    test = 0
                if test:
                    last_i, last_j = 0, 0
                else:
                    last_i, last_j = int(row['i']), int(row['j'])

        else :
            with open(f_name, 'w') as csvfile :
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                last_i, last_j = 0, 0

        with open(f_name, 'a') as csvfile :
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            if not ((last_i == nsteps_x-1) and (last_i == nsteps_y-1)):
                for i in range(last_i,nsteps_x):
                    for j in range(last_j+1 if i == last_i else 0, nsteps_y):
                        print i,j
                        r, d, s_r, s_d = tab[i][j].hypergeometric_lyap_exp(verbose = False, nb_iterations=nb_iterations, nb_experiments=nb_experiments,return_error=True)
                        writer.writerow({'Exp' : tab[i][j], 'i' : i, 'j' : j, 
                                         'x_plot' : tab[i][j].x_plot, 'y_plot' : tab[i][j].y_plot,
                                         'res' : r, 'std' : d, 'sum' : s_r, 'sum_std' : s_d})

        return True

    def test_monodromy(self, nsteps_x, nsteps_y, save=False):
        import os.path, csv
        from csv import DictReader, DictWriter
        import sage.all 
        from sage.matrix.special import identity_matrix
        from sage.rings.complex_double import CDF
        from sage.symbolic.constants import e

        x_list = lin_space(self._xmin, self._xmax, nsteps_x)
        y_list = lin_space(self._ymin, self._ymax, nsteps_y)

        tab = [[self.section_f(x,y) for y in y_list] for x in x_list]
        fieldnames = ['Exp', 'i', 'j', 'x_plot', 'y_plot', 'dist_id', 'dist_ev']

        f_name = self.f_name(nsteps_x, nsteps_y, -1, -1)

        if os.path.isfile(f_name) :
            with open(f_name, 'r') as csvfile :
                reader= csv.DictReader(csvfile, delimiter=',')
                test = 1
                for row in reader :
                    test = 0
                if test:
                    last_i, last_j = 0, 0
                else:
                    last_i, last_j = int(row['i']), int(row['j'])

        else :
            with open(f_name, 'w') as csvfile :
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                last_i, last_j = 0, 0

        with open(f_name, 'a') as csvfile :
            from cmath import exp, pi

            from __builtin__ import set
            from sage.misc.functional import numerical_approx
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            if not ((last_i == nsteps_x-1) and (last_i == nsteps_y-1)):
                for i in range(last_i,nsteps_x):
                    for j in range(last_j+1 if i == last_i else 0, nsteps_y):
                        print i,j
                        M0, M1, MInfty = tab[i][j].monodromy_matrices()
                        n = tab[i][j]._dimension
                        beta = tab[i][j]._beta
                        dist_id = (M0*M1*MInfty - identity_matrix(CDF,n)).norm()
                        ev = set(MInfty.eigenvalues())
                        d = []
                        for b in beta:
                            val = exp(-2*1j*pi*b)
                            c_val = min(ev, key=lambda z: (z-val).norm())
                            d.append((c_val-val).norm())
                            ev.remove(c_val)
                        dist_ev = max(d)
                        print max(dist_ev, dist_id)
                        writer.writerow({'Exp' : tab[i][j], 'i' : i, 'j' : j, 
                                         'x_plot' : tab[i][j].x_plot, 'y_plot' : tab[i][j].y_plot,
                                         'dist_id' : dist_id, 'dist_ev' : dist_ev})

        return True


#    def slope(self, n_zone, nb_tests=5, nb_points, **args): 


class TorSecZone(TorusPlanarSection):
    def __copy__(self):
        from copy import copy
        q = self.__class__()
        q.__hash__ = self.__hash__
        q.dim = self.dim
        q.section_f = self.section_f
        q._xmin = self._xmin
        q._ymin = self._ymin
        q._xmax = self._xmax
        q._ymax = self._ymax
        q._zone_test= self._zone_test
        q._zone_name = self._zone_name

    def __init__(self, dim, section_f, section_name, zone_test=(lambda x,y: True), zone_name='all', xmin=None, xmax=None, ymin=None, ymax=None):
        self.dim = dim
        self.section_f = section_f
        self.section_name = section_name
        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax
        self._zone_test = zone_test
        self._zone_name = zone_name

    def rand(self):
        if not (self._xmin and self._xmax and self._ymin and self._ymax):
            raise NameError('Please bound the domain to define a rand experiment')
        from random import uniform
        for i in range(1,1000):
            x = uniform(self._xmin, self._xmax)
            y = uniform(self._ymin, self._ymax)
            if self._zone_test(x,y): return self.section_f(x,y)
        raise NameError('Test Zone too small, cannot find random values')

    def random_points(self, nb_points=5, nb_iterations=10**4, nb_experiments=10, plot=False, reg=None):
        from os import path 
        from csv import DictReader, DictWriter
        print self._zone_name
        f_name = 'res/random_d' + str(self.dim) + '_' + self._zone_name + '_' + str(nb_iterations) + '_' +str(nb_experiments)
        fieldnames = ['Exp', 'i', 'x_plot', 'y_plot', 'res', 'std', 'sum', 'sum_std']

        if path.isfile(f_name) :
            with open(f_name, 'r') as csvfile :
                reader= DictReader(csvfile, delimiter=',')
                test = 1
                for row in reader :
                    test = 0
                if test:
                    last_i = 0
                else:
                    last_i = int(row['i'])

        else :
            with open(f_name, 'w') as csvfile :
                writer = DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                last_i = 0

        with open(f_name, 'a') as csvfile :
            writer = DictWriter(csvfile, fieldnames=fieldnames)
            for i in range(last_i, nb_points) :
                print i
                exp = self.rand()
                r, d, s_r, s_d = exp.hypergeometric_lyap_exp(verbose = False, nb_iterations=nb_iterations, nb_experiments=nb_experiments)
                
                writer.writerow({'Exp' : exp, 'i' : i+1,
                                 'x_plot' : exp.x_plot, 'y_plot' : exp.y_plot,
                                 'res' : r, 'std' : d, 'sum' : s_r, 'sum_std' : s_d})

        res = []
        with open(f_name, 'r') as csvfile :
            reader= DictReader(csvfile, delimiter=',')
            for row in reader :
                res.append((float(row['x_plot']), float(row['y_plot']), float(row['sum']), float(row['sum_std'])))

        if reg or plot:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = Axes3D(fig)

            aux = [(x,y,z-(x*reg[0] + y*reg[1]) if reg else z,ds) for (x,y,z,ds) in res]

            for [x, y, z, err] in aux:
                ax.plot([x, x],[y, y],zs=[z+err, z-err], color='lightsalmon')
            ax.scatter([val[0] for val in aux], [val[1] for val in aux], zs=[val[2] for val in aux], color='brown')


            plt.show()
            f_name = 'fig/random_d' + str(self.dim) + '_' + self._zone_name + '_' + str(nb_iterations) + '_' +str(nb_experiments) + '.png'
            plt.savefig(f_name)

        return res

    def plot_discretized(self, nsteps_x, nsteps_y, nb_iterations=10**6, nb_experiments=10, plot=False, color='brown', save=False, reg=[]):
        r"""
        Return a table of points with experiments.

        The computation need to be done with an other function compute_discretized.
        The result of this function is stored in a file which is read when plotting.
        The parameters are use to call the file, since the file name is choosen depending
        of them.

        INPUT::
            - nsteps_x -- table of parameters with format Experiments.
            - nsteps_y --
            - nb_iterations --
            - nb_experiments --
            - plot -- option for the type of plotting we want. Should be in the list - 'scatter', '2d', 'animate'.
            - color -- color of points
            - save -- whether save the plotting result or show it.
            - reg -- enter linear regression parameter to plot distance from the closest plane.

        OUTPUT::
            - Table of points
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        import os.path, csv
        from csv import DictReader, DictWriter

        f_name = self.f_name(nsteps_x, nsteps_y, nb_iterations, nb_experiments)

        if os.path.isfile(f_name) :
            with open(f_name, 'r') as csvfile :
                reader= csv.DictReader(csvfile, delimiter=',')
                test = 1
                for row in reader :
                    test = 0
                if test:
                    last_i, last_j = 0, 0
                else:
                    last_i, last_j = int(row['i']), int(row['j'])

        else :
            last_i, last_j = 0, 0

        if last_i < nsteps_x-1 or last_j < nsteps_y-1 :
            print f_name
            raise NameError('You should compute the section before drawing a domain in it')

        res = []
        with open(f_name, 'r') as csvfile :
            reader= DictReader(csvfile, delimiter=',')
            for row in reader :
                x, y  = float(row['x_plot']), float(row['y_plot'])
                i, j = int(row['i']), int(row['j'])
                z, sd = float(row['sum']), float(row['sum_std'])
                # Select the zone we are ploting
                if self._zone_test(x,y):
                    res.append((x,y,i,j,z,sd))

        # Applying regression, choosing parameters we need
        if len(reg) == 2:
            res = [(x,y,i,j,z-(x*reg[0] + y*reg[1]) if reg else z,ds) for (x,y,i,j,z,ds) in res]
        elif len(reg) == 3:
            res = [(x,y,i,j,z-(x*reg[0] + y*reg[1] + reg[2]) if reg else z,ds) for (x,y,i,j,z,ds) in res]

        aux = [(x,y,z,ds) for (x,y,i,j,z,ds) in res]


        if plot == 'scatter':
            fig = plt.figure()
            ax = Axes3D(fig)
            for [x, y, z, err] in aux:
                ax.plot([x, x],[y, y],zs=[z+err, z-err], color='lightsalmon')
            ax.scatter([val[0] for val in aux], [val[1] for val in aux], zs=[val[2] for val in aux], color=color)

            if save: plt.savefig('fig/' + f_name[4:] + '.png')
            else: plt.show()

        if plot == '2d':            
            x_list = lin_space(self._xmin, self._xmax, nsteps_x)
            y_list = lin_space(self._ymin, self._ymax, nsteps_y)
            X, Y = meshgrid(x_list, y_list)
            Z = [[0. for _ in y_list] for _ in x_list]
            for (x,y,i,j,z,sd) in res:
                Z[i][j] = z

            A = array(Z)
            fig, ax = plt.subplots()
            p = ax.pcolor(X, Y, A.T, cmap=cm.RdBu)
            cb = fig.colorbar(p, ax=ax)

            plt.xlabel('r')
            plt.ylabel('x')

            if save: plt.savefig('fig/'  + f_name[4:] + '_flat' + '_' + self._zone_name + '.png')
            else:
                plt.axis([x_list[0], x_list[-1], y_list[0], y_list[-1]])
                plt.show()

        if plot == 'animate':
            from matplotlib import animation
            def init():
                for [x, y, z, err] in aux:
                    ax.plot([x, x],[y, y],zs=[z+err, z-err], color='lightsalmon')
                    #ax.scatter([val[0] for val in aux], [val[1] for val in aux], zs=[val[2] for val in aux], color='brown')
            def animate(i):
                ax.view_init(elev=10., azim=i)

            #Animate
            anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
            anim.save('video/' + f_name[4:] + '_' + self._zone_name + '_anim.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        return res

def lyap_exp_CY(C, d, nb_vectors=None, nb_experiments=10, nb_iterations=10**4, verbose=False, output_file=None, return_error=False):
        r"""
        Compute the Lyapunov exponents of the geodesic flow in the hypergeometric function
        space.
        
        INPUT:

        - ``nb_vectors`` -- the number of vectors to use
        
        - ``nb_experiments`` -- number of experimets
        
        - ``nb_iterations`` -- the number of iterations of the induction to perform
        
        - ``output_file`` -- place where we print the results
        
        - ``verbose`` -- do we print the result with an extensive number of information or not
    
        """
        import time
        import lyapunov_exponents    # the cython bindings
        from math import sqrt
        
        if output_file is None:
            from sys import stdout
            output_file = stdout
            closed = True
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")
            closed = False

        if nb_vectors <> None and nb_vectors <= 0:
            raise ValueError("the number of vectors must be positive")
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        #recall that the lyapunov exponents are symmetric
        if nb_vectors == None:
            nb_vectors = 4

        t0 = time.time()
        res = lyapunov_exponents.lyapunov_exponents([0]*4, [0]*4, 4, nb_vectors, nb_experiments, nb_iterations, [C, d])
        t1 = time.time()

        res_final = []
        std_final = []
        s_m, s_d = 0, 0

        if verbose:
            from math import floor, log
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2**%d)\n"%(nb_iterations, floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
        for i in xrange(nb_vectors):
            m,d = mean_and_std_dev(res[i])
            s_m += m
            s_d += d**2
            if verbose:
                output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    i,m,d, 2.576*d/sqrt(nb_experiments)))
            res_final.append(m)
            std_final.append(2.576*d/sqrt(nb_experiments))

        s_d = sqrt(s_d)
        s_d_final = 2.576*s_d/sqrt(nb_experiments)
        if verbose:
            output_file.write("sum_theta        : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n\n"%(
                s_m,s_d, 2.576*s_d/sqrt(nb_experiments)))
        
        if not closed :
            output_file.close()
            print "file closed"

        if return_error:
            return (res_final, std_final, s_m, s_d_final)
        else:
            return res_final
