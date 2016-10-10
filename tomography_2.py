from template import Experiment, TorusPlanarSection, lin_space

from pylab import *
from matplotlib import *

#############################################################################################
#                                    EXPERIMENT 1                                           #
#############################################################################################

# Dimension 2
# Tomography

n_step = 100
d_min = 0.01
d_max = 0.99

def shift(s):
    def section_fun(rel, dist):
        return Experiment([rel, s + 2*rel], [0., dist], rel, dist)
    return(section_fun)

for s in np.arange(0,1,.1):
    t = TorusPlanarSection(2, shift(s), '0,x_r,2r+' + str(s), d_min, d_max/2, d_min, d_max)
    print t.section_name
    t.compute_discretized(100, 100, nb_iterations=10**5)

for s in np.arange(.5,1,.02):
    t = TorusPlanarSection(2, shift(s), '0,x_r,2r+%.2f'%s, d_min, d_max/2, d_min, d_max)
    print t.section_name
    t.compute_discretized(100, 100, nb_iterations=10**4)
