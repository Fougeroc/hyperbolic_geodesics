from plot_exp import big_plot, plot_show, plot_show_, compute_and_plot, random_points
from template import Experiment, TorusPlanarSection, lin_space
from random import uniform

from pylab import *
from matplotlib import *



#############################################################################################
#                                     BABY CASE                                             #
#############################################################################################

# n_step = 10
# d_min = 1./(10*n_step)
# d_max = 1./4 - d_min

# def exp(rel, dist):
#     return Experiment([0., -dist, -2*dist], [rel, dist + rel, dist + 2*rel], rel, dist)

#slope(zone_1, 'r>2x', rand, exp, 5, nb_iterations=100)

#big_plot(e, 'res/baby_case', nb_iterations=10**4)
#plot_show_('res/baby_case')


#############################################################################################
#                                    EXPERIMENT 0                                           #
#############################################################################################

# Dimension 2

n_step = 100
d_min = 0.01
d_max = 0.99

def section_fun(rel, dist):
    return Experiment([rel, 2*rel], [0., dist], rel, dist)

t = TorusPlanarSection(2, section_fun, '0,x_r,2r', d_min, d_max/2, d_min, d_max)
#t.zone().discretized_values(100,100, nb_iterations=10**6, plot=True, plot_2d=True)
#t.zone().discretized_values(100,100, nb_iterations=10**7, plot=True, plot_2d=True)
#t.zone().discretized_values(100,100, nb_iterations=10**8, nb_experiments=5, plot=True, plot_2d=True)

#print E.x_plot
#print E.y_plot

#print E._alpha, E._beta
#print E.hypergeometric_lyap_exp()

#print Experiment([0,.3],[.5, .7]).monodromy_matrices()
#t.compute_discretized(100, 100, nb_iterations=10**5
#t.compute_discretized(100, 100, nb_iterations=10**8, nb_experiments=5)
#t.compute_discretized(100, 100, nb_iterations=10**7)
#t.test_monodromy(100,100)

# zone_0 = (lambda r, x: r>x and 3*r<1+x, 'r>x, 3r<1+x')
# zone_1 = (lambda r, x: r>x and 3*r>1+x, 'r>x, 3r>1+x')
# zone_2 = (lambda r, x: r<x and 2*r>x, 'alternate')
# zone_3 = (lambda r, x: 2*r<x and 3*r>x, '2r<x,3r>x')
# zone_4 = (lambda r, x: 3*r<x, '3r<x')

# zones = t.zone_list([zone_0, zone_1, zone_2, zone_3, zone_4])
# reg = [(2,-2,0), (-4,0,2), (0,0,0), (-4,2,0), (2,0,0)]

# all = t.zone()
# p_1 = all.plot_discretized(100,100, nb_iterations=10**5, plot='2d', save=True)
# #p_2 = all.plot_discretized(100,100, nb_iterations=10**7, plot='2d', save=True)


# i = 0
# for i in range(5):
#    zones[i].plot_discretized(100,100, nb_iterations=10**5, reg=reg[i], plot='2d', save=True)

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# x_list = lin_space(t._xmin, t._xmax, 100)
# y_list = lin_space(t._ymin, t._ymax, 100)
# X, Y = meshgrid(x_list, y_list)
# Z = [[0. for _ in y_list] for _ in x_list]
# for k in range(5):
#     for (x,y,i,j,z,sd) in zones[k].plot_discretized(100,100, nb_iterations=10**5, reg=reg[k]):
#         Z[i][j] = z

# A = array(Z)
# fig, ax = plt.subplots()
# p = ax.pcolor(X, Y, A.T, cmap=cm.RdBu)
# cb = fig.colorbar(p, ax=ax)
# plt.xlabel('r')
# plt.ylabel('x')

# plt.savefig('fig/all_reg_e6.png')


    
# from pylab import *
# from matplotlib import *

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = Axes3D(fig)


# def lin_space(a, b, s):
#     delta = (b-a)/(s-1.)
#     return [a + delta*k for k in range(s)]

# x_list = lin_space(t._xmin, t._xmax, 100)
# y_list = lin_space(t._ymin, t._ymax, 100)
# X, Y = meshgrid(x_list, y_list)
# Z = [[0. for _ in y_list] for _ in x_list]
# # for (x,y,i,j,z,sd) in p_1:
# #     Z[i][j] = z
# # for (x,y,i,j,z,sd) in p_2:
# #     Z[i][j] -= z

# A = array(Z)
# fig, ax = plt.subplots()
# p = ax.pcolor(X, Y, A.T, cmap=cm.RdBu)
# cb = fig.colorbar(p, ax=ax)

# plt.xlabel('r')
# plt.ylabel('x')

# plt.savefig('fig/test_diff.png')

#ax.scatter([val[0] for val in p_1], [val[1] for val in p_1], zs=[val[2] for val in p_1], color='brown')
#ax.scatter([val[0] for val in p_2], [val[1] for val in p_2], zs=[val[2] for val in p_2], color='red')

#plt.show()

#all = t.zone()

#for i in range(len(zones)):
    #print zones[i].section_name
    #print zones[i]._zone_name
    #zones[i].discretized_values(100, 100, nb_iterations=10**7, plot=True)
    #all.plot_discretized(100,100,plot=True,nb_iterations=10**7)
    #zones[i].plot_discretized(100,100, nb_iterations=10**7, reg=reg[i])

#############################################################################################
#                                    EXPERIMENT 1                                           #
#############################################################################################

#Dimension 3

n_step = 100
d_min = 1./(10*n_step)
d_max = 1./2 - d_min

def section_fun(rel, dist):
    return Experiment([rel, dist+rel, dist+2*rel], [0., dist, 2*dist], rel, dist)

# zone_0 = (lambda r, x: r>2*x, 'r>2x')
# zone_1 = (lambda r, x: r<2*x and r>x, 'r<2x,r>x')
# #zone_2 is alternate, then zero lyapunov exponents
# zone_2 = (lambda r, x: r<x and 2*r>x, 'alternate')
# zone_3 = (lambda r, x: 2*r<x and 4*r>x, '2r<x,4r>x')
# zone_4 = (lambda r, x: 4*r<x, '4r<x')

T = TorusPlanarSection(3, section_fun, 'r,d+r,d+2r_0,d,2d', d_min, d_max/2, d_min, d_max)
t.compute_discretized(100, 100, nb_iterations=10**5)
#t.zone().discretized_values(100,100, nb_iterations=10**7, reg=reg[i])
#t.zone().discretized_values(100,100, nb_iterations=10**7,plot=True, plot_2d=True)
zones = t.zone_list([zone_0, zone_1, zone_2, zone_3, zone_4])

all = t.zone()
all.plot_discretized(100,100,nb_iterations=10**5, plot='2d', save=True)

#for i in range(len(zones)):
#    print zones[i]._zone_name
#    if i != 2: zones[i].plot_discretized(100,100, plot=True, nb_iterations=10**7, reg=reg[i])


# print str(e).__hash__()
# big_plot(e, 'res/k3_0,x,2x_r,x+r,x+2r')
# compute_and_plot(e, linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), f_name='k3_0,x,2x_r,x+r,x+2r', dim_rep=3, rep_type='points')

# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', suffixe='all', mode='show')
# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', suffixe='all')

# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', zone_1, reg=(3.94,-4.446), suffixe='r>2x_line')
# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', zone_2, reg=(3.446,-3.452), suffixe='r<2x,r>x_line')
# #plot_show_('res/k3_0,x,2x_r,x+r,x+2r', zone_3, reg=(0,0), suffixe='r<x,2r>x_line')
# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', zone_4, reg=(-3.5, 1.75), suffixe='2r<x,3r>x_line')
# plot_show_('res/k3_0,x,2x_r,x+r,x+2r', zone_5, reg=(3.5, 0), suffixe='4r<x_line')


#############################################################################################
#                                    EXPERIMENT 2                                           #
#############################################################################################

# Dimension 4
# First type of alternance

# n_step = 100
# d_min = 1./(10*n_step)
# d_max = 1./3 - d_min

#e = [[Experiment([0., -dist, -2*dist, -3*dist],[rel, dist + rel, dist + 2*rel, dist + 3*rel]) 

#print str(e).__hash__()
#big_plot(e)
#compute_and_plot(e, linspace(d_min, (d_max-d_min)/2, n_step), linspace(d_min, d_max, n_step), f_name='k4_0,x,2x,3x_r,x+r,x+2r,x+3r', dim_rep=3, rep_type='points')

#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: r>x)
#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: r<x and 3*r>2*x)
#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 3*r<2*x and 2*r>x)
#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 2*r<x and 7*r>3*x)
#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 7*r<3*x and 3*r>x)
#plot_show_('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 3*r<x)
                

#############################################################################################
#                                    EXPERIMENT 2'                                          #
#############################################################################################

# Dimension 4
# Second type of alternance


# n_step = 100
# d_min = 1./(10*n_step)
# d_max = 1./3 - d_min

#e = [[Experiment([0., -dist, -2*dist, -3*dist],[rel, dist + rel, 2*dist + rel, 2*dist + 2*rel]) 
#      for dist in linspace(d_min, d_max, n_step)] 
#     for rel in linspace(d_min, (d_max-d_min)/2, n_step)]

#print str(e).__hash__()
#big_plot(e)
#plot_show(e, linspace(d_min, (d_max-d_min)/2, n_step), linspace(d_min, d_max, n_step), dim_rep=3)

#compute_and_plot(e, linspace(d_min, (d_max-d_min)/2, n_step), linspace(d_min, d_max, n_step), f_name='k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', dim_rep=2, rep_type='points')

# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: r>x)
# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: r<x and 2*r>x)
# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 2*r<x and 5*r>x)
# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 5*r<x)
# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 7*r<3*x and 3*r>x)
# plot_show_('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), lambda r, x: 3*r<x)
