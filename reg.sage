from template import Experiment, TorusPlanarSection, mean_and_std_dev

def reg():
    data = [(i, 1.2 * sin(0.5*i-0.2) + 0.1 * normalvariate(0, 1)) for i in xsrange(0, 4*pi, 0.2)]
    var('a, b, c, x')
    model(x) = a * sin(b * x - c)
    f = find_fit(data, model)
    return f

def linspace(a, b, s):
    delta = (b-a)/(s-1.)
    return [a + delta*k for k in range(s)]

n_step = 100
d_min = 1./(10*n_step)
d_max = 1./2 - d_min


def load(f_name, X_init, Y_init, zone_test):
    import csv
    from matplotlib import animation
    from mpl_toolkits.mplot3d import Axes3D

    data = []

    with open(f_name, 'r') as csvfile :
        reader= csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            x_data = X_init[int(row['i'])]
            y_data = Y_init[int(row['j'])]
            if zone_test(x_data, y_data):
                data.append((x_data, y_data, float(row['sum'])))
        var('x,y,a,b')
        model(x,y) = a*x + b*y
        f = find_fit(data, model, solution_dict=True)

        dz = 0
        for (x_data, y_data, z) in data:
            dz += (f[a]*x_data + f[b]*y_data - z)**2

    return (f, sqrt(dz/len(data))*2, dz, len(data))



def load_slope(tor_sec_zone, nb_tests, nb_points=5, **args):
    import csv
    from math import sqrt

    data = [d[0:3] for d in tor_sec_zone.random_points(nb_points*nb_tests, **args)]
    a_res, b_res = [], []

    for i in range(nb_tests):
        var('x,y,a,b')
        model(x,y) = a*x + b*y
        f = find_fit(data[nb_points*i: nb_points*(i+1)], model, solution_dict=True)
	a_res.append(f[a]), b_res.append(f[b])

    m_a, s_a = mean_and_std_dev(a_res)
    m_b, s_b = mean_and_std_dev(b_res)
    return (m_a, 2.576*s_a/sqrt(nb_tests), m_b, 2.576*s_a/sqrt(nb_tests))

def reg_slope(tor_sec_zone, nsteps_x, nsteps_y, nb_iterations, nb_experiments):
    import csv
    from math import sqrt

    data = []
    f_name = tor_sec_zone.f_name(nsteps_x, nsteps_y, nb_iterations, nb_experiments)


    with open(f_name, 'r') as csvfile :
        reader= csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            x_data = float(row['x_plot'])
            y_data = float(row['y_plot'])
            if tor_sec_zone._zone_test(x_data, y_data):
                data.append((x_data, y_data, float(row['sum'])))
        var('x,y,a,b')
        model(x,y) = a*x + b*y
        f = find_fit(data, model, solution_dict=True)

        dz = 0
        for (x_data, y_data, z) in data:
            dz += (f[a]*x_data + f[b]*y_data - z)**2

    return (f, sqrt(dz/len(data))*2, dz, len(data))


def reg_slope_aff(tor_sec_zone, nsteps_x, nsteps_y, nb_iterations, nb_experiments):
    r"""
    Look for the closest plane approximating the function in the given zone.
    
    RETURN::
        - A dictionary with three entries: 'a', 'b', 'c'
        where the plane is of equation a*x + b*y + c*x
    """
    import csv
    from math import sqrt

    data = []
    f_name = tor_sec_zone.f_name(nsteps_x, nsteps_y, nb_iterations, nb_experiments)


    with open(f_name, 'r') as csvfile :
        reader= csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            x_data = float(row['x_plot'])
            y_data = float(row['y_plot'])
            if tor_sec_zone._zone_test(x_data, y_data):
                data.append((x_data, y_data, float(row['sum'])))
        var('x,y,a,b,c')
        model(x,y) = a*x + b*y + c
        f = find_fit(data, model, solution_dict=True)

        dz = 0
        for (x_data, y_data, z) in data:
            dz += (f[a]*x_data + f[b]*y_data + f[c] - z)**2

    return (f, sqrt(dz/len(data))*2, dz, len(data))



n_step = 100
d_min = 1./(10*n_step)
d_max = 1./2 - d_min

def exp(rel, dist):
    return Experiment([0., -dist, -2*dist], [rel, dist + rel, dist + 2*rel], rel, dist)

zone_1 = (lambda r, x: r>2*x, 'r>2x')
zone_2 = (lambda r, x: r<2*x and r>x, 'r<2x,r>x')
#zone_3 is alternate, then zero lyapunov exponents
zone_3 = (lambda r, x: r<x and 2*r>x, 'alternate')
zone_4 = (lambda r, x: 2*r<x and 4*r>x, '2r<x,4r>x')
zone_5 = (lambda r, x: 4*r<x, '4r<x')

zones = [zone_1, zone_2, zone_3, zone_4, zone_5]

t = TorusPlanarSection(3, exp, '0,x,2x_r,x+r,x+2r', d_min, d_max/2, d_min, d_max)
zones = t.zone_list([zone_1, zone_2, zone_3, zone_4, zone_5])

#print load_slope(zones[1], nb_tests=10, nb_points=5, nb_iterations=10**8, nb_experiments=10)

# print "k3"
# print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r,x: r>2*x)
# print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r,x: r<2*x and r>x)
# print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<x and 2*r>x)
# #print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
# #           lambda r, x: 2*r<x and 3*r>x)
# #print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
# #           lambda r, x: 3*r<x and 4*r>x)
# print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 2*r<x and 4*r>x)
# print load('res/k3_0,x,2x_r,x+r,x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r,x: 4*r<x)


# print "k4 un"
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r>3*x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<3*x and r>2*x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<2*x and r<x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<x and 3*r>2*x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 3*r<2*x and 2*r>x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 2*r<x and 7*r>3*x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 7*r<3*x and 3*r>x)
# print load('res/k4_0,x,2x,3x_r,x+r,x+2r,x+3r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 3*r<x)

# print "k4 deux"
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r>3*x)
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<3*x and r>2*x)
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<2*x and r>x)
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: r<x and 2*r>x)
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 2*r<x and 5*r>x)
# print load('res/k4_0,x,2x,3x_r,x+r,2x+r,2x+2r', linspace(d_min,(d_max-d_min)/2,n_step), linspace(d_min,d_max,n_step), 
#            lambda r, x: 5*r<x)



n_step = 100
d_min = 0.01
d_max = 0.99

def exp(rel, dist):
    return Experiment([0., -dist], [rel, 2*rel], rel, dist)

t = TorusPlanarSection(2, exp, '0,x_r,2r', d_min, d_max/2, d_min, d_max)
#t.zone().discretized_values(100,100, nb_iterations=10**7, plot=True, plot_2d=True)

zone_0 = (lambda r, x: r>x and 3*r<1+x, 'r>x, 3r<1+x')
zone_1 = (lambda r, x: r>x and 3*r>1+x, 'r>x, 3r>1+x')
zone_2 = (lambda r, x: r<x and 2*r>x, 'alternate')
zone_3 = (lambda r, x: 2*r<x and 3*r>x, '2r<x,3r>x')
zone_4 = (lambda r, x: 3*r<x, '3r<x')

zones = t.zone_list([zone_0, zone_1, zone_2, zone_3, zone_4])
for z in zones:
    print z._zone_name
    #print reg_slope(z, 100, 100, nb_iterations=10**7, nb_experiments=10)

for i in range(5):
    print i
    print reg_slope_aff(zones[i],100, 100, nb_iterations=10**5, nb_experiments=10)
