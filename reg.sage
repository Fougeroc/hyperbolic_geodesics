from template import Experiment, TorusPlanarSection, lin_space

#############################################################################################
#                                    EXPERIMENT 1                                           #
#############################################################################################

# Dimension 2
# Tomography

n_step = 100
d_min = 0.01
d_max = 0.99

def linspace(a, b, s):
    delta = (b-a)/(s-1.)
    return [a + delta*k for k in range(s)]

def f_range(a,b,s):
    res = []
    aux = a
    while aux < b:
        res.append(aux)
	aux+=s
    return(res)

def shift(s):
    def section_fun(rel, dist):
        return Experiment([rel, s + 2*rel], [0., dist], rel, dist)
    return(section_fun)


reg_tot=[]
for s in f_range(.0,.2,.1):
    t = TorusPlanarSection(2, shift(s), '0,x_r,2r+%.1f'%(s), d_min, d_max/2, d_min, d_max)
    #t.compute_discretized(100, 100, nb_iterations=10**5)

    all = t.zone()

    zone = []
    zone.append((lambda r, x: r>x and 3*r+s<1+x, 'r>x, 3r+s<1+x'))
    zone.append((lambda r, x: 3*r+s>1+x and 2*r<1-s, '3r+s>1+x, 2r>1-s'))
    zone.append((lambda r, x: x<2*r+s-1, 'x<2r+s-1'))
    zone.append((lambda r, x: r>x and 2*r<1-s and x<2*r+s-1, 'r>x, 2r<1-s, x<2r+s-1'))
    zone.append((lambda r, x: 3*r+s>1+x and x>r, '3r+s>1+x, x>r'))
    zone.append((lambda r, x: 2*r>1-s and 3*r+s<1+x, '2r>1-s, 3r+s<1+x'))
    zone.append((lambda r, x: x>2*r+s and 3*r+s>x, 'x>2x+s, 3r+s>x'))
    zone.append((lambda r, x: 3*r+s<x, '3r+s<x'))

    zones = t.zone_list(zone)

    reg = []
	
    for z in zones:
        print z.section_name
        print z._zone_name
        aux = z.plot_discretized(100,100, nb_iterations=10**5, plot=False)
	res = [[float(x),float(y),float(z)] for (x,y,i,j,z,ds) in aux]
	var('a,b,c')
        model(x,y) = a*x + b*y + c
	if res:
	    f = find_fit(res, model, solution_dict=True)
	    print f
	    reg.append([f[a],f[b],f[c]])
	else:
	    print "EMPTY"
	    reg.append([0,0,0])

    reg_tot.append(reg)

print reg_tot