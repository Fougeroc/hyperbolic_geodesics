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
reg_s = {'a': [[] for _ in range(8)], 'b': [[] for _ in range(8)], 'c': [[] for _ in range(8)]}
print reg_s['a']

for s in f_range(.5,.9,.02):
    t = TorusPlanarSection(2, shift(s), '0,x_r,2r+%.2f'%(s), d_min, d_max/2, d_min, d_max)

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

    for k in [2,4,5]:
    	z = zones[k]
        print z.section_name
        print z._zone_name
        aux = z.plot_discretized(100,100, nb_iterations=10**4, plot=False)
	res = [[float(x),float(y),float(z)] for (x,y,i,j,z,ds) in aux]
	var('a,b,c')
        model(x,y) = a*x + b*y + c
	if res:
	    f = find_fit(res, model, solution_dict=True)
	    for x in [a,b,c]:
	        (reg_s[str(x)][k]).append((f[x],s))
	else:
	    print "EMPTY"

var('s')
model(s) = a*s + c
for i in range(3):
    for x in [a,b,c]:
    	print i, x
	print reg_s[str(x)][i]
	if reg_s[str(x)][i]:
            f = find_fit(reg_s[str(x)][i], model, solution_dict=True)
	else:
	    print "EMPTY"
	print f
