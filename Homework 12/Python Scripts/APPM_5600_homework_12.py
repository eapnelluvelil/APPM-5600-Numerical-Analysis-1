import numpy as np

def main():
	g = lambda x: x*np.log(x)
	a_g = 1e-16;
	b_g = 1;
	n = 512;

	comp_simp_g = -4*comp_simp(a_g, b_g, g, n);
	print(comp_simp_g)


def comp_simp(a, b, f, n):
	nodes = np.linspace(a, b, n+1);
	h = (b-a)/n;
	f_nodes = f(nodes);
	f_nodes[1:2:n-1] = 4*f_nodes[1:2:n-1];
	f_nodes[2:2:n-2] = 2*f_nodes[2:2:n-2];
	return (h/3)*np.sum(f_nodes)

main()