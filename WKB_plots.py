import matplotlib.pyplot as plt

from WKB_Disc import WKB_Disc
def sort_plotting_lst(r, k):
	r_new, k_new = [], []

	while len(k) != 0:
		i = k.index(min(k))
		r_new.append(r.pop(i))
		k_new.append(k.pop(i))

	return r_new, k_new
def plot_Dispersion_relations():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	Q = [1.01, 1.20, 1.50, 2.00]
	colors =  ["navy","cornflowerblue","lightcoral","firebrick"]

	for q, c in zip(Q, colors):
		disc = WKB_Disc(0.366, epsilon = 0)
		disc.set_sigma_with_Q(q)
		r, k, CR = *disc.k_func_r(1, 0.99999999), disc.CR(1)
		r, k = sort_plotting_lst(r, k)

		k_star = disc.k_turning(1)
		fr = disc.forbidden_radius(1)


		plt.plot([k[i]/disc.k_crit(r[i]) for i in range(len(r))], [i/CR for i in r], color = c, label = f"{q:.2f}")
		plt.plot([0, k_star], [fr,fr] , color = c, linestyle = '--')
		plt.plot([k_star, k_star], [0, fr], linestyle = 'dotted', color = c)

	disc = WKB_Disc(0.366, epsilon = 0)
	plt.axhline(disc.ILR(1)/disc.CR(1), linestyle = '--', color = 'black')
	plt.axhline(1, linestyle = '--', color = 'black')
	plt.ylim([0.95*disc.ILR(1)/disc.CR(1), 1.05])
	plt.xlim([0,6])
	plt.legend(title = "Q", title_fontsize = 15, fontsize = 15)
	plt.ylabel(r"$R/R_{CR}$", fontsize = 15)
	plt.xlabel(r"$|k|/k_{crit}$", fontsize = 15)


	plt.show()

## Can I write a function that generaters 


plot_Dispersion_relations()