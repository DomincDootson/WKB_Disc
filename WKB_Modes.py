from WKB_Disc import *

def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []
		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data, dtype=np.float64)

def deRijkeMode():
	disc = WKB_Disc(1/sqrt(12.4))

	print(disc.modeFinder([0.4, 0.9], 1)) 




def varyingSoftening(innerPosition, filename, softeningValues):
	omega0, rInner = [], []
	f = open(filename, 'w')
	for rIn in innerPosition:
		print(f"Finding mode for {rIn:.2f}")
		densityFile = f"Disc_Density/Tapered_R_{str(int(rIn*10))}_W_25_D_-95_G.csv"
		print("Density file: " + densityFile)
		disc = WKB_Disc(1/sqrt(12.4), epsilon = softeningValues, densityFile = "Tapered_Density.csv")
		omega0.append(disc.modeFinder([0.2, 0.9], rScar = rIn))
		rInner.append(disc.forbidden_radius(omega0[-1].real) * disc.CR(omega0[-1].real))
		f.write(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")
		print(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")

	f.close()


	
	
def varyingInnerPositon(innerPosition, filename = "Cavity_Modes/VaryingInnerPosition.csv", softeningValues = 0):
	omega0, rInner = [], []
	f = open(filename, 'w')
	for rIn in innerPosition:
		print(f"Finding mode for {rIn:.2f}")
		
		densityFile = f"Disc_Density/Tapered_R_{str(int(rIn*10))}_W_25_D_-95_G.csv"
		print("Density file: " + densityFile)
		

		disc = WKB_Disc(1/sqrt(12.4), epsilon = softeningValues, densityFile = densityFile)
		omega0.append(disc.modeFinder([0.4, 0.9], rScar = rIn))
		rInner.append(disc.forbidden_radius(omega0[-1].real) * disc.CR(omega0[-1].real))
		
		f.write(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")
		print(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")

	f.close()

	
def plottingInnerPosition(files = ["Cavity_Modes/VaryingInnerPosition.csv", "Cavity_Modes/VaryingInnerPosition_Tapered.csv", "Softening_Data/Tapered_Unscared_Soft_1.csv", "Cavity_Modes/VaryingInnerPosition_Tapered_Scarred.csv"], scatterFiles = None):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	fig, axs = plt.subplots(ncols = 2)
	labels, colors, ls = ["Untapered","Tapered", "Tapered Softened"], ["royalblue", "firebrick", "firebrick"], ['-', '-', '--']
	for file, l,c,s in zip(files, labels, colors, ls):
		data = readingInRealCSV(file)
		axs[0].plot(data[:,0], data[:,1], label = l, color = c, linestyle = s) 
		axs[1].plot(data[:,0], data[:,2]*(2/0.693), color = c, linestyle = s) # 2/ln(2) to account for the wrong factor we were using perviously 

	if scatterFiles is not None:
		for file in scatterFiles:
			data = readingInRealCSV(file)
			axs[0].scatter(data[:,0], data[:, 1], color = 'firebrick')
			axs[1].scatter(data[:,0], (data[:, 2]), color = 'firebrick', label = "Volterra Method") 		

	axs[0].legend(fontsize = 12, title = "Density Model For WKB", title_fontsize=12)
	axs[0].set_xlabel(r"Inner Scar Radius, $R_{in}$", fontsize = 12)
	axs[0].set_ylabel(r"$\omega_{0}$", fontsize = 12)
	axs[0].set_title(r"Spatial Frequency", fontsize = 15)

	
	axs[1].legend(fontsize = 12)
	axs[1].set_xlabel(r"Inner Scar Radius, $R_{in}$", fontsize = 12)
	axs[1].set_ylabel(r"$\eta$")
	axs[1].set_title(r"Growth Rate", fontsize = 15)

	plt.show() 

def varyingOuterPositon(outerPosition, innerPosition = 1.2):
	omega0, rOuter = [], []

	for rO in outerPosition:
		print(f"Finding mode for {rO:.2f}")
		disc = WKB_Disc(1/sqrt(12.4), epsilon = 0)
		omega0.append(disc.modeFinder([0.2, 2], rScar = innerPosition, rUpperScar = rO))
		rOuter.append(disc.forbidden_radius(omega0[-1].real) * disc.CR(omega0[-1].real))
		print(omega0[-1])

	fig, axs = plt.subplots(ncols = 3)

	axs[0].plot(outerPosition, [each.real for each in omega0])
	axs[0].set_xlabel(r"Outer Scar Radius")
	axs[0].set_ylabel(r"$\omega_{0}$")

	axs[1].plot(outerPosition, [each.imag for each in omega0]) 
	axs[1].set_xlabel(r"Outer Scar Radius")
	axs[1].set_ylabel(r"$\eta$")

	axs[2].plot(outerPosition, rOuter)
	axs[2].plot(outerPosition, [rOuter[i] / disc.CR(omega0[i].real) for i in range(len(omega0))]) 
	axs[2].set_xlabel(r"Outer Scar Radius")
	axs[2].set_ylabel(r"Outer Position")

	plt.show()

## Softened Data ##
## ------------- ##

def softenedComparison():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	files = ["Softening_Data/Tapered_Scared_Hard.csv", "Softening_Data/Tapered_Scared_Soft.csv", "Softening_Data/Tapered_Uncared_Hard.csv", "Softening_Data/Tapered_Unscared_Soft.csv","Softening_Data/Untapered_Unscared_Hard.csv","Softening_Data/Untapered_Unscared_Soft.csv"]
	data = [readingInRealCSV(file) for file in files]
	fig, axs = plt.subplots(ncols = 2, sharey = True)

	# Hard #

	axs[0].set_title(r"$\epsilon = 0$", fontsize = 12)
	axs[0].plot(data[4][:,0]+0.125, data[4][:,1], color = 'royalblue', label = "Untapered")
	axs[0].plot(data[2][:,0]+0.125, data[2][:,1], color = 'firebrick', linestyle = '--', label = "Tapered")
	axs[0].plot(data[0][:,0]+0.125, data[0][:,1], color = 'firebrick', label = "Tapered Scarred")
	axs[0].set_xlabel(r"$R_{s}$", fontsize = 12)
	axs[0].set_ylabel(r"$\omega_{0}$", fontsize = 12)


	# Soft #
	axs[1].set_title(r"$\epsilon = 1/8$", fontsize = 12)
	axs[1].plot(data[5][:,0]+0.125, data[5][:,1], color = 'royalblue', label = "Untapered")
	axs[1].plot(data[3][:,0]+0.125, data[3][:,1], color = 'firebrick', linestyle = '--', label = "Tapered")
	axs[1].plot(data[1][:,0]+0.125, data[1][:,1], color = 'firebrick', label = "Tapered Scarred")
	axs[1].set_xlabel(r"$R_{s}$", fontsize = 12)
	
	axs[1].legend()

	plt.show()
	




def WKBapproximation():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots()
	disc = WKB_Disc(1/sqrt(12.4), densityFile = "Disc_Density/Tapered_R_20_W_25_D_-95_G.csv", epsilon = 0)
	innerScar = disc.k_func_r(0.62)
	outerScar = disc.k_func_r(0.39)

	innerR = [innerScar[0][i] for i in range(0,len(innerScar[0]), 2)] + [innerScar[0][i] for i in range(len(innerScar[0])-1, 0, -2)]
	innerk = [innerScar[1][i] for i in range(0,len(innerScar[1]), 2)] + [innerScar[1][i] for i in range(len(innerScar[1])-1, 0, -2)]

	axs.plot([i -1 for i in innerR], [k*r for k, r in zip(innerR, innerk)], color = "royalblue", label = r"$1.0$")
	

	outerR = [outerScar[0][i] for i in range(0,len(outerScar[0]), 2)] + [outerScar[0][i] for i in range(len(outerScar[0])-1, 0, -2)]
	outerk = [outerScar[1][i] for i in range(0,len(outerScar[1]), 2)] + [outerScar[1][i] for i in range(len(outerScar[1])-1, 0, -2)]

	axs.plot([i-2 for i in outerR], [k*r for k, r in zip(outerR, outerk)], color ="firebrick", label = r"$2.0$")
	

	axs.legend(title = r"Scar Radius, $R_{i}$", fontsize = 12)

	axs.set_xlabel(r"$R-R_{i}$", fontsize = 12)
	axs.set_ylabel(r"$kR$", fontsize = 12)
	axs.set_ylim([0, 17.5])
	axs.axvline(0, linestyle ='--', color ="black")
	
	plt.show()

def forbiddenRadii(rScar, omega0 = []):
	
	for o, rIn in zip(omega0, rScar):
		densityFile = f"Disc_Density/Tapered_R_{str(int(rIn*10))}_W_25_D_-95_G.csv"
		disc =  WKB_Disc(1/sqrt(12.4), epsilon = 0, densityFile = densityFile)
		print(disc.CR(o) * disc.forbidden_radius(o))



## Varying Depth of Scar ##
## --------------------- ##

def varyingScarDepth(file = "Cavity_Modes/VaryingScarDepth.csv", r = "15"):
	depths = np.arange(0, 1+0.1, 0.1)
	densityFile = lambda depth : f"Disc_Density/Tapered_R_{r}_W_25_D_-{depth*100:.0f}_G.csv"
	discs = [WKB_Disc(1/sqrt(12.4), densityFile = densityFile(d), epsilon = 0) for d in depths]
	
	
	with open(file, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		for d, disc in zip(depths, discs):	
			print(f'Depth of Scar: {d}')
			writer.writerow([d, disc.modeFinder(range_omega0 = [0.7, 0.3], rScar = float(r)/10).real]) # Please be careful when setting the initial guess
	
def plottingVaryingScarDepth():


	fig, axs = plt.subplots()

	data = readingInRealCSV("Cavity_Modes/VaryingScarDepth_10.csv")
	axs.plot(data[:,0], data[:,1])

	data = readingInRealCSV("Cavity_Modes/VaryingScarDepth_15.csv")
	axs.plot(data[:,0], data[:,1])

	data = readingInRealCSV("Cavity_Modes/VaryingScarDepth_20.csv")
	axs.plot(data[:,0], data[:,1])
	axs.set_ylabel(r"$\omega_{0}$", fontsize = 15)
	axs.set_xlabel(r"$A_{s}$", fontsize = 15)

	plt.show()



#varyingSoftening(np.linspace(0.08, 0.10, 3))
#varyingInnerPositon(np.linspace(1, 2, 11), "Cavity_Modes/VaryingInnerPosition_Tapered_Scarred.csv")
#plottingInnerPosition()
#varyingOuterPositon(np.linspace(1.5, 2, 1))
# disc = WKB_Disc(1/sqrt(12.4))
# print(disc.forbidden_radius(0.371186)*disc.CR(0.371186))
#varyingSoftening(np.linspace(1, 2, 11), "Softening_Data/Tapered_Unscared_Soft_1.csv", 0.125)
plottingInnerPosition(scatterFiles = ["Cavity_Modes/VaryingInnerPositionResponse_25_-95.csv"])
#forbiddenRadii(rScar = [1.2, 1.4, 1.6, 1.8, 2.0], omega0 = [0.564407, 0.510169, 0.469492, 0.428814, 0.394915])
#deRijkeMode()

'''disc = WKB_Disc(1/sqrt(12.4), epsilon = 0)
disc.plotting_k_vs_r(0.40453948974609377, scars = [2])'''


# disc = WKB_Disc(1/sqrt(12.4), densityFile = "Disc_Density/Tapered_R_20_W_25_D_-95_G.csv", epsilon = 0)
# print(disc.forbidden_radius(0.9288)*disc.CR(0.9288), disc.ILR(0.9288), disc.CR(0.9288), disc.OLR(0.9288))


# disc.modeFinder([0.35, 0.55], rScar = 2)


#varyingSoftening(np.linspace(1, 2, 11), "Softening_Data/Tapered_Uncared_Hard.csv", 0)
#varyingSoftening(np.linspace(1, 2, 11), "Softening_Data/Tapered_Unscared_Soft_1.csv", 0.125)




# disc = WKB_Disc(0.377, activeFraction = 0.5, densityFile = "Tapered_Density.csv")

# r = np.linspace(0.5, 5)
# rho = [disc.Sigma(r_e) for r_e in r]
# plt.plot(r, rho)

# disc = WKB_Disc(0.377, activeFraction = 0.5, densityFile = "Disc_Density/Tapered_R_20_W_25_D_-95_G.csv")

# r = np.linspace(0.5, 5)
# rho = [disc.Sigma(r_e) for r_e in r]
# plt.plot(r, rho)
# plt.show()


#varyingScarDepth()
# varyingScarDepth(file = "Cavity_Modes/VaryingScarDepth_20.csv", r = "20")
# varyingScarDepth(file = "Cavity_Modes/VaryingScarDepth_10.csv", r = "10")
# plottingVaryingScarDepth()