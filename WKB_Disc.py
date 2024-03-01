from math import *
import csv

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

class WKB_Disc():

	def __init__(self, sigmaR, densityFile = None, epsilon = 0.125, activeFraction = 0.5):

		self.sigmaR = sigmaR

		self.vc = 1
		self.activeFraction = activeFraction
		self.Sigma0 = self.activeFraction * 1/(2*pi)
		self.epsilon = epsilon # For the time being we will ignore softening 
		self.m = 2
		self.taperedDensity = False
		if densityFile != None:
			self.taperedDensity = True
			self.readin_density(densityFile)

	def __repr__(self):
		return f"WKB_Disc({self.sigmaR!r}, {self.taperedDensity}, {self.epsilon}, {self.activeFraction})"

	## General Functions ##
	## ----------------- ##

	def Omega(self, radius):
		return self.vc/radius

	def kappa(self, radius):
		return self.Omega(radius) * sqrt(2)

	def Sigma(self, radius):
		if self.taperedDensity == False:
			return self.Sigma0 * (1/radius)
		else:
			return ((11.6873/self.discMass)) *self.tapered_Sigma(radius) ## ((11.6873/self.discMass)) * is the mass ratio

	def k_crit(self, radius):
		return self.kappa(radius)**2 /(2*pi*self.Sigma(radius))

	def CR(self, omega, m = 2):
		return m / omega

	def ILR(self, omega, m = 2):
		return self.CR(omega, m) * (2-sqrt(2))*0.5

	def OLR(self, omega, m =2):
		return self.CR(omega, m) * (2+sqrt(2))*0.5

	def Q(self):
		return self.sigmaR / (self.activeFraction * (3.36/(sqrt(2)*2*pi)))

	def set_sigma_with_Q(self, Q):
		self.sigmaR = Q * self.activeFraction * (3.36/(sqrt(2)*2*pi))
	## Tapered Density ## 
	## --------------- ##

	def readin_density(self, filename):
		with open(filename) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ',')
			data = []
			for row in csv_reader:
				lst = [float(i) for i in row]
				data.append(lst)

		data = np.asarray(data, dtype=np.float64)
		self.spacing, self.upper_val, self.discMass = data[0,1], data[-1, 0], data[1,0]
		self.radii_values, self.v_tapered_density = data[2:, 0], data[2:, 1]

	def tapered_Sigma(self, radius):
		if (radius > self.upper_val):
			return self.Sigma0 * (1/radius)
		else:

			lowerIndex, trueIndex, upperIndex = (int(floor(radius/self.spacing))), radius/self.spacing, (int(ceil(radius/self.spacing)))
			
			if lowerIndex == upperIndex:
				return self.v_tapered_density[lowerIndex] * self.activeFraction
			return self.activeFraction * ((upperIndex-trueIndex) * self.v_tapered_density[lowerIndex] + (trueIndex-lowerIndex) * self.v_tapered_density[upperIndex])


	def check_interploation_density(self):
		plt.scatter(self.radii_values, ((11.6873/self.discMass)) *0.5*self.v_tapered_density, s = 0.7)
		plt.plot(np.linspace(0.5, 8), ((11.6873/self.discMass)) *(0.5/(2 * pi)) * 1/ np.linspace(0.5, 8))
		plt.plot(np.linspace(0, 8, 1000), list(map(self.tapered_Sigma, np.linspace(0, 8, 1000))))
		#plt.plot(disc.radii_values, list(map(disc.tapered_Sigma, disc.radii_values)))

		plt.xlabel("Radius")
		plt.ylabel("Density")
		plt.show()


	## Reduction Factor ##
	## ---------------- ##

	def reduction_factor(self, s, chi, stepSize = 0.1):
		if chi ==0:
			return 1
		if s==1:
			s=0.9999999 # The two lines tend to each other 
    
		tau = np.arange(0, pi, stepSize)
	    
		integrand = np.exp(-chi *(1 + np.cos(tau)))*np.sin(s * tau)*np.sin(tau)
		return stepSize*((1-s**2)/sin(pi * s)) * np.sum(integrand)

	def partialRFpartialS(self, s, chi, deltaS = 0.0001):
		return (1/(2*deltaS)) * (self.reduction_factor(s+deltaS, chi) - self.reduction_factor(s-deltaS, chi))

	def partialRFpartialChi(self, s, chi, deltaChi = 0.0001):
		return (1/(2*deltaChi)) * (self.reduction_factor(s, chi+deltaChi) - self.reduction_factor(s, chi-deltaChi))

	def chi(self, k, radius):
		return ((self.sigmaR*k)/self.kappa(radius))**2

	def s(self, omega0, radius):
		return (omega0 - self.m *self.Omega(radius))/self.kappa(radius)

	def omegaFromS(self, s, radius):
		return s * self.kappa(radius) + self.m * self.Omega(radius)

	## Dispersion relation ## 
	## ------------------- ##
		

	def LHS_disp(self, s, k, radius):
		chi = self.chi(k, radius)
		return s**2 + 2 * pi * self.Sigma(radius)/(self.kappa(radius)**2) * self.reduction_factor(s, chi)*abs(k) * exp(-self.epsilon * abs(k))-1

	def s_from_k(self, k, radius, nstep = 100): # returns s
	    upperVal, middleVal, lowerVal = .99999, 0.5, 0.000001
	    upper, lower = self.LHS_disp(upperVal, k, radius), self.LHS_disp(lowerVal, k, radius)
	    if k ==0: # by inspection this must be true
	        return 1
	    for n in range(nstep):
	        middle = self.LHS_disp(middleVal, k, radius)
	        if (upper*middle < 0): 
	            lowerVal = middleVal
	            middleVal = (lowerVal + upperVal)*0.5
	            lower = middle
	        
	        else:
	            upperVal = middleVal
	            middleVal = (lowerVal + upperVal)*0.5
	            upper = middle
	    
	    return middleVal

	def k_from_omega(self, omega, r, nstep = 300, upper = 20):
	    sVal, kc = self.s(omega, r), self.k_crit(r)
	    k_list = np.linspace(0,upper*kc, nstep)
	    
	    
	    RHS = lambda k : sVal**2 + 2 * pi * self.Sigma(r)/(self.kappa(r)**2) * self.reduction_factor(sVal, self.chi(k, r))*abs(k)*exp(-abs(k) *self.epsilon)-1
	    lst = list(map(RHS, k_list))
	    
	    k_index = [i for i in range(1, len(lst)) if lst[i]*lst[i-1] < 0]	    

	    if (len(k_index) ==1): # If we find one, we haven't looked at large enough k, therefore repeat this function with larger k range
	    	self.k_from_omega(omega, r, nstep, upper * 2) 
	    
	    interpolate = lambda i : k_list[i-1] + ((k_list[i] - k_list[i-1])/(1 + abs(lst[i]/lst[i-1])))
	    return [interpolate(i) for i in k_index]

	def forbidden_radius(self, omega0, nstep=100, rUpperScar = 10000, lowerLim = 1.01, upperLim = 0.99): # Returns the inner forbidden radius/CR 
	    cr = self.CR(omega0)
	    radii = np.linspace(lowerLim*self.ILR(omega0), upperLim*cr, nstep)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii))
	    r_plot = [radii[i]/cr for i, klst in enumerate(k) for void in klst]
	    
	    if rUpperScar/cr < max(r_plot):
	    	return rUpperScar/cr

	    if upperLim == 0.99: 	## THIS WILL RUN INTO PROBLEMS IF THERE IS AN UPPER SCAR (WILL GET TWO FACTORS OF 1/CR)
	    	r_plot = [self.forbidden_radius(omega0, nstep, rUpperScar, lowerLim =(cr* max(r_plot)) /self.ILR(omega0), upperLim = max(r_plot) * 1.05)] #]

	    return max(r_plot)

	def outer_forbidden_radius(self, omega0):
		return self.CR(omega0)*(2  - self.forbidden_radius(omega0))


	def k_func_r(self, omega0, forbidden_radius = 0.99): # Use this for plotting r vs k
	    cr = self.CR(omega0)
	    fr = self.forbidden_radius(omega0, 500)
	    radii = np.linspace(1.01*self.ILR(omega0), fr*cr, 100)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii, [1000 for e in radii]))
	    r_plot, k_plot = [radii[i] for i, klst in enumerate(k) for void in klst], [kvalue for klst in k for kvalue in klst]
	    
	    return r_plot, k_plot

	def k_turning(self, omega0):	 
		fr, CR = self.forbidden_radius(omega0, 500), self.CR(omega0)
		r, k = self.k_func_r(omega0, fr)

		index = r.index(max(r))
		return k[index]/self.k_crit(r[index])

	def plotting_k_vs_r(self, omega0 = 1, scars = []):
		fr, CR = self.forbidden_radius(omega0, 500), self.CR(omega0)
		r, k = self.k_func_r(omega0, fr)
		plt.scatter([i/CR for i in r], [k[i]/self.k_crit(r[i]) for i in range(len(r))], s =1)

		for scar in scars:
			plt.axvline(scar/self.CR(omega0), linestyle = '--', color = 'firebrick')
		plt.axvline(self.ILR(omega0)/CR, linestyle = '--', color = 'firebrick')
		plt.axhline(self.k_turning(omega0), linestyle = ':', color = 'k')
		

		plt.title(f"{omega0}")
		plt.xlabel(r"$R/R_{crit}$")
		plt.ylabel(r"$k/k_{crit}$")
		plt.show()

	## Group Velocity ##
	## -------------- ##

	def dRFds(self, s, chi):
		return (1/(2*0.001)) * (self.reduction_factor(s+0.001, chi) - self.reduction_factor(s-0.001, chi))

	def dRFdChi(self, s, chi):
		return (1/(2*0.001)) * (self.reduction_factor(s, chi+0.001) - self.reduction_factor(s, chi-0.001))


	def stellar_vg(self, k, radius, deltaX = 0.001, omega0 = 1):
	    kc = self.k_crit(radius)
	    x0, x1, x2 = abs(k/kc)- deltaX, abs(k/kc), abs(k/kc) + deltaX
	    s0, s1, s2 = self.s_from_k(x0*kc, radius), self.s_from_k(x1*kc, radius), self.s_from_k(x2*kc, radius)
	    grad = -((1/(x2-x0)) * (abs(s2)-abs(s0)))
	    
	    return (self.kappa(radius)/kc)* np.sign(k * s1) * grad

	

	## Groove Search for omega_0 ##
	## ------------------------- ##

	def integrate_k(self, omega0, r_inf, nstep = 500, rUpperScar = 0):
	    r_sup = self.forbidden_radius(omega0, rUpperScar  = rUpperScar)
	    radii = self.CR(omega0) * np.linspace(r_inf, r_sup, nstep)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii))
	    k_unpacked = [kvalue for klst in k for kvalue in klst]

	    return (sum(k_unpacked))*(radii[1] - radii[0]) # We need to integrate over both branches 
        

	def find_mode_omega0(self, inital_guess, r_inf, nstep = 10, rUpperScar = 10000): 
	    
	    omegaL, omegaM, omegaU = inital_guess[0], sum(inital_guess)/len(inital_guess), inital_guess[1]
	    integralL = self.integrate_k(omegaL, r_inf/self.CR(omegaL), rUpperScar = rUpperScar) - pi
	    integralU = self.integrate_k(omegaU, r_inf/self.CR(omegaU), rUpperScar = rUpperScar) - pi
	    
	    for n in range(nstep):
	        integralM = self.integrate_k(omegaM, r_inf/self.CR(omegaM), rUpperScar = rUpperScar) - pi
		# integralL = self.integrate_k(omegaL, 0.5*(2-sqrt(2)), rUpperScar = rUpperScar) - pi
		# integralU = self.integrate_k(omegaU, 0.5*(2-sqrt(2)), rUpperScar = rUpperScar) - pi
	    
		# for n in range(nstep):
		# 	integralM = self.integrate_k(omegaM, 0.5*(2-sqrt(2)), rUpperScar = rUpperScar) - pi
	        #print(omegaL, omegaU)
	        if (integralM*integralL < 0):
	            omegaU, integralU = omegaM, integralM
	            omegaM = 0.5 * (omegaL + omegaU)
	        else:
	            omegaL, integralL = omegaM, integralM
	            omegaM = 0.5 * (omegaL + omegaU)
	    
	    return omegaM;

	## Groove Search for eta ##
	## --------------------- ## 

	def find_mode_eta(self, omega0, r_inf, nstep = 500, rUpperScar = 10000): ## Can we make this better
		r_sup = self.forbidden_radius(omega0, rUpperScar = rUpperScar, nstep = nstep)
		radii = np.linspace(r_inf, self.CR(omega0) * r_sup, nstep)
		k = list(map(self.k_from_omega, [omega0 for e in radii], radii, [nstep for e in radii]))
		r, k_plot   = [radii[i] for i, klst in enumerate(k) for void in klst], [kvalue for klst in k for kvalue in klst]

		  
		vg_2_invert = list(map(self.stellar_vg, k_plot, r, [0.001 for e in r], [omega0 for e in r]))
		vg_inverted = [1/abs(v_g) for v_g in vg_2_invert]

		
		integral = sum(vg_inverted) - 0.5*(vg_inverted[0] + vg_inverted[1] + vg_inverted[-1] + vg_inverted[-2])

		return 0.25/(integral * (radii[1] - radii[0]))


	## Find Mode ##
	## ------------ ##

	def modeFinder(self, range_omega0 = [0.5, 0.6], rScar = 1.2, rUpperScar = 100000):
		omega = self.find_mode_omega0(range_omega0, rScar, rUpperScar = rUpperScar)
		eta = self.find_mode_eta(omega, rScar, rUpperScar = rUpperScar) ## This can sometimes produce an indexing error. 

		return omega + 1j* eta


    ## Wave Motion in WKB Disc ##
    ## ----------------------- ##

	def correct_k(self, omega0, radius, bool_lst): # bool_lst = [isLong, isTrailing]
		sgn = 1 if bool_lst[1] == True else -1
		k_lst = self.k_from_omega(omega0, radius, 300, 20)
		#print(radius/self.ILR(omega0), k_lst, bool_lst)

		return sgn * min(k_lst) if bool_lst[0] == True else sgn * max(k_lst)

	def update_region(self, omega0, radius, bool_lst, fr):
		if radius > 0.995 * fr *self.CR(omega0): 
			if ((bool_lst[0] == bool_lst[1])):
				bool_lst[0] = not (bool_lst[0])
			

		if ((bool_lst[1] == False) & (bool_lst[0] == True) & (radius < 1.05 * self.ILR(omega0))): 
			bool_lst[1] = not bool_lst[1]

		return bool_lst




	def motion_of_wave(self, omega0 = 0.4, rad = 8): 
		time, region, fr = np.linspace(0, 57, 800), [False, False], self.forbidden_radius(omega0)
		radius, k, vg, deltaT = np.zeros_like(time), np.zeros_like(time), np.zeros_like(time), time[1]
		radius[0] = 1.3*self.ILR(omega0) 
		k[0]  = self.correct_k(omega0, radius[0], region)
		vg[0] = self.stellar_vg(k[0], radius[0])
		
		cmap = ScalarMappable(cmap = 'viridis')
		
		
		for t in range(1, np.shape(time)[0]):
			radius[t] = radius[t-1] + vg[t-1] * deltaT 
			region = self.update_region(omega0, radius[t], region, fr)
			k[t]  = self.correct_k(omega0, radius[t], region)
			vg[t] = self.stellar_vg(k[t], radius[t])

		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')			
		fig, axs = plt.subplots()
		for i in range(np.shape(time)[0]-1):
			axs.plot([k[i]/self.k_crit(radius[i]) for i in [i,i+1]], [radius[i]/self.CR(omega0) for i in [i, i+1]], color = cmap.to_rgba(time * omega0/(2*pi))[i])
			axs.plot([k[i]/self.k_crit(radius[i]) for i in [i,i+1]], [2-(radius[i])/self.CR(omega0) for i in [i, i+1]], color = cmap.to_rgba(time * omega0/(2*pi))[i])

		axs.axhline(fr, linestyle = ':',  color = 'navy')
		axs.axhline(2-fr, linestyle = ':', color = 'navy')
		#axs.fill_between([k[i]/self.k_crit(radius[i]) for i in range(np.shape(time)[0])], [fr for _ in range(np.shape(time)[0])], [2-fr for _ in range(np.shape(time)[0])])

		axs.axhline(self.ILR(omega0)/self.CR(omega0), color = 'firebrick', linestyle = '--')
		axs.axhline(1, color = 'firebrick', linestyle = '--')
		axs.axhline(self.OLR(omega0)/self.CR(omega0), color = 'firebrick', linestyle = '--')

		axs.set_ylabel(r"$R/R_{CR}$", fontsize = 15)
		axs.set_xlabel(r"$k/k_{crit}$", fontsize = 15)

		axs.text(-3, 0.40, "A", fontsize = 12)
		axs.text(-1.2, 0.58, "B", fontsize = 12)
		axs.text(-0.15, 0.35, "C", fontsize = 12)
		axs.text(1.0, 0.58, "D", fontsize = 12)
		axs.text(3, 0.40, "E", fontsize = 12)

		
		k = disc.k_from_omega(0.4, rad)
		axs.scatter([-k[1]/self.k_crit(rad)],  [rad/self.CR(omega0)], color = 'firebrick')

		fig.colorbar(cmap, ax=axs, label = r"$t/(2\pi/\omega_{0})$").set_label(label = r"$t/(2\pi/\omega_{0})$", fontsize = 15)
		plt.show()


# disc = WKB_Disc(0.2835, epsilon = 0)
# disc.motion_of_wave(0.4)
# print(disc.Q(), 

# 	disc.outer_forbidden_radius(0.4),
# 	0.98*disc.OLR(0.4), 
# 	disc.k_from_omega(0.4, 1.5))




# disc.set_sigma_with_Q(1.2)
# print(disc.sigmaR, disc.CR(1))
# disc.plotting_k_vs_r(1)
# disc = WKB_Disc(0.366, epsilon = 0, activeFraction = 0.5)
# disc.set_sigma_with_Q(1.3)

# print(disc.k_turning(1))

# disc = WKB_Disc(sigmaR = 1/sqrt(12.4), activeFraction = 0.5, densityFile = "Disc_Density/Tapered_R_10_W_25_D_-0_G.csv")
# print(disc.modeFinder())
