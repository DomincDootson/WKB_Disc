from WKB_Disc import *

import matplotlib.animation as animation

class WKB_Wave():
	def __init__(self, omega0):

		self.disc = WKB_Disc(0.40, epsilon =  0, activeFraction = 1)
		self.omega0 = omega0

		self.n_points =100
		self.radii_limits = [1.001, 0.999]
		self.FR = self.disc.forbidden_radius(self.omega0)
		#self.vg_long, self.vg_short = self.vg()

		self.deltaT = 0.05
		self.end_time = 50
		self.n_step = 800

		self.grid_radii = self.grid_radii()
		self.grid_step = self.grid_radii[1] - self.grid_radii[0]
		self.grid_envelope = np.zeros((self.n_step, np.shape(self.grid_radii)[0]))

		self.grid_k = self.grid_k()
		self.grid_vg = self.grid_vg()


	def __repr__(self):
		return f"WKB_Wave({self.disc})"


	
	def radii(self):
		return np.linspace(self.radii_limits[0]*self.disc.ILR(self.omega0), 
							self.radii_limits[-1]*self.disc.CR(self.omega0)*self.disc.forbidden_radius(self.omega0),
							self.n_points)

	def vg(self):
		radii = self.radii()

		k_long   = list(map(self.disc.correct_k, [self.omega0 for _ in radii], radii, [[True,  True] for _ in radii])) 
		vg_long  = list(map(self.disc.stellar_vg, k_long, radii))
		
		k_short  = list(map(self.disc.correct_k, [self.omega0 for _ in radii], radii, [[False, True] for _ in radii])) 
		vg_short = list(map(self.disc.stellar_vg, k_short, radii))

		return vg_long, vg_short

	def k_func_r(self, isLong = True, isTrailing = True):
		return list(map(self.disc.correct_k, [self.omega0 for _ in self.radii()], self.radii(), [[isLong, isTrailing] for _ in self.radii()]))


	def initial_envelope(self):
		rad = self.radii() 
		envelope = rad * 0
		
		amp, mu, sigma = 1, 1.25*self.disc.ILR(self.omega0), 0.07 * 1.1*self.disc.ILR(self.omega0)

		envelope_func = lambda r : amp * exp(-0.5 * ((r-mu)/sigma)**2) if (abs((r-mu)/sigma < 4)) else 0 # Use Guassian profile

		return rad, list(map(envelope_func, rad))
		

	def initial_wave(self):
		rad, envelope = self.initial_envelope()
		k_list = self.k_func_r(isLong = False, isTrailing = False)

		return [env*cos(k *r) for env, k, r in zip(envelope, k_list, rad)]



	## Evolution Functions ##
	## ------------------- ##

	def grid_radii(self):
		radii = list(self.radii())
		return np.asarray(radii + radii[::-1] + radii + radii[::-1])

	def grid_k(self):
		k_short, k_long = self.k_func_r(isLong = False, isTrailing = True), self.k_func_r(isLong = True, isTrailing = True) # Want all velocities > 0
		return 	[-k for k in k_short] + [-k for k in k_long[::-1]] + k_long + k_short[::-1]

	def grid_vg(self):
		return [abs(vg) for vg in map(self.disc.stellar_vg, self.grid_k, self.grid_radii)]


	def envelope_step(self, time):
		new_envelope = np.zeros_like(self.grid_envelope[time, :]) 

		for index in range(np.shape(new_envelope)[0]):
			mass2move = self.grid_radii[index] * self.grid_envelope[time-1, index] # Conservation of mass, not density 
			newIndex = index + (self.grid_vg[index]*self.deltaT) * (1/self.grid_step)
			
			if (ceil(newIndex) >= np.shape(new_envelope)[0]): # Deal with edge case
				new_envelope[-1] += mass2move/(self.grid_radii[-1])
			
			else:
				new_envelope[int(newIndex//1)]     += (ceil(newIndex) - newIndex)  * mass2move * (1/self.grid_radii[int(newIndex//1)])
				new_envelope[int(newIndex//1) + 1] += (newIndex - floor(newIndex)) * mass2move * (1/self.grid_radii[int(newIndex//1) + 1])
				
		return new_envelope

	def envelope_evolution(self):
		self.grid_envelope[0, 0:self.n_points] = self.initial_envelope()[1]
		
		for time in range(1, self.n_step):
			self.grid_envelope[time, :] = self.envelope_step(time)
			#print(sum(self.grid_radii * self.grid_envelope[time, :])) # Use to check that mass is conserved 

	def split_envelope(self, time):
		return self.grid_envelope[time,:self.n_points], self.grid_envelope[time,self.n_points:(2*self.n_points)], self.grid_envelope[time,(2*self.n_points):(3*self.n_points)], self.grid_envelope[time,(3*self.n_points):]

	def split_radii(self):
		return self.grid_radii[:self.n_points], self.grid_radii[self.n_points:(2*self.n_points)], self.grid_radii[(2*self.n_points):(3*self.n_points)], self.grid_radii[(3*self.n_points):]


	## Wave Evolution ## 
	## -------------- ## 


	def grid_integrate_k(self):
		rad_short, rad_long, *_ = self.split_radii()
		rad_long = rad_long[::-1]

		k_int_long = [self.grid_step*rad_long[0]*self.grid_k[2*self.n_points]]#*r* 
		for i in range(1,self.n_points):
			k_int_long.append(k_int_long[-1] + self.grid_step*rad_long[i] * self.grid_k[2*self.n_points + i])

		k_int_short = [k_int_long[-1]]
		for i in range(1,self.n_points):
			k_int_short.append(k_int_short[-1] + self.grid_step*rad_long[i] * self.grid_k[3*self.n_points+i])

		self.k_grid_integrated =  np.asarray([-k for k in k_int_short][::-1] + [-k for k in k_int_long][::-1] +k_int_long + k_int_short)


	def density_at_radius(self, time):
		density_radius_unwrapped = self.grid_envelope[time, :] * np.exp(1j * self.k_grid_integrated)
		return density_radius_unwrapped[:self.n_points] + density_radius_unwrapped[2*self.n_points-1:self.n_points-1:-1] + density_radius_unwrapped[2*self.n_points:3*self.n_points] + density_radius_unwrapped[2*self.n_points-1:1*self.n_points-1:-1]
				
	def density_at_point(self, den_at_rad, R, phi):
		if (R < self.radii_limits[0]*self.disc.ILR(self.omega0)) or (R > self.FR*self.radii_limits[-1]*self.disc.CR(self.omega0)): 
			return 0

		# First interpolate to find correct index
		index = (R - self.radii_limits[0]*self.disc.ILR(self.omega0))/self.grid_step
		return np.exp(1j * 2 * phi)*complex((index - floor(index)) * den_at_rad[int(index//1 +1)] + (ceil(index) - index) * den_at_rad[int(index)])
		
		


	def density_grid(self, time):
		den_at_rad = self.density_at_radius(time)

		n_grid =int((2*self.radii_limits[-1]*self.disc.CR(self.omega0))/self.grid_step)
		centre = (self.radii_limits[-1]*self.disc.CR(self.omega0))/self.grid_step
		grid = np.zeros((n_grid, n_grid), dtype = np.csingle)


		for i in range(n_grid):
			for j in range(n_grid):
				x, y = (j - centre) * self.grid_step, (i - centre) * self.grid_step
				R, phi = sqrt(x**2+y**2), np.arctan2(y, x)

				grid[i,j] = self.density_at_point(den_at_rad, R, phi)

		return 2*np.real(grid)



	## Some General Plotting Functions ##
	## ------------------------------- ## 

	def group_velocity_plot(self):
		radii = self.radii()

		plt.plot(radii, self.vg_long)
		plt.plot(radii, self.vg_short)
		plt.show()

	def envelope_evolution_animation(self):
		fig, axs = plt.subplots()

		self.envelope_evolution()
		radii = self.split_radii()

		ims = []
		for time in range(self.n_step):
			line0, = axs.plot(radii[0], self.split_envelope(time)[0], animated = True, color = 'firebrick')
			line1, = axs.plot(radii[1], self.split_envelope(time)[1], animated = True, color = 'royalblue')
			line2, = axs.plot(radii[2], self.split_envelope(time)[2], animated = True, color = 'royalblue')
			line3, = axs.plot(radii[3], self.split_envelope(time)[3], animated = True, color = 'firebrick')
		
			ims.append([line0,line1, line2, line3])

		ani = animation.ArtistAnimation(fig, ims, interval=30)
		plt.show()


	def density_evolution_animation(self, skip = 5):
		fig, axs = plt.subplots()

		self.envelope_evolution()
		self.grid_integrate_k()

		ims = []

		for time in range(0, self.n_step, skip):
			print(time)
			plot = axs.imshow(self.density_grid(time), animated = True)
			ims.append([plot])

		ani = animation.ArtistAnimation(fig, ims, interval=30)
		plt.show()


wave = WKB_Wave(0.6)
#plt.plot(*wave.initial_envelope())
#plt.plot(wave.radii(), wave.initial_wave())
#print(len(wave.grid_radii))
wave.envelope_evolution_animation()

plt.show()
