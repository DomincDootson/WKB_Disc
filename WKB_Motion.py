from WKB_Disc import *

from matplotlib.cm import *
#import matplotlib as mpl

def waveMotionPlot():
	disc = WKB_Disc(0.40, epsilon =  0, activeFraction = 1)
	print(disc.forbidden_radius(0.6))
	disc.motion_of_wave(0.6)


waveMotionPlot()

