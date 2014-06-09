import sys
import os
import subprocess as sp
import numpy as np
import struct
from array import array
import matplotlib as mpl
from matplotlib  import pyplot as plt
from matplotlib import animation as animation
from matplotlib.colors import LogNorm

"""
Python code to read tso file based on matlab code
function [zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename )
%--------------------------------------------------------------------------
%[zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename ) 
% Reads XNet binary output file.
% Inputs>  ts_filename: name of binary file 
% Outputs< zz: proton number for each isotope in the network.
%          aa: total nucleon number for each isotope in the network
%          xmf: time series of mass fractions 
%          time: discrete temporal evolution
%          temperature:  time series of temperature
%          density: time series of density
%          timestep:  time series of discrete timesteps 
%          edot: time series of energy generation 
%          flux_end, the starting and ending points for each reaction flux.
%          flux: time seris of the reaction fluxes in the network.
%--------------------------------------------------------------------------
"""

gith = 0
desc1 = []
desc2 = []
desc3 = []
data_desc = []
abund_file = []
abund_desc = []
thermo_file = []
thermo_desc = []
zz = []
aa = []
nn = []
xmf = []
flx = [[3]]
flx_end = []
time = []
temperature = []
density = []
timestep = []
edot = []
ims = []

with open("tso25",'rb') as file_id:
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	for ith in range (0,80):
		desc1.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc2.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc3.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		data_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('<i',file_id.read(4))[0]
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	kstmx    =struct.unpack('<i',file_id.read(4))[0]
	kitmx    =struct.unpack('<i',file_id.read(4))[0]
	iweak    =struct.unpack('<i',file_id.read(4))[0]
	iscrn    =struct.unpack('<i',file_id.read(4))[0]
	iconvc   =struct.unpack('<i',file_id.read(4))[0]
	changemx =struct.unpack('<d',file_id.read(8))[0]
	tolm     =struct.unpack('<d',file_id.read(8))[0]
	tolc     =struct.unpack('<d',file_id.read(8))[0]
	yacc     =struct.unpack('<d',file_id.read(8))[0]
	ymin     =struct.unpack('<d',file_id.read(8))[0]
	tdel_mm  =struct.unpack('<d',file_id.read(8))[0]
	record_length2=struct.unpack('<i',file_id.read(4))[0]

	# Read Abundance Info
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	for ith in range (0,80):
		abund_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		abund_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('<i',file_id.read(4))

	# Read Thermodynamic Info
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	for ith in range (0,80):
		thermo_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		thermo_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('<i',file_id.read(4))[0]

	# Read Nuclear Info
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	ny            =struct.unpack('<i',file_id.read(4))[0]
	print "ny=",ny
	for ith in range (0,ny):
		zz.append(struct.unpack('<d',file_id.read(8))[0])
	for ith in range (0,ny):
		aa.append(struct.unpack('<d',file_id.read(8))[0])
	record_length2=struct.unpack('<i',file_id.read(4))[0]
	#print zz,aa

	# Read Flux Info
	record_length1=struct.unpack('<i',file_id.read(4))[0]
	nflx          =struct.unpack('<i',file_id.read(4))[0]
	if nflx > 0:
		for ith in range(0,nflx):
			flx_end1=struct.unpack('<i',file_id.read(4))[0]
			flx_end2=struct.unpack('<i',file_id.read(4))[0]
			#flx_end[ith].append(struct.unpack('<i',file_id.read(4))[0])
			#flx_end[ith].append(struct.unpack('<i',file_id.read(4))[0])
	record_length2=struct.unpack('<i',file_id.read(4))[0]
	print record_length2

	# Read data from each timestep
	print "kstmx=",kstmx
	#kstmx=2000
	for k in range(1, kstmx):
		try:
			record_length1 =struct.unpack('<i',file_id.read(4))[0]
		except:
			break
		# If end of file, exit 
		if not record_length1:
			break
		# Otherwise read data
		kstep          =struct.unpack('<i',file_id.read(4))[0]
		#time.append(struct.unpack('<d',file_id.read(8))[0])
		#temperature.append(struct.unpack('<d',file_id.read(8))[0])
		#density.append(struct.unpack('<d',file_id.read(8))[0])
		#timestep.append(struct.unpack('<d',file_id.read(8))[0])
		#edot.append(struct.unpack('<d',file_id.read(8))[0])
		timel = struct.unpack('<d',file_id.read(8))[0]
		temperaturel = struct.unpack('<d',file_id.read(8))[0]
		densityl = struct.unpack('<d',file_id.read(8))[0]
		timestepl = struct.unpack('<d',file_id.read(8))[0]
		edotl = struct.unpack('<d',file_id.read(8))[0]
		time.append(timel)
		xmf.append([])
		for ith in range(0,ny):
			xmfl = struct.unpack('<d',file_id.read(8))[0]
			xmf[k-1].append(aa[ith]*xmfl)
		if nflx > 0: 
			for ith in range(0,nflx):
				flxl = struct.unpack('<d',file_id.read(8))[0]
				#print "flxl=",flxl
		record_length2 =struct.unpack('<i',file_id.read(4))[0]

nn[:] = np.array(aa) - np.array(zz)

#Set up the chart. Import the A,Z,N of known nuclei

chart = np.genfromtxt('zna.dat',delimiter=" ", dtype = float) 
stable = np.genfromtxt('Stable_Nuclides.txt',delimiter=" ",dtype=float)
Abig = [row[0] for row in chart]
Zbig = [row[1] for row in chart]
Nbig = [row[2] for row in chart]

Abig = np.asarray(Abig)
Zbig = np.asarray(Zbig)
Nbig = np.asarray(Nbig) 

Astable = np.asarray([row[0] for row in stable]) 
Zstable =  np.asarray([row[1] for row in stable]) 
Nstable = Astable - Zstable 

color_data = np.asarray(abundance) 

#transpose the graph

#color_data=zip(*abundance) 
color_data=np.transpose(color_data) 
#print type(color_data)
color_data = np.asarray(color_data) 

#Make it log
color_data=color_data+0.000000000000000000000000000000001
color_data = np.log10(color_data)


fig=plt.figure()
ax = fig.add_subplot(111)
plot, = ax.plot([], [])
print "start"
def init():
	plot=fig.clf()
	return plot
def animate(k):
	global gith
	#while time[i]<k:
		#i=i+1
	print "animate", gith, time[gith]
	plt.clf()
	plt.ylim([0,150])
	plt.xlim([0,250])
	map=np.zeros((150,250))
	
	for l in xrange(0,ny):
		map[int(zz[l])][int(nn[l])]=xmf[gith][l]
	plot=plt.pcolor(map,norm=LogNorm(vmin=1.0e-25, vmax=1.0e+0))
	plt.title(sys.argv[1]+" t="+'%12e s'%time[gith])
	plt.colorbar()
	gith = gith+100
	return plot
#writer = animation.MovieWriter()
#ani = animation.FuncAnimation(fig,animate,np.arange(0, int(time[-1]),100000),blit=False,init_func=init)
ani = animation.FuncAnimation(fig,animate,blit=False,init_func=init)
#ani.save("ab_"+sys.argv[1]+".mp4", writer=writer,fps=15)
plt.show()
