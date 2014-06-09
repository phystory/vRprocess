import os
import numpy as np
import struct
import binascii
import array
from matplotlib  import pyplot as plt
from matplotlib import animation as animation

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
xmin = 0
xmax = 200
ymin = 1E-20
ymax = 1

fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(ymin,ymax))
data, = ax.plot([], [], lw=2)
plt.yscale('log')

def init():
	data.set_data([], [])
	time_text.set_text('')
	return data,time_text

def animate(i):
	x = np.array(aa)
	y = np.array(xmf[i])
	data.set_xdata(x)
	data.set_ydata(y)
	time_text.set_text(time_template%(time[i]))
	return data, time_text

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

with open("tso1",'rb') as file_id:
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		desc1.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc2.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc3.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		data_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]
	record_length1=struct.unpack('i',file_id.read(4))[0]
	kstmx    =struct.unpack('i',file_id.read(4))[0]
	kitmx    =struct.unpack('i',file_id.read(4))[0]
	iweak    =struct.unpack('i',file_id.read(4))[0]
	iscrn    =struct.unpack('i',file_id.read(4))[0]
	iconvc   =struct.unpack('i',file_id.read(4))[0]
	changemx =struct.unpack('d',file_id.read(8))[0]
	tolm     =struct.unpack('d',file_id.read(8))[0]
	tolc     =struct.unpack('d',file_id.read(8))[0]
	yacc     =struct.unpack('d',file_id.read(8))[0]
	ymin     =struct.unpack('d',file_id.read(8))[0]
	tdel_mm  =struct.unpack('d',file_id.read(8))[0]
	record_length2=struct.unpack('i',file_id.read(4))[0]

	# Read Abundance Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		abund_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		abund_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))

	# Read Thermodynamic Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		thermo_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		thermo_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]

	# Read Nuclear Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	ny            =struct.unpack('i',file_id.read(4))[0]
	print "ny=",ny
	for ith in range (0,ny):
		zz.append(struct.unpack('d',file_id.read(8))[0])
	for ith in range (0,ny):
		aa.append(struct.unpack('d',file_id.read(8))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]
	#print zz,aa

	# Read Flux Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	nflx          =struct.unpack('i',file_id.read(4))[0]
	if nflx > 0:
		for ith in range(0,nflx):
			flx_end1=struct.unpack('i',file_id.read(4))[0]
			flx_end2=struct.unpack('i',file_id.read(4))[0]
			#print "flx_end1=",flx_end1,"flx_end2=",flx_end2
			#flx_end[ith].append(struct.unpack('i',file_id.read(4))[0])
			#flx_end[ith].append(struct.unpack('i',file_id.read(4))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]
	print record_length2

	# Read data from each timestep
	print "kstmx=",kstmx
	kstmx = 2000
	for k in range(1, kstmx):
		record_length1 =struct.unpack('i',file_id.read(4))[0]
		# If end of file, exit 
		if not record_length1:
			break
		# Otherwise read data
		kstep          =struct.unpack('i',file_id.read(4))[0]
		#time.append(struct.unpack('d',file_id.read(8))[0])
		#temperature.append(struct.unpack('d',file_id.read(8))[0])
		#density.append(struct.unpack('d',file_id.read(8))[0])
		#timestep.append(struct.unpack('d',file_id.read(8))[0])
		#edot.append(struct.unpack('d',file_id.read(8))[0])
		timel = struct.unpack('d',file_id.read(8))[0]
		temperaturel = struct.unpack('d',file_id.read(8))[0]
		densityl = struct.unpack('d',file_id.read(8))[0]
		timestepl = struct.unpack('d',file_id.read(8))[0]
		edotl = struct.unpack('d',file_id.read(8))[0]
		time.append(timel)
		xmf.append([])
		for ith in range(0,ny):
			xmfl = struct.unpack('d',file_id.read(8))[0]
			xmf[k-1].append(xmfl)
		if nflx > 0: 
			for ith in range(0,nflx):
				flxl = struct.unpack('d',file_id.read(8))[0]
				#print "flxl=",flxl
		record_length2 =struct.unpack('i',file_id.read(4))[0]
		#update_line(hl, aa, xmf)
		#print k,timel,xmf

	#print k,kstep,time,temperature,density
	# Convert abundance to mass fraction
	#for n in range(1, ny):
		#xmf[n] = aa[n]*xmf[n]

	# Calculate # of timesteps
	#nstep = len(time,2)
nn[:] = np.array(aa) - np.array(zz)
time_template = 'time = %.12e s'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
"""
#print nn,zz
for k in range(1, kstmx):
	if k == 1:
		#points, = ax.plot(nn, xmf[k-1], marker='s', ms=5, linestyle='None')
		points, = ax.plot(nn, xmf[k-1], linestyle='-')
	else:
		new_x = nn
		new_y = xmf[k-1]
		print new_y
		points.set_data(new_x, new_y)
	plt.pause(0.001)
#print aa,zz
"""
#ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=300,interval=100, blit=True)
#ani.save('rprocess1.mp4', fps=15, clear_temp=True)
plt.show()

