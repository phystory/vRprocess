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

gstep = 1
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

cm = mpl.cm.get_cmap('Reds')
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

Nnetwork = np.asarray(nn)
Znetwork = np.asarray(zz)
color_data = np.asarray(xmf) 

#transpose the graph

#color_data=zip(*abundance) 
#color_data=np.transpose(color_data) 
#print type(color_data)
#color_data = np.asarray(color_data) 

#Make it log
color_data=color_data+0.000000000000000000000000000000001
color_data = np.log10(color_data)

#---------Plotting begins here------

#Plot a white set of boxes
fig = plt.figure()
ax = fig.add_subplot(111) 
big=ax.scatter(Nbig,Zbig,c='w',marker='s', lw=0.1, s=40)
#initialize animation 
scatt = ax.scatter(Nnetwork,Znetwork,c='w',marker ='s', lw=0.1, s=40, cmap=cm, norm=LogNorm(vmin=1.0e-25, vmax=1.0e+0))
plt.title("Nuclear Chart for Abundance change in the r-process over time")
ax.set_xlabel("Number of Neutrons",fontsize=12)
ax.set_ylabel("Number of Protons",fontsize=12)
ax.set_xlim(-2,250)
ax.set_ylim(-2,150)
time_template = 'time = %.12e s'
time_text = ax.text(0.1, 0.95, '', transform=ax.transAxes)
time_text.set_text('')

def update_plot(k, data, scatt):
    global time_text, gstep
    gith = k*gstep
    if gith > (kstep-2) :
        gith = kstep-2
    scatt.set_array(data[gith])
    scatt.set_cmap(cm)
    #print gith, data[gith]
    time_text.set_text(time_template%(time[gith]))
    return scatt

norstep = int(round(float(kstep)/float(gstep) + 1.49))
ani=animation.FuncAnimation(fig,update_plot, frames=norstep, interval=1, fargs=(color_data,scatt))

stablep=ax.scatter(Nstable,Zstable,c='k',marker='x', s=20)

cbar = plt.colorbar(scatt)
cbar.set_label('Abundance', rotation=270)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('np150chart.mp4', fps=70, extra_args=['-vcodec', 'libx264'])
#ani.save('np150chart.mp4', fps=70)


plt.show()
