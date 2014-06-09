#A python script for plotting abundances on the nuclear chart as a function of time
#KJC 2/6/2014


import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter 
import matplotlib.animation as animation 
#from itertools import islice 

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

#print Nstable 
#Import 150 isotopes from netwinv 

#with open('netwinv') as fin, open('output','w') as fout: 
#    fout.writelines(islice(fin,None,None,4))

#print fout

network = np.genfromtxt('150network.dat',delimiter=" ", dtype = float) 
Znetwork = [row[2] for row in network]
Nnetwork = [row[3] for row in network]

#print Znetwork
#print Nnetwork

#Import the data from XNet ev1 outputs

data = np.genfromtxt('ev150np1.dat', dtype = float) 

abundance=[]
for i in range (6,156):
    a = [row[i] for row in data]
    abundance.append(a) 

time = []
time = [row[1] for row in data]

print type(abundance)

color_data = np.asarray(abundance) 

#transpose the graph

#color_data=zip(*abundance) 
color_data=np.transpose(color_data) 
#print type(color_data)
color_data = np.asarray(color_data) 

#Make it log
color_data=color_data+0.000000000000000000000000000000001
color_data = np.log10(color_data)


#---------Plotting begins here------

#Plot a white set of boxes
fig = plt.figure()
ax = fig.add_subplot(111) 
big=ax.scatter(Nbig,Zbig,c='w',marker='s', s=40)
#big1=ax.scatter(Nbig,Nbig)
#initialise animation 
scatt = ax.scatter(Nnetwork,Znetwork,c='w',marker ='s', s=40)
#ax.text(0,25,'time')
ax.set_xlabel("Number of Neutrons",fontsize=12)
ax.set_ylabel("Number of Protons",fontsize=12)
ax.set_xlim(-2,40)
ax.set_ylim(-2,31)



def update_plot(i, data, scatt):
    scatt.set_array(data[i])
    return scatt,

ani=animation.FuncAnimation(fig,update_plot,frames=len(color_data), interval=1, fargs=(color_data,scatt))


cbar = plt.colorbar(scatt)

stablep=ax.scatter(Nstable,Zstable,c='k',marker='x', s=20)
cbar.set_label('Log10(Abundance)', rotation=270)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('np150chart.mp4', fps=70, extra_args=['-vcodec', 'libx264'])
ani.save('np150chart.mp4', fps=70)


plt.show()


