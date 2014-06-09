import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from array import array
from matplotlib.colors import LogNorm
# two parameters, input binary file and time (in seconds) to look for the nearest time step after.
# output is abundance of all isotopes at the time step found in the form "Z A ab" in file ab_$source_$time
chemsym=["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"]

record_length=0

if(len(sys.argv) > 1):
	f=open(sys.argv[1],'rb')
	record_length1=struct.unpack('<i', f.read(4))[0]
	desc1 = f.read(80)
	print desc1
	desc2 = f.read(80)
	print desc2
	desc3 = f.read(80)
	print desc3
	data_desc = f.read(80)
	print data_desc
	record_length2=struct.unpack('<i',f.read(4))[0]
	
	record_length1=struct.unpack('<i', f.read(4))[0]
	kstmx=struct.unpack('<i',f.read(4))[0]
	kitmx=struct.unpack('<i',f.read(4))[0]
	iweak=struct.unpack('<i',f.read(4))[0]
	iscrn=struct.unpack('<i',f.read(4))[0]
	iconvc=struct.unpack('<i',f.read(4))[0]
	changemx=struct.unpack('<d',f.read(8))[0]
	tolm=struct.unpack('<d',f.read(8))[0]
	tolc=struct.unpack('<d',f.read(8))[0]
	yacc=struct.unpack('<d',f.read(8))[0]
	ymin=struct.unpack('<d',f.read(8))[0]
	tdel_mm=struct.unpack('<d',f.read(8))[0]
	#print kstmk, " ", kitmx, " ", iweak
	record_length2=struct.unpack('<i',f.read(4))[0]
	
	record_length1=struct.unpack('<i',f.read(4))[0]
	abund_file  = f.read(80)
	print abund_file
	abund_desc = f.read(80)
	print abund_desc
	record_length2=struct.unpack('<i',f.read(4))[0]
	print record_length1, " ", record_length2
	
	record_length1=struct.unpack('<i',f.read(4))[0]
	thermo_file  = f.read(80)
	print thermo_file
	thermo_desc = f.read(80)
	print thermo_desc
	record_length2=struct.unpack('<i',f.read(4))[0]
	print record_length1, " ", record_length2
	
	record_length1=struct.unpack('<i',f.read(4))[0]
	print record_length1
	ny=struct.unpack('<i',f.read(4))[0]
	print ny
	zz=[0.]*(ny)
	aa=[0.]*(ny)
	for n in range(0,ny):
		zz[n]=struct.unpack('<d',f.read(8))[0]
	for n in range(0,ny):
		aa[n]=struct.unpack('<d',f.read(8))[0]
	record_length2=struct.unpack('<i',f.read(4))[0]
	
	record_length1=struct.unpack('<i',f.read(4))[0]
	print record_length1
	nflx=struct.unpack('<i',f.read(4))[0]
	print "nflx", nflx
	flx_end=[]
	if(nflx>0):
		for n in range(0,nflx):
			flx_end.append([])
			flx_end[n].append(struct.unpack('<i',f.read(4)))
			flx_end[n].append(struct.unpack('<i',f.read(4)))
	record_length2=struct.unpack('<i',f.read(4))[0]
	time=[]
	temperature=[]
	density=[]
	timestep=[]
	edot=[]
	xmf=[]
	flx=[]
	fxmf=[]*ny
	for k in range(0,kstmx):
		#print k
		try:
			record_length1=struct.unpack('<i',f.read(4))[0]
		except:
			break
		#print record_length1, "size"
		if record_length1==0:
			break
		kstep=struct.unpack('<i',f.read(4))[0]
		x = struct.unpack('<d',f.read(8))[0]
		#print x
		time.append(x)
		temperature.append(struct.unpack('<d',f.read(8))[0])
		density.append(struct.unpack('<d',f.read(8))[0])
		timestep.append(struct.unpack('<d',f.read(8))[0])
		edot.append(struct.unpack('<d',f.read(8))[0])
		xmf.append([])
		for i in xrange(0,ny):
			xmf[k].append(struct.unpack('<d',f.read(8))[0])
		fxmf=xmf[k]
		flx.append([])
		for j in xrange(0,nflx):
			flx[k].append(struct.unpack('<d',f.read(8))[0])
		record_length2=struct.unpack('<i',f.read(4))[0]
		#print "length 2", record_length2
	#for l in xrange(0,ny):
	#	print xmf[0][l]," ",zz[l], " ", aa[l], " ", fxmf[l]
	#print time
	f.close();
	if(len(sys.argv)>2):
		t=float(sys.argv[2])
		i=0
		while t>time[i] and i<len(time)-1:
			i=i+1
		w=open("ab_"+sys.argv[1]+"_"+str(time[i]),"w")
		for j in xrange(0,ny):
			w.write(str(int(zz[j]))+" "+str(int(aa[j]))+" "+str(xmf[i][j])+"\n")
		w.close()
fig=plt.figure()
ax = fig.add_subplot(111)
plot, = ax.plot([], [])
print "start"
def init():
	plot=fig.clf()
	return plot
def animate(k):
	i=0
	while time[i]<k:
		i=i+1
	print "animate", time[i]
	plt.clf()
	plt.ylim([0,70])
	plt.xlim([0,110])
	map=np.zeros((84,180))
	
	for l in xrange(0,ny):
		map[int(zz[l])][int(aa[l]-zz[l])]=xmf[i][l]
	plot=plt.pcolor(map,norm=LogNorm(vmin=1.0e-25, vmax=1.0e+0))
	plt.title(sys.argv[1]+" t="+str(time[i]))
	plt.colorbar()
	return plot
#writer = animation.MovieWriter()
ani = animation.FuncAnimation(fig,animate,np.arange(0, int(time[-1]),100000),blit=False,init_func=init)
#ani.save("ab_"+sys.argv[1]+".mp4", writer=writer,fps=15)
plt.show()