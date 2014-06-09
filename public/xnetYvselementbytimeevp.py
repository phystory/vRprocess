import numpy as np
from matplotlib  import pyplot as plt
from matplotlib import animation as animation

xmin = 0
xmax = 20
ymin = 1E-30
ymax = 1
ithmax = 0
args = []
indexes = []
legends = []


fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(ymin,ymax))
data, = ax.plot([], [], lw=2)

def init():
	data.set_data([], [])
	return data,

def animate(i):
	x = indexes
	y = args[i]
	data.set_xdata(np.arange(len(y)))
	data.set_ydata(y)
	return data,

with open('evp','r') as ifile:
	for ith,line in enumerate(ifile):
		linemember=line.split()
		if ith == 0: 
			ithmax = len(linemember)
			for jth,value in enumerate(linemember):
				if jth > 5 and jth < ithmax:
					legends.append(value)
					indexes.append(float(jth-5))
		else :
			args.append([])
			index = int(linemember[0]) - 1
			for jth,value in enumerate(linemember):
				valuef = float(value)
				if jth > 5 and jth < ithmax:
					args[index].append(valuef)
plt.xlabel('Mass Index')
plt.ylabel('Mass fraction')
#plt.xlim(xmin, xmax)
#plt.ylim(ymin, ymax)
plt.title('Mass fraction along time')
plt.yscale('log')
line_ani = animation.FuncAnimation(fig, animate, init_func=init, frames=124,interval=1, blit=True)
#line_ani.save('lines.mp4', fps=20)

#plt.legend(legends, loc = 'upper right', prop = {'size':10})
plt.show()
