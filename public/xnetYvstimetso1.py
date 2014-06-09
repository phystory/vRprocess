import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_line(num, data, line):
    line.set_data(data[...,:num])
    print line
    return line,

fig1 = plt.figure()

args = []
legends = []
xmin = 0
xmax = 100
ymin = 1E-9
ymax = 1
ithmax = 0
with open('ev1','r') as ifile:
	for ith,line in enumerate(ifile):
		linemember=line.split()
		if ith == 0: 
			ithmax = len(linemember)
			for ith,value in enumerate(linemember):
				if ith > 5 and ith < ithmax: legends.append(value)
		else :
			index = float(linemember[0])
			for ith,value in enumerate(linemember):
				valuef = float(value)
				if index == 0: args.append([])
				args[ith].append(valuef)
plt.xlabel('time')
plt.ylabel('Mass fraction')
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.title('Mass fraction along time')
plt.yscale('log')
for i in range(6,ithmax-1):
	plt.plot(args[1], args[i], '-')
#line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l), interval=50, blit=True)
#line_ani.save('lines.mp4')

plt.legend(legends, loc = 'upper right', prop = {'size':10})
plt.show()
