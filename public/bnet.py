from numpy import *
from scipy import *
from pylab import *
from read_ascii import *
import matplotlib 
matplotlib.rc('text', usetex=True)
    
# files
file_path  = './'
file_plots = './' 

 
# axes
axes_left   = 0.15
axes_bottom = 0.17
axes_all  = [ axes_left, axes_bottom, 0.97-axes_left, 0.95-axes_bottom ]

# colors
cols = ['g',  'r', 'k', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y', \
        'k',  'g', 'r', 'b',  'orange','c','m', 'y' ]

linest = ['-','-','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']
# size
sizesmall = (6,4)
sizesmall = (5.25,3.5)
#sizesmall = (5,3) #talk
sizebig   = (12,8)
figformat = '.png'



#def plot_finabund(runs,figname):
def plot_finabund(runs,figname):

    figname = file_plots + figname
    for i in range(len(runs)):
        # read files
        filename = file_path# + runs[i]
        file     = runs[i]
#        file     = 'solar_mfrac_anders_grevesse.txt'
        print(file)
        i_z, i_a, y = read_ascii( file, [0,1,2])

        # build a matrix ZxN for the abundances
        a = r_[ 0:300. ] 
        z = r_[ 0:110. ] 
        abun = ones( ( 110, 300 ) ) * (0.)
        n    = ones( ( 110, 300 ) ) * (0.)
        for k in range(len(i_a)):
            abun[ int(i_z[k]), int(i_a[k]) ] = y[k]
            n   [ int(i_z[k]), int(i_a[k]) ] = i_a[k]-i_z[k]
        abun_a = abun.sum(axis=0)   # sum in A
        abun_z = abun.sum(axis=1)   # sum in Z
        
        
        
        # Abundance vs Z      
        figure(12,figsize=sizesmall)
        if i == 0:    # initial setup for the plot
            clf()
            axes(axes_all)

        semilogy( z, abun_z, c=cols[i], ls=linest[i])  
        if i == len(runs) - 1:
            xticks(r_[0:83:5])     # draw ticks in the x axis
#            axis([29,83,1e-10,20.])
            axis([0,83,1e-18,1.])
            xlabel('Z')
            ylabel('abundance')
            legend(loc=1)
        show()
        savefig( figname + '_YvsZ'+figformat )

        # Abundance vs A      
        figure(13,figsize=sizesmall)
        if i == 0:
            clf()
            axes(axes_all)
        semilogy(a,abun_a,c=cols[i],ls=linest[i])
        if i == len(runs) - 1:
            xticks(r_[0:203:20])
#            axis([63,203,1e-10,20.])
            axis([1,203,1e-18,1.])
            xlabel('A')
            ylabel('abundance')
            legend(loc=1)
        show()
        savefig( figname + '_YvsA'+figformat )

