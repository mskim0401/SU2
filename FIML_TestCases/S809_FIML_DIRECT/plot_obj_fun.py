#import csv
import numpy as np
import sys
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec



hist_data = genfromtxt('history_project.dat', delimiter =',', skip_header=1)

# hist_data = hist_data[0:11,:]

#sys.stdout.write(str(hist_data)+'\n')

num_iters = len(hist_data[:,0])

sys.stdout.write('Number of Iterations Loaded: ' +str(num_iters)+'\n')

font = {'size'   : 14}
plt.rc('font', **font)

fig, ax = plt.subplots()

ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

#fig.set_figheight(5)
#fig.set_figwidth(13)

#gs = gridspec.GridSpec(1, 1, width_ratios=[4,10], height_ratios=[4])
#ax0 = fig.add_subplot(gs[0])
#ax2 = fig.add_subplot(gs[1])

def update(i) :
	#plt.figure(1)
	#ax0 = plt.subplot(1,2,1)
	#plt.plot(hist_data[:,0],hist_data[:,1],'-x',label="Lift")
	ax0.plot(hist_data[i-1,0],hist_data[i-1,1],'or',label="Design")
	ax0.set_ylabel('Lift Coefficient', color='r')
	ax0.tick_params('y', colors='r')
	ax0.set_xlabel('Evaluation',color = 'k')
	ax0.set_ylim(1.0,1.17)
	

	
	#ax2 = plt.subplot(1,2,2)
	if i < 10 : img = mpimg.imread('Beta_Plots/DSN_00'+str(i)+'.png')
	else : img = mpimg.imread('Beta_Plots/DSN_0'+str(i)+'.png')
	imgplot = ax2.imshow(img)
	ax2.axes.get_xaxis().set_visible(False)
	ax2.axes.get_yaxis().set_visible(False)
	ax2.set_xlim(182,1231)
	ax2.set_ylim(900,250)
	
	#plt.xlim(182,1231)
	#plt.ylim(900,250)
	
	return plt

plt.rcParams.update({'font.size': 14})

plt.figure(1)
plt.subplot(1,2,1)
plt.plot(hist_data[:,0],hist_data[:,1],'-*r',linewidth=2,markersize=8,label=r"$C_L$")
plt.plot(hist_data[:,0],np.ones((num_iters,))*1.0546,'-r',label=r"$C_{L_{exp}}$")
plt.xlabel('Evaluation')
plt.ylabel(r'Lift Coefficient ($C_L$)')
plt.legend()
#plt.plot(hist_data[0,0],hist_data[0,1],'or',label="Design")
plt.grid()

plt.subplot(1,2,2)
plt.semilogy(hist_data[:,0],hist_data[:,13],'-*',linewidth=2,markersize=8,label=r"$J_c$")
plt.semilogy(hist_data[:,0],hist_data[:,15],'--o',linewidth=2,markersize=8,label=r"$J_c'$")
#ax1 = ax0.twinx()
#plt.semilogy(hist_data[:,0],hist_data[:,13],'ob',label="OF")
#ax1.plot(hist_data[:,0],hist_data[:,13],'ob',label="Design_OF")
plt.ylabel('Objective Function')
plt.xlabel('Evaluation')
plt.grid()
#plt.tick_params('y', colors='b')
plt.legend()
plt.show()

#img = mpimg.imread('Beta_Plots/DSN_0'+str(num_iters)+'.png')
#imgplot = ax2.imshow(img)

# FuncAnimation will call the 'update' function for each frame; here
# animating over 10 frames, with an interval of 200ms between frames.
#anim = FuncAnimation(fig, update, frames=np.arange(1, num_iters+1), interval=600)
	
#anim.save('Beta_Plots/conv_beta.gif', dpi=150, writer='imagemagick')

