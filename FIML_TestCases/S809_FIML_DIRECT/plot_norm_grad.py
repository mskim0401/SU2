#import csv
import numpy as np
import sys
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import os



hist_data = genfromtxt('history_project.dat', delimiter =',', skip_header=1)

# hist_data = hist_data[0:11,:]

#sys.stdout.write(str(hist_data)+'\n')

num_iters = len(hist_data[:,0])

sys.stdout.write('Number of Iterations Loaded: ' +str(num_iters)+'\n')

for i in range(0,num_iters) :
	
	if i < 9 :
		exists = os.path.isfile('DESIGNS/DSN_00'+str(i+1)+'/ADJOINT_INVERSE_DESIGN_LIFT_FIML/beta_fiml_grad.dat')
		if exists:
			curr_beta = genfromtxt('DESIGNS/DSN_00'+str(i+1)+'/ADJOINT_INVERSE_DESIGN_LIFT_FIML/beta_fiml_grad.dat')
		else :
			curr_beta = np.float('nan')
	else :
		exists = os.path.isfile('DESIGNS/DSN_0'+str(i+1)+'/ADJOINT_INVERSE_DESIGN_LIFT_FIML/beta_fiml_grad.dat')
		if exists :
			curr_beta = genfromtxt('DESIGNS/DSN_0'+str(i+1)+'/ADJOINT_INVERSE_DESIGN_LIFT_FIML/beta_fiml_grad.dat')
		else :
			curr_beta = np.float('nan')
		
		
	if i == 0 : norm_beta = np.linalg.norm(curr_beta,2)
	else :
		if exists :
			norm_beta = np.append(norm_beta,np.linalg.norm(curr_beta,2))
		else :
			norm_beta = np.append(norm_beta,np.float('nan'))

font = {'size'   : 14}
plt.rc('font', **font)


#fig.set_figheight(5)
#fig.set_figwidth(13)

#gs = gridspec.GridSpec(1, 1, width_ratios=[4,10], height_ratios=[4])
#ax0 = fig.add_subplot(gs[0])
#ax2 = fig.add_subplot(gs[1])



plt.rcParams.update({'font.size': 14})

sys.stdout.write(str(hist_data[:,0])+'\n')

sys.stdout.write(str(norm_beta)+'\n')

plt.figure(1)
#plt.subplot(1,2,1)
plt.semilogy(hist_data[:,0],norm_beta,'-*r',linewidth=2,markersize=8,label=r"$\norm \nabla J_c \norm _2$")
plt.xlabel('Evaluation')
plt.ylabel(r'$\Vert \nabla J_c \Vert _2$')
#plt.legend()
#plt.plot(hist_data[0,0],hist_data[0,1],'or',label="Design")
plt.grid()
plt.show()



