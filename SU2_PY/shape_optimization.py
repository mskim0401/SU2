#!/usr/bin/env python 

## \file shape_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author T. Economon, T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os, sys, shutil, copy
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2
import numpy as np

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-r", "--name", dest="projectname", default='',
                      help="try to restart from project file NAME", metavar="NAME")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="CONTINUOUS_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-o", "--optimization", dest="optimization", default="SLSQP",
                      help="OPTIMIZATION techique (SLSQP, CG, BFGS, POWELL)", metavar="OPTIMIZATION")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    parser.add_option("-c", "--n_configs", dest="num_configs", default=1,
                      help="Number of config files (for weights NN FIML case)", metavar="NUM_CONFIGS")
    parser.add_option("-b", "--baseline_nn", dest="weights_file", default=None,
                      help="Baseline Weights File for Neural Network (If none defaults to random init)", metavar="BASE_NETWORK")
    
    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    
    sys.stdout.write('\n-------------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )   Release 5.0.0 \"Raven\"                             |\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                      |\n')
    sys.stdout.write('|   |___/\\___//___|   Aerodynamic Shape Optimization Script             |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Lead Dev.: Dr. Francisco Palacios, Francisco.D.Palacios@boeing.com|\n')
    sys.stdout.write('|                Dr. Thomas D. Economon, economon@stanford.edu          |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Developers:                                                       |\n')
    sys.stdout.write('| - Prof. Juan J. Alonso\'s group at Stanford University.                |\n')
    sys.stdout.write('| - Prof. Piero Colonna\'s group at Delft University of Technology.      |\n')
    sys.stdout.write('| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. |\n')
    sys.stdout.write('| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  |\n')
    sys.stdout.write('| - Prof. Rafael Palacios\' group at Imperial College London.            |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| Copyright (C) 2012-2017 SU2, the open-source CFD code.                |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is free software; you can redistribute it and/or                  |\n')
    sys.stdout.write('| modify it under the terms of the GNU Lesser General Public            |\n')
    sys.stdout.write('| License as published by the Free Software Foundation; either          |\n')
    sys.stdout.write('| version 2.1 of the License, or (at your option) any later version.    |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is distributed in the hope that it will be useful,                |\n')
    sys.stdout.write('| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n')
    sys.stdout.write('| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n')
    sys.stdout.write('| Lesser General Public License for more details.                       |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| You should have received a copy of the GNU Lesser General Public      |\n')
    sys.stdout.write('| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')

    shape_optimization( options.filename     ,
                        options.projectname  ,
                        options.partitions   ,
                        options.gradient     ,
                        options.optimization ,
                        options.quiet        ,
                        options.num_configs  ,
                        options.weights_file)
    
#: main()

def shape_optimization( filename                           ,
                        projectname = ''                   ,
                        partitions  = 0                    ,
                        gradient    = 'CONTINUOUS_ADJOINT' ,
                        optimization = 'SLSQP'             ,
                        quiet       = False                ,
                        num_configs = 1                    ,
                        weights_file = None):
  
    # Config
    
    config = SU2.io.Config(filename)
    
    #if num_configs > 1 :
#	config_vec = config #Use first configuration as default for rest of script JRH 09272018
#	for i in range(2,num_configs) :
#	    config_vec.append(SU2.io.Config(str(i)+'.cfg'))
#	    config_vec[i].NUMBER_PART = partitions
#	    config_vec[i].GRADIENT_METHOD = 'gradient'
#	    if quiet: config_vec[i].CONSOLE = 'CONCISE'
    
    config.NUMBER_PART = partitions
    
    if quiet: config.CONSOLE = 'CONCISE'
    
    config.GRADIENT_METHOD = gradient
    
    its         = int ( config.OPT_ITERATIONS )
    accu        = float ( config.OPT_ACCURACY )
    bound_upper = float ( config.OPT_BOUND_UPPER )
    bound_lower = float ( config.OPT_BOUND_LOWER )
    def_dv      = config.DEFINITION_DV
    kind_dv     = config.DV_KIND
    if config.NPOIN == 0:
        n_dv        = sum(def_dv['SIZE'])
        x0          = [0.0]*n_dv
    else :
	if config.KIND_TRAIN_NN == 'WEIGHTS':
	    #ANY CHANGES NEED TO BE MIMICKED IN PROJECTION.PY
	    if weights_file == None : 
		layers = int(config.N_HIDDEN_LAYERS)+2
		num_nodes = [5] #5 if using bias nodes
		for iLayer in range(1,layers-1) :
			num_nodes.append(int(config.N_NEURONS))
		num_nodes.append(1) #2 if using bias nodes
		num_inputs = [5] #5 if using bias nodes
		#num_inputs.append(int(5))
		for iLayer in range(1,layers) :
		    num_inputs.append(int(config.N_NEURONS))
		n_dv = 0
		x0 = []
		for iLayer in range(1,layers) :
		    for iNode in range(num_nodes[iLayer]) :
			for iInput in range(num_inputs[iLayer-1]) :
			    #x0.append(np.random.randn()/np.sqrt(5.0)/10.0+0.1)
			    #x0.append(0.0)
			    if iLayer > 0 :  #JRH 09232018 - Removing input layer weights from costly computation
				x0.append(np.random.randn()*np.sqrt(1.0/num_inputs[iLayer-1]))
				n_dv+=1
		sys.stdout.write('Total Number of Weights Set to ' + str(n_dv) + ' in shape_optimization.py\n')
	    else :
		sys.stdout.write('Reading weights from' + str(weights_file) + '\n')
		weights = open(weights_file)
		x0 = []
		lines = weights.readlines()     #i am reading lines here
		n_dv = 0
		for line in lines:            #taking each line
		    x0.append(float(line))         #converting string to float	
		    n_dv = n_dv+1
		    #sys.stdout.write('n_dv ' + str(x0[n_dv-1]) +'\n')
		sys.stdout.write('Total Number of Weights Set to ' + str(n_dv) + ' in shape_optimization.py \n')
	else:
	    if weights_file == None :
		n_dv = int (config.NPOIN)
		x0          = [float(config.DV_VALUE[0])-1.0]*n_dv # initial design
	    else :
		sys.stdout.write('Reading betas from' + str(weights_file) + '\n')
		weights = open(weights_file)
		x0 = []
		lines = weights.readlines()     #i am reading lines here
		n_dv = 0
		for line in lines:            #taking each line
		    x0.append(float(line))         #converting string to float	
		    n_dv = n_dv+1
		    #sys.stdout.write('n_dv ' + str(x0[n_dv-1]) +'\n')
		sys.stdout.write('Total Number of Betas Set to ' + str(n_dv) + ' in shape_optimization.py \n')
		
    xb_low      = [float(bound_lower)]*n_dv # lower dv bound
    xb_up       = [float(bound_upper)]*n_dv # upper dv bound
    xb          = zip(xb_low,xb_up) # design bounds
    
    sys.stdout.write('NPOIN = %s\n' % config.NPOIN)
    
    # State
    state = SU2.io.State()
    state.find_files(config)
    
    # Project
    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project.config = config
    else:
	#if num_configs > 1 :
	    #project = SU2.opt.Project(config,state, config_vector = config_vec)
	#else :
	project = SU2.opt.Project(config,state)
    
    # Optimize
    if optimization == 'SLSQP':
      SU2.opt.SLSQP(project,x0,xb,its,accu)
    if optimization == 'CG':
      SU2.opt.CG(project,x0,xb,its,accu)
    if optimization == 'BFGS':
      SU2.opt.BFGS(project,x0,xb,its,accu)
    if optimization == 'BFGSJ':
      SU2.opt.BFGSJ(project,x0,xb,its,accu)  
    if optimization == 'STEEPJ':
      SU2.opt.STEEPJ(project,x0,xb,its,accu) 
    if optimization == 'POWELL':
      SU2.opt.POWELL(project,x0,xb,its,accu)
    if optimization == 'BFGSG':
      SU2.opt.BFGSG(project,x0,xb,its,accu) 
    if optimization == 'STEEP':
      SU2.opt.STEEP(project,x0,xb,its,accu)
    if optimization == 'TNC':
      SU2.opt.TNC(project,x0,xb,its,accu)      

    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)
    
    return project

#: shape_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

