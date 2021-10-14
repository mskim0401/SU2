-----------------------------------------------------------
SU2 FIML Code Axisymmety Modifided Version by M. S. Kim
-----------------------------------------------------------

In progeress:: bug fix and validation

This Branch is from SU2 FIML branch by J. Holland. <br />
To ensure axisymmetric problem in SU2 v5.0 that has a viscous source term issue, source code was modified from SU2 v7. Some bugs in SU2 v5.0 are also fixed.
To clear output log, some debugging outputs are eliminated.

See below pull requests.<br />
Viscous axisymmetric source term: https://github.com/su2code/SU2/pull/1106/<br />
Axisymmetric source term bug in energy equation: https://github.com/su2code/SU2/pull/1366<br />
Axisymmetric source term sign bug: https://github.com/su2code/SU2/pull/1366/<br />
SST 2D axisymmetric source term: https://github.com/su2code/SU2/pull/1195<br />
SST Jacobian bug: https://github.com/su2code/SU2/pull/491<br />
Axisymmetric source term implicit decision: https://github.com/su2code/SU2/pull/456<br />
Viscous Jacobian sign error: https://github.com/su2code/SU2/pull/612<br />
Symmetric and Euler BC issues: https://github.com/su2code/SU2/pull/657, https://github.com/su2code/SU2/pull/740<br />
::These modifications may make another problem at symmetric axis, so these are manually available only when AXISYMMETRIC= NO.<br />
Turbulent viscosity calculation in SST: https://github.com/su2code/SU2/pull/905<br />

Add missing term to strain magnitude for 2D: https://github.com/su2code/SU2/pull/670<br />
Periodic Green Gauss bug: https://github.com/su2code/SU2/pull/496<br />
Incompressible outlet BC and 'Set_MPI_AxiAuxVar_Gradient': SU2 v6<br />
Fixed Coord_j for boundary viscous numerics: https://github.com/su2code/SU2/pull/1189<br />

CFL reduction factors has been modified to be available for Discrete Adjoint not only Continuous Adjoint.<br />
Low Mach Pre-conditioning is available and modified as version 6. (not recommemded to use)

2D Axisymmetric SA equation is now conducted with the help of S. Heo.
Jacobian of axisymmetric source term for adjoint calculation in compressible N-S and SA equations are added to see how affects adjoint solution. (need to validate)

Currently, only compressible flow is guaranteed.

The following is readme provided by J. Holland:

-----------------------------------------------------------
  SU2 Field Inversion and Machine Learning (FIML) Development Branch
-----------------------------------------------------------

Branch of SU2 v5.0 to support research on Field Inversion and Machine Learning. FIML modifications performed by myself (Jon Holland), and results are documented on this ResearchGate Page: https://www.researchgate.net/profile/Jonathan_Holland5

This branch is not stable, development of FIML capabilities in SU2 is ongoing. Development was performed to support my dissertation research. I intend to include test cases for FIML as well as a short tutorial in the near future. If interested in using SU2 for FIML, please message me.

My researchgate page containing publications documetning the numerical methodology and results is here: https://www.researchgate.net/profile/Jonathan_Holland5

A presentation about this effort was presented at the 110th NIA CFD seminar, recording and slides available here: http://ossanworld.com/hiroakinishikawa/niacfds/

Many thanks to all of the SU2 developers!

-Jon Holland


-----------------------------------------------------------
Overview of Modifications to the Configuration File for FIML Problems
-----------------------------------------------------------

The problem definition is largely the same as for shape optimization problems with the autodifferentiated discrete adjoint solver. In the FIML problem, the design variables have been re-defined. For FIML-Classic, the design variables are a correction to the production term on the turbulence model. For FIML-Embedded, the design variables are a target correction field that a neural network is then trained to match. For FIML-Direct, the design variables are the weights of the neural network. The problem definition is still controled by the *.cfg file with some additional options, and the solution is started and controlled by the shape_optimization.py script. The user should be familiar with SU2 shape optimization problems before attempting a FIML case. The primary changes to the SU2 config file options are listed below. This is not a comprehensive list:

The "SA_FIML" turbulence model option specifies that we are performing a FIML problem (Classic, Embedded, or Direct). This is a modified "SA" model with a correction to the right hand side of the turbulence equation.

% Specify turbulent model (NONE, SA, SA_NEG, SST, SA_FIML) <br />
KIND_TURB_MODEL= SA_FIML <br />

There are several added objective functions including: INVERSE_DESIGN_PRESSURE_FIML, INVERSE_DESIGN_LIFT, INVERSE_DESIGN_LIFT_FIML, INVERSE_DESIGN_DRAG, INVERSE_DESIGN_DRAG_FIML. The INVERSE_* indicates that the user will specify a TARGET_C* value which is the higher fidelity data we want our turbulence model to be corrected to match. The *_FIML objective functions include a regularization parameter that penalizes overly intrusive corrections to the turbulence model. <br />
OBJECTIVE_FUNCTION= INVERSE_DESIGN_LIFT_FIML

FIML cases are restricted to SU2 for now. FIML-Classic design variable initialization requires SU2 format grid.
% Mesh input file format (SU2, CGNS, NETCDF_ASCII)<br />
MESH_FORMAT= SU2<br />


This is another flag to indicate that we are doing a FIML problem, so set DV_KIND to "FIML" for all FIML types (Classic, Embedded, and Direct). For all FIML cases the meaning of DV_PARAM and DV_VALUE has changed as it would be impossible to write the starting value of all the design variables by hand. (FIML problems are very high dimensional / thousands of design variables). For FIML-Classic this sets the value of the correction at every point, so DV_VALUE= 0.0 would be no turbulent production (not recommended), and DV_VALUE=1.0 would be the baseline SA model (no correction, recommended starting point).<br />
% Kind of deformation (FFD_SETTING, FFD_CONTROL_POINT_2D, FFD_CAMBER_2D, FFD_THICKNESS_2D,<br />
%                      HICKS_HENNE, PARABOLIC,<br />
%                      NACA_4DIGITS, DISPLACEMENT, ROTATION, FFD_CONTROL_POINT,<br />
%                      FFD_NACELLE, FFD_TWIST, FFD_ROTATION,<br />
%                      FFD_CAMBER, FFD_THICKNESS, SURFACE_FILE, FIML)<br />
DV_KIND= FIML<br />
DV_PARAM= ( 1.0 )<br />
DV_VALUE= 1.0<br />

Objective function definition is the same as standard SU2 with a scaling factor. Note that for FIML-Classic and Embedded the gradients are extremely small (becuase changing the production term at any single point in the domain only has a small influence on the objective function). Some of the SciPy optimizers (L-BFGS-B and BFGS for example) take a unit step in the descent direction. For very small gradients this is not a substantial step and will result in no change of the objective function unless the objective function is scaled properly. This will be problem dependent, but scale factors of 1.0E16 was not uncommon for FIML-Classic cases. Essentially the scale factor should be chosen to have a small, but significant first step in the descent direction. In practice, if the largest gradients were of the order 1.0E-17, then an objective function scale factor of 1.0E16 typically worked well. The same is true for FIML-Embedded. For FIML-Direct the scale factors were much smaller, on the order of 1.0E-2 depending on the chosen neural network structure and the test case.<br />
OPT_OBJECTIVE= INVERSE_DESIGN_LIFT_FIML * 1.0E-2<br />

Again, for FIML cases these parameters set limits for all of the design variables with a single number:<br />
OPT_CONSTRAINT= NONE<br />
OPT_BOUND_UPPER = 2.0<br />
OPT_BOUND_LOWER = -2.0<br />
DEFINITION_DV= ( 103, 1.0 | airfoil | 1.0 )<br />

This parameter has to be set for FIML problems. It is the value defined in the SU2 grid file, and should be copied into the config file. This is a hack that I'm not really happy about, but it works. It is required to tell the shape_optimization.py script how many design variables we are using so that it can initialize the starting value of the deisgn variable vector. Note that NPOIN= 59400 means there are 59400 design variables for the FIML-Classic and FIML-Embedded problems. For FIML-Direct the number of design variables is determined by the chosen structure of the neural network. Again, NPOIN must be copied from the SU2 grid file.<br />
NPOIN= 59400<br />

This sets the scale factor for the regularization portion of the objective function, only needs to be set when using the *_FIML objective functions (that have the regularization term included)<br />
LAMBDA_FIML = 1.0E-4<br />

This sets the target lift coefficient (lift coefficient we want our augmented turbulence model to match). Also have TARGET_INVERSE_CD and TARGET_INVERSE_CP available<br />
TARGET_INVERSE_CL = 1.0546<br />

This tells SU2 that we are doing either FIML-Embedded or FIML-Direct and therefore will be using the neural network portion of the code in solver_direct_turbulent.cpp. Defaults to "NO" for either FIML-Classic or no FIML cases.<br />
TRAIN_NN= YES<br />

Set to "BACKPROP" for FIML-Embedded problems, or "WEIGHTS" for FIML-Direct.<br />
KIND_TRAIN_NN= WEIGHTS<br />

Required for FIML Embedded and Direct problems, sets the number of neurons in the hidden layers of the network<br />
N_NEURONS= 20<br />

Sets the scaling for the inputs to the neural network. Options are NO_SCALE, MIN_MAX, Z_SCALE, Q_TRANSFORM, BOX_COX. BOX_COX is recommended, but parameters are hard-coded for now. Q_TRANSFORM applies the inverse CDF of a normal distribution and may be the most robust but is substantially more involved than the Box-Cox scaling. Note that the Box-Cox scaling parameters were calibrated for the S809 test case and may not be optimal for other cases. Input scaling can dramatically affect the performance for FIML-Embedded and FIML-Direct as there are HUGE outliers in the (unscaled) chosen neural network inputs that must be properly mitigated or the convergence of the algorithms is affected.<br />
KIND_NN_SCALING= BOX_COX<br />

Only required for Q_TRANSFORM scaling, sets the number of bins for the CDF estimate. In practice results do not appear to be senstive to this value, typically used 10 or 20 bins.<br />
N_BINS= 10<br />

Only used in FIML-Embedded, and it sets the learning rate for the backpropagation algorithm.<br />
LEARNING_RATE= 0.01<br />

This is very confusing I know, but LEARN_RATE sets the step size for the STEEPJ optimizer option. STEEPJ option (set by shape_optimization.py) will take a step in the steepest descent direction of size LEARN_RATE. This is a very inefficient optimizer, but was useful for testing.
%This one is for STEEPJ<br />
LEARN_RATE= 0.01<br />

This parameter sets the number of iterations that should be performed by the flow solver prior to using the neural network for either FIML-Embedded or FIML-Direct. In practice, this should be the smallest number of iterations that doesn't blow up the flow solver. If it is too low the neural network can output large corrections in the initial iterations which can cause the flow solver to not converge. 1000 worked well for the S809 cases, but could be lower or higher depending on the problem.<br />
ITER_START_NN= 1000<br />

This sets the number of backpropagation epochs (iterations) that should be done per flow solver. This feature is almost completely untested so use with caution. All cases so far have used NUM_EPOCH= 1<br />
NUM_EPOCH= 1<br />

This parameter sets the number of hidden layers in the neural network. Not tested for networks with more than 3 hidden layers (N_HIDDEN_LAYERS= 3)<br />
N_HIDDEN_LAYERS= 3<br />

This sets whether we should filter out points in the domain that we don't think are relevant to the problem. This limits the dimensionality of the problem and potentially can improve convergence of FIML-Embedded and FIML-Direct. Filter parameters are hard-coded and calibrated to the S809 airfoil case.<br />
FILTER_SHIELD= YES<br />

Set to YES for multi-case FIML-Direct problems with different meshes. This feature still in development, but has worked.<br />
MULTI_MESH= NO<br />



Calling the modified shape_optimization.py script for FIML cases:<br />
Example command line to run a FIML-Classic Case with 4 processors. Note -g DISCRETE_ADJOINT must be used. Also note that depending on the problem, the limited memory BFGS option (L-BFGS-G) must be used specified by -o BFGSG. For FIML-CLassic the number of design variables is huge, and constructing the full Hessian estimate in BFGS can require huge amounts of memory and is very slow:<br />
shape_optimization.py -n 4 -g DISCRETE_ADJOINT -f config.cfg -o BFGSG<br />

A case can be restarted in the usual way by adding -r project.pkl but the value of the design variable vector must also be provided by -b beta_fiml.dat. This can be found by examining the objective function convergence and selecting the beta_fiml.dat file from the best Direct solution in the convergence history. So for a restart we would enter:<br />
shape_optimization.py -n 4 -g DISCRETE_ADJOINT -f config.cfg -o BFGSG -r project.pkl -b beta_fiml_best.dat<br />

The -b option is also useful for initializing any FIML method design variables to a set starting value. For FIML-Classic this could be the output of a machine learned model (completing the "ML" part offline from SU2), or for FIML-Direct this could be used to initialize the neural network to a previously trained value. This is particularly useful for multi-case FIML-Direct applications as new cases can be added and training continued using the previously trained network weights. Also, sometimes the optimizers take extremely large line search steps that can cause convergence issues. If this occurs the solution can be restarted using a previous step and typically the optimizer won't follow the same path as the Hessian estimate is reset by a restart.<br />


-----------------------------------------------------------
  Field Inversion and Machine Learning SU2 Tutorials
-----------------------------------------------------------

Three tutorials are provided for FIML-Classic, single case FIML-Direct, and multi-case FIML-Direct. A test case for FIML-Embedded is in the works. The multi-case FIML-Direct is still in active development and will be updated in the future. THe cases were all run on UMD's Deepthought2 HPC with the submission script "job.sub". The last line on the job.sub file is the command for the shape_optimization.py script for that case. The solution files are not provided, except for the "history_project.dat" which should give an idea of the convergence history expeced by the optimizer. Note that quite a bit of disk space is required (similar to shape_optimization problems), especially for the AoA_Seven case which is running 7 DIRECT solution and 7 ADJOINT solutions for each optimizer iteration. A multi-mesh FIML-Direct case is also in development. All of the test cases use experimentally obtained lift coefficient as the higher fidelity data (TARGET_INVERSE_CL). The AoA_Seven case considers 7 separate angles of attack simultaneously to train the neural network augmentation.<br />


The following is readme provided by the SU2 developers:<br />

-----------------------------------------------------------
  SU2 (ver. 5.0.0 "Raven"): The Open-Source CFD Code
-----------------------------------------------------------

Computational analysis tools have revolutionized the way we design aerospace systems, but most established codes are proprietary, unavailable, or prohibitively expensive for many users. The SU2 team is changing this, making computational analysis and design freely available as open-source software and involving everyone in its creation and development.

For an overview of the technical details in SU2, please see the following AIAA Journal article:

"SU2: An open-source suite for multiphysics simulation and design," AIAA Journal, 54(3):828-846, 2016. http://arc.aiaa.org/doi/10.2514/1.J053813

[![Build Status](https://travis-ci.org/su2code/SU2.svg?branch=develop)](https://travis-ci.org/su2code/SU2)

----------------------------------------------------------
  SU2 INTRODUCTION
----------------------------------------------------------

SU2 is a suite of open-source software tools written in C++ for the numerical solution of partial differential equations (PDE) and performing PDE constrained optimization.

The primary applications are computational fluid dynamics and aerodynamic shape optimization, but has been extended to treat more general equations such as electrodynamics and chemically reacting flows.

You will find more information and the latest news in:
   - GitHub:    https://github.com/su2code
   - CFD-online http://www.cfd-online.com/Forums/su2/
   - Twitter:   https://twitter.com/su2code

---------------------------------------------------
  SU2 INSTALLATION
---------------------------------------------------

To build SU2 from the source code, just open a terminal and run the './configure', 'make', and 'make install' commands in the root directory of the source distribution. You can provide an install location using the prefix option to configure. If there are issues with autotool version requirements, run the ./bootstrap script provided in the root directory in order to build local versions of the tools and reset the makefiles (before trying configure/make/make install again). Please note that more detailed instructions on the configure and build processes can be found within the INSTALL file.

----------------------------------------------------------
  SU2 PATH SETUP
----------------------------------------------------------

SU2 is built using a typical configure/make/make install process. When make install is complete, please be sure to add the $SU2_HOME and $SU2_RUN environment variables, and update your $PATH with $SU2_RUN.

For example, add these lines to your .bashrc file:

- export SU2_RUN="your_prefix/bin"
- export SU2_HOME="/path/to/SU2vX.X.X/"
- export PATH=$PATH:$SU2_RUN
- export PYTHONPATH=$SU2_RUN:$PYTHONPATH

$SU2_RUN should point to the folder where all binaries and python scripts were installed. This is the prefix you set with the --prefix option to configure. Note that the bin/ directory is automatically added to your prefix path.

$SU2_HOME should point to the root directory of the source code distribution, i.e., /path/to/SU2vX.X.X/.

Thanks for building, and happy optimizing!

- The SU2 Development Team

----------------------------------------------------------
  SU2 DEVELOPERS
----------------------------------------------------------

SU2 is being developed by individuals and organized teams all around the world.

The SU2 Lead Developers are:

   - Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com)
   - Dr. Thomas D. Economon (economon@stanford.edu)

and the most active groups developing SU2 are:

   - Prof. Juan J. Alonso's group at Stanford University.
   - Prof. Piero Colonna's group at Delft University of Technology.
   - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
   - Prof. Alberto Guardone's group at Polytechnic University of Milan.
   - Prof. Rafael Palacios' group at Imperial College London.
   - Prof. Edwin van der Weide's group at the University of Twente.
   - Prof. Vincent Terrapon's group at the University of Liege.
