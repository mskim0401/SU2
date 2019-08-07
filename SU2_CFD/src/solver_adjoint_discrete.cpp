/*!
 * \file solver_adjoint_discrete.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/solver_structure.hpp"

CDiscAdjSolver::CDiscAdjSolver(void) : CSolver () {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config)  : CSolver() {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  unsigned short iVar, iMarker, iDim;

  bool restart = config->GetRestart();
  train_NN = false;
  jrh_debug = false;
  unsigned long iVertex, iPoint, index;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;
  su2double dull_val;
  bool done_before;

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  int rank = MASTER_NODE;
  int size;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  unsigned long nDV = config->GetnDV(); //JRH - 05082017
  nDV_total = nDV;
  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  Sens_Geo   = NULL;
  Sens_Mach  = NULL;
  Sens_AoA   = NULL;
  Sens_Press = NULL;
  Sens_Temp  = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;


  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Allocate the node variables ---*/

  node = new CVariable*[nPoint];


  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CSensitivity[iMarker]        = new su2double [geometry->nVertex[iMarker]];
  }

  Sens_Geo  = new su2double[nMarker];
  Sens_Mach = new su2double[nMarker];
  Sens_AoA  = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Sens_Geo[iMarker]  = 0.0;
      Sens_Mach[iMarker] = 0.0;
      Sens_AoA[iMarker]  = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          CSensitivity[iMarker][iVertex] = 0.0;
      }
  }

  //JRH - 04262017 Not sure where this should be, but I have a hunch it shouldn't be in RegisterVariables...
  if (KindDirect_Solver == RUNTIME_TURB_SYS) {
	  //Val_Beta_Fiml = new su2double[nPoint+1]; //This worked in serial
	  //Total_Sens_Beta_Fiml = new su2double[nPoint+1]; //This worked in serial
	  Val_Beta_Fiml = new su2double[nDV]; //JRH 05082017
	  Total_Sens_Beta_Fiml = new su2double[nDV]; //JRH 05082017
	  Local_Sens_Beta_Fiml = new su2double[nDV]; //JRH 09122017 - Made global variable for speed...
  }

  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  Local2Global = new unsigned long[nPointDomain];//JRH 06122018
  if (!restart || (iMesh != MESH_0)) {
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {//JRH 06122018
	      Local2Global[iPoint] = geometry->node[iPoint]->GetGlobalIndex();
	}
    /*--- Restart the solution from zero ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);
  }
  else {

    /*--- Restart the solution from file information ---*/
    mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no file ---*/
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }

    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    map<unsigned long,unsigned long> Global2Local; //Made protected variables of class so they could be accessed below - JRH 05102017
    map<unsigned long,unsigned long>::const_iterator MI;
    
    /*--- Now fill array with the transform values only for local points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
      Local2Global[iPoint] = geometry->node[iPoint]->GetGlobalIndex();
    }

    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- Skip coordinates ---*/
    unsigned short skipVars = nDim;

    /*--- Skip flow adjoint variables ---*/
    if (Kind_Solver == RUNTIME_TURB_SYS) {
      if (compressible) {
        skipVars += nDim + 2;
      }
      if (incompressible) {
        skipVars += nDim + 1;
      }
    }

    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPoint(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
        
        point_line >> index;
        for (iVar = 0; iVar < skipVars; iVar++) { point_line >> dull_val;}
        for (iVar = 0; iVar < nVar; iVar++) { point_line >> Solution[iVar];}
        node[iPoint_Local] = new CDiscAdjVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
      }
      
    }


    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }

    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);
    }

    /*--- Close the restart file ---*/
    restart_file.close();

  }

  /*--- Store the direct solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution_Direct(direct_solver->node[iPoint]->GetSolution());

    //Required if doing NN Training JRH 05022018
    if(config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
    	//Changed to not do this if == WEIGHTS, 07292019, not sure why I did this... Necessary? Bug? JRH
    	node[iPoint]->SetBetaFiml(direct_solver->node[iPoint]->GetBetaFiml()); //05022018
    	node[iPoint]->SetBetaFimlGrad(1.0e-16);
    }
  }

 //##################################################################################################
  //JRH - Broadcast nPointDomainProcs to all procs so that we know how to map the FIML design variables
  //Below code counts number of points in domain on each proc and sends it to all other procs
  //Variable is then created on local proc that corresponds to the starting index for that procs DVs - 05102017
  Fiml_Skip_Index = 0;
#ifdef HAVE_MPI
  //if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_TURB_SYS)) {
  if (KindDirect_Solver == RUNTIME_TURB_SYS) {
  //if (config->GetKind_Turb_Model() == SA_FIML) {
	  unsigned long *nPointFiml_Local;
	  unsigned long *nPointFiml;
	  nPointFiml_Local = new unsigned long[size];
	  nPointFiml = new unsigned long[size];
	  for (unsigned short iRank = 0; iRank < size ; iRank++) nPointFiml_Local[iRank] = 0;
	  for (iPoint = 0; iPoint < nPoint; iPoint++) {
		  if (geometry->node[iPoint]->GetDomain()) nPointFiml_Local[rank]++;
	  }
	  //SU2_MPI::Allreduce(&Local_Sens_Beta_Fiml[iDV],  &Total_Sens_Beta_Fiml[iDV],  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  for (unsigned short iRank = 0; iRank < size ; iRank++) {
		  SU2_MPI::Allreduce(&nPointFiml_Local[iRank],  &nPointFiml[iRank],  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	  }
	  for (unsigned short iRank = 0 ; iRank < rank ; iRank++) Fiml_Skip_Index += nPointFiml[iRank];
	  unsigned long Total_Fiml_Nodes = 0;
	  for (unsigned short iRank = 0 ; iRank < size ; iRank++) Total_Fiml_Nodes += nPointFiml[iRank];
	  //Fiml_Skip_Index += rank;
	  //cout << "JRH Debugging: MPI Process " << rank << " nPointDomain = " << nPointDomain << " nPoint = " << nPoint << " and begins with fiml variable " << Fiml_Skip_Index << endl;
	  if (rank == MASTER_NODE) cout << "Total_Fiml_Nodes = " << Total_Fiml_Nodes << endl;

	  //JRH - Every "new" needs a "delete"!! JRH 12182017
	  delete [] nPointFiml_Local;
	  delete [] nPointFiml;
  }
#else
#endif


  //################################################NEW BETA MAPPING CODE#########################################################
//  map<unsigned long,unsigned long> Global2Local; //Made protected variables of class so they could be accessed below - JRH 05102017
//  map<unsigned long,unsigned long>::const_iterator MI;
//  long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
//  /*--- Now fill array with the transform values only for local points ---*/
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
//    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
//  }
//  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPoint(); iPoint_Global++ ) {
//
////    /getline (restart_file, text_line);
//
//    //istringstream point_line(text_line);
//
//    /*--- Retrieve local index. If this node from the restart file lives
//     on the current processor, we will load and instantiate the vars. ---*/
//
//    MI = Global2Local.find(iPoint_Global);
//    if (MI != Global2Local.end()) {
//
//      iPoint_Local = Global2Local[iPoint_Global];
//
//     // point_line >> index;
//      //for (iVar = 0; iVar < skipVars; iVar++) { point_line >> dull_val;}
//      //for (iVar = 0; iVar < nVar; iVar++) { point_line >> Solution[iVar];}
//      //node[iPoint_Local] = new CDiscAdjVariable(Solution, nDim, nVar, config);
//      if(!config->GetTrainNN()) direct_solver->node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(iPoint_Global,0));
//      direct_solver->node[iPoint_Local]->SetBetaFimlTrain(config->GetDV_Value(iPoint_Global,0));
//      iPoint_Global_Local++;
//
//    }
//
//  }
  //################################################NEW BETA MAPPING CODE###############################################

  //JRH 04302018 - Initialize NN variables so we can store derivative (adjoint) information every Discrete Adjoint Iteration
  //Extract Adjoint in ExtractAdjoint_Solution()
  //Set Adjoint in SetAdjoint_Output
  //Initialize Weight Values to 1.0 (SetDerivative() called before ComputeAdjoint())
  //Initializing all variables for now in case we need later...could be cleaned up a bunch!!
  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
	  train_NN = true;
	  if (jrh_debug) cout << "JRH Debug: In solver_adjoint_discrete.cpp Beginning to set neural network vars" << endl;
 	  num_epoch = config->GetNumEpoch();
 	  learn_rate = config->GetLearningRate();
 	  num_nn_inputs = 4; //Hard-coded number of neural network inputs
 	  nLayers = config->GetNHiddenLayers()+2;//0 - Input Layer, 2<->nLayers-2 - Hidden Layers, nLayer-1 - Output Layer
 	  nNeurons = config->GetNNeurons();


 	  //initialize num_nodes[]
 	  num_nodes = new unsigned long [nLayers];
 	  num_nodes[0] = num_nn_inputs+1; //+1 if using bias nodes
 	  num_nodes[nLayers-1] = 1; //+1 if using bias nodes
 	  for(unsigned short iLayer = 1; iLayer < nLayers-1; iLayer++) num_nodes[iLayer] = nNeurons;

 	  //initialize num_inputs[]
 	  num_inputs = new unsigned long [nLayers];
 	  //num_inputs[0] = 0; //Input layer has no inputs
 	  num_inputs[0] = num_nn_inputs+1; //+1 if using bias nodes
 	  //num_inputs[1] = num_nn_inputs+1;  //JRH 09232018 - Removing input layer weights from costly computation
	  for (unsigned short iLayer = 1; iLayer<nLayers;iLayer++) num_inputs[iLayer] = nNeurons;  //JRH 09232018 - Removing input layer weights from costly computation


 	  feat_send = new su2double [num_nn_inputs];
 	  feat_recv = new su2double [num_nn_inputs];

 	  //Initialize variables to store intermediate network values
 	  ai = new su2double * [nLayers]; //"Activations"
 	  for (unsigned short iLayer = 0; iLayer < nLayers; iLayer++) ai[iLayer] = new su2double[nNeurons];

 	  inputs = new su2double * [nLayers]; //(o)Array for storing outputs of each layer (hidden 1 outputs at index 0) Does not include output layer
 	  //inputs[0] = new su2double[num_nn_inputs];
 	  //inputs[1] = new su2double[num_nn_inputs];
 	  for (unsigned short iLayer = 0; iLayer < nLayers; iLayer++) inputs[iLayer] = new su2double[nNeurons];

 	  //Initilize variables to hold delta (gradient intermediates) function
 	  deltas = new su2double * [nLayers];
 	  deltas[nLayers-1] = new su2double[1];
 	  deltas[0] = new su2double[1];
 	  for (unsigned short iLayer = nLayers-1; iLayer > 0; iLayer--) deltas[iLayer] = new su2double[nNeurons];

 	  //Initialize weights and similarly dimensioned array to store derivative of error w.r.t each weight
 	  weights = new su2double ** [nLayers];
 	  lweights = new su2double ** [nLayers];
 //	  weight_send = new su2double ** [nLayers-2];
 //	  weight_recv = new su2double ** [nLayers-2];
 	  Ep = new su2double ** [nLayers];
 	  su2double rand_temp = 0.0;
 	  //initialize hidden layers
 	  su2double scale_range = 1.0/sqrt(num_nn_inputs);
 	  num_weights = 0;
 	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
 		  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
 			  weights[iLayer] = new su2double * [num_inputs[iLayer]];
 			  lweights[iLayer] = new su2double * [num_inputs[iLayer]];
 			  Ep[iLayer] = new su2double * [num_inputs[iLayer]];
 //			  if (iLayer>0 && iLayer<nLayers-1) {
 //				  weight_send[iLayer-1] = new su2double * [num_inputs[iLayer]];
 //				  weight_recv[iLayer-1] = new su2double * [num_inputs[iLayer]];
 //			  }
 			  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
 				  weights[iLayer][iInput] = new su2double [num_nodes[iLayer]];
 				  lweights[iLayer][iInput] = new su2double [num_nodes[iLayer]];
 				  Ep[iLayer][iInput] = new su2double [num_nodes[iLayer]];
 //				  if (iLayer>0 && iLayer<nLayers-1) {
 //					  weight_send[iLayer-1][iInput] = new su2double [num_inputs[iLayer]];
 //					  weight_recv[iLayer-1][iInput] = new su2double [num_inputs[iLayer]];
 //				  }
 				  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
 					  //rand_temp = su2double(rand()*scale_range*2.0/RAND_MAX-scale_range);
 					  weights[iLayer][iInput][iNode] = 1.0e-16;
 					  lweights[iLayer][iInput][iNode] = direct_solver->GetWeight(iLayer,iInput,iNode);
 					  Ep[iLayer][iInput][iNode] = 0.0;
 					  num_weights++;  //JRH 09232018 - Removing input layer weights from costly computation
 					  //if (jrh_debug) cout << "iLayer " << iLayer << " iInput " << iInput << " iNode " << " weight init to -> " << rand_temp << endl;
 				  }
 			  }
 		  }
 	  }
 	  if (jrh_debug) cout << "num_weights set to " << num_weights << " in solver_adjoint_discrete.cpp " << endl;
 	  //lloss = direct_solver->GetTotal_Loss();

 	  weight_send = new su2double [num_weights];
 	  weight_recv = new su2double [num_weights];

 	  //Allocate regardless - makes delete[] easier...
       restart_f1 = new su2double[nPointDomain];
       restart_f2 = new su2double[nPointDomain];
       restart_f3 = new su2double[nPointDomain];
       restart_f4 = new su2double[nPointDomain];
       sse = 0.0;
       if (jrh_debug) cout << "JRH Debug: In solver_adjoint_discrete.cpp Done setting neural network vars" << endl;
  }

}

CDiscAdjSolver::~CDiscAdjSolver(void) { 

  unsigned short iMarker;
  unsigned short iDV; //JRH 05172017

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

  if (Sens_Geo   != NULL) delete [] Sens_Geo;
  if (Sens_Mach  != NULL) delete [] Sens_Mach;
  if (Sens_AoA   != NULL) delete [] Sens_AoA;
  if (Sens_Press != NULL) delete [] Sens_Press;
  if (Sens_Temp  != NULL) delete [] Sens_Temp;

  if (Local2Global != NULL) delete [] Local2Global;

  //Delete variables needed for fiml if they are defined - JRH 05172017
  if (KindDirect_Solver == RUNTIME_TURB_SYS) {
	  if (Total_Sens_Beta_Fiml != NULL) delete [] Total_Sens_Beta_Fiml;
	  if (Val_Beta_Fiml != NULL) delete [] Val_Beta_Fiml;
	  if (Local_Sens_Beta_Fiml != NULL) delete [] Local_Sens_Beta_Fiml;



	  if (train_NN) {
		  if (jrh_debug) cout << "JRH Debug: In solver_adjoint_discrete.cpp Beginning to delete NN vars" << endl;
		  if (weights != NULL) { //JRH 04172018 - Delete neural network weights in deconstructor
			  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
				for (unsigned short iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
					delete [] weights[iLayer][iInput];
					delete [] lweights[iLayer][iInput];
					delete [] Ep[iLayer][iInput];
		//			if (iLayer>0 && iLayer<nLayers-1) {
		//				  delete [] weight_send[iLayer-1][iInput];
		//				  delete [] weight_recv[iLayer-1][iInput];
		//			}
				}
				delete [] weights[iLayer];
				delete [] lweights[iLayer];
				delete [] Ep[iLayer];
		//		if (iLayer>0 && iLayer<nLayers-1) {
		//			  delete [] weight_send[iLayer-1];
		//			  delete [] weight_recv[iLayer-1];
		//		}
			  }
			  delete [] weights;
			  delete [] lweights;
			  delete [] Ep;
			  delete [] weight_send;
			  delete [] weight_recv;
			  delete [] feat_send;
			  delete [] feat_recv;
		  }
		  for (unsigned short iLayer = 0; iLayer < nLayers; iLayer++) {
			  delete [] inputs[iLayer];
			  delete [] deltas[iLayer];
			  delete [] ai[iLayer];
		  }
		  delete [] inputs;
		  delete [] deltas;
		  delete [] ai;
		  if (num_inputs != NULL) delete [] num_inputs;
		  if (num_nodes != NULL) delete [] num_nodes;
		  if (restart_f1 != NULL) {
			  delete [] restart_f1;
			  delete [] restart_f2;
			  delete [] restart_f3;
			  delete [] restart_f4;
		  }
	  } //<- train_NN
  }
}

void CDiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config, unsigned short kind_recording) {


  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->SetSolution(node[iPoint]->GetSolution_Direct());

    //Need to reset beta at each node if NN training 05022018
    //if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) direct_solver->node[iPoint]->SetBetaFiml(node[iPoint]->GetBetaFiml());
  }

  //JRH 04312018
  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolver::SetRecording about to reset weights with SetWeight()"<<endl;
	  unsigned long iDV = 0;
 	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
			  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
				  direct_solver->SetWeight(iLayer,iInput,iNode,lweights[iLayer][iInput][iNode]);
				  iDV++;
			  }
		  }
	  }
 	  //if (config->GetLambdaLossFiml() > 1.0e-16) direct_solver->SetTotal_SSE(lloss);
  }

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the meanflow iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjSolver::RegisterNNSolution(CGeometry *geometry, CConfig *config) {
//	  JRH 04302018
	  bool input = true;
	  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//	  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
		  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterNNSolution() about to call direct_solver->RegisterWeights();" << endl;
		  direct_solver->RegisterWeights(input);
	//	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterSolution() about to call direct_solver->ForwardPropagate();" << endl;
	//	  direct_solver->ForwardPropagate();
	  }
}

void CDiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, nPoint = geometry->GetnPoint();

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND,
  input = true;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);

//    if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) direct_solver->node[iPoint]->RegisterBeta(input);//JRH 05022018
  }
  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n1();
    }
  }
//  JRH 04302018
  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterSolution() about to call direct_solver->RegisterWeights();" << endl;
	  direct_solver->RegisterWeights(input);
//	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterSolution() about to call direct_solver->ForwardPropagate();" << endl;
//	  direct_solver->ForwardPropagate();
  }
}

void CDiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {
	bool fiml; //JRH - 04192017
	unsigned long nDV = config->GetnDV(); //JRH - 04192017
	unsigned long iDV = 0; //JRH - 04192017
	unsigned long nPoint = geometry->GetnPoint(); //JRH - 04192017
	unsigned long GlobalIndex; // JRH - 05082017
	unsigned short KindTrainNN = config->GetKindTrainNN();//JRH 06282018
	//JRH 05082017 - Needed for debugging comments only
	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	/*--- Register farfield values as input ---*/
  if((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)) {



    su2double Velocity_Ref = config->GetVelocity_Ref();
    Alpha                  = config->GetAoA()*PI_NUMBER/180.0;
    Beta                   = config->GetAoS()*PI_NUMBER/180.0;
    Mach                   = config->GetMach();
    Pressure               = config->GetPressure_FreeStreamND();
    Temperature            = config->GetTemperature_FreeStreamND();

    su2double SoundSpeed = 0.0;
    
    if (nDim == 2) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*Mach); }
    if (nDim == 3) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*cos(Beta)*Mach); }

    if (!reset) {
      AD::RegisterInput(Mach);
      AD::RegisterInput(Alpha);
      AD::RegisterInput(Temperature);
      AD::RegisterInput(Pressure);
    }

    /*--- Recompute the free stream velocity ---*/

    if (nDim == 2) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[2] = sin(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
    }

    config->SetTemperature_FreeStreamND(Temperature);
    direct_solver->SetTemperature_Inf(Temperature);
    config->SetPressure_FreeStreamND(Pressure);
    direct_solver->SetPressure_Inf(Pressure);




    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/


  }
  //JRH Added below blocks to extract beta_fiml sensitivities. 04192017
  //JRH - All code below (in this routine) is by JRH, was blank before he touched it.
  //Find any fiml variables (if any)
  //if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_TURB_SYS)) { //JRH 05012017 - I think we should only do this for the turbulent solver...
  if (KindDirect_Solver == RUNTIME_TURB_SYS) { //JRH - Trying to run compressible, also changed in ExtractAdjoint... 09212017
  	  for (unsigned long count = 0; count < nDV ; count++) {
		  if (config->GetDesign_Variable(count) == FIML) {
			  if (!fiml) {
				  iDV = count;
				 // cout << "In solver_adjoint_discrete.cpp - Setting fiml == " << fiml << " and iDV = " << iDV << endl;
			  }
			  fiml = true;
		  }
	  }

	  if (fiml) {
		  //Val_Beta_Fiml = new su2double[nDV-iDV+1];
		  //cout << "JRH Debugging: Process " << rank << " in solver_adjoint_discrete.cpp - Starting to set Val_Beta_Fiml from config class values" << endl;
		  //cout << "nDV = " << nDV << endl;

		  for (iDV = 0; iDV < nDV; iDV++) {
			  if (KindTrainNN != WEIGHTS) Val_Beta_Fiml[iDV] = config->GetDV_Value(iDV, 0)+1.0;
			  else Val_Beta_Fiml[iDV] = config->GetDV_Value(iDV, 0); //JRH 06282018
		  }
		  if (jrh_debug) cout << endl;
	  }
	  //Find total sensitivity to beta at each point (each DV) 04192017
	  if (fiml && !reset) {
	  //if (fiml) {
		  //cout << "JRH Debugging: Process "<< rank << " starting to register Val_Beta_Fiml values as input for AD" << endl;
		  //Val_Beta_Fiml = new su2double[nDV-iDV+1];
		  su2double beta_fiml_val;
		  for (iDV = 0; iDV < nDV; iDV++) {
			  AD::RegisterInput(Val_Beta_Fiml[iDV]);
			  //beta_fiml_val = config->GetDV_Value(count,0)+1.0; //This worked in serial
			  //AD::RegisterInput(beta_fiml_val); //This worked in serial
			  //Val_Beta_Fiml[count] = beta_fiml_val; //This worked in serial
		  }
		  //if (config->GetTrainNN()) direct_solver->ForwardPropagate();
		  //cout << "JRH Debugging: Done registering Val_Beta_Fiml as input for AD" << endl;
	  }
	  if (fiml && KindTrainNN != WEIGHTS) {
		  //cout << "JRH Debugging: In solver_adjoint_discrete.cpp - Starting to set beta via SetDV_Value()" << endl;
		  //cout << "Length of Val_Beta_Fiml = " << nDV-iDV << " and nPoint = " << nPoint << endl;
//#ifdef HAVE_MPI
				  for (iDV = 0; iDV < nDV ; iDV++){
					  config->SetDV_Value(iDV, 0, Val_Beta_Fiml[iDV]-1.0);
				  }
				  unsigned long nDV_Local = 0;
				  //cout << "JRH Debugging: proc " << rank << " starting with beta_fiml index " << nDV_Local+Fiml_Skip_Index << endl;
				  for (unsigned long iPoint = 0; iPoint < nPoint ; iPoint++) {
					  if (geometry->node[iPoint]->GetDomain()) {
						  //cout << "JRH Debugging: Setting Val_Beta_Fiml[" << nDV_Local+Fiml_Skip_Index << "] on process " << rank << endl;
						  //config->SetDV_Value(nDV_Local+Fiml_Skip_Index, 0, Val_Beta_Fiml[nDV_Local+Fiml_Skip_Index]-1.0);

					      //THIS WORKED SINGLE NODE 06122018
//						  if (!config->GetTrainNN()) direct_solver->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[nDV_Local+Fiml_Skip_Index]);
//						  direct_solver->node[iPoint]->SetBetaFimlTrain(Val_Beta_Fiml[nDV_Local+Fiml_Skip_Index]);
//						  geometry->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[nDV_Local+Fiml_Skip_Index]);

						  if (!config->GetTrainNN()) direct_solver->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[Local2Global[iPoint]]);
						  direct_solver->node[iPoint]->SetBetaFimlTrain(Val_Beta_Fiml[Local2Global[iPoint]]);
						  geometry->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[Local2Global[iPoint]]);

						  nDV_Local++;
					  }
				  }
				  //cout << "JRH Debugging: proc " << rank << " ending with beta_fiml index " << nDV_Local+Fiml_Skip_Index-1 << endl;
				/* This didn't work...not sure why, trying to rewrite using different function calls, can't seem to go back and forth between
				 * global and local as proven by the comments - JRH 05102017
				 */
			      /*long iPoint_Local;
			      cout << "JRH Debugging: In solver_adjoint_discrete RegisterVariables() nDV = " << nDV << " Global_nPoint_Domain = " << geometry->GetGlobal_nPointDomain() << endl;
				  for (unsigned long iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPoint(); iPoint_Global++ ) {
					  iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);
					  //if (iPoint_Local == -1 && rank == MASTER_NODE) cout << "JRH Debugging: iPoint_Global " << iPoint_Global << " not in master node partition" << endl;
					  if ( iPoint_Local != -1 && geometry->node[iPoint_Local]->GetDomain()) {
						 //if (rank == MASTER_NODE) cout << "JRH Debugging: Setting iPoint_Local " << iPoint_Local << " iPoint_Global " << iPoint_Global << " geometry->node[iPoint_Local]->GetGlobalIndex() = " << geometry->node[iPoint_Local]->GetGlobalIndex() << " Process " << rank << endl;
						 direct_solver->node[iPoint_Local]->SetBetaFiml(Val_Beta_Fiml[iPoint_Global]);
					  }
				  }*/

//#else
//		  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
//			  //cout << "JRH Debugging: Setting DV " << count << " to value Val_Beta_Fiml[" << count << "]= " << Val_Beta_Fiml[count] << endl;
//
//			  //THIS WORKED SINGLE NODE 06122018
////			  config->SetDV_Value(iPoint, 0, Val_Beta_Fiml[iPoint]-1.0);
////			  //cout << "JRH Debugging: Done setting in config class, now setting in solver class..." << endl;
////			  if (!config->GetTrainNN()) direct_solver->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[iPoint]);
////			  direct_solver->node[iPoint]->SetBetaFimlTrain(Val_Beta_Fiml[iPoint]);
//
//			  config->SetDV_Value(iPoint, 0, Val_Beta_Fiml[Local2Global[iPoint]]-1.0);
//			  //cout << "JRH Debugging: Done setting in config class, now setting in solver class..." << endl;
//			  if (!config->GetTrainNN()) direct_solver->node[iPoint]->SetBetaFiml(Val_Beta_Fiml[Local2Global[iPoint]]);
//			  direct_solver->node[iPoint]->SetBetaFimlTrain(Val_Beta_Fiml[Local2Global[iPoint]]);
//		  }
//		  //cout << "JRH Debugging: In solver_adjoint_discrete.cpp - Done setting beta via SetBetaFiml()" << endl;
//#endif
	  }
	  else if (KindTrainNN == WEIGHTS) { //JRH 06282018
		  unsigned long iDV = 0;
	 	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
			  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
				  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
					  direct_solver->SetWeight(iLayer,iInput,iNode,Val_Beta_Fiml[iDV]);
					  iDV++;
				  }
			  }
		  }
	  }
  }

}

void CDiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, nPoint = geometry->GetnPoint();

  /*--- Register variables as output of the solver iteration ---*/

  bool input = false;

  /*--- Register output variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);
//    if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) direct_solver->node[iPoint]->RegisterBeta(input); //JRH 05022015
  }

  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterOutput() about to call direct_solver->RegisterWeights();" << endl;
	  direct_solver->RegisterWeights(input);
	  //if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::RegisterSolution() about to call direct_solver->ForwardPropagate();" << endl;
	  //direct_solver->ForwardPropagate();
  }
}

void CDiscAdjSolver::RegisterObj_Func(CConfig *config) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Here we can add new (scalar) objective functions ---*/
  if (config->GetnObj()==1) {
    switch (config->GetKind_ObjFunc()) {
    case DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CD();
      if (config->GetFixed_CL_Mode()) ObjFunc_Value -= config->GetdCD_dCL() * direct_solver->GetTotal_CL();
      if (config->GetFixed_CM_Mode()) ObjFunc_Value -= config->GetdCD_dCM() * direct_solver->GetTotal_CMy();
      break;
    case LIFT_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CL();
      break;
    case AERO_DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_AeroCD();
      break;
    case RADIAL_DISTORTION:
      ObjFunc_Value = direct_solver->GetTotal_RadialDistortion();
      break;
    case CIRCUMFERENTIAL_DISTORTION:
      ObjFunc_Value = direct_solver->GetTotal_CircumferentialDistortion();
      break;
    case SIDEFORCE_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CSF();
      break;
    case EFFICIENCY:
      ObjFunc_Value = direct_solver->GetTotal_CEff();
      break;
    case MOMENT_X_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMx();
      break;
    case MOMENT_Y_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMy();
      break;
    case MOMENT_Z_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMz();
      break;
    case EQUIVALENT_AREA:
      ObjFunc_Value = direct_solver->GetTotal_CEquivArea();
      break;
    case AVG_TOTAL_PRESSURE:
      ObjFunc_Value = direct_solver->GetOneD_TotalPress();
      break;
    case AVG_OUTLET_PRESSURE:
      ObjFunc_Value = direct_solver->GetOneD_FluxAvgPress();
      break;
    case MASS_FLOW_RATE:
      ObjFunc_Value = direct_solver->GetOneD_MassFlowRate();
      break;
    case INVERSE_DESIGN_LIFT:
      ObjFunc_Value = direct_solver->GetTotal_ClDiff();
      break;
    case INVERSE_DESIGN_LIFT_FIML:
      ObjFunc_Value = direct_solver->GetTotal_ClDiff_FIML();
      break;
    case INVERSE_DESIGN_DRAG:
      ObjFunc_Value = direct_solver->GetTotal_CdDiff();
      break;
    case INVERSE_DESIGN_DRAG_FIML:
      ObjFunc_Value = direct_solver->GetTotal_CdDiff_FIML();
      break;
    case INVERSE_DESIGN_PRESSURE_FIML: //TESTING 02242019 - WASN'T HERE PRIOR
    	ObjFunc_Value = direct_solver->GetTotal_CpDiff_FIML();
      break;
    case INVERSE_DESIGN_PRESSURE: 	//TESTING 02242019 - WASN'T HERE PRIOR
    	ObjFunc_Value = direct_solver->GetTotal_CpDiff();
      break;
    //case INVERSE_DESIGN_PRESSURE: //NOT in original SU2 - testing whether nObj set == 1 if INVERSE_DESIGN_PRESSURE and objective not being computed correctly JRH 09272017
    //  ObjFunc_Value = direct_solver->GetTotal_ComboObj();
    //  break;
    }

    /*--- Template for new objective functions where TemplateObjFunction()
     *  is the routine that returns the obj. function value. The computation
     * must be done while the tape is active, i.e. between AD::StartRecording() and
     * AD::StopRecording() in DiscAdjMeanFlowIteration::Iterate(). The best place is somewhere
     * inside MeanFlowIteration::Iterate().
     *
     * case TEMPLATE_OBJECTIVE:
     *    ObjFunc_Value = TemplateObjFunction();
     *    break;
     * ---*/
  }
  else{
    ObjFunc_Value = direct_solver->GetTotal_ComboObj();
  }
  if (rank == MASTER_NODE) {
	//CEulerSolver::Compute_ComboObj computes INVERSE_DESIGN_PRESSURE case, this should work as is - JRH 09272017
    AD::RegisterOutput(ObjFunc_Value);
  }
}

void CDiscAdjSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) {
  int rank = MASTER_NODE;

  bool time_stepping = config->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  unsigned long ExtIter = config->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping) {
    if (ExtIter < IterAvg_Obj) {
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else {
      seeding = 0.0;
    }
  }

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CDiscAdjSolver::ExtractNNAdjoint(CGeometry *geometry, CConfig *config) {
//	  JRH 04302018 - Get Derivative of all weights in NN and store in weights variable in Adjoint Solver
	  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//	  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
		  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::ExtractNNAdjoint() about to call direct_solver->GetNNGradient();"<<endl;
		  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
			  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
					  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
						  weights[iLayer][iInput][iNode] = direct_solver->GetNNGradient(iLayer,iInput,iNode);
						  if (jrh_debug) cout << weights[iLayer][iInput][iNode] << " ";
					  }
			  }
		  }
		  if (jrh_debug) cout << endl;
//		  if (config->GetLambdaLossFiml() > 1.0e-16) {
//			  sse = direct_solver->GetNNLossGradient();
//			  if (jrh_debug) cout << "JRH Debugging: in CDiscAdjSolution::ExtractNNAdjoint(), gradient of sse = " << sse << endl;
//		  }
	  }
}

void CDiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config) {

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }
  //if (jrh_debug) cout << "JRH Debuggin: In CDiscAdjSolution::ExtractAdjoint_Solution() about to call direct_solver->node->GetAdjointBeta()" << endl;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Set the old solution ---*/

    node[iPoint]->Set_OldSolution();

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

//    if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
//    	node[iPoint]->SetBetaFimlGrad(direct_solver->node[iPoint]->GetAdjointBeta());
//    	if (jrh_debug) cout << node[iPoint]->GetBetaFimlGrad() << " ";
//    }

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }
//  if (jrh_debug) cout << endl;

  //JRH 04302018 - Get Derivative of all weights in NN and store in weights variable in Adjoint Solver
  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
	  if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::ExtractAdjoint_Solution() about to call direct_solver->GetNNGradient();"<<endl;
	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
		  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
				  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
					  weights[iLayer][iInput][iNode] = direct_solver->GetNNGradient(iLayer,iInput,iNode);
					  if (jrh_debug) cout << weights[iLayer][iInput][iNode] << " ";
				  }
		  }
	  }
	  if (jrh_debug) cout << endl;
//	  if (config->GetLambdaLossFiml() > 1.0e-16) {
//		  sse = direct_solver->GetNNLossGradient();
//		  if (jrh_debug) cout << "JRH Debugging: in CDiscAdjSolution::ExtractNNAdjoint(), gradient of sse = " << sse << endl;
//	  }
  }


  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);

      /*--- Store the adjoint solution at time n ---*/

      node[iPoint]->Set_Solution_time_n(Solution);
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n-1 ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);

      /*--- Store the adjoint solution at time n-1 ---*/

      node[iPoint]->Set_Solution_time_n1(Solution);
    }
  }

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->GetSolution_Old(iVar);

          AddRes_RMS(iVar,residual*residual);
          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_RMS(geometry, config);
}

void CDiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) {
  su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;
  //su2double *Local_Sens_Beta_Fiml; //JRH 04192017
  bool fiml = false;
  bool l2_reg = config->GetL2Reg();
  unsigned long nDV = config->GetnDV();
  unsigned long iDV = 0;
  unsigned long Global_Index;
  bool fiml_of = false;
  unsigned short kind_of = config->GetKind_ObjFunc();
  su2double lambda_fiml = config->GetLambdaFiml();
  if (kind_of == INVERSE_DESIGN_LIFT_FIML || kind_of == INVERSE_DESIGN_DRAG_FIML || kind_of == INVERSE_DESIGN_PRESSURE_FIML) fiml_of = true;

  //JRH - For MPI Debugging only - 05082017
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)) {
    Local_Sens_Mach  = SU2_TYPE::GetDerivative(Mach);
    Local_Sens_AoA   = SU2_TYPE::GetDerivative(Alpha);
    Local_Sens_Temp  = SU2_TYPE::GetDerivative(Temperature);
    Local_Sens_Press = SU2_TYPE::GetDerivative(Pressure);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_Mach,  &Total_Sens_Mach,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_AoA,   &Total_Sens_AoA,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temp,  &Total_Sens_Temp,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    Total_Sens_Mach  = Local_Sens_Mach;
    Total_Sens_AoA   = Local_Sens_AoA;
    Total_Sens_Temp  = Local_Sens_Temp;
    Total_Sens_Press = Local_Sens_Press;
#endif


  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/


  }

  //JRH Added below blocks to extract beta_fiml sensitivities. 04192017
  //Find any fiml variables (if any)
  //cout << "JRH Debugging: Process " << rank << " in CDiscAdjSolver::ExtractAdjointVariables()" << endl;

  //Find total sensitivity to beta at each point (each DV) 04192017
  //cout << "JRH Debugging: Process " << rank << " starting to extract FIML sensitivities nPoint = " << nPoint << " and nDV = " << nDV << endl;
  //if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_TURB_SYS)) { //JRH 05012017 - I think we should only do this for the turbulent solver... RUNTIME_TURB_SYS never called?
  if (KindDirect_Solver == RUNTIME_TURB_SYS) {
  // for (unsigned long count = 0; count < nDV ; count++) {
	 //	  if (config->GetDesign_Variable(count) == FIML) {
	 //		  if (!fiml) iDV = count;
	 //		  fiml = true;
	 //	  }
	  //}
	  if (config->GetKind_Turb_Model()==SA_FIML) fiml = true;
	  if (fiml) {
		su2double beta_fiml_val;
		//Local_Sens_Beta_Fiml = new su2double[nDV];
		//su2double weight = config->GetWeight_ObjFunc(0);
		su2double weight = 1.0;
#ifdef HAVE_MPI
		  //cout << "rank = " << rank << " nDV = " << nDV << " nPoint = " << nPoint << " nPointDomain = " << nPointDomain << endl;
		  for (iDV = 0; iDV < nDV; iDV++) {
			  //beta_fiml_val = config->GetDV_Value(iDV,0)+1.0;
			  //Local_Sens_Beta_Fiml[iDV] = SU2_TYPE::GetDerivative(beta_fiml_val)*weight;
			  Local_Sens_Beta_Fiml[iDV] = SU2_TYPE::GetDerivative(Val_Beta_Fiml[iDV]);
			  //cout << "JRH Debugging: Starting to do SU2_MPI::Allreduce() on element count = " << count << endl;
			  //SU2_MPI::Reduce(&Local_Sens_Beta_Fiml[iDV],  &Total_Sens_Beta_Fiml[iDV],  1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

			  //SU2_MPI::Allreduce(&Local_Sens_Beta_Fiml[iDV],  &Total_Sens_Beta_Fiml[iDV],  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			  //cout << "JRH Debugging: In CDiscAdjSolver::ExtractAdjointVariables() - Total_Sens_Beta_Fiml[" << count << "]= " << Total_Sens_Beta_Fiml[count] << endl;
			  if (jrh_debug) cout << " Total_Sens_Beta_Fiml[" << iDV << "]= " << Total_Sens_Beta_Fiml[iDV] << " ";
		  }

		  SU2_MPI::Allreduce(Local_Sens_Beta_Fiml,  Total_Sens_Beta_Fiml,  nDV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		  //Add in Regularization penalty
		  if (config->GetKindTrainNN() != WEIGHTS && fiml_of) {
			  for (iDV = 0; iDV < nDV; iDV++) {
				  //THIS IS WRONG!! 02232019
				  //Total_Sens_Beta_Fiml[iDV] += lambda_fiml*Val_Beta_Fiml[iDV]*(Val_Beta_Fiml[iDV]-1.0);
				  //THIS IS RIGHT!!
				  Total_Sens_Beta_Fiml[iDV] += lambda_fiml*(Val_Beta_Fiml[iDV]-1.0);
			  }
		  }

		  unsigned long nDV_Local = 0;
		  //if (rank == MASTER_NODE) { //JRH 12192017
			  for (unsigned long iPoint = 0; iPoint < nPoint ; iPoint++) {
				  if (geometry->node[iPoint]->GetDomain()) {
					  //cout << "proc " << rank << " Setting node[" << iPoint << "] with beta = " << Val_Beta_Fiml[nDV_Local+Fiml_Skip_Index] << " Total_Sens_Beta_Fiml[" << nDV_Local+Fiml_Skip_Index << "] = " << Total_Sens_Beta_Fiml[nDV_Local+Fiml_Skip_Index] << endl;
					  //JRH 03232018 - ADDED VOLUME NORMALIZATION - JRH 03232018
					  //Total_Sens_Beta_Fiml[nDV_Local+Fiml_Skip_Index] = Total_Sens_Beta_Fiml[nDV_Local+Fiml_Skip_Index]/geometry->node[iPoint]->GetVolume();

					  //THIS WORKED SINGLE NODE 06122018
//					  direct_solver->node[iPoint]->SetBetaFimlGrad(Total_Sens_Beta_Fiml[nDV_Local+Fiml_Skip_Index]);

					  if (config->GetKindTrainNN() != WEIGHTS) direct_solver->node[iPoint]->SetBetaFimlGrad(Total_Sens_Beta_Fiml[Local2Global[iPoint]]);

					  nDV_Local++;
				  }
				  else {
					  //direct_solver->node[iPoint]->SetBetaFimlGrad(-10.0);
					  //direct_solver->node[iPoint]->SetBetaFiml(-10.0);
				  }
			  }
		 // }
		  if (jrh_debug) cout << endl;
		  	/*unsigned long iPoint_Local;
			for (unsigned long iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPoint(); iPoint_Global++ ) {
			  iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);
			  if ( iPoint_Local != -1 && geometry->node[iPoint_Local]->GetDomain()) {
				 direct_solver->node[iPoint_Local]->SetBetaFimlGrad(Total_Sens_Beta_Fiml[iPoint_Global]);
			  }
			}*/
#else
		  for (unsigned long count = 0; count < nDV; count++) {
			  //THIS WORKED SINGLE NODE 06122018
//			  Local_Sens_Beta_Fiml[count] = SU2_TYPE::GetDerivative(Val_Beta_Fiml[count]);//Should this be iDV vs count?? Changed to count vs count JRH 03232018
//			  //Total_Sens_Beta_Fiml[count] = Local_Sens_Beta_Fiml[count]/geometry->node[count]->GetVolume(); //JRH 03232018 - ADDED VOLUME NORMALIZATION - JRH 03232018
//			  Total_Sens_Beta_Fiml[count] = Local_Sens_Beta_Fiml[count];
//			  //cout << "JRH Debugging: In CDiscAdjSolver::ExtractAdjointVariables() - Beta = " << Val_Beta_Fiml[count] << " Total_Sens_Beta_Fiml[" << count << "]= " << Total_Sens_Beta_Fiml[count] << endl;
//			  //cout << " Total_Sens_Beta_Fiml[" << count << "]= " << Total_Sens_Beta_Fiml[count];
//			  direct_solver->node[count]->SetBetaFimlGrad(Total_Sens_Beta_Fiml[count]);

			  Local_Sens_Beta_Fiml[count] = SU2_TYPE::GetDerivative(Val_Beta_Fiml[count]);//Should this be iDV vs count?? Changed to count vs count JRH 03232018
			  Total_Sens_Beta_Fiml[count] = Local_Sens_Beta_Fiml[count];
			  //Add in Regularization penalty
			  //THIS IS WRONG 02231019
			  //Total_Sens_Beta_Fiml[count] += lambda_fiml*Val_Beta_Fiml[count]*(Val_Beta_Fiml[count]-1.0);
			  //THIS IS RIGHT 02232019
			  Total_Sens_Beta_Fiml[count] += lambda_fiml*(Val_Beta_Fiml[count]-1.0);
			  //cout << "JRH Debugging: In CDiscAdjSolver::ExtractAdjointVariables() - Beta = " << Val_Beta_Fiml[count] << " Total_Sens_Beta_Fiml[" << count << "]= " << Total_Sens_Beta_Fiml[count] << endl;
			  //cout << " Total_Sens_Beta_Fiml[" << count << "]= " << Total_Sens_Beta_Fiml[count];
			  if (config->GetKindTrainNN() != WEIGHTS) direct_solver->node[count]->SetBetaFimlGrad(Total_Sens_Beta_Fiml[Local2Global[count]]);

		  }
		  //cout << endl;
		  //cout << "JRH Debugging: Process " << rank << " done in CDiscAdjSolver::ExtractAdjointVariables()" << endl;
#endif

			  // Tested on 04192017 and was getting non-zero Total_Sens_Beta_Fiml[iDV] at every iteration.
		  //delete [] Local_Sens_Beta_Fiml; // JRH 05192017 - Delete to prevent memory leak
		  //cout<< "fiml_of = " << fiml_of << endl;
		  //cout<< "kind_of = " << kind_of << endl;
		  if (config->GetKindTrainNN() == WEIGHTS && fiml_of) {

			  //bool jrh_debug = true;
			  unsigned long iDV = 0;
		 	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
				  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
					  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
						  //direct_solver->SetWeight(iLayer,iInput,iNode,lweights[iLayer][iInput][iNode]);
						  if (l2_reg) {
							Total_Sens_Beta_Fiml[iDV] += lambda_fiml*weights[iLayer][iInput][iNode];
						  }
						  else {
							  Total_Sens_Beta_Fiml[iDV] += lambda_fiml*direct_solver->GetEp(iLayer,iInput,iNode);
						  }
						  if (jrh_debug && rank == MASTER_NODE) cout << "JRH Debug: Total_Sens_Beta_Fiml[iDV]=" << Total_Sens_Beta_Fiml[iDV] << " + lambda_fiml*direct_solver->GetEp(iLayer,iInput,iNode)=" << lambda_fiml*direct_solver->GetEp(iLayer,iInput,iNode) << endl;
						  iDV++;
					  }
				  }
			  }
		  }
  }

  }

}

void CDiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config) {

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
      config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  unsigned short iVar;
  unsigned long iPoint;
  //cout << "JRH Debugging: In CDiscAdjSolver::SetAdjoint_Output(), Starting to call GetSolution() and SetAdjointSolution() for every node" << endl;
//  if (jrh_debug) {cout << "In SetAdjoint_Output() Setting beta_fiml gradients to: ";}
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution[iVar] = node[iPoint]->GetSolution(iVar);
    }
    if (dual_time) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += node[iPoint]->GetDual_Time_Derivative(iVar);
      }
    }
    direct_solver->node[iPoint]->SetAdjointSolution(Solution);
  }
//    if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
//    	node[iPoint]->SetBetaFimlGrad(direct_solver->node[iPoint]->GetAdjointBeta());
//    	direct_solver->node[iPoint]->SetAdjointBeta(node[iPoint]->GetBetaFimlGrad());
//    	if (jrh_debug) cout << iPoint << " " << direct_solver->node[iPoint]->GetAdjointBeta() << " | " ;
//    }
//  }
//  if (jrh_debug) cout << endl;
  //JRH 04302018 - Get Derivative of all weights in NN and store in weights variable in Adjoint Solver
//  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
//  	if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::SetAdjointOutput() about to call direct_solver->GetNNGradient();"<<endl;
//	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
//		  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
//				  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
//					  weights[iLayer][iInput][iNode] = direct_solver->GetNNGradient(iLayer,iInput,iNode);
//					  if (jrh_debug) cout << weights[iLayer][iInput][iNode] << " ";
//				  }
//		  }
//	  }
//	  if (jrh_debug) cout << endl;
    if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS && config->GetKindTrainNN() != WEIGHTS) {
//  	  if (config->GetTrainNN() && KindDirect_Solver == RUNTIME_TURB_SYS) {
		if (jrh_debug) cout << "JRH Debugging: In CDiscAdjSolution::SetAdjoint_Output() about to call direct_solver->SetAdjointNNSolution();"<<endl;
		  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
			  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
					  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
						  if (jrh_debug) cout << weights[iLayer][iInput][iNode] << " ";
						  direct_solver->SetAdjointNNSolution(iLayer,iInput,iNode,weights[iLayer][iInput][iNode]);
					  }
			  }
		  }
		  if (jrh_debug) cout << endl;
//		  if (config->GetLambdaLossFiml() > 1.0e-16) {
//			  direct_solver->SetAdjointNNLoss(sse);
//			  if (jrh_debug) cout << "JRH Debugging: in CDiscAdjSolution::SetAdjoint_Output(), gradient of sse = " << sse << endl;
//		  }
	  }

  //cout << "JRH Debugging:  Done in CDiscAdjSolver::SetAdjoint_Output()" << endl;
}

void CDiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;
  //bool fiml = (config->GetKind_Turb_Model()==SA_FIML); //JRH 04282017

  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++) {

      Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);

      /*--- Set the index manually to zero. ---*/

     AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetLimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
          Sensitivity = 0.0;
      }
      //if (!fiml) { //JRH - 04282017 - Didn't work, commenting out in case I want to try something similar later.
		  if (!time_stepping) {
			node[iPoint]->SetSensitivity(iDim, Sensitivity);
		  } else {
			node[iPoint]->SetSensitivity(iDim, node[iPoint]->GetSensitivity(iDim) + Sensitivity);
		  }
      //} else { //JRH 04282017 - Testing theory that geometry sensitivities are wiping out adjoints of beta_fiml;
      //	  node[iPoint]->SetSensitivity(iDim, 1.0);
      //}
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker,iDim;
  unsigned long iVertex, iPoint;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;
  su2double Total_Sens_Geo_local = 0.0;
  Total_Sens_Geo = 0.0;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Prod = 0.0;
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          /*--- retrieve the gradient calculated with AD -- */
          SensDim = node[iPoint]->GetSensitivity(iDim);

          /*--- calculate scalar product for projection onto the normal vector ---*/
          Prod += Normal[iDim]*SensDim;

          Area += Normal[iDim]*Normal[iDim];
        }

        Area = sqrt(Area);

        /*--- projection of the gradient
         *     calculated with AD onto the normal
         *     vector of the surface ---*/
        Sens = Prod/Area;

        /*--- Compute sensitivity for each surface point ---*/
        CSensitivity[iMarker][iVertex] = -Sens;
        if (geometry->node[iPoint]->GetDomain()) {
          Sens_Geo[iMarker] += Sens*Sens;
        }
      }
      Total_Sens_Geo_local += sqrt(Sens_Geo[iMarker]);
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Total_Sens_Geo_local,&Total_Sens_Geo,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  Total_Sens_Geo = Total_Sens_Geo_local;
#endif
}

void CDiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  bool dual_time_1st = (config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  su2double *solution_n, *solution_n1;
  unsigned long iPoint;
  unsigned short iVar;
  if (dual_time) {
      for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++) {
          solution_n = node[iPoint]->GetSolution_time_n();
          solution_n1 = node[iPoint]->GetSolution_time_n1();

          for (iVar=0; iVar < nVar; iVar++) {
              node[iPoint]->SetDual_Time_Derivative(iVar, solution_n[iVar]+node[iPoint]->GetDual_Time_Derivative_n(iVar));
              node[iPoint]->SetDual_Time_Derivative_n(iVar, solution_n1[iVar]);

            }

        }

    }
}

su2double CDiscAdjSolver::GetBetaFiml(unsigned long iPoint) {
	su2double val_beta_fiml;
	val_beta_fiml = direct_solver->node[iPoint]->GetBetaFiml();
	return val_beta_fiml;
}

su2double CDiscAdjSolver::GetBetaFimlGrad(unsigned long iPoint) {
	su2double val_beta_fiml_grad;
	val_beta_fiml_grad = direct_solver->node[iPoint]->GetBetaFimlGrad();
	//val_beta_fiml_grad = Total_Sens_Beta_Fiml[iPoint];
	return val_beta_fiml_grad;
}

void CDiscAdjSolver::WriteBetaFimlGrad(void) {
	//Write out beta_fiml_grad.dat which contains derivatives in the same sequential order of the design variables
	//JRH 07202017
	ofstream restart_file;
	string filename = "beta_fiml_grad.dat";
    restart_file.open(filename.c_str(), ios::out);
    restart_file.precision(15);
	for (unsigned long iDV = 0; iDV < nDV_total; iDV++) {
		restart_file << Total_Sens_Beta_Fiml[iDV] << "\n";
	}
	restart_file.close();
}
