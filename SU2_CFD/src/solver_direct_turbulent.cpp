/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios, A. Bueno
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

CTurbSolver::CTurbSolver(void) : CSolver() {
  train_NN = false; //JRH 05202018
  filter_shield = false;
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;
  
}

CTurbSolver::CTurbSolver(CConfig *config) : CSolver() {
  train_NN = false; //JRH 05202018
  filter_shield = false;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;

}

CTurbSolver::~CTurbSolver(void) {
  if (jrh_debug) cout << "In CTurbSolver() Destructor" << endl;
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit != NULL) delete [] lowerlimit;
  if (upperlimit != NULL) delete [] upperlimit;
  if (Local2Global!=NULL) delete [] Local2Global;
  if (train_NN) {
	  cout << "JRH Debugging - Beginning to Delete Neural Network Arrays" << endl;
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
		  if (weight_send !=NULL) delete [] weight_send;
		  if (weight_recv !=NULL) delete [] weight_recv;
		  if (feat_send !=NULL) delete [] feat_send;
		  if (feat_recv !=NULL) delete [] feat_recv;
		  if (isHoldout !=NULL) delete [] isHoldout;
	  }
	  for (unsigned short iLayer = 0; iLayer < nLayers; iLayer++) {
		  delete [] inputs[iLayer];
		  delete [] deltas[iLayer];
		  delete [] ai[iLayer];
	  }
	  if (inputs !=NULL) delete [] inputs;
	  if (deltas !=NULL) delete [] deltas;
	  if (deltas !=NULL) delete [] ai;

	  if (num_inputs != NULL) delete [] num_inputs;
	  if (num_nodes != NULL) delete [] num_nodes;
	  if (restart_f1 != NULL) delete [] restart_f1;
	  if (restart_f2 != NULL)  delete [] restart_f2;
	  if (restart_f3 != NULL) delete [] restart_f3;
	  if (restart_f4 != NULL) delete [] restart_f4;

	  if (min_max_send != NULL) delete [] min_max_send;
	  if (min_max_recv != NULL) delete [] min_max_recv;
	  if (edf1 != NULL) delete [] edf1;
	  if (edf2 != NULL) delete [] edf2;
	  if (edf3 != NULL) delete [] edf3;
	  if (edf4 != NULL) delete [] edf4;
	  if (ledf1 != NULL) delete [] ledf1;
	  if (ledf2 != NULL) delete [] ledf2;
	  if (ledf3 != NULL) delete [] ledf3;
	  if (ledf4 != NULL) delete [] ledf4;
	  if (jrh_debug) cout << "Done With ~CTurbSolver() Destructor" << endl;
  }
}

void CTurbSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_muT = NULL, *Buffer_Send_muT = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
  
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      Buffer_Receive_muT = new su2double [nBufferR_Scalar];
      Buffer_Send_muT = new su2double[nBufferS_Scalar];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_muT, nBufferS_Scalar, MPI_DOUBLE, send_to, 1,
                   Buffer_Receive_muT, nBufferR_Scalar, MPI_DOUBLE, receive_from, 1, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      delete [] Buffer_Send_muT;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy conservative variables. ---*/
        node[iPoint]->SetmuT(Buffer_Receive_muT[iVertex]);
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_muT;
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CTurbSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CTurbSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
 
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CTurbSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
  
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}


void CTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
  su2double *Turb_i, *Turb_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool second_order  = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
  bool limiter       = (config->GetSpatialOrder() == SECOND_ORDER_LIMITER);
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Turbulent variables w/o reconstruction ---*/
    
    Turb_i = node[iPoint]->GetSolution();
    Turb_j = node[jPoint]->GetSolution();
    numerics->SetTurbVar(Turb_i, Turb_j);
    
    /*--- Grid Movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    if (second_order) {

      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->node[iPoint]->GetLimiter_Primitive();
        Limiter_j = solver_container[FLOW_SOL]->node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Turbulent variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter();
        Limiter_j = node[jPoint]->GetLimiter();
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          Solution_i[iVar] = Turb_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Solution_j[iVar] = Turb_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetTurbVar(Solution_i, Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Menter's first blending function (only SST)---*/
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetF1blending(node[iPoint]->GetF1blending(), node[jPoint]->GetF1blending());
    
    /*--- Compute residual, and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Convective and viscous fluxes across symmetry plane are equal to zero. ---*/

}

// mskim
void CTurbSolver::BC_Euler_Wall(CGeometry      *geometry,
								CSolver        **solver_container,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
								CConfig        *config,
								unsigned short val_marker) {
  
  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}

void CTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, density_old = 0.0, density = 0.0;
  
  bool adjoint = config->GetContinuous_Adjoint();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta = Vol / (config->GetCFLRedCoeff_Turb()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    
    /*--- Update and clip trubulent solution ---*/
    
    switch (config->GetKind_Turb_Model()) {
        
      case SA:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddClippedSolution(0, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint], lowerlimit[0], upperlimit[0]);
        }
        
        break;
        
      case SA_NEG:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddSolution(0, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint]);
        }
        
        break;

      case SST:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          
          if (compressible) {
            density_old = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_Old(0);
            density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          }
          if (incompressible) {
            density_old = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
            density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          }
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->AddConservativeSolution(iVar, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint*nVar+iVar], density, density_old, lowerlimit[iVar], upperlimit[iVar]);
          }
          
        }
        
        break;
      case SA_FIML:
        //JRH 04262017 - Same as case SA:
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddClippedSolution(0, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint], lowerlimit[0], upperlimit[0]);
        }

        break;
        
    }
  }
  
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CTurbSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double Density_nM1, Density_n, Density_nP1;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        
        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
        Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                              +1.0*Density_nM1*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
        }
        
      } else {
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                              +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
        }
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  } else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      
      /*--- Multiply by density at node i for the SST model ---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      } else {
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      }
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      
      /*--- Multiply by density at node j for the SST model ---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        Density_n = solver_container[FLOW_SOL]->node[jPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      } else {
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      }
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---  Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- Get the index for node i plus the boundary face normal ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Grid velocities stored at boundary node i ---*/
        
        GridVel_i = geometry->node[iPoint]->GetGridVel();
        
        /*--- Compute the GCL term by dotting the grid velocity with the face
         normal. The normal is negated to match the boundary convention. ---*/
        
        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
        
        /*--- Compute the GCL component of the source term for node i ---*/
        
        U_time_n = node[iPoint]->GetSolution_time_n();
        
        /*--- Multiply by density at node i for the SST model ---*/
        
        if (config->GetKind_Turb_Model() == SST) {
          Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
        } else {
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        }
        LinSysRes.AddBlock(iPoint, Residual);
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
    of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        
        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
        Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
            + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
        
      } else {
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
            + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}


void CTurbSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) {

  /*--- Restart the solution from file information ---*/
  
  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double dull_val, Area_Children, Area_Parent, *Solution_Fine;
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  string UnstExt, text_line;
  ifstream restart_file;
  string restart_filename = config->GetSolution_FlowFileName();
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time|| time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Open the restart file, throw an error if this fails. ---*/
  
  restart_file.open(restart_filename.data(), ios::in);
  if (restart_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
    exit(EXIT_FAILURE);
  }

  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  map<unsigned long,unsigned long> Global2Local;
  map<unsigned long,unsigned long>::const_iterator MI;

  /*--- Now fill array with the transform values only for local points ---*/
  
  for (iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[MESH_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
  }

  /*--- Read all lines in the restart file ---*/
  
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Skip flow variables ---*/
  
  unsigned short skipVars = 0;

  if (compressible) {
    if (nDim == 2) skipVars += 6;
    if (nDim == 3) skipVars += 8;
  }
  if (incompressible) {
    if (nDim == 2) skipVars += 5;
    if (nDim == 3) skipVars += 7;
  }

  /*--- The first line is the header ---*/
  
  getline (restart_file, text_line);

  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
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
      node[iPoint_Local]->SetSolution(Solution);

    }

  }

  /*--- Close the restart file ---*/
  
  restart_file.close();

  /*--- MPI solution and compute the eddy viscosity ---*/
  
  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

}

CTurbSASolver::CTurbSASolver(void) : CTurbSolver() { }

CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel) : CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint, index;
  su2double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Factor_nu_Engine, Factor_nu_ActDisk, dull_val;
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  restart_gate = false; //JRH 05022018
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  int rank = MASTER_NODE;
  int size;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  
  nVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nPoint_Global = nPointDomain;
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointDomain, &nPoint_Global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  inv_nPoint_Global = 1.0/su2double(nPoint_Global);
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  
  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {
    
    /*--- Define some auxiliar vector related with the residual ---*/
    
    Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    
    /*--- Define some structures for locating max residuals ---*/
    
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
    
    /*--- Define some auxiliar vector related with the solution ---*/
    
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
    
    /*--- Define some auxiliar vector related with the flow solution ---*/
    
    FlowPrimVar_i = new su2double [nDim+7]; FlowPrimVar_j = new su2double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    if (config->GetExtraOutput()) {
      if (nDim == 2) { nOutputVariables = 13; }
      else if (nDim == 3) { nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }
    
    /*--- Computation of gradients by least squares ---*/
    
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new su2double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new su2double [nDim];
      
      /*--- c vector := transpose(WA)*(Wb) ---*/
      
      Cvector = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Cvector[iVar] = new su2double [nDim];
    }
    
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0;
  

  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/

  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Inf  = 0.005*Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_Engine ---*/
  Factor_nu_Engine = config->GetNuFactor_Engine();
  nu_tilde_Engine  = Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Engine  = 0.005*Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_ActDisk ---*/
  Factor_nu_ActDisk = config->GetNuFactor_Engine();
  nu_tilde_ActDisk  = Factor_nu_ActDisk*Viscosity_Inf/Density_Inf;

  /*--- Eddy viscosity at infinity ---*/
  su2double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
  su2double muT_Inf;
  Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
  Ji_3 = Ji*Ji*Ji;
  fv1 = Ji_3/(Ji_3+cv1_3);
  muT_Inf = Density_Inf*fv1*nu_tilde_Inf;
  
  Local2Global = new unsigned long [nPointDomain];

  //Need to do this here to read values from text files to map to correct index on local grid
  if (restart || config->GetKind_Turb_Model() == SA_FIML) {
	  /*--- In case this is a parallel simulation, we need to perform the
	   Global2Local index transformation first. ---*/

	  map<unsigned long,unsigned long> Global2Local;
	  map<unsigned long,unsigned long>::const_iterator MI;

	  /*--- Now fill array with the transform values only for local points ---*/
	  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
		Local2Global[iPoint] = geometry->node[iPoint]->GetGlobalIndex();
	  }
  }
  /*--- Read all lines in the restart file ---*/
 // if (config->GetKind_Turb_Model()==SA_FIML) {
  unsigned long Fiml_Skip_Index = 0;
#ifdef HAVE_MPI
  if (config->GetKind_Turb_Model() == SA_FIML) {
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
	  //cout << "JRH Debugging: In Turb Solver MPI Process " << rank << " nPointDomain = " << nPointDomain << " nPoint = " << nPoint << " and begins with fiml variable " << Fiml_Skip_Index << endl;
	  //if (rank == MASTER_NODE) cout << "Total_Fiml_Nodes = " << Total_Fiml_Nodes << endl;
  }
#else
#endif
//}
  long iPoint_Local; unsigned long iPoint_Global = 0; string text_line; unsigned long iPoint_Global_Local = 0;
  /*--- Restart the solution from file information ---*/
  if (!restart || (iMesh != MESH_0)) {
	  //cout << "JRH Debugging Message: Number of Points (nPoint): " << nPoint << endl;
	  //cout << "JRH Debugging Message: Number of Design Variables (nDV): " << config->GetnDV() << endl;
	  //nPoints are the number of points
	  //nElem are the total number of elements regardless of shape of the element
	  if (config->GetKind_Turb_Model() == SA_FIML) {
		  /*--- In case this is a parallel simulation, we need to perform the
		   Global2Local index transformation first. ---*/

		  map<unsigned long,unsigned long> Global2Local;
		  map<unsigned long,unsigned long>::const_iterator MI;

		  /*--- Now fill array with the transform values only for local points ---*/
		  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
			Local2Global[iPoint] = geometry->node[iPoint]->GetGlobalIndex();
		  }
	  }
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
    	//JRH - Aliased CTurbSAVariable and CTurbVariable to accept iPoint so beta_fiml can be initialized 04122017
      if (config->GetKind_Turb_Model() != SA_FIML) {
    	  node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, config);
      }
      else {
    	  //node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, iPoint, config);
    	  node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, config);
      }
    }

    if (config->GetKind_Turb_Model() == SA_FIML && config->GetKindTrainNN() != WEIGHTS) {
		/*--- In case this is a parallel simulation, we need to perform the
		 Global2Local index transformation first. ---*/

		//map<unsigned long,unsigned long> Global2Local;
		//map<unsigned long,unsigned long>::const_iterator MI;

    	//JRH - Mimicking code in solver_adjoint_discrete.cpp
    	unsigned long nDV_Local = 0;
        for (iPoint_Local = 0; iPoint_Local < nPoint; iPoint_Local++ ) {
          /*--- Retrieve local index. If this node from the restart file lives
           on the current processor, we will load and instantiate the vars. ---*/

          //MI = Global2Local.find(iPoint_Global);
          //if (MI != Global2Local.end()) {
        	if (geometry->node[iPoint_Local]->GetDomain()) {
            //iPoint_Local = Global2Local[iPoint_Global];

        		//THIS WORKED SINGLE NODE 06122018
//        		node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
//        		node[iPoint_Local]->SetBetaFimlTrain(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
//        		geometry->node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);

        		node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);
        		node[iPoint_Local]->SetBetaFimlTrain(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);
        		geometry->node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);

        		//if (rank == MASTER_NODE) cout << "JRH Debugging: In turb solver setting node " << iPoint_Local << " to iDV =  " << nDV_Local+Fiml_Skip_Index << " = " << config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0 << endl;
        		nDV_Local++;
          }
        	else {
        		node[iPoint_Local]->SetBetaFiml(1.0);
        		node[iPoint_Local]->SetBetaFimlTrain(1.0);
        		geometry->node[iPoint_Local]->SetBetaFiml(1.0);
        	}
        }
    }

  }
  else {
#ifdef HAVE_MPI
	  if (rank == MASTER_NODE) cout << "JRH Debugging: Instantiating Turb SA Variables in CTurbSASolver constructor from restart file" << endl;
#else
#endif
    /*--- Restart the solution from file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    su2double Density, StaticEnergy, Laminar_Viscosity, nu, nu_hat, muT = 0.0, U[5];
    int Unst_RestartIter;

    /*--- Modify file name for multizone problems ---*/
    if (nZone >1)
      filename= config->GetMultizone_FileName(filename, iZone);
    
    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time) {
      if (adjoint) {
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

    /*--- Modify file name for a simple unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) {
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter()) - 1;
      } else {
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      }
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no turbulent restart file!!" << endl;
      exit(EXIT_FAILURE);
    }

    //Moved commented section to above because need to do this for beta_fiml regardless
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/

    map<unsigned long,unsigned long> Global2Local;
    map<unsigned long,unsigned long>::const_iterator MI;

    /*--- Now fill array with the transform values only for local points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }

    /*--- Read all lines in the restart file ---*/

    long iPoint_Local; unsigned long iPoint_Global = 0; string text_line; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
    
    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      su2double beta_temp;
      su2double Fw_temp;
      su2double delta_temp;
      su2double Prod_temp;
      su2double Dest_temp;
      su2double Chi_temp;
      su2double gam_temp;
      su2double fd_temp;
      su2double S_temp;
      su2double O_temp;
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
//        "PointID"	"x"	"y"	"Conservative_1"	"Conservative_2"	"Conservative_3"	"Conservative_4"	"Conservative_5"	"Pressure"
//"Temperature"	"C<sub>p</sub>"	"Mach"	"<greek>m</greek>"	"C<sub>f</sub>_x"	"C<sub>f</sub>_y"	"h"	"y<sup>+</sup>"	"<greek>m</greek><sub>t</sub>"
//"Beta_Fiml"	"Production"	"Destruction"	"S<sub>hat</sub>"	"<greek>X</greek>"	"<greek>d</greek>"	"f<sub>w</sub>"	"r"	"Strain_Mag"	"Vort_Mag"	"wall_dist"

        if (compressible) {
			  if (nDim == 2) {
				  if (!config->GetTrainNN()) point_line >> index >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];
	           	  else if (config->GetKind_Trans_Model() == BC) {
	           		  	  	  	  	  	  	 point_line >> index >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0] >> //1-8
	            		          		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> //9-16
	            		          		  dull_val>> dull_val >> beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> //17-24
	            		          		  Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp >> gam_temp >> dull_val;

	            	  }
	            	  else {
		  	  	  	  	  	  point_line >> index >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0] >> //1-8
		          		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> //9-16
		          		  dull_val>> dull_val >> beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> //17-24
		          		  Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp;
	            	  }
	          }
          if (nDim == 3) {
        	  if (!config->GetTrainNN()) point_line >> index >> dull_val >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0];
           	  else if (config->GetKind_Trans_Model() == BC) {
           		point_line >> index >> dull_val >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0] >>
            		          		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>>
            		          		  beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp >> gam_temp >> dull_val;

            	  }
            	  else {
            		  point_line >> index >> dull_val >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0] >>
            		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>>
            		  beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp;
            	  }
          }
          Density = U[0];
          if (nDim == 2)
            StaticEnergy = U[3]/U[0] - (U[1]*U[1] + U[2]*U[2])/(2.0*U[0]*U[0]);
          else
            StaticEnergy = U[4]/U[0] - (U[1]*U[1] + U[2]*U[2] + U[3]*U[3] )/(2.0*U[0]*U[0]);

          FluidModel->SetTDState_rhoe(Density, StaticEnergy);
          Laminar_Viscosity = FluidModel->GetLaminarViscosity();
          nu     = Laminar_Viscosity/Density;
          nu_hat = Solution[0];
          Ji     = nu_hat/nu;
          Ji_3   = Ji*Ji*Ji;
          fv1    = Ji_3/(Ji_3+cv1_3);
          muT    = Density*fv1*nu_hat;
          
        }
//        "PointID"	"x"	"y"	"Conservative_1"	"Conservative_2"	"Conservative_3"	"Conservative_4"	"Pressure"
//        "Temperature"	"C<sub>p</sub>"	"Mach"	"<greek>m</greek>"	"C<sub>f</sub>_x"	"C<sub>f</sub>_y"	"h"	"y<sup>+</sup>"	"<greek>m</greek><sub>t</sub>"
//        "Beta_Fiml"	"Production"	"Destruction"	"S<sub>hat</sub>"	"<greek>X</greek>"	"<greek>d</greek>"	"f<sub>w</sub>"	"r"	"Strain_Mag"	"Vort_Mag"	"wall_dist"
        if (incompressible) {
          if (nDim == 2) {
        	  if (!config->GetTrainNN()) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        	  else if (config->GetKind_Trans_Model() == BC) {
        		  	  	  	  	  	  	 point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> dull_val >> dull_val>> //1-9
        		  	  	  	  	  	  	 dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val >> //10-17
        		  	  	  	  	  	  	 beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> Fw_temp >> dull_val >>//18-25
        		  	  	  	  	  	  	 S_temp >> O_temp >> dull_val >> fd_temp >> gam_temp;

        	  }
        	  else {
        		  point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >>
        		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>>
        		  beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp;
        	  }
          }
          if (nDim == 3) {
        	  if (!config->GetTrainNN()) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
           	  else if (config->GetKind_Trans_Model() == BC) {
           		point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >>
            		          		  dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>>
            		          		  beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> delta_temp >> Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp >> gam_temp >> dull_val;

            	  }
            	  else {
            		  point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> //8
            		  Solution[0] >> dull_val >> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>> dull_val>>  //9-16
            		  dull_val>> dull_val>> dull_val>> beta_temp >> Prod_temp >> Dest_temp >> dull_val >> Chi_temp >> //17-24
            		  delta_temp >> Fw_temp >> dull_val >> S_temp >> O_temp >> dull_val >> fd_temp >> gam_temp >> dull_val;
            	  }
              }
      	  muT = muT_Inf;
        }
        
        /*--- Instantiate the solution at this node, note that the eddy viscosity should be recomputed ---*/
        node[iPoint_Local] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
        //node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(iPoint_Global, 0)+1.0);
        if (config->GetTrainNN()){
			node[iPoint_Local]->SetBetaFiml(beta_temp);
			node[iPoint_Local]->SetProduction(Prod_temp);
			node[iPoint_Local]->SetDestruction(Dest_temp);
			node[iPoint_Local]->SetChiSA(Chi_temp);
			node[iPoint_Local]->SetFwSA(Fw_temp);
			node[iPoint_Local]->SetDeltaCriterion(delta_temp);
			node[iPoint_Local]->SetGammaTrans(gam_temp);
			node[iPoint_Local]->SetDES_fd(fd_temp);
			node[iPoint_Local]->SetStrainMagnitude(S_temp);
			node[iPoint_Local]->SetVorticityMagnitude(O_temp);
        }
        iPoint_Global_Local++;
      }

    }
    
    unsigned long nDV_Local = 0;
    for (iPoint_Local = 0; iPoint_Local < nPoint; iPoint_Local++ ) {
        /*--- Retrieve local index. If this node from the restart file lives
         on the current processor, we will load and instantiate the vars. ---*/

        //MI = Global2Local.find(iPoint_Global);
        //if (MI != Global2Local.end()) {
      	if (geometry->node[iPoint_Local]->GetDomain()) {
          //iPoint_Local = Global2Local[iPoint_Global];

      		//THIS WORKED SINGLE NODE 06122018
//      		if (!config->GetTrainNN()) {
//      			node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
//      		}
//      		node[iPoint_Local]->SetBetaFimlTrain(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
//      		geometry->node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
      		if (config->GetKindTrainNN() != WEIGHTS) { //Don't do this if we're using weights as design variables
				if (!config->GetTrainNN()) {
					node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);
				}
				node[iPoint_Local]->SetBetaFimlTrain(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);
				geometry->node[iPoint_Local]->SetBetaFiml(config->GetDV_Value(Local2Global[iPoint_Local],0)+1.0);

				//if (rank == MASTER_NODE) cout << "JRH Debugging: In turb solver setting node " << iPoint_Local << " to iDV =  " << nDV_Local+Fiml_Skip_Index << " = " << config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0 << endl;
				nDV_Local++;
      		}
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
      node[iPoint] = new CTurbSAVariable(Solution[0], muT_Inf, nDim, nVar, config);
    }
    
    /*--- Close the restart file ---*/
    restart_file.close();

  }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  


  //INITIALIZE NEURAL NETWORK MATRICES
  //Initialize Weights for Neural Network Training - JRH 04172018
  //nBins = 5; //For some reason this also needs to be defined in ComputeEDF();
  jrh_debug = false;
  if (jrh_debug) cout << "JRH Debugging - Outside of NN Training in CTurbSASolver Constructor" << endl;
  if (config->GetTrainNN()) {
	  train_NN = true;
	  if (config->GetFilterShield()) filter_shield = true;
	  if (restart) restart_gate = true;
	  num_epoch = config->GetNumEpoch();
	  learn_rate = config->GetLearningRate();
	  num_nn_inputs = 4; //Hard-coded number of neural network inputs
	  nLayers = config->GetNHiddenLayers()+2;//0 - Input Layer, 2<->nLayers-2 - Hidden Layers, nLayer-1 - Output Layer
	  nNeurons = config->GetNNeurons();
	  unsigned short nBins = config->GetnBins();
	  kind_scale = config->GetKind_NN_Scaling();

	  edf1 = new su2double[nBins+1];
	  edf2 = new su2double[nBins+1];
	  edf3 = new su2double[nBins+1];
	  edf4 = new su2double[nBins+1];

	  ledf1 = new su2double[nBins+1];
	  ledf2 = new su2double[nBins+1];
	  ledf3 = new su2double[nBins+1];
	  ledf4 = new su2double[nBins+1];

	  isHoldout = new unsigned long [nPointDomain];

	  //initialize num_nodes[]
	  num_nodes = new unsigned long [nLayers];
	  num_nodes[0] = num_nn_inputs+1; //+1 if using bias nodes
	  num_nodes[nLayers-1] = 1; //2 if using bias nodes
	  for(unsigned short iLayer = 1; iLayer < nLayers-1; iLayer++) num_nodes[iLayer] = nNeurons;

	  //initialize num_inputs[]
	  num_inputs = new unsigned long [nLayers];
	  //num_inputs[0] = 0; //Input layer has no inputs
	  num_inputs[0] = num_nn_inputs+1; //Added +1 07042018 if using bias nodes
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

	  su2double frac_holdout = config->GetPercentHoldout()/100.0;
	  unsigned long nTrainSamples_Local = 0;
	  nTrainSamples = 0;
	  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		  if (restart_gate) {
			  isHoldout[iPoint] = 0; //Don't hold out any samples for adjoint
			  nTrainSamples_Local++;
		  }
		  else {
			  if (su2double(rand())/su2double(RAND_MAX) > frac_holdout) {
				  isHoldout[iPoint] = 0;
				  nTrainSamples_Local++;
			  }
			  else isHoldout[iPoint] = 1;
		  }
	  }
#ifdef HAVE_MPI
	  nTrainSamples = 0;
	  SU2_MPI::Allreduce(&nTrainSamples_Local,&nTrainSamples,1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
	  nTrainSamples = nTrainSamples_Local;
#endif

	  min_max_send = new su2double[size];
	  min_max_recv = new su2double[size];

	  if (rank == MASTER_NODE) cout << "Number of Training Samples: " << nTrainSamples << " out of " << nPoint_Global << " nodes --> " << (1.0-su2double(nTrainSamples)/su2double(nPoint_Global))*100.0 << " % Holdout" << endl;

	  su2double rand_temp = 0.0;
	  //initialize hidden layers
	  su2double scale_range = 1.0/su2double(sqrt(num_nn_inputs));
	  num_weights = 0;
	  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
		//  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
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
					  rand_temp = su2double(rand())*scale_range*2.0/su2double(RAND_MAX)-scale_range;
					  weights[iLayer][iInput][iNode] = rand_temp;
					  //Actually, get weight from design variable if we're using weights as the design variables  JRH 06282018
					  if (jrh_debug && rank==MASTER_NODE) cout << "iLayer = " << iLayer << ", iInput = " << iInput << ", iNode = " << iNode << endl;
					  if (jrh_debug && rank==MASTER_NODE && config->GetKindTrainNN() == WEIGHTS) cout << "Setting Weight " << num_weights << " using GetDV_Value to " << config->GetDV_Value(num_weights,0) << endl;
					  if (config->GetKindTrainNN() == WEIGHTS) weights[iLayer][iInput][iNode] = config->GetDV_Value(num_weights,0);
					  Ep[iLayer][iInput][iNode] = 0.0;
					  num_weights++; //JRH 09232018 - Removing input layer weights from costly computation
					  //if (jrh_debug) cout << "iLayer " << iLayer << " iInput " << iInput << " iNode " << " weight init to -> " << rand_temp << endl;
				  }
			  }
		  //}
	  }
	  weight_send = new su2double [num_weights];
	  weight_recv = new su2double [num_weights];

      /*************************************THIS CODE TO AVERAGE NN WEIGHTS ACROSS ALL PROCS***********************************/
      if (rank == MASTER_NODE && jrh_debug) cout << rank << " JRH Debugging - NN SSE is currently: " << sse << endl;
      //if (jrh_debug && rank == MASTER_NODE) cout << rank << " JRH Debugging - NN SSE is currently: " << sse << endl;
//      solver_container[FLOW_SOL]->SetTotal_Loss(sse);
#ifdef HAVE_MPI //Perform weight averaging over all nodes (all nodes have same weights after averaging (reduction
      //Unwrap weight array so reduction can be performed
      su2double invSize = 1.0/su2double(size);
      unsigned long iWeight = 0;
      for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) { //JRH 09232018 - Removing input layer weights from costly computation
    	  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
    		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
    			  weight_send[iWeight] = weights[iLayer][iInput][iNode]*invSize;
    			  iWeight++;
    		  }
    	  }
      }

      //Reduce unwrapped weights and broadcast to all nodes (all procs will have same weights after this op
      SU2_MPI::Allreduce(weight_send,weight_recv,num_weights, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      //Re-Wrap received weights into weight array
      iWeight = 0;
      for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
    	  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
    		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
//    			  cout << iLayer << " " << iNode << " " << iInput << weight_recv[iWeight] << endl;
    			  weights[iLayer][iInput][iNode] = weight_recv[iWeight];
    			  iWeight++;
    		  }
    	  }
      }
#endif

	  //Allocate regardless - makes delete[] easier...
      restart_f1 = new su2double[nPointDomain];
      restart_f2 = new su2double[nPointDomain];
      restart_f3 = new su2double[nPointDomain];
      restart_f4 = new su2double[nPointDomain];

	  //If restarting (like for adjoint solutions, need to set weights to restart values
	  if (restart && config->GetKindTrainNN() != WEIGHTS) {
		  string nn_file = "nn_weights.dat";
		  ifstream nn_restart_file;

		    /*--- Open the restart file, throw an error if this fails. ---*/
		   nn_restart_file.open(nn_file.data(), ios::in);
		   if (nn_restart_file.fail()) {
		     cout << "Restart detected with NN Training Requested but no weights file found!!" << endl;
		     exit(EXIT_FAILURE);
		   }
		   else {
			   if (rank == MASTER_NODE) cout << "JRH: Found nn_weights.dat: Loading neural network weights from file" << endl;
		   }


		  //First line is header
		  getline (nn_restart_file, text_line);

	      su2double f1,f2,f3,f4;
		  //next are scaling
//		  for (unsigned short iLine = 0; iLine < num_nn_inputs*2; iLine++) {
		  getline (nn_restart_file, text_line);
		  istringstream point_line(text_line);
		  point_line >> mean_f1 >> mean_f2 >> mean_f3 >> mean_f4 >> std_f1 >> std_f2 >> std_f3 >> std_f4;



	      //Re-Wrap received weights into weight array
	      unsigned long iWeight = 0;
	      for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) { //JRH 09232018 - Removing input layer weights from costly computation
	    	  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
	    		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
	//    			  cout << iLayer << " " << iNode << " " << iInput << weight_recv[iWeight] << endl;
	    			  getline (nn_restart_file, text_line);
	    			  istringstream point_line(text_line);
	    			  point_line >> weight_send[iWeight];
	    			  weights[iLayer][iInput][iNode] = weight_send[iWeight];
	    			  iWeight++;
	    		  }
	    	  }
	      }

	      nn_restart_file.close();
	  }

	  if (jrh_debug) cout << "JRH Debugging - Done Declaring NN vectors and arrays and initializing weights" << endl;
  }

}

CTurbSASolver::~CTurbSASolver(void) {
  
}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;

  unsigned long ExtIter      = config->GetExtIter();
  bool limiter_flow          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool restart = (config->GetRestart() || config->GetRestart_Flow()); //JRH 04252018

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  /*--- Upwind second order reconstruction ---*/

  if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) SetSolution_Limiter(geometry, config);

  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);


  //Do one epoch of NN training. This registers these ops on the tape??
  if (config->GetTrainNN() && restart) {
	  //ForwardPropagate();//JRH 04242018
	  //SetDES_LengthScale(solver_container, geometry, config);
	  solver_container[FLOW_SOL]->SetTotal_Loss(sse);
  }
}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  
  su2double rho = 0.0, mu = 0.0, nu, *nu_hat, muT, Ji, Ji_3, fv1;
  su2double cv1_3 = 7.1*7.1*7.1;
  unsigned long iPoint;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool neg_spalart_allmaras = (config->GetKind_Turb_Model() == SA_NEG);
  
  /*--- Compute eddy viscosity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    if (compressible) {
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    if (incompressible) {
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    
    nu  = mu/rho;
    nu_hat = node[iPoint]->GetSolution();
    
    Ji   = nu_hat[0]/nu;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    
    muT = rho*fv1*nu_hat[0];
    
    if (neg_spalart_allmaras && (muT < 0.0)) muT = 0.0;
    
    node[iPoint]->SetmuT(muT);
    
  }
  //if (config->GetTrainNN()) ForwardPropagate();//JRH 04242018
  
}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                    CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
//seoyeon-0915. mskim: To support SA axisymmetric
  bool axisymmetric = config->GetAxisymmetric();

  bool harmonic_balance = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool transition    = (config->GetKind_Trans_Model() == LM);
  bool beta_fiml	=	(config->GetKind_Turb_Model() == SA_FIML); //JRH - 04032017

  //START OF JRH NEURAL NETWORK CODE!!!!

  //Now Prepare Features
  su2double f1,f2,f3,f4,lmean_f1,lmean_f2,lmean_f3,lmean_f4,
  	  lstd_f1,lstd_f2,lstd_f3,lstd_f4;
  su2double inv_nPointDomain = 1.0/su2double(nPointDomain);
//  int nPoint_Global;
  int rank = MASTER_NODE;
  int size;

  /* Compute Current Scaling for Features */
  //Note this must be done across all domains so need MPI
  //COMPUTE MEANS
//  if(jrh_debug) cout << "JRH Debugging - Computing mean and stdDev for scale factors"<<endl;
  sse = 0.0;

  //Need to make sure training path happens within the first 2 iterations if building AD graph.
  //But in direct computation need to converge flow vars a little for stability before training on them
  //Logic to delay training in direct mode but run after first iteration if we have restart variables - JRH 05022018
//  bool train_NN = ((config->GetTrainNN() && config->GetExtIter() >= config->GetIterStartNN()) || (restart_gate && config->GetExtIter() >= 0));
  bool dont_stop = (config->GetExtIter() < config->GetIterStopNNScaling() || config->GetIterStopNNScaling() == 0); //New logic to stop updating beta
  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 //at user-defined point, JRH 07052019
  bool train_NN = ((config->GetTrainNN() && config->GetExtIter() >= config->GetIterStartNN() && dont_stop) || (restart_gate && config->GetExtIter() >= 0));
    //bool train_NN = (config->GetTrainNN() && config->GetExtIter() >= config->GetIterStartNN() && !restart_gate);
//  if (train_NN && config->GetExtIter() == 0) {
//	  node[iPoint]->SetStrainMagnitude(numerics->GetStrainMagnitude());
//	  node[iPoint]->SetVorticityMagnitude(numerics->GetVorticityMagnitude());
//  }
  if (train_NN) {
	  ForwardPropagate(config, solver_container, geometry);
  }
  solver_container[FLOW_SOL]->SetTotal_Loss(sse);//end of NN code

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Conservative variables w/o reconstruction ---*/

      numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);

      /*--- Gradient of the primitive and conservative variables ---*/

      numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

      /*--- Set vorticity and strain rate magnitude ---*/

      numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);

      numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);

      /*--- Set intermittency ---*/

      if (transition) {
        numerics->SetIntermittency(solver_container[TRANS_SOL]->node[iPoint]->GetIntermittency());
      }

      /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

      numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
      numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);


      /*--- JRH - Set FIML Correction Term ---*/
      if (beta_fiml) {
      	numerics->SetBetaFiml(node[iPoint]->GetBetaFiml()); //JRH - 04032017
      	geometry->node[iPoint]->SetBetaFiml(node[iPoint]->GetBetaFiml()); //JRH -04132018
      	solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(node[iPoint]->GetBetaFiml());
      }


      /*--- Set volume ---*/

      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

//seoyeon's contribution -0915
// mskim: To support SA axisymmetric. 
      if (axisymmetric){
      /*--- Set y coordinate ---*/
        numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      }


      /*--- Compute the source term ---*/

      numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

      /*--- Subtract residual and the Jacobian ---*/

      LinSysRes.SubtractBlock(iPoint, Residual);

      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

      //Since we have computed residual, store features for learning at each node in CVariable class - JRH 02062018
      if (beta_fiml) {
      	node[iPoint]->SetProduction(numerics->GetProduction());
      	node[iPoint]->SetDestruction(numerics->GetDestruction());
      	node[iPoint]->SetSTildeSA(numerics->GetSTildeSA());
      	node[iPoint]->SetChiSA(numerics->GetChiSA());
      	node[iPoint]->SetDeltaCriterion(numerics->GetDeltaCriterion());
      	node[iPoint]->SetFwSA(numerics->GetFwSA());
      	node[iPoint]->SetRSA(numerics->GetRSA());
      	node[iPoint]->SetStrainMagnitude(numerics->GetStrainMagnitude());
      	node[iPoint]->SetVorticityMagnitude(numerics->GetVorticityMagnitude());
      	node[iPoint]->SetGammaTrans(numerics->GetGammaTrans());
      	node[iPoint]->SetWallDist(numerics->GetWallDist());
      	node[iPoint]->SetkSALSA(numerics->GetkSALSA());
      }

    }
//    SetDES_LengthScale(solver_container, geometry, config);

  if (harmonic_balance) {
    
    su2double Volume, Source;
    unsigned short nVar_Turb = solver_container[TURB_SOL]->GetnVar();
    
    /*--- Loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Access stored harmonic balance source term ---*/
      
      for (unsigned short iVar = 0; iVar < nVar_Turb; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
}

void CTurbSASolver::WriteNNWeights() {
	ofstream restart_file;
	string filename = "nn_weights.dat";
	restart_file.open(filename.c_str(), ios::out);
	restart_file.precision(15);
	restart_file << num_nn_inputs << " Inputs," << nLayers-2 << " Hidden Layers," << num_nodes[2] << " Hidden Nodes \n";
	restart_file << mean_f1 << " " << mean_f2 << " " << mean_f3 << " " << mean_f4 << " " << std_f1 << " " << std_f2 << " " << std_f3 << " " << std_f4 << "\n";
    unsigned long iWeight = 0;
    for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) { //JRH 09232018 - Removing input layer weights from costly computation
  	  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
  		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
  			  weight_send[iWeight] = weights[iLayer][iInput][iNode];
  			  restart_file << weight_send[iWeight] << "\n";
  			  iWeight++;
  		  }
  	  }
    }

//    su2double f1,f2,f3,f4;
//    for (unsigned long iPoint = 0; iPoint < nPointDomain ; iPoint++) {
    	  //su2double gam = node[iPoint]->GetGammaTrans();
//    	  su2double gam = 1.0;
////		  f1 = node[iPoint]->GetFwSA()*gam;
//    	  f1 = node[iPoint]->GetDES_fd();
//		  f2 = node[iPoint]->GetChiSA()*gam;
//		  f3 = node[iPoint]->GetDeltaCriterion()*gam;
//		  f4 = node[iPoint]->GetProduction()*gam/(node[iPoint]->GetDestruction()+1.0);

    	//f1 = node[iPoint]->GetStrainMag();


//		  f1 = (f1-mean_f1)/std_f1;
//		  f2 = (f2-mean_f2)/std_f2;
//		  f3 = (f3-mean_f3)/std_f3;
//		  f4 = (f4-mean_f4)/std_f4;

		  //restart_file << node[iPoint]->GetBetaFiml() << " " << f1 << " " << f2 << " " << f3 << " " << f4 << "\n";
//    }
    restart_file.close();

//	ofstream restart_file;
//	string filename = "beta_fiml_grad.dat";
//    restart_file.open(filename.c_str(), ios::out);
//    restart_file.precision(15);
//	for (unsigned long iDV = 0; iDV < nDV_total; iDV++) {
//		restart_file << Total_Sens_Beta_Fiml[iDV] << "\n";
//	}
//	restart_file.close();
}
void CTurbSASolver::ForwardPropagate(CConfig *config, CSolver **solver_container,CGeometry *geometry) {

	su2double f1,f2,f3,f4,lmean_f1,lmean_f2,lmean_f3,lmean_f4,lstd_f1,lstd_f2,lstd_f3,lstd_f4,local_sse,gam;
	su2double inv_nPointDomain = 1.0/su2double(nPointDomain);
	su2double inv_nTrainSamples = 1.0/su2double(nTrainSamples);
	unsigned short nBins = config->GetnBins();
	su2double dBin = 1.0/su2double(nBins);
	su2double invMaxmMinf1, invMaxmMinf2, invMaxmMinf3, invMaxmMinf4;
	su2double inv_dBin = 1.0/dBin;
	unsigned short iBin;
	su2double bin_max;
	unsigned long iPoint;
	kind_scale = config->GetKind_NN_Scaling();
	unsigned long iter_stop_rescale = config->GetIterStopNNScaling();
//	bool stop_rescale = (config->GetExtIter() >= iter_stop_rescale && iter_stop_rescale != 0);
	bool stop_rescale = false;
	bool filter_point = false; //JRH 07182019
	bool first_point = true; //JRH - Sentry to detect first iPoint that passes filter. 07192019, set to false at first iPoint that passes filter
	su2double node_count = 0;
	su2double l1,l2,l3,l4,inv_l1,inv_l2,inv_l3,inv_l4; //Lambdas for Box-Cox Scaling JRH 07182019
	l1 = 0;
	l2 = 0;
	l3 = 0;
	l4 = 0;

	int rank = MASTER_NODE;
	int size;
#ifdef HAVE_MPI
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

	  if (jrh_debug && rank == MASTER_NODE) cout << "stop_rescale = " << stop_rescale << " Iteration = " << config->GetExtIter() << " iter_stop_rescale = " << iter_stop_rescale << endl;

	  su2double inv_std_f1 = 0.0;
	  su2double inv_std_f2 = 0.0;
	  su2double inv_std_f3 = 0.0;
	  su2double inv_std_f4 = 0.0;
	  if (kind_scale == Z_SCALE) {
		  if (jrh_debug && rank == MASTER_NODE) cout << " means " << mean_f1 << " " << mean_f2 << " " << mean_f3 << " " << mean_f4 << endl;
		  if (jrh_debug && rank == MASTER_NODE) cout << " stddevs " << std_f1 << " " << std_f2 << " " << std_f3 << " " << std_f4 << endl;
		  if (stop_rescale == false) {
			  //if (jrh_debug) cout << " means " << mean_f1 << " " << mean_f2 << " " << mean_f3 << " " << mean_f4 << endl;
			  //if (jrh_debug) cout << " stddevs " << std_f1 << " " << std_f2 << " " << std_f3 << " " << std_f4 << endl;
			  lmean_f1 = 0.0;
			  lmean_f2 = 0.0;
			  lmean_f3 = 0.0;
			  lmean_f4 = 0.0;

			  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
					  filter_point = false;
					  if (filter_shield == true) {
						  f3 = node[iPoint]->GetDeltaCriterion();
						  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
						  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
					  }

					  if ((filter_point == false && filter_shield == true) || filter_shield == false) {
//				  if (1) {
				  //if (isHoldout[iPoint] == 0) {
//				  gam = node[iPoint]->GetGammaTrans();
//				  gam = 1.0;
////				  f1 = node[iPoint]->GetFwSA()*gam;
//				  f1 = node[iPoint]->GetDES_fd();
//				  f2 = node[iPoint]->GetChiSA()*gam;
//				  f3 = node[iPoint]->GetDeltaCriterion()*gam;
//				  f4 = node[iPoint]->GetProduction()*gam/(node[iPoint]->GetDestruction()+1.0);

//					  su2double wall_dist = geometry->node[iPoint]->GetWall_Distance();
//					  su2double strain = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//					  su2double density = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
//					  su2double lam_visc = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//					  f1 = 1.0; //Bias Node
//					  f2 = node[iPoint]->GetDeltaCriterion();
//					  f3 = wall_dist;
//					  f4 = density*strain*pow(wall_dist,2.0)/lam_visc;

//			    	  f1 = node[iPoint]->GetDES_fd();
//					  f1 = 1.0;
					  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
			    	  f2 = node[iPoint]->GetChiSA();
					  f3 = node[iPoint]->GetDeltaCriterion();
//					  f4 = node[iPoint]->GetFwSA();
					  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));

					  lmean_f1 += f1;
					  lmean_f2 += f2;
					  lmean_f3 += f3;
					  lmean_f4 += f4;
				  } // if (node[iPoint]->GetDES_fd() > 0.94)
				  node_count += 1.0;
			  }
#ifdef HAVE_MPI
			  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			  MPI_Comm_size(MPI_COMM_WORLD, &size);
			  feat_send[0] = lmean_f1;
			  feat_send[1] = lmean_f2;
			  feat_send[2] = lmean_f3;
			  feat_send[3] = lmean_f4;
		  //	  SU2_MPI::Allreduce(&lmean_f1, &mean_f1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lmean_f2, &mean_f2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lmean_f3, &mean_f3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lmean_f4, &mean_f4, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			  SU2_MPI::Allreduce(feat_send, feat_recv, num_nn_inputs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			  mean_f1 = feat_recv[0];
			  mean_f2 = feat_recv[1];
			  mean_f3 = feat_recv[2];
			  mean_f4 = feat_recv[3];
			  //SU2_MPI::Allreduce(&nPointDomain, &nPoint_Global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); //Should do this in the constructor!!
#else
		   mean_f1 = lmean_f1;
		   mean_f2 = lmean_f2;
		   mean_f3 = lmean_f3;
		   mean_f4 = lmean_f4;
		  //nPoint_Global = nPointDomain;
#endif
//					mean_f1 = mean_f1*inv_nTrainSamples;
//					mean_f2 = mean_f2*inv_nTrainSamples;
//					mean_f3 = mean_f3*inv_nTrainSamples;
//					mean_f4 = mean_f4*inv_nTrainSamples;
		   	   	    mean_f1 = mean_f1/node_count;
		   	   		mean_f2 = mean_f2/node_count;
		   	   		mean_f3 = mean_f3/node_count;
		   	   		mean_f4 = mean_f4/node_count;
					//COMPUTE STD DEVs
					lstd_f1 = 0.0;
					lstd_f2 = 0.0;
					lstd_f3 = 0.0;
					lstd_f4 = 0.0;

				  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
						  filter_point = false;
						  if (filter_shield == true) {
							  f3 = node[iPoint]->GetDeltaCriterion();
							  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
							  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
						  }

						  if ((filter_point == false && filter_shield == true) || filter_shield == false) {//					  if ((node[iPoint]->GetDES_fd() < 0.999 && filter_shield == true) || filter_shield == false) {
//					  if (1) {
						  //if (isHoldout[iPoint] == 0) {
	////					  gam = node[iPoint]->GetGammaTrans();
	//					  gam = 1.0;
	////					  f1 = node[iPoint]->GetFwSA()*gam;
	//					  f1 = node[iPoint]->GetDES_fd();
	//					  f2 = node[iPoint]->GetChiSA()*gam;
	//					  f3 = node[iPoint]->GetDeltaCriterion()*gam;
	//					  f4 = node[iPoint]->GetProduction()*gam/(node[iPoint]->GetDestruction()+1.0);

//						  su2double wall_dist = geometry->node[iPoint]->GetWall_Distance();
//						  su2double strain = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//						  su2double density = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
//						  su2double lam_visc = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//						  f1 = 1.0; //Bias Node
//						  f2 = node[iPoint]->GetDeltaCriterion();
//						  f3 = wall_dist;
//						  f4 = density*strain*pow(wall_dist,2.0)/lam_visc;

						  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
				    	  f2 = node[iPoint]->GetChiSA();
						  f3 = node[iPoint]->GetDeltaCriterion();
//						  f4 = node[iPoint]->GetFwSA();
						  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));


						  lstd_f1 += pow(f1-mean_f1,2);
						  lstd_f2 += pow(f2-mean_f2,2);
						  lstd_f3 += pow(f3-mean_f3,2);
						  lstd_f4 += pow(f4-mean_f4,2);
					  } //if (node[iPoint]->GetDES_fd() > 0.94) {
			  }
#ifdef HAVE_MPI
			  feat_send[0] = lstd_f1;
			  feat_send[1] = lstd_f2;
			  feat_send[2] = lstd_f3;
			  feat_send[3] = lstd_f4;
			  SU2_MPI::Allreduce(feat_send, feat_recv, num_nn_inputs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			  std_f1 = feat_recv[0];
			  std_f2 = feat_recv[1];
			  std_f3 = feat_recv[2];
			  std_f4 = feat_recv[3];
		  //	  SU2_MPI::Allreduce(&lstd_f1, &std_f1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lstd_f2, &std_f2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lstd_f3, &std_f3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  //	  SU2_MPI::Allreduce(&lstd_f4, &std_f4, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
				std_f1 = lstd_f1;
				std_f2 = lstd_f2;
				std_f3 = lstd_f3;
				std_f4 = lstd_f4;

#endif
//					std_f1 = sqrt(std_f1*inv_nTrainSamples);
//					std_f2 = sqrt(std_f2*inv_nTrainSamples);
//					std_f3 = sqrt(std_f3*inv_nTrainSamples);
//					std_f4 = sqrt(std_f4*inv_nTrainSamples);
				std_f1 = sqrt(std_f1/node_count);
				std_f2 = sqrt(std_f2/node_count);
				std_f3 = sqrt(std_f3/node_count);
				std_f4 = sqrt(std_f4/node_count);
		  } //<-- if (stop_rescale == false)
				inv_std_f1 = 1.0/std_f1;
				inv_std_f2 = 1.0/std_f2;
				inv_std_f3 = 1.0/std_f3;
				inv_std_f4 = 1.0/std_f4;
	  }
	  else if (kind_scale == Q_TRANSFORM) {
		  if (stop_rescale == false) {
			  if (jrh_debug && rank == MASTER_NODE) cout << rank << "JRH Debugging: About to compute min and max of all features" << endl;
			  FindFeaturesMinMax();

			  if (jrh_debug && rank == MASTER_NODE) cout << rank << ": min_f1: " << min_f1 << ": min_f2: " << min_f2 << ": min_f3: " << min_f3 << ": min_f4: " << min_f4 << endl;
			  if (jrh_debug && rank == MASTER_NODE) cout << rank << ": max_f1: " << max_f1 << ": max_f2: " << max_f2 << ": max_f3: " << max_f3 << ": max_f4: " << max_f4 << endl;


			  if (jrh_debug && rank == MASTER_NODE) cout << rank << "JRH Debugging: About to compute EDF of all features" << endl;
			  ComputeEDF(nBins);
			  if (jrh_debug && rank == MASTER_NODE) cout << rank << "JRH Debugging: Back from computing EDF of all features" << endl;
		  }
    	  invMaxmMinf1 = 1.0/(max_f1-min_f1);
    	  invMaxmMinf2 = 1.0/(max_f2-min_f2);
    	  invMaxmMinf3 = 1.0/(max_f3-min_f3);
    	  invMaxmMinf4 = 1.0/(max_f4-min_f4);
       }
      else if (kind_scale == MIN_MAX) {
    	  if (stop_rescale == false) FindFeaturesMinMax();
    	  invMaxmMinf1 = 1.0/(max_f1-min_f1);
    	  invMaxmMinf2 = 1.0/(max_f2-min_f2);
    	  invMaxmMinf3 = 1.0/(max_f3-min_f3);
    	  invMaxmMinf4 = 1.0/(max_f4-min_f4);
      }
      else if (kind_scale == MAN_Z_SCALE) {
    	  //Hard coded from RAE2822_fine baseline with backprop NN_Weights file for C6 - JRH 07152018
    	  //Testing to see if this helps convergence (scale factor not chaning every iteration
    	  //Also scaling required for initialization when using weights as design variables to avoid very large initial beta's
//    	  0.265498448168419 43.3673570094839 0.102649025382011 8.78702527792828e-07 0.436832600884541 124.934481484434 0.199031871975153 4.22143821731118e-06
    	  mean_f1 = 0.265498448168419;
    	  mean_f2 = 43.3673570094839;
    	  mean_f3 = 0.102649025382011;
    	  mean_f4 = 8.78702527792828e-07;

    	  std_f1 = 0.436832600884541;
    	  std_f2 = 124.934481484434;
    	  std_f3 = 0.199031871975153;
    	  std_f4 = 4.22143821731118e-06;

		  inv_std_f1 = 1.0/std_f1;
		  inv_std_f2 = 1.0/std_f2;
		  inv_std_f3 = 1.0/std_f3;
		  inv_std_f4 = 1.0/std_f4;
      }
      else if (kind_scale == BOX_COX) {
    	  //Init lambdas for box-cox JRH 07182019
    	  if (filter_shield) {
			  l1 = -0.09146872915503056;
			  l2 = -0.049611844047549786;
			  l3 = -0.3085491501297445;
			  l4 = 0.3506666390043312;
    	  }
		  else {
			  l1 = -0.0350928110267145;
			  l2 = -0.40273798876730166;
			  l3 = -0.11597545566454645;
			  l4 = -0.0077125590994702965;
		  }

    	  inv_l1 = 1.0/l1;
    	  inv_l2 = 1.0/l2;
    	  inv_l3 = 1.0/l3;
    	  inv_l4 = 1.0/l4;
      }
	  //EPOCH LOOP START!!
	  if (config->GetKindTrainNN() != WEIGHTS || restart_gate) {
//	  if (config->GetKindTrainNN() != WEIGHTS) {
      //su2double sse = 0.0;
      for (unsigned long iEpoch = 0; iEpoch < num_epoch; iEpoch++) {

    	  sse = 0.0;
    	  //First re-initialize Ep to zero everywhere:
    	  // (Done More Efficiently Now Inside iPoint Loop when iPoint == 0;
	  // Not sure the more efficient code was working, definitely had a bug at some point, uncommented JRH  07202019
	  if (jrh_debug) cout << "Initializing Ep to zero in ForwardPropagate()" << endl;
	  //for (unsigned short iLayer = 0; iLayer < nLayers; iLayer++) {
    	//	  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
    	//		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer]; iInput++) {
    	//			  Ep[iLayer][iInput][iNode] = 0.0;
    	//		  }
    	//	  }
    	  //}
	  for (unsigned short iLayer = nLayers-1; iLayer > 0 ; iLayer--) {
				  //iNode --> j
				  //iInput --> i
		for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
			 for (unsigned long idInput = 0; idInput < num_nodes[iLayer-1]; idInput++) {
				Ep[iLayer][idInput][iNode] = 0.0;
			 }
		 }
	  }
	  if (jrh_debug) cout << "Done Initializing Ep" << endl;

      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    		  filter_point = false;
    		  if (filter_shield == true) {
    			  f3 = node[iPoint]->GetDeltaCriterion();
    			  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
    			  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
    		  }

    		  if ((filter_point == false && filter_shield == true) || filter_shield == false) {    	  
			
//    	  if ((node[iPoint]->GetDES_fd() < 0.999 && filter_shield == true) || filter_shield == false) {
//		  f1 = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//		  f2 = node[iPoint]->GetChiSA();
//		  f3 = node[iPoint]->GetDeltaCriterion();

////    	  gam = node[iPoint]->GetGammaTrans();
//    	  gam = 1.0;
////		  f1 = node[iPoint]->GetFwSA()*gam;
//    	  f1 = node[iPoint]->GetDES_fd();
//    	  f2 = node[iPoint]->GetChiSA()*gam;
//		  f3 = node[iPoint]->GetDeltaCriterion()*gam;
//		  f4 = node[iPoint]->GetProduction()*gam/(node[iPoint]->GetDestruction()+1.0);

//		  su2double wall_dist = geometry->node[iPoint]->GetWall_Distance();
//		  su2double strain = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//		  su2double density = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
//		  su2double lam_visc = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//		  f1 = 1.0; //Bias Node
//		  f2 = node[iPoint]->GetDeltaCriterion();
//		  f3 = wall_dist;
//		  f4 = density*strain*pow(wall_dist,2.0)/lam_visc;

		  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
		  f2 = node[iPoint]->GetChiSA();
		  f3 = node[iPoint]->GetDeltaCriterion();
//		  f4 = node[iPoint]->GetFwSA();
		  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));

		  if (kind_scale == Z_SCALE || kind_scale == MAN_Z_SCALE) {
//			  f1 = (f1-mean_f1)*inv_std_f1;
			  f2 = (f2-mean_f2)*inv_std_f2;
			  f3 = (f3-mean_f3)*inv_std_f3;
			  f4 = (f4-mean_f4)*inv_std_f4;
		  }
		  else if (kind_scale == Q_TRANSFORM) {
			  	//if (jrh_debug) cout << "iPoint = " << iPoint << " applying Q_TRANSFORM" << endl;
				f1 = (f1-min_f1)*invMaxmMinf1*0.998+0.001;
				f2 = (f2-min_f2)*invMaxmMinf2*0.998+0.001;
				f3 = (f3-min_f3)*invMaxmMinf3*0.998+0.001;
				f4 = (f4-min_f4)*invMaxmMinf4*0.998+0.001;
			  //interpoplate feature onto edf and apply logit function
			  //NOTE: Logit approximates the true inverse cdf function for Gaussian distribution, but is not exact! JRH 05132018


				//if (jrh_debug) cout << "f1 min_max scaled = " << f1 << endl;
				//FEATURE 1 INTERPOLATION
				iBin = 1;
				bin_max = dBin;
				while (bin_max < f1) {
					iBin++;
					bin_max += dBin;
				}
				iBin--; //Now at
				f1 = edf1[iBin]+(f1-(bin_max-dBin))*(edf1[iBin+1]-edf1[iBin])*inv_dBin;


				//if (jrh_debug) cout << "f2 min_max scaled = " << f2 << endl;
				//FEATURE 2 INTERPOLATION
				iBin = 1;
				bin_max = dBin;
				while (bin_max < f2) {
					iBin++;
					bin_max += dBin;
				}
				iBin--; //Now at
				f2 = edf2[iBin]+(f2-(bin_max-dBin))*(edf2[iBin+1]-edf2[iBin])*inv_dBin;

				//if (jrh_debug) cout << "f3 min_max scaled = " << f3 << endl;
				//FEATURE 3 INTERPOLATION
				iBin = 1;
				bin_max = dBin;
				while (bin_max < f3) {
					iBin++;
					bin_max += dBin;
				}
				iBin--; //Now at
				f3 = edf3[iBin]+(f3-(bin_max-dBin))*(edf3[iBin+1]-edf3[iBin])*inv_dBin;


				//if (jrh_debug) cout << "f4 min_max scaled = " << f4 << endl;
				//FEATURE 4 INTERPOLATION
				iBin = 1;
				bin_max = dBin;
				while (bin_max < f4) {
					iBin++;
					bin_max += dBin;
				}
				iBin--; //Now at
				f4 = edf4[iBin]+(f4-(bin_max-dBin))*(edf4[iBin+1]-edf4[iBin])*inv_dBin;


				f1 = log(f1/(1.0-f1));
				f2 = log(f2/(1.0-f2));
				f3 = log(f3/(1.0-f3));
				f4 = log(f4/(1.0-f4));

				if (jrh_debug && rank==MASTER_NODE) cout << "Features Q-Scaled = " << f1 << " " << f2 << " " << f3 << " " << f4 << endl;

		  }
		  else if (kind_scale == MIN_MAX) {
				f1 = (f1-min_f1)*invMaxmMinf1*0.998+0.001;
				f2 = (f2-min_f2)*invMaxmMinf2*0.998+0.001;
				f3 = (f3-min_f3)*invMaxmMinf3*0.998+0.001;
				f4 = (f4-min_f4)*invMaxmMinf4*0.998+0.001;
		  }
		  else if (kind_scale == BOX_COX) {
			  if (jrh_debug) cout << "Lambdas = " << l1 <<" "<< l2<<" "<< l3<<" " << l4<<" " << " inv Lambdas = " <<inv_l1 <<" " << inv_l2<<" " << inv_l3<<" " << inv_l4<<" "<<"\n";
			  if (filter_shield) {
				  f1 = (pow(f1,l1)-1.0)*inv_l1;
				  f2 = (pow(f2,l2)-1.0)*inv_l2;
				  f3 = (pow(f3,l3)-1.0)*inv_l3;
				  f4 = (pow(f4,l4)-1.0)*inv_l4;
			  }
			  else {
				  f1 = (pow(f1+pow(10.0,-16.0),l1)-1.0)*inv_l1;
				  f2 = (pow(f2+pow(10.0,-16.0),l2)-1.0)*inv_l2;
				  f3 = (pow(f3+pow(10.0,-16.0),l3)-1.0)*inv_l3;
				  f4 = (pow(f4+pow(10.0,-16.0),l4)-1.0)*inv_l4;
			  }
			  if (jrh_debug) cout << "BOX_COX Features = " << f1<<" " <<f2<<" "<<f3<<" "<<f4<<"\n";
		  }
		  else if (kind_scale != NO_SCALE) {
			  cout << "JRH: ERROR, ONLY TWO IMPLEMENTED NN SCALINGS ARE Z_SCALE, MIN_MAX, AND Q_TRANSFORM" << endl;
		  }

		  //if (jrh_debug) cout << "iPoint " << iPoint << " Features " << f1 << "   " <<  f2 << "  " << f3 << endl;
		  su2double temp = 0.0;
		  ai[0][0] = 1.0;
		  ai[0][1] = f1;
		  ai[0][2] = f2;
		  ai[0][3] = f3;
		  ai[0][4] = f4;
		  inputs[0][0] = 1.0;
		  inputs[0][1] = f1;
		  inputs[0][2] = f2;
		  inputs[0][3] = f3;
		  inputs[0][4] = f4;
		  //cout << "Starting Forward Prop" << endl;
		  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
			  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
				  ai[iLayer][iNode] = 0.0;
				  //Apply Weights to Inputs
				  //if (iNode>0) {
					  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
						  //if (jrh_debug) cout << " iLayer " << iLayer << " iNode " << iNode << " iInput " << iInput << endl;
						  ai[iLayer][iNode] += weights[iLayer][iInput][iNode]*inputs[iLayer-1][iInput];
					  }
				  //}
					  //Apply Activation Function

					  if (iLayer < nLayers-1){
						  //tanh Activation Function
						  inputs[iLayer][iNode] = (exp(2.0*ai[iLayer][iNode])-1.0)/(exp(2.0*ai[iLayer][iNode])+1.0);

						  //Sigmoid Activation Function
						  //inputs[iLayer][iNode] = 1.0/(1.0+exp(-ai[iLayer][iNode]));

						  //Leaky RELU Activation Function - https://towardsdatascience.com/activation-functions-neural-networks-1cbd9f8d91d6
						  //if (ai[iLayer][iNode] <= 0.0) inputs[iLayer][iNode] = 0.01*ai[iLayer][iNode];
						  //else inputs[iLayer][iNode] = ai[iLayer][iNode];

						  //TESTING LINEAR ACTIVATION FUNCTION for WEIGHTS AS DVs...CHANGE TO BACKPROP NOT IMPLEMENTED!!!!!
//						  inputs[iLayer][iNode] = ai[iLayer][iNode];
					  }
				  else inputs[iLayer][iNode] = ai[iLayer][iNode];
			  }
			  //BIAS NODES
			  if (iLayer < nLayers-1) inputs[iLayer][0] = 1.0; //Bias Node (no bias node on output layer (-1))...or hidden layer (-2)
		  }
//		  inputs[nLayers-1][0] += 1.0;//Bias node in output layer
		  //Set predicted beta_fiml to beta_fiml in CVariable Class
		  //if (jrh_debug) cout << "Predicted beta = " << inputs[nLayers-1][0] << endl;
		  node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0);
		  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0);
		  sse += pow(node[iPoint]->GetBetaFiml()-node[iPoint]->GetBetaFimlTrain(),2)*0.5;


      //cout << "JRH Debugging - Starting Backprop (when implemented)"<<endl;
  /* --- JRH - Perform Backpropagation Algorithm on Neural Network if Required ---*/ //JRH 04182018
		  //BUT ONLY IF THIS iPoint is in our training set
  /*********************************************************************************/
		  //if (isHoldout[iPoint] == 0) {
		  //Compute deltas for output layer
		  //cout << "Starting Back Prop Deltas" << endl;
		  if (config->GetKindTrainNN() == BACKPROP) deltas[nLayers-1][0] = -(node[iPoint]->GetBetaFimlTrain()-node[iPoint]->GetBetaFiml());

		  else deltas[nLayers-1][0] = -(1.0-node[iPoint]->GetBetaFiml());//To compute change in regularization term of objective function (not cost!!)
		  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 //FLIPPED SIGN OF ABOVE to positive 06102019 JRH
		  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 //FLIPPED BACK AFTER I THOUGHT ABOUT IT.
		  //Compute deltas for hidden layers (output layer not necessary)
		  //for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			  for (unsigned short iLayer = nLayers-2; iLayer > 0 ; iLayer--) {
					  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
						  deltas[iLayer][iNode] = 0.0;
						  //idInput --> l -> Starts at 1 because of bias node... REMOVING BIAS NODE, setting to 0
						  //iNode --> j
						  for (unsigned long idInput = 0; idInput < num_nodes[iLayer+1]; idInput++) {
							deltas[iLayer][iNode] += deltas[iLayer+1][idInput]*weights[iLayer+1][iNode][idInput]; // note flip on weights index!!
						  }
						  //Multiply by derivative of activation function for hidden nodes
						  //ELU
//						  if (ai[iLayer][iNode] <= 0.0) deltas[iLayer][iNode] = exp(ai[iLayer][iNode]-1.0)*deltas[iLayer][iNode];
//						  else deltas[iLayer][iNode] = ai[iLayer][iNode]*deltas[iLayer][iNode];
						  //Hyperbolic Tangent (applying caching trick: https://theclevermachine.wordpress.com/2014/09/08/derivation-derivatives-for-common-neural-network-activation-functions/
						  if (iNode>0) deltas[iLayer][iNode] = (1.0-pow(inputs[iLayer][iNode],2.0))*deltas[iLayer][iNode];
						  //else deltas[iLayer][iNode] = deltas[iLayer][iNode]; //Bias Nodes https://stackoverflow.com/questions/3775032/how-to-update-the-bias-in-neural-network-backpropagation
					  }

			  }

			  //cout << "Starting Back Prop Ep" << endl;
			  //idInput --> i
			  //iNode --> j
			  for (unsigned short iLayer = nLayers-1; iLayer > 0 ; iLayer--) {
				  //iNode --> j
				  //iInput --> i
				  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
					  for (unsigned long idInput = 0; idInput < num_nodes[iLayer-1]; idInput++) {
						  //if (iPoint == 0) Ep[iLayer][idInput][iNode] = 0.0; //BUG HERE, DON'T ALWAYS START WITH ipoint==1 if we are filtering!! JRH 07192019 - THIS ALSO AFFECTS ALL f_d FILTERING RUNS PRIOR TO THIS DATE...
						  //if (first_point) Ep[iLayer][idInput][iNode] = 0.0;
//						  Ep[iLayer][idInput][iNode] += inputs[iLayer-1][idInput]*deltas[iLayer][iNode]*inv_nPointDomain;
						  Ep[iLayer][idInput][iNode] += inputs[iLayer-1][idInput]*deltas[iLayer][iNode];
						  //if (iNode > 0) Ep[iLayer][idInput][iNode] += inputs[iLayer-1][idInput]*deltas[iLayer][iNode];
						  //else Ep[iLayer][idInput][iNode] += deltas[iLayer][iNode]*inv_nPointDomain;
					  }
				  }
			  }
			  // if (first_point == true) first_point == false; //Deactivate sentry for first point JRH 07182019
		  //} //end of isHoldout if statement
    	  } //end of fd filtering if statement
    	  else {
    		  node[iPoint]->SetBetaFiml(1.0);
    		  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(1.0);
    	  }
     }//end of nPoint loop
      /*************************************THIS CODE TO SUM GRADIENTS ACROSS ALL PROCS***********************************/
      //Need to initialize weights to same values otherwise this won't work...
      if (rank == MASTER_NODE && jrh_debug) cout << rank << " JRH Debugging - NN SSE is currently: " << sse << endl;
      //if (jrh_debug && rank == MASTER_NODE) cout << rank << " JRH Debugging - NN SSE is currently: " << sse << endl;
//      solver_container[FLOW_SOL]->SetTotal_Loss(sse);
      //cout << "Starting MPI Share of Ep" << endl;
#ifdef HAVE_MPI //Perform weight averaging over all nodes (all nodes have same weights after averaging (reduction
      //Unwrap weight array so reduction can be performed
//      su2double invSize = 1.0/su2double(size);
      unsigned long iWeight = 0;
      for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
    	  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
    		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
    			  weight_send[iWeight] = Ep[iLayer][iInput][iNode];
    			  iWeight++;
    		  }
    	  }
      }

      //Reduce unwrapped weights and broadcast to all nodes (all procs will have same weights after this op
      SU2_MPI::Allreduce(weight_send,weight_recv,num_weights, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      //Re-Wrap received weights into weight array
      iWeight = 0;
      for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
    	  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
    		  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
//    			  cout << iLayer << " " << iNode << " " << iInput << weight_recv[iWeight] << endl;
    			  Ep[iLayer][iInput][iNode] = weight_recv[iWeight];
    			  iWeight++;
    		  }
    	  }
      }
#endif


  /*********************************************************************************/
      //Now Increment Weights
      //cout << "Starting Back Prop Weight Increment" << endl;
      if (config->GetKindTrainNN() != WEIGHTS) { //Don't increment weights if WEIGHTS and in Discrete Adjoint solver
		  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {  //JRH 09232018 - Removing input layer weights from costly computation
			  for (unsigned short iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
				  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
					  weights[iLayer][iInput][iNode] = weights[iLayer][iInput][iNode]-learn_rate*Ep[iLayer][iInput][iNode];
				  }
			  }
		  }
      }
      sse = sse*inv_nPointDomain;

      } //End of Epoch loop

} // <<-----if (config->GetKindTrainNN() != WEIGHTS || config->GetKind_Solver() == RANS)
  //end of NN code
      if (restart_gate == false || config->GetKindTrainNN() == BACKPROP) {
//	  if (1) {
	  //Forward propagate one more time to get update
      sse = 0.0;
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		  filter_point = false;
		  if (filter_shield == true) {
			  f3 = node[iPoint]->GetDeltaCriterion();
			  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
//			  if (jrh_debug) cout << "f3 = " << f3 << " f4 = " << f4 << endl;
			  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
//			  if (jrh_debug) cout << "filter_shield = " << filter_shield << " filter_point = " << filter_point << endl;
		  }

		  if ((filter_point == false && filter_shield == true) || filter_shield == false) {
    	  //    	  if ((node[iPoint]->GetDES_fd() < 0.999 && filter_shield == true) || filter_shield == false) {
//    	  if (1) { //commenting out filtering for now
	//		  f1 = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
	//		  f2 = node[iPoint]->GetChiSA();
	//		  f3 = node[iPoint]->GetDeltaCriterion();
	//    	  gam = node[iPoint]->GetGammaTrans();
	//    	  gam = 1.0;
	//		  f1 = node[iPoint]->GetFwSA()*gam;

//			  su2double wall_dist = geometry->node[iPoint]->GetWall_Distance();
//			  su2double strain = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//			  su2double density = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
//			  su2double lam_visc = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//			  f1 = 1.0; //Bias Node
//			  f2 = node[iPoint]->GetDeltaCriterion();
//			  f3 = wall_dist;
//			  f4 = density*strain*pow(wall_dist,2.0)/lam_visc;

			  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
	    	  f2 = node[iPoint]->GetChiSA();
			  f3 = node[iPoint]->GetDeltaCriterion();
//			  f4 = node[iPoint]->GetFwSA();
			  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));

			  if (kind_scale == Z_SCALE || kind_scale == MAN_Z_SCALE) {
//				  f1 = (f1-mean_f1)*inv_std_f1;
				  f2 = (f2-mean_f2)*inv_std_f2;
				  f3 = (f3-mean_f3)*inv_std_f3;
				  f4 = (f4-mean_f4)*inv_std_f4;
			  }
			  else if (kind_scale == Q_TRANSFORM) {
					//if (jrh_debug) cout << "iPoint = " << iPoint << " applying Q_TRANSFORM" << endl;
					f1 = (f1-min_f1)*invMaxmMinf1*0.998+0.001;
					f2 = (f2-min_f2)*invMaxmMinf2*0.998+0.001;
					f3 = (f3-min_f3)*invMaxmMinf3*0.998+0.001;
					f4 = (f4-min_f4)*invMaxmMinf4*0.998+0.001;
				  //interpoplate feature onto edf and apply logit function
				  //NOTE: Logit approximates the true inverse cdf function for Gaussian distribution, but is not exact! JRH 05132018


					//if (jrh_debug) cout << "f1 min_max scaled = " << f1 << endl;
					//FEATURE 1 INTERPOLATION
					iBin = 1;
					bin_max = dBin;
					while (bin_max < f1) {
						iBin++;
						bin_max += dBin;
					}
					iBin--; //Now at
					f1 = edf1[iBin]+(f1-(bin_max-dBin))*(edf1[iBin+1]-edf1[iBin])*inv_dBin;


					//if (jrh_debug) cout << "f2 min_max scaled = " << f2 << endl;
					//FEATURE 2 INTERPOLATION
					iBin = 1;
					bin_max = dBin;
					while (bin_max < f2) {
						iBin++;
						bin_max += dBin;
					}
					iBin--; //Now at
					f2 = edf2[iBin]+(f2-(bin_max-dBin))*(edf2[iBin+1]-edf2[iBin])*inv_dBin;

					//if (jrh_debug) cout << "f3 min_max scaled = " << f3 << endl;
					//FEATURE 3 INTERPOLATION
					iBin = 1;
					bin_max = dBin;
					while (bin_max < f3) {
						iBin++;
						bin_max += dBin;
					}
					iBin--; //Now at
					f3 = edf3[iBin]+(f3-(bin_max-dBin))*(edf3[iBin+1]-edf3[iBin])*inv_dBin;


					//if (jrh_debug) cout << "f4 min_max scaled = " << f4 << endl;
					//FEATURE 4 INTERPOLATION
					iBin = 1;
					bin_max = dBin;
					while (bin_max < f4) {
						iBin++;
						bin_max += dBin;
					}
					iBin--; //Now at
					f4 = edf4[iBin]+(f4-(bin_max-dBin))*(edf4[iBin+1]-edf4[iBin])*inv_dBin;


					f1 = log(f1/(1.0-f1));
					f2 = log(f2/(1.0-f2));
					f3 = log(f3/(1.0-f3));
					f4 = log(f4/(1.0-f4));

					if (jrh_debug && rank==MASTER_NODE) cout << "Features Q-Scaled = " << f1 << " " << f2 << " " << f3 << " " << f4 << endl;

			  }
			  else if (kind_scale == MIN_MAX) {
					f1 = (f1-min_f1)*invMaxmMinf1*0.998+0.001;
					f2 = (f2-min_f2)*invMaxmMinf2*0.998+0.001;
					f3 = (f3-min_f3)*invMaxmMinf3*0.998+0.001;
					f4 = (f4-min_f4)*invMaxmMinf4*0.998+0.001;
			  }
			  else if (kind_scale == BOX_COX) {
				  if (jrh_debug) cout << "Lambdas = " << l1 <<" "<< l2<<" "<< l3<<" " << l4<<" " << " inv Lambdas = " <<inv_l1 <<" " << inv_l2<<" " << inv_l3<<" " << inv_l4<<" "<<"\n";
				  if (jrh_debug) cout << "Unscaled Features = " << f1<<" " <<f2<<" "<<f3<<" "<<f4<<"\n";
				  if (filter_shield) {
					  f1 = (pow(f1,l1)-1.0)*inv_l1;
					  f2 = (pow(f2,l2)-1.0)*inv_l2;
					  f3 = (pow(f3,l3)-1.0)*inv_l3;
					  f4 = (pow(f4,l4)-1.0)*inv_l4;
				  }
				  else {
					  f1 = (pow(f1+pow(10.0,-16.0),l1)-1.0)*inv_l1;
					  f2 = (pow(f2+pow(10.0,-16.0),l2)-1.0)*inv_l2;
					  f3 = (pow(f3+pow(10.0,-16.0),l3)-1.0)*inv_l3;
					  f4 = (pow(f4+pow(10.0,-16.0),l4)-1.0)*inv_l4;
				  }
				  if (jrh_debug) cout << "BOX_COX Scaled Features = " << f1<<" " <<f2<<" "<<f3<<" "<<f4<<"\n";

			  }
			  else if (kind_scale != NO_SCALE) {
				  cout << "JRH: ERROR, ONLY TWO IMPLEMENTED NN SCALINGS ARE Z_SCALE, MIN_MAX, AND Q_TRANSFORM" << endl;
			  }
			  //if (jrh_debug) cout << "iPoint " << iPoint << " Features " << f1 << "   " <<  f2 << "  " << f3 << endl;


			  su2double temp = 0.0;
			  ai[0][0] = 1.0;
			  ai[0][1] = f1;
			  ai[0][2] = f2;
			  ai[0][3] = f3;
			  ai[0][4] = f4;
			  inputs[0][0] = 1.0;
			  inputs[0][1] = f1;
			  inputs[0][2] = f2;
			  inputs[0][3] = f3;
			  inputs[0][4] = f4;
			  for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
				  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
					  ai[iLayer][iNode] = 0.0;
					  //Apply Weights to Inputs
					  //if (iNode>0) {
						  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
							  //if (jrh_debug) cout << " iLayer " << iLayer << " iNode " << iNode << " iInput " << iInput << endl;
							  ai[iLayer][iNode] += weights[iLayer][iInput][iNode]*inputs[iLayer-1][iInput];
						  }
					  //}
					  //Apply Activation Function

					  if (iLayer < nLayers-1){
						  //tanh Activation Function
						  inputs[iLayer][iNode] = (exp(2.0*ai[iLayer][iNode])-1.0)/(exp(2.0*ai[iLayer][iNode])+1.0);

						  //Sigmoid Activation Function
						  //inputs[iLayer][iNode] = 1.0/(1.0+exp(-ai[iLayer][iNode]));

						  //Leaky RELU Activation Function - https://towardsdatascience.com/activation-functions-neural-networks-1cbd9f8d91d6
						  //if (ai[iLayer][iNode] <= 0.0) inputs[iLayer][iNode] = 0.01*ai[iLayer][iNode];
						  //else inputs[iLayer][iNode] = ai[iLayer][iNode];

						  //TESTING LINEAR ACTIVATION FUNCTION for WEIGHTS AS DVs...CHANGE TO BACKPROP NOT IMPLEMENTED!!!!!
//						  inputs[iLayer][iNode] = ai[iLayer][iNode];
					  }
					  else inputs[iLayer][iNode] = ai[iLayer][iNode];
				  }
				  //BIAS NODES
				  if (iLayer < nLayers-1) inputs[iLayer][0] = 1.0; //Bias Node (no bias node on output layer (-1))...or hidden layer (-2)
			  }
	//		  inputs[nLayers-1][0] += 1.0;//Bias node in output layer
			  //Set predicted beta_fiml to beta_fiml in CVariable Class
//			  if (jrh_debug) cout << "Predicted beta = " << inputs[nLayers-1][0]+1 << endl;
			  if (config->GetKindTrainNN() == WEIGHTS) {
				  node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0); //JRH 07212018
				  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0);
			  }
			  else {
				  //JRH 01032018 - LIMITING FOR FIML-CLASSIC TESTING
				  //if (inputs[nLayers-1][0] > 1.0) inputs[nLayers-1][0] = 1.0;
				  //else if (inputs[nLayers-1][0] < -1.0) inputs[nLayers-1][0] = -1.0;
				  //JRH 01032018 - LIMITING FOR FIML-CLASSIC TESTING

				  node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0);
				  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(inputs[nLayers-1][0]+1.0);
				  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFimlTrain(node[iPoint]->GetBetaFimlTrain());
			  }
			  sse += pow(node[iPoint]->GetBetaFiml()-node[iPoint]->GetBetaFimlTrain(),2)*0.5;
	//		  solver_container[FLOW_SOL]->SetTotal_Loss(sse);
    	  } //if (node[iPoint]->GetDES_fd() > 0.94) {
    	  else {
    		  node[iPoint]->SetBetaFiml(1.0);
    		  solver_container[FLOW_SOL]->node[iPoint]->SetBetaFiml(1.0);
    	  }
    	  sse += pow(node[iPoint]->GetBetaFiml()-node[iPoint]->GetBetaFimlTrain(),2)*0.5;
      } //end of forward propagate npoint loop

      local_sse = sse*inv_nPointDomain;
//if (config->GetKindTrainNN() != WEIGHTS) { //No need to waste overhead computing exact sse if we're not using it (not doing backprop) JRH 07152018
#ifdef HAVE_MPI
      //Sum error across all procs to get total estimate JRH 05082018
      local_sse = sse*inv_nPointDomain/size;
//      sse = 0.0;
      SU2_MPI::Allreduce(&local_sse,&sse,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sse = local_sse;
#endif
//}
//else sse = local_sse;
      solver_container[FLOW_SOL]->SetTotal_Loss(sse);
      } //If statement to not run this forward prop if WEIGHTS and DISC_ADJ solver running
}

su2double CTurbSASolver::GetEp(unsigned short iLayer, unsigned short iInput, unsigned short iNode) {
	//cout << "JRH Debug: Ep = " << Ep[iLayer][iInput][iNode] << endl;
	return Ep[iLayer][iInput][iNode];
}

void CTurbSASolver::SetDES_LengthScale(CSolver **solver, CGeometry *geometry, CConfig *config){

//  unsigned short kindHybridRANSLES = config->GetKind_HybridRANSLES();
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short iDim = 0, jDim = 0, iNeigh = 0, nNeigh = 0;

//  su2double constDES = config->GetConst_DES();

  su2double density = 0.0, laminarViscosity = 0.0, kinematicViscosity = 0.0,
      eddyViscosity = 0.0, kinematicViscosityTurb = 0.0, wallDistance = 0.0, lengthScale = 0.0;

  su2double maxDelta = 0.0, deltaAux = 0.0, distDES = 0.0, uijuij = 0.0, k2 = 0.0, r_d = 0.0, f_d = 0.0,
      deltaDDES = 0.0, deltaAuxDDES = 0.0, omega = 0.0, ln_max = 0.0, ln[3] = {0.0, 0.0, 0.0},
      aux_ln = 0.0, f_kh = 0.0;

  su2double nu_hat, fw_star = 0.424, cv1_3 = pow(7.1, 3.0); k2 = pow(0.41, 2.0);
  su2double cb1   = 0.1355, ct3 = 1.2, ct4   = 0.5;
  su2double sigma = 2./3., cb2 = 0.622, f_max=1.0, f_min=0.1, a1=0.15, a2=0.3;
  su2double cw1 = 0.0, Ji = 0.0, Ji_2 = 0.0, Ji_3 = 0.0, fv1 = 0.0, fv2 = 0.0, ft2 = 0.0, psi_2 = 0.0;
  su2double *coord_i = NULL, *coord_j = NULL, **primVarGrad = NULL, *vorticity = NULL, delta[3] = {0.0,0.0,0.0},
      ratioOmega[3] = {0.0, 0.0, 0.0}, vortexTiltingMeasure = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    coord_i                 = geometry->node[iPoint]->GetCoord();
    nNeigh                  = geometry->node[iPoint]->GetnPoint();
    wallDistance            = geometry->node[iPoint]->GetWall_Distance();
    primVarGrad             = solver[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
    vorticity               = solver[FLOW_SOL]->node[iPoint]->GetVorticity();
    density                 = solver[FLOW_SOL]->node[iPoint]->GetDensity();
    laminarViscosity        = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    eddyViscosity           = solver[TURB_SOL]->node[iPoint]->GetmuT();
    kinematicViscosity      = laminarViscosity/density;
    kinematicViscosityTurb  = eddyViscosity/density;

    uijuij = 0.0;
    for(iDim = 0; iDim < nDim; iDim++){
      for(jDim = 0; jDim < nDim; jDim++){
        uijuij += primVarGrad[1+iDim][jDim]*primVarGrad[1+iDim][jDim];
      }
    }
    uijuij = sqrt(fabs(uijuij));
    uijuij = max(uijuij,1e-10);

    /*--- Low Reynolds number correction term ---*/

    nu_hat = node[iPoint]->GetSolution()[0];
    Ji   = nu_hat/kinematicViscosity;
    Ji_2 = Ji * Ji;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    cw1 = cb1/k2+(1.0+cb2)/sigma;

    psi_2 = (1.0 - (cb1/(cw1*k2*fw_star))*(ft2 + (1.0 - ft2)*fv2))/(fv1 * max(1.0e-10,1.0-ft2));
    psi_2 = min(100.0,psi_2);

//    switch(kindHybridRANSLES){
//      case SA_DES:
//        /*--- Original Detached Eddy Simulation (DES97)
//        Spalart
//        1997
//        ---*/
//        maxDelta=0.;
//        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
//          jPoint  = geometry->node[iPoint]->GetPoint(iNeigh);
//          coord_j = geometry->node[jPoint]->GetCoord();
//
//          deltaAux = 0.;
//          for (iDim = 0;iDim < nDim; iDim++){
//            deltaAux += pow((coord_j[iDim]-coord_i[iDim]),2.);
//          }
//
//          maxDelta = max(maxDelta,sqrt(deltaAux));
//        }
//
//        distDES         = constDES * maxDelta;
//        lengthScale = min(distDES,wallDistance);
//
//        break;
//
//      case SA_DDES:
        /*--- A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
         Spalart et al.
         Theoretical and Computational Fluid Dynamics - 2006
         ---*/

        maxDelta = 0.0;
        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
          jPoint  = geometry->node[iPoint]->GetPoint(iNeigh);
          coord_j = geometry->node[jPoint]->GetCoord();

          deltaAux = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            deltaAux += pow((coord_j[iDim]-coord_i[iDim]),2.);
          }

          maxDelta = max(maxDelta,sqrt(deltaAux));
        }

        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));

//        distDES = constDES * maxDelta;
//        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));

//        break;
//      case SA_ZDES:
//        /*--- Recent improvements in the Zonal Detached Eddy Simulation (ZDES) formulation.
//         Deck
//         Theoretical and Computational Fluid Dynamics - 2012
//         ---*/
//
//        deltaDDES = 0.0;
//        for (iNeigh = 0; iNeigh < nNeigh; iNeigh++){
//            jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
//            coord_j = geometry->node[jPoint]->GetCoord();
//            deltaAuxDDES = 0.0;
//            for ( iDim = 0; iDim < nDim; iDim++){
//              deltaAux       = abs(coord_j[iDim] - coord_i[iDim]);
//              delta[iDim]     = max(delta[iDim], deltaAux);
//              deltaAuxDDES += pow((coord_j[iDim]-coord_i[iDim]),2.);
//            }
//            deltaDDES = max(deltaDDES,sqrt(deltaAuxDDES));
//        }
//
//        omega = sqrt(vorticity[0]*vorticity[0] +
//                     vorticity[1]*vorticity[1] +
//                     vorticity[2]*vorticity[2]);
//
//        for (iDim = 0; iDim < 3; iDim++){
//          ratioOmega[iDim] = vorticity[iDim]/omega;
//        }
//
//        maxDelta = sqrt(pow(ratioOmega[0],2.0)*delta[1]*delta[2] +
//                        pow(ratioOmega[1],2.0)*delta[0]*delta[2] +
//                        pow(ratioOmega[2],2.0)*delta[0]*delta[1]);
//
//        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
//        f_d = 1.0-tanh(pow(8.0*r_d,3.0));
//
//        if (f_d < 0.99){
//          maxDelta = deltaDDES;
//        }
//
//        distDES = constDES * maxDelta;
//        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));
//
//        break;
//
//      case SA_EDDES:
//
//        /*--- An Enhanced Version of DES with Rapid Transition from RANS to LES in Separated Flows.
//         Shur et al.
//         Flow Turbulence Combust - 2015
//         ---*/
//
//        vortexTiltingMeasure = node[iPoint]->GetVortex_Tilting();
//
//        omega = sqrt(vorticity[0]*vorticity[0] +
//                     vorticity[1]*vorticity[1] +
//                     vorticity[2]*vorticity[2]);
//
//        for (iDim = 0; iDim < 3; iDim++){
//          ratioOmega[iDim] = vorticity[iDim]/omega;
//        }
//
//        ln_max = 0.0;
//        deltaDDES = 0.0;
//        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
//          jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
//          coord_j = geometry->node[jPoint]->GetCoord();
//          deltaAuxDDES = 0.0;
//          for (iDim = 0; iDim < nDim; iDim++){
//            delta[iDim] = fabs(coord_j[iDim] - coord_i[iDim]);
//            deltaAuxDDES += pow((coord_j[iDim]-coord_i[iDim]),2.);
//          }
//          deltaDDES=max(deltaDDES,sqrt(deltaAuxDDES));
//          ln[0] = delta[1]*ratioOmega[2] - delta[2]*ratioOmega[1];
//          ln[1] = delta[2]*ratioOmega[0] - delta[0]*ratioOmega[2];
//          ln[2] = delta[0]*ratioOmega[1] - delta[1]*ratioOmega[0];
//          aux_ln = sqrt(ln[0]*ln[0] + ln[1]*ln[1] + ln[2]*ln[2]);
//          ln_max = max(ln_max,aux_ln);
//          vortexTiltingMeasure += node[jPoint]->GetVortex_Tilting();
//        }
//
//        vortexTiltingMeasure = (vortexTiltingMeasure/fabs(nNeigh + 1.0));
//
//        f_kh = max(f_min, min(f_max, f_min + ((f_max - f_min)/(a2 - a1)) * (vortexTiltingMeasure - a1)));
//
//        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
//        f_d = 1.0-tanh(pow(8.0*r_d,3.0));
//
//        maxDelta = (ln_max/sqrt(3.0)) * f_kh;
//        if (f_d < 0.999){
//          maxDelta = deltaDDES;
//        }
//
//        distDES = constDES * maxDelta;
//        lengthScale=wallDistance-f_d*max(0.0,(wallDistance-distDES));
//
//        break;
//
//    }

//    node[iPoint]->SetDES_LengthScale(lengthScale); //COMMENTED OUT JRH 07182019 - COSTLY AND NO LONGER USED....
      node[iPoint]->SetDES_fd(f_d);
  }
}

//JRH 04302018
su2double CTurbSASolver::GetNNGradient(unsigned short iLayer, unsigned long iInput, unsigned long iNode) {
	//Called By Discrete Adjoint Solver in InitializeAdjoint()
	return SU2_TYPE::GetDerivative(weights[iLayer][iInput][iNode]);
}
//JRH 04302018
su2double CTurbSASolver::GetNNLossGradient() {
	//Called By Discrete Adjoint Solver in InitializeAdjoint()
	return SU2_TYPE::GetDerivative(sse);
}
void CTurbSASolver::SetAdjointNNSolution(unsigned short iLayer, unsigned long iInput, unsigned long iNode, su2double val_gradient) {
	//Called By Discrete Adjoint Solver in InitializeAdjoint()
	SU2_TYPE::SetDerivative(weights[iLayer][iInput][iNode], SU2_TYPE::GetValue(val_gradient));
}
void CTurbSASolver::SetAdjointNNLoss(su2double val_gradient) {
	//Called By Discrete Adjoint Solver in InitializeAdjoint()
	SU2_TYPE::SetDerivative(sse, SU2_TYPE::GetValue(val_gradient));
}
void CTurbSASolver::RegisterWeights(bool input) {
	if (input) {
		for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
					  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
							  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
								  AD::RegisterInput(weights[iLayer][iInput][iNode]);
							  }
					  }
		}
//		AD::RegisterInput(sse);
	}
	else {
		for (unsigned short iLayer = 1; iLayer < nLayers; iLayer++) {
						  for (unsigned long iNode = 0; iNode < num_nodes[iLayer]; iNode++) {
								  for (unsigned long iInput = 0; iInput < num_inputs[iLayer-1]; iInput++) {
									  AD::RegisterOutput(weights[iLayer][iInput][iNode]);
								  }
						  }
		}
//		AD::RegisterOutput(sse);
	}
}

void CTurbSASolver::SetWeight(unsigned short iLayer, unsigned long iInput, unsigned long iNode, su2double val_gradient) {
	weights[iLayer][iInput][iNode] = val_gradient;
}

su2double CTurbSASolver::GetWeight(unsigned short iLayer, unsigned long iInput, unsigned long iNode) {
	return weights[iLayer][iInput][iNode];
}

void CTurbSASolver::FindFeaturesMinMax() {
    su2double f1,f2,f3,f4,gam;
    su2double lmin_f1,lmin_f2,lmin_f3,lmin_f4,lmax_f1,lmax_f2,lmax_f3,lmax_f4;
	int rank = MASTER_NODE;
	int size;
	bool filter_point;

#ifdef HAVE_MPI
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	bool started = false;
	for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
			  filter_point = false;
			  if (filter_shield == true) {
				  f3 = node[iPoint]->GetDeltaCriterion();
				  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
				  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
			  }

			  if ((filter_point == false && filter_shield == true) || filter_shield == false) {//		if ((node[iPoint]->GetDES_fd() < 0.999 && filter_shield == true) || filter_shield == false) { //JRH 11232018
//		  f1 = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
//		  f2 = node[iPoint]->GetChiSA();
//		  f3 = node[iPoint]->GetDeltaCriterion();
//		  gam = node[iPoint]->GetGammaTrans();
		gam = 1.0;
//		  f1 = node[iPoint]->GetFwSA()*gam;
//		  f1 = node[iPoint]->GetDES_fd();
//		  f2 = node[iPoint]->GetChiSA();
//		  f3 = node[iPoint]->GetDeltaCriterion();
//		  f4 = node[iPoint]->GetChiSA();

//		f1 = node[iPoint]->GetDES_fd();
//		f2 = node[iPoint]->GetChiSA();
//		f3 = node[iPoint]->GetDeltaCriterion();
//		f4 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction());
		  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
  	  f2 = node[iPoint]->GetChiSA();
		  f3 = node[iPoint]->GetDeltaCriterion();
//		  f4 = node[iPoint]->GetFwSA();
		  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));




		  if (started == false) {
			  lmin_f1 = f1;
			  lmin_f2 = f2;
			  lmin_f3 = f3;
			  lmin_f4 = f4;

			  lmax_f1 = f1;
			  lmax_f2 = f2;
			  lmax_f3 = f3;
			  lmax_f4 = f4;

			  started = true;
		  }
		  if (f1>lmax_f1) lmax_f1 = f1;
		  if (f2>lmax_f2) lmax_f2 = f2;
		  if (f3>lmax_f3) lmax_f3 = f3;
		  if (f4>lmax_f4) lmax_f4 = f4;

		  if (f1<lmin_f1) lmin_f1 = f1;
		  if (f2<lmin_f2) lmin_f2 = f2;
		  if (f3<lmin_f3) lmin_f3 = f3;
		  if (f4<lmin_f4) lmin_f4 = f4;
		} //fd filtering
    } //npoint loop
#ifdef HAVE_MPI
	//MIN FEATURE 1
	SU2_MPI::Allreduce(&lmin_f1, &min_f1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	//MIN FEATURE 2
	SU2_MPI::Allreduce(&lmin_f2, &min_f2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	//MIN FEATURE 3
	SU2_MPI::Allreduce(&lmin_f3, &min_f3, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	//MIN FEATURE 4
	SU2_MPI::Allreduce(&lmin_f4, &min_f4, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);



	//MAX FEATURE 1
	SU2_MPI::Allreduce(&lmax_f1, &max_f1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	//MAX FEATURE 2
	SU2_MPI::Allreduce(&lmax_f2, &max_f2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	//MAX FEATURE 3
	SU2_MPI::Allreduce(&lmax_f3, &max_f3, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	//MAX FEATURE 4
	SU2_MPI::Allreduce(&lmax_f4, &max_f4, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else
	min_f1 = lmin_f1;
	min_f2 = lmin_f2;
	min_f3 = lmin_f3;
	min_f4 = lmin_f4;

	max_f1 = lmax_f1;
	max_f2 = lmax_f2;
	max_f3 = lmax_f3;
	max_f4 = lmax_f4;
#endif
}

void CTurbSASolver::ComputeEDF(unsigned short nBins) {
	//compute EDFs
	//nBins = 5;
	su2double dBin = 1.0/su2double(nBins);
	su2double f1, f2, f3, f4,gam;
	su2double invMaxmMinf1, invMaxmMinf2, invMaxmMinf3, invMaxmMinf4;
	su2double bin_max;
	bool filter_point;
	unsigned short iBin;
	int rank = MASTER_NODE;
	int size;

#ifdef HAVE_MPI
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

	invMaxmMinf1 = 1.0/(max_f1-min_f1);
	invMaxmMinf2 = 1.0/(max_f2-min_f2);
	invMaxmMinf3 = 1.0/(max_f3-min_f3);
	invMaxmMinf4 = 1.0/(max_f4-min_f4);
	if (jrh_debug) cout << "in ComputeEDF about to zero out local edf vectors" << endl;


	//Zero out edf
	for (iBin = 0; iBin < nBins+1; iBin++) {
		if (jrh_debug) cout << rank << " iBin=" << iBin << " nBins=" << nBins << endl;
		ledf1[iBin] = 0.0;
		ledf2[iBin] = 0.0;
		ledf3[iBin] = 0.0;
		ledf4[iBin] = 0.0;

//		if (jrh_debug) cout << rank << " done with ledf's" << endl;
//		edf1[iBin] = 0.0;
//		edf2[iBin] = 0.0;
//		edf3[iBin] = 0.0;
//		edf4[iBin] = 0.0;
	}
	if (jrh_debug) cout << "In ComputeEDF() About to start nPointDomain loop" << endl;
	for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
//		  gam = node[iPoint]->GetGammaTrans();
			  filter_point = false;
			  if (filter_shield == true) {
				  f3 = node[iPoint]->GetDeltaCriterion();
				  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));
				  if ( f3 < 0.25 || f4 > 2.0 ) filter_point = true; //True means don't this point JRH 07182019
			  }

			  if ((filter_point == false && filter_shield == true) || filter_shield == false) {//		if ((node[iPoint]->GetDES_fd() < 0.999 && filter_shield == true) || filter_shield == false) {
//		if (1) {
			gam = 1.0;
	//		  f1 = node[iPoint]->GetFwSA()*gam;
	//		  f1 = node[iPoint]->GetDES_fd();
	//		  f2 = node[iPoint]->GetChiSA();
	//		  f3 = node[iPoint]->GetDeltaCriterion();
	//		  f4 = node[iPoint]->GetChiSA();

	//		f1 = node[iPoint]->GetDES_fd();
	//		f2 = node[iPoint]->GetChiSA();
	//		f3 = node[iPoint]->GetDeltaCriterion();
	//		f4 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction());
			  f1 = node[iPoint]->GetProduction()/(node[iPoint]->GetDestruction()+pow(10.0,-16.0));
		  f2 = node[iPoint]->GetChiSA();
			  f3 = node[iPoint]->GetDeltaCriterion();
//			  f4 = node[iPoint]->GetFwSA();
			  f4 = node[iPoint]->GetStrainMagnitude()/(node[iPoint]->GetVorticityMagnitude()+pow(10.0,-16.0));



			f1 = (f1-min_f1)*invMaxmMinf1*0.998+0.001;
			f2 = (f2-min_f2)*invMaxmMinf2*0.998+0.001;
			f3 = (f3-min_f3)*invMaxmMinf3*0.998+0.001;
			f4 = (f4-min_f4)*invMaxmMinf4*0.998+0.001;
			bin_max = 0;
			for (unsigned short iBin = 0; iBin < nBins+1; iBin++) {
				if (f1 < bin_max) ledf1[iBin] += 1.0;
				if (f2 < bin_max) ledf2[iBin] += 1.0;
				if (f3 < bin_max) ledf3[iBin] += 1.0;
				if (f4 < bin_max) ledf4[iBin] += 1.0;
				//if (jrh_debug) cout << "iBin = " << iBin << " " << " bin_max = " << bin_max << endl;
				bin_max += dBin;

			} //if (node[iPoint]->GetDES_fd() < 0.94) {
		}
	}
	if (jrh_debug) cout << "In ComputeEDF() at start of MPI #ifdef statement" << endl;

#ifdef HAVE_MPI
//	SU2_MPI::Allreduce(Local_Sens_Beta_Fiml,  Total_Sens_Beta_Fiml,  nDV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (jrh_debug) cout << rank << " Inside MPI #ifdef in ComputeEDF()" << endl;
	SU2_MPI::Allreduce(ledf1, edf1, nBins+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(ledf2, edf2, nBins+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(ledf3, edf3, nBins+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(ledf4, edf4, nBins+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	for (unsigned short iBin = 0; iBin<nBins+1; iBin++) {
		edf1[iBin] = ledf1[iBin];
		edf2[iBin] = ledf2[iBin];
		edf3[iBin] = ledf3[iBin];
		edf4[iBin] = ledf4[iBin];
	}
#endif

	for (unsigned short iBin = 0; iBin < nBins+1; iBin++) {
		edf1[iBin] = edf1[iBin]/su2double(nTrainSamples);
		edf2[iBin] = edf2[iBin]/su2double(nTrainSamples);
		edf3[iBin] = edf3[iBin]/su2double(nTrainSamples);
		edf4[iBin] = edf4[iBin]/su2double(nTrainSamples);
		if (jrh_debug) cout << "iBin: " << iBin << " edf1: " << edf1[iBin] << " edf2: " << edf2[iBin]<< " edf3: " << edf3[iBin]<< " edf4: " << edf4[iBin] << endl;
	}
}

void CTurbSASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
  
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[0] = nu_tilde_Inf;
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      su2double Coord_Reflected[nDim];
      geometry->PointPointReflect(nDim, geometry->node[Point_Normal]->GetCoord(),
                                        geometry->node[iPoint]->GetCoord(), Coord_Reflected);
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), Coord_Reflected);


      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim	  
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      su2double Coord_Reflected[nDim];
      geometry->PointPointReflect(nDim, geometry->node[Point_Normal]->GetCoord(),
                                        geometry->node[iPoint]->GetCoord(), Coord_Reflected);
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), Coord_Reflected);


      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iDim;
  su2double *V_inflow, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_inflow = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/
      
      conv_numerics->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim. Is Point_Normal needed?
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_exhaust, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the infinity ---*/
      
      V_exhaust = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Engine;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim. Is Point_Normal needed?      
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics,
             config,  val_marker, true);
  
}

void CTurbSASolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics,
             config,  val_marker, false);
  
}

void CTurbSASolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker, bool inlet_surface) {
  
  unsigned long iPoint, iVertex, GlobalIndex_donor, GlobalIndex, iPoint_Normal;
  su2double *V_outlet, *V_inlet, *V_domain, *Normal, *UnitNormal, Area, Vn;
  bool ReverseFlow;
  unsigned short iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    iPoint_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    GlobalIndex_donor = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex != GlobalIndex_donor)) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/
      
      Vn = 0.0; ReverseFlow = false;
      for (iDim = 0; iDim < nDim; iDim++) {  Vn += V_domain[iDim+1]*UnitNormal[iDim]; }
      
      if ((inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
      if ((!inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }
      
      /*--- Do not anything if there is a
       reverse flow, Euler b.c. for the direct problem ---*/
      
      if (!ReverseFlow) {
        
        /*--- Allocate the value at the infinity ---*/
        
        if (inlet_surface) {
          V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_outlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_inlet);
        }
        else {
          V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_inlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_outlet);
        }
        
        /*--- Set the turb. variable solution
         set  the turbulent variables. Here we use a Neumann BC such
         that the turbulent variable is copied from the interior of the
         domain to the outlet before computing the residual.
         or set the turbulent variable states (prescribed for an inflow)  ----*/
        
        Solution_i[0] = node[iPoint]->GetSolution(0);
        
        //      if (inlet_surface) Solution_j[0] = 0.5*(node[iPoint]->GetSolution(0)+V_outlet [nDim+9]);
        //      else Solution_j[0] = 0.5*(node[iPoint]->GetSolution(0)+V_inlet [nDim+9]);
        
        //      /*--- Inflow analysis (interior extrapolation) ---*/
        //      if (((inlet_surface) && (!ReverseFlow)) || ((!inlet_surface) && (ReverseFlow))) {
        //        Solution_j[0] = 2.0*node[iPoint]->GetSolution(0) - node[iPoint_Normal]->GetSolution(0);
        //      }
        
        //      /*--- Outflow analysis ---*/
        //      else {
        //        if (inlet_surface) Solution_j[0] = Factor_nu_ActDisk*V_outlet [nDim+9];
        //        else { Solution_j[0] = Factor_nu_ActDisk*V_inlet [nDim+9]; }
        //      }
        
        /*--- Inflow analysis (interior extrapolation) ---*/
        if (((inlet_surface) && (!ReverseFlow)) || ((!inlet_surface) && (ReverseFlow))) {
          Solution_j[0] = node[iPoint]->GetSolution(0);
        }
        
        /*--- Outflow analysis ---*/
        else {
          Solution_j[0] = nu_tilde_ActDisk;
        }
        
        conv_numerics->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Grid Movement ---*/
        
        if (grid_movement)
          conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
        
        /*--- Compute the residual using an upwind scheme ---*/
        
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.AddBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
        /*--- Viscous contribution ---*/
        
        visc_numerics->SetNormal(Normal);
// mskim		
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());
        su2double Coord_Reflected[nDim];
        geometry->PointPointReflect(nDim, geometry->node[iPoint_Normal]->GetCoord(),
                                          geometry->node[iPoint]->GetCoord(), Coord_Reflected);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), Coord_Reflected);



        /*--- Conservative variables w/o reconstruction ---*/
        
        if (inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
        else visc_numerics->SetPrimitive(V_domain, V_outlet);
        
        /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
        
        visc_numerics->SetTurbVar(Solution_i, Solution_j);
        
        visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
        
        /*--- Compute residual, and Jacobians ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        
        /*--- Subtract residual, and update Jacobians ---*/
        
        //        LinSysRes.SubtractBlock(iPoint, Residual);
        //        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  delete[] UnitNormal;
  
}

void CTurbSASolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {
//
//  unsigned long iVertex, iPoint, jPoint;
//  unsigned short iVar, iDim;
//
//  su2double *Vector = new su2double[nDim];
//
//#ifndef HAVE_MPI
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
//
//      if (iPoint != jPoint) {
//
//        /*--- Store the solution for both points ---*/
//        for (iVar = 0; iVar < nVar; iVar++) {
//          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//          Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
//        }
//
//        /*--- Set Conservative Variables ---*/
//        numerics->SetTurbVar(Solution_i, Solution_j);
//
//        /*--- Retrieve flow solution for both points ---*/
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
//          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
//        }
//
//        /*--- Set Flow Variables ---*/
//        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
//
//        /*--- Set the normal vector ---*/
//        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
//        for (iDim = 0; iDim < nDim; iDim++)
//          Vector[iDim] = -Vector[iDim];
//        numerics->SetNormal(Vector);
//
//        /*--- Add Residuals and Jacobians ---*/
//        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.AddBlock(iPoint, Residual);
//        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
//      }
//    }
//  }
//
//#else
//
//  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
//  su2double *Conserv_Var, *Flow_Var;
//  bool compute;
//
//  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
//  su2double *Buffer_Send_U = new su2double [Buffer_Size];
//  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
//
//  /*--- Do the send process, by the moment we are sending each
//   node individually, this must be changed ---*/
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
//      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
//
//      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
//      else compute = true;
//
//      /*--- We only send the information that belong to other boundary ---*/
//      if ((jProcessor != rank) && compute) {
//
//        Conserv_Var = node[iPoint]->GetSolution();
//        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//
//        for (iVar = 0; iVar < nVar; iVar++)
//          Buffer_Send_U[iVar] = Conserv_Var[iVar];
//
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
//          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
//
//        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
//
//      }
//    }
//  }
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
//      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
//
//      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
//      else compute = true;
//
//      if (compute) {
//
//        /*--- We only receive the information that belong to other boundary ---*/
//        if (jProcessor != rank) {
//          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
//        }
//        else {
//
//          for (iVar = 0; iVar < nVar; iVar++)
//            Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
//
//          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
//            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
//
//        }
//
//        /*--- Store the solution for both points ---*/
//        for (iVar = 0; iVar < nVar; iVar++) {
//          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//          Solution_j[iVar] = Buffer_Receive_U[iVar];
//        }
//
//        /*--- Set Turbulent Variables ---*/
//        numerics->SetTurbVar(Solution_i, Solution_j);
//
//        /*--- Retrieve flow solution for both points ---*/
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
//          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
//        }
//
//        /*--- Set Flow Variables ---*/
//        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
//
//        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
//        for (iDim = 0; iDim < nDim; iDim++)
//          Vector[iDim] = -Vector[iDim];
//        numerics->SetNormal(Vector);
//
//        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.AddBlock(iPoint, Residual);
//        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
//      }
//    }
//  }
//
//  delete[] Buffer_Send_U;
//  delete[] Buffer_Receive_U;
//
//#endif
//
//  delete[] Vector;
//
}

void CTurbSASolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {
//
//  unsigned long iVertex, iPoint, jPoint;
//  unsigned short iVar, iDim;
//
//  su2double *Vector = new su2double[nDim];
//
//#ifndef HAVE_MPI
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
//
//      if (iPoint != jPoint) {
//
//        /*--- Store the solution for both points ---*/
//        for (iVar = 0; iVar < nVar; iVar++) {
//          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//          Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
//        }
//
//        /*--- Set Conservative Variables ---*/
//        numerics->SetTurbVar(Solution_i, Solution_j);
//
//        /*--- Retrieve flow solution for both points ---*/
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
//          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
//        }
//
//        /*--- Set Flow Variables ---*/
//        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
//
//        /*--- Set the normal vector ---*/
//        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
//        for (iDim = 0; iDim < nDim; iDim++)
//          Vector[iDim] = -Vector[iDim];
//        numerics->SetNormal(Vector);
//
//        /*--- Add Residuals and Jacobians ---*/
//        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.AddBlock(iPoint, Residual);
//        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
//      }
//    }
//  }
//
//#else
//
//  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
//  su2double *Conserv_Var, *Flow_Var;
//  bool compute;
//
//  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
//  su2double *Buffer_Send_U = new su2double [Buffer_Size];
//  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
//
//  /*--- Do the send process, by the moment we are sending each
//   node individually, this must be changed ---*/
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
//      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
//
//      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
//      else compute = true;
//
//      /*--- We only send the information that belong to other boundary ---*/
//      if ((jProcessor != rank) && compute) {
//
//        Conserv_Var = node[iPoint]->GetSolution();
//        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//
//        for (iVar = 0; iVar < nVar; iVar++)
//          Buffer_Send_U[iVar] = Conserv_Var[iVar];
//
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
//          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
//
//        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
//
//      }
//    }
//  }
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Find the associate pair to the original node ---*/
//      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
//      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
//
//      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
//      else compute = true;
//
//      if (compute) {
//
//        /*--- We only receive the information that belong to other boundary ---*/
//        if (jProcessor != rank) {
//          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
//        }
//        else {
//
//          for (iVar = 0; iVar < nVar; iVar++)
//            Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
//
//          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
//            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
//
//        }
//
//        /*--- Store the solution for both points ---*/
//        for (iVar = 0; iVar < nVar; iVar++) {
//          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//          Solution_j[iVar] = Buffer_Receive_U[iVar];
//        }
//
//        /*--- Set Turbulent Variables ---*/
//        numerics->SetTurbVar(Solution_i, Solution_j);
//
//        /*--- Retrieve flow solution for both points ---*/
//        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
//          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
//        }
//
//        /*--- Set Flow Variables ---*/
//        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
//
//        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
//        for (iDim = 0; iDim < nDim; iDim++)
//          Vector[iDim] = -Vector[iDim];
//        numerics->SetNormal(Vector);
//
//        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.AddBlock(iPoint, Residual);
//        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
//      }
//    }
//  }
//
//  delete[] Buffer_Send_U;
//  delete[] Buffer_Receive_U;
//
//#endif
//
//  delete[] Vector;
//
}

CTurbSSTSolver::CTurbSSTSolver(void) : CTurbSolver() {
  
  /*--- Array initialization ---*/
  constants = NULL;
  
}

CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint, index;
  su2double dull_val;
  ifstream restart_file;
  string text_line;
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Array initialization ---*/
  
  constants = NULL;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  
  nVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  
  if (iMesh == MESH_0) {
    
    /*--- Define some auxiliary vector related with the residual ---*/
    
    Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    
    /*--- Define some structures for locating max residuals ---*/
    
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
    
    /*--- Define some auxiliary vector related with the solution ---*/
    
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    
    /*--- Define some auxiliary vector related with the geometry ---*/
    
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
    
    /*--- Define some auxiliary vector related with the flow solution ---*/
    
    FlowPrimVar_i = new su2double [nDim+7]; FlowPrimVar_j = new su2double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize value for model constants ---*/
  constants = new su2double[10];
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  
  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e10;
  
  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;
  
  /*--- Flow infinity initialization stuff ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;
  
  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();
  
  su2double VelMag = 0;
  for (iDim = 0; iDim < nDim; iDim++)
  VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);
  
  kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
  omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);
  
  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;
  
  /*--- Restart the solution from file information ---*/
  if (!restart || (iMesh != MESH_0)) {
    for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nDim, nVar, constants, config);
  }
  else {
    
    /*--- Restart the solution from file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    
    /*--- Modify file name for multizone problems ---*/
    if (nZone >1)
      filename= config->GetMultizone_FileName(filename, iZone);

    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time || time_stepping) {
      int Unst_RestartIter;
      if (adjoint) {
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no turbulent restart file!!" << endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/

    map<unsigned long,unsigned long> Global2Local;
    map<unsigned long,unsigned long>::const_iterator MI;
    
    /*--- Now fill array with the transform values only for local points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }
    
    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; string text_line; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
        
        if (compressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
        }
        if (incompressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
        }
        
        /*--- Instantiate the solution at this node, note that the muT_Inf should recomputed ---*/
        node[iPoint_Local] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
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
      node[iPoint] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
    }
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
  }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CTurbSSTSolver::~CTurbSSTSolver(void) {
  
  if (constants != NULL) delete [] constants;
  
}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;

  unsigned long ExtIter      = config->GetExtIter();
  bool limiter_flow          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  su2double rho = 0.0, mu = 0.0, dist, omega, kine, strMag, F2, muT, zeta;
  su2double a1 = constants[7];
  unsigned long iPoint;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Compute mean flow and turbulence gradients ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_GG(geometry, config);
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_LS(geometry, config);
    SetSolution_Gradient_LS(geometry, config);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Compute blending functions and cross diffusion ---*/
    
    if (compressible) {
      rho  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu   = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    if (incompressible) {
      rho  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu   = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    
    dist = geometry->node[iPoint]->GetWall_Distance();
    
    strMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
// mskim. PR#905. S -> W
    su2double *Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    su2double VorticityMag = sqrt(Vorticity[0]*Vorticity[0] +
                                  Vorticity[1]*Vorticity[1] +
                                  Vorticity[2]*Vorticity[2]);


    node[iPoint]->SetBlendingFunc(mu, dist, rho);
    
    F2 = node[iPoint]->GetF2blending();
    
    /*--- Compute the eddy viscosity ---*/
    
    kine  = node[iPoint]->GetSolution(0);
    omega = node[iPoint]->GetSolution(1);
// mskim. PR#905. S -> W
//    zeta = min(1.0/omega, a1/(strMag*F2));
    zeta = min(1.0/omega, a1/(VorticityMag*F2));
    muT = min(max(rho*kine*zeta,0.0),1.0);
    node[iPoint]->SetmuT(muT);
    
  }
  
}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {

// mskim
  bool axisymmetric = config->GetAxisymmetric();

  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set distance to the surface ---*/
    
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Menter's first blending function ---*/
    
    numerics->SetF1blending(node[iPoint]->GetF1blending(),0.0);
    
    /*--- Menter's second blending function ---*/
    
    numerics->SetF2blending(node[iPoint]->GetF2blending(),0.0);
    
    /*--- Set vorticity and strain rate magnitude ---*/
    
    numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);
    
    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);
    
    /*--- Cross diffusion ---*/
    
    numerics->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);

// mskim
    if (axisymmetric){
    /*--- Set y coordinate ---*/
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
    }

    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
  
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      if (compressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      if (incompressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
  
}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      if (compressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      if (incompressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
  
}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0] = kine_Inf;
      Solution_j[1] = omega_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0]= kine_Inf;
      Solution_j[1]= omega_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim	  
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      su2double Coord_Reflected[nDim];
      geometry->PointPointReflect(nDim, geometry->node[Point_Normal]->GetCoord(),
                                        geometry->node[iPoint]->GetCoord(), Coord_Reflected);
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), Coord_Reflected);


      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
// mskim	  
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      su2double Coord_Reflected[nDim];
      geometry->PointPointReflect(nDim, geometry->node[Point_Normal]->GetCoord(),
                                        geometry->node[iPoint]->GetCoord(), Coord_Reflected);
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), Coord_Reflected);


      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

su2double* CTurbSSTSolver::GetConstants() {
  return constants;
}
