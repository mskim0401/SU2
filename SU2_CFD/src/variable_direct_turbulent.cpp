/*!
 * \file variable_direct_turbulent.cpp
 * \brief Definition of the solution fields.
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

#include "../include/variable_structure.hpp"

CTurbVariable::CTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  HB_Source = NULL;
  
}

CTurbVariable::CTurbVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
  
  unsigned short iVar;

  /*--- Array initialization ---*/

  HB_Source = NULL;

  /*--- Allocate space for the harmonic balance source terms ---*/

  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    HB_Source = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      HB_Source[iVar] = 0.0;
  }

  /*--- Initialize beta if FIML corrections required ---*/ //JRH - 03312017

  if (config->GetKind_Turb_Model() == SA_FIML) {
	  beta_fiml = 1.0; //JRH - 03312017
  }


  /*--- Allocate space for the limiter ---*/

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nVar];
  Solution_Min = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

}

//Alias of CTurbVariable that includes additional input to initialize beta from FIML DVs set in CConfig - JRH 04122017
CTurbVariable::CTurbVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_iPoint, CConfig *config)
: CVariable(val_nDim, val_nvar,config) {

  unsigned short iVar;
  //unsigned long nDV; //JRH 04122017
  //unsigned long iDV; //JRH 04122017
  /*--- Array initialization ---*/
  
  HB_Source = NULL;
  
  /*--- Allocate space for the harmonic balance source terms ---*/
  
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    HB_Source = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      HB_Source[iVar] = 0.0;
  }
  
  /*--- Initialize beta if FIML corrections required ---*/ //JRH - 03312017

  if (config->GetKind_Turb_Model() == SA_FIML) {
	  //beta_fiml = 1.0; //JRH - 03312017
	  //Initialize beta to
	  //if (config->GetDesign_Variable(0) != NULL) { //JRH 04122017
	  	  beta_fiml_grad = 0.0; //Default JRH 05032017
		  //iDV = 0;
		  //JRH - Find the first occurrence of FIML DV (if any)
		  //nDV = config->GetnDV();
		  //while (config->GetDesign_Variable(iDV) != FIML && iDV < nDV) iDV++;
		  beta_fiml = 1.0;
		  fd = 0.0;
		  Production = 0.0;
		  Destruction = 0.0;
		  STildeSA = 0.0;
		  ChiSA = 0.0;
		  Delta_Criterion = 0.0;
		  FwSA = 0.0;
		  RSA = 0.0;
		  Strain_Magnitude = 0.0;
		  Vorticity_Magnitude = 0.0;
		  k_SALSA = 0.0;
	  //}
  }


  /*--- Allocate space for the limiter ---*/
  
  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nVar];
  Solution_Min = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
}

CTurbVariable::~CTurbVariable(void) { }

su2double CTurbVariable::GetmuT() { return muT; }

//JRH - Added routine to retrieve beta FIML correction factor - 04012017
su2double CTurbVariable::GetBetaFiml() { return beta_fiml;}

su2double CTurbVariable::GetDES_fd() {return fd;}

void CTurbVariable::SetDES_fd(su2double val_fd) {fd = val_fd;}

su2double CTurbVariable::GetBetaFimlTrain() { return beta_fiml_train; }
void CTurbVariable::SetBetaFimlTrain(su2double val_beta_fiml_train) { beta_fiml_train = val_beta_fiml_train; }

void CTurbVariable::RegisterBeta(bool input) {
	if (input) {
		AD::RegisterInput(beta_fiml); //JRH 05022018
	}
	else {
		//AD::RegisterOutput(beta_fiml);
	}
}
void CTurbVariable::SetAdjointBeta(su2double val_adjoint_beta) {
	SU2_TYPE::SetDerivative(beta_fiml,SU2_TYPE::GetValue(val_adjoint_beta));//JRH 05022018
}
su2double CTurbVariable::GetAdjointBeta(void) {
	return SU2_TYPE::GetDerivative(beta_fiml);//JRH 05022018
}

//JRH - Added routine to set beta FIML correction factor - Called after registering DV as input in solver_adjoint_discrete 04242017
void CTurbVariable::SetBetaFiml(su2double val_beta_fiml) { beta_fiml = val_beta_fiml; }

//JRH - Added routine to retrieve beta FIML correction factor gradient - 05032017
su2double CTurbVariable::GetBetaFimlGrad() { return beta_fiml_grad;}

//JRH - Added routine to set beta FIML correction factor gradient - Called after registering DV as input in solver_adjoint_discrete 05032017
void CTurbVariable::SetBetaFimlGrad(su2double val_beta_fiml_grad) { beta_fiml_grad = val_beta_fiml_grad; }

//JRH - Added routines to get ML features from turbulence models
su2double CTurbVariable::GetProduction() {return Production;}
su2double CTurbVariable::GetDestruction() {return Destruction;}
su2double CTurbVariable::GetSTildeSA() {return STildeSA;}
su2double CTurbVariable::GetChiSA() {return ChiSA;}
su2double CTurbVariable::GetDeltaCriterion() {return Delta_Criterion;}
su2double CTurbVariable::GetFwSA() {return FwSA;}
su2double CTurbVariable::GetRSA() {return Production;}
su2double CTurbVariable::GetStrainMagnitude() {return Strain_Magnitude;}
su2double CTurbVariable::GetVorticityMagnitude() {return Vorticity_Magnitude;}
su2double CTurbVariable::GetGammaTrans() {return gamma_trans;}
su2double CTurbVariable::GetWallDist() {return wall_dist;}
su2double CTurbVariable::GetkSALSA() {return k_SALSA;}

void CTurbVariable::SetProduction(su2double val_Production) {Production = val_Production;}
void CTurbVariable::SetDestruction(su2double val_Destruction) {Destruction = val_Destruction;}
void CTurbVariable::SetSTildeSA(su2double val_STildeSA) {STildeSA = val_STildeSA;}
void CTurbVariable::SetChiSA(su2double val_ChiSA) {ChiSA = val_ChiSA;}
void CTurbVariable::SetDeltaCriterion(su2double val_Delta_Criterion) {Delta_Criterion = val_Delta_Criterion;}
void CTurbVariable::SetFwSA(su2double val_FwSA) {FwSA = val_FwSA;}
void CTurbVariable::SetRSA(su2double val_RSA) {RSA = val_RSA;}
void CTurbVariable::SetStrainMagnitude(su2double val_StrainMag_i) {Strain_Magnitude = val_StrainMag_i;}
void CTurbVariable::SetVorticityMagnitude(su2double val_Omega) {Vorticity_Magnitude = val_Omega;}
void CTurbVariable::SetGammaTrans(su2double val_Gamma_Trans) {gamma_trans = val_Gamma_Trans;}
void CTurbVariable::SetWallDist(su2double val_Wall_Dist) {wall_dist = val_Wall_Dist;}
void CTurbVariable::SetkSALSA(su2double val_k_SALSA) {k_SALSA = val_k_SALSA;}

void CTurbVariable::SetmuT(su2double val_muT) { muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Initialization of S-A variables ---*/
  Solution[0] = val_nu_tilde;    Solution_Old[0] = val_nu_tilde;
  
  /*--- Initialization of the eddy viscosity ---*/
  muT = val_muT;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_nu_tilde;
    Solution_time_n1[0] = val_nu_tilde;
  }

}

//JRH - Overload to inherit iPoint from new alias of CTurbVariable
CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar,unsigned short val_iPoint, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, val_iPoint, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of S-A variables ---*/
  Solution[0] = val_nu_tilde;    Solution_Old[0] = val_nu_tilde;

  /*--- Initialization of the eddy viscosity ---*/
  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_nu_tilde;
    Solution_time_n1[0] = val_nu_tilde;
  }

}

CTurbSAVariable::~CTurbSAVariable(void) {
  
  if (HB_Source != NULL) delete [] HB_Source;
  
}

CTurbMLVariable::CTurbMLVariable(void) : CTurbVariable() { }

CTurbMLVariable::CTurbMLVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Initialization of S-A variables ---*/
  Solution[0] = val_nu_tilde;    Solution_Old[0] = val_nu_tilde;
  
  /*--- Initialization of the eddy viscosity ---*/
  muT = val_muT;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_nu_tilde;
    Solution_time_n1[0] = val_nu_tilde;
  }
  
}

CTurbMLVariable::~CTurbMLVariable(void) {
  
  if (HB_Source != NULL) delete [] HB_Source;
  
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

CTurbSSTVariable::CTurbSSTVariable(su2double val_kine, su2double val_omega, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar,
                                   su2double *constants, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Initialization of variables ---*/
  
  Solution[0] = val_kine;     Solution_Old[0] = val_kine;
  Solution[1] = val_omega;  Solution_Old[1] = val_omega;
  
  sigma_om2 = constants[3];
  beta_star = constants[6];
  
  F1   = 1.0;
  F2   = 0.0;
  CDkw = 0.0;
  
  /*--- Initialization of eddy viscosity ---*/
  
  muT = val_muT;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
  
  if (dual_time) {
    Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_omega;
    Solution_time_n1[0]  = val_kine; Solution_time_n1[1]  = val_omega;
  }
    
}

CTurbSSTVariable::~CTurbSSTVariable(void) {

  if (HB_Source != NULL) delete [] HB_Source;
  
}

void CTurbSSTVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) {
  unsigned short iDim;
  su2double arg2, arg2A, arg2B, arg1;
  
  /*--- Cross diffusion ---*/
  
  CDkw = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    CDkw += Gradient[0][iDim]*Gradient[1][iDim];
  CDkw *= 2.0*val_density*sigma_om2/Solution[1];
  CDkw = max(CDkw, pow(10.0, -20.0));
  
  /*--- F1 ---*/
  
  arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist+EPS*EPS);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution[1]+EPS*EPS);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution[0] / (CDkw*val_dist*val_dist+EPS*EPS));
  F1 = tanh(pow(arg1, 4.0));
  
  /*--- F2 ---*/
  
  arg2 = max(2.0*arg2A, arg2B);
  F2 = tanh(pow(arg2, 2.0));
  
}
