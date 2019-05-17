/*!
 * \file output_elasticity.cpp
 * \brief Main subroutines for FEA output
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/output/output_elasticity.hpp"

CElasticityOutput::CElasticityOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  dynamic = (config->GetTime_Domain() || (config->GetDynamic_Analysis() == DYNAMIC));  // Dynamic analysis.
  
  nDim = geometry->GetnDim();

  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;

  /*--- Default fields for screen output ---*/
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  /*--- Default fields for screen output ---*/
  if (nRequestedScreenFields == 0){
    if (dynamic) RequestedScreenFields.push_back("TIME_ITER");
    if (multizone) RequestedScreenFields.push_back("OUTER_ITER");
    RequestedScreenFields.push_back("INNER_ITER");
    if(linear_analysis){
      RequestedScreenFields.push_back("RMS_DISP_X");
      RequestedScreenFields.push_back("RMS_DISP_Y");
      RequestedScreenFields.push_back("RMS_DISP_Z");
    }
    if(nonlinear_analysis){
      RequestedScreenFields.push_back("RMS_UTOL");
      RequestedScreenFields.push_back("RMS_RTOL");
      RequestedScreenFields.push_back("RMS_ETOL");
    }
    RequestedScreenFields.push_back("VMS");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  /*--- Default fields for volume output ---*/
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("SOLUTION");
    RequestedVolumeFields.push_back("STRESS");    
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Structure)";
  MultiZoneHeaderString = ss.str();
  
  /*--- Set the volume filename --- */
  
  VolumeFilename = config->GetVolume_FileName();
  
  /*--- Set the surface filename --- */
  
  SurfaceFilename = config->GetSurfCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  RestartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (Conv_Field.size() == 0 ) Conv_Field = "RMS_DISP_X";

}

CElasticityOutput::~CElasticityOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();

  }

}

void CElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* fea_solver = solver[FEA_SOL];
  
  SetHistoryOutputValue("TIME_ITER",  curr_TimeIter);  
  SetHistoryOutputValue("INNER_ITER", curr_InnerIter);
  SetHistoryOutputValue("OUTER_ITER", curr_OuterIter); 

  /*--- Residuals: ---*/
  /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
  /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

  
  if (linear_analysis){
    SetHistoryOutputValue("RMS_DISP_X", log10(fea_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISP_Y", log10(fea_solver->GetRes_RMS(1)));
    if (nDim == 3){
      SetHistoryOutputValue("RMS_DISP_Z", log10(fea_solver->GetRes_RMS(2)));
    }
  } else if (nonlinear_analysis){
    SetHistoryOutputValue("RMS_UTOL", log10(fea_solver->GetRes_FEM(0)));
    SetHistoryOutputValue("RMS_RTOL", log10(fea_solver->GetRes_FEM(1)));
    if (nDim == 3){
      SetHistoryOutputValue("RMS_ETOL", log10(fea_solver->GetRes_FEM(2)));
    }
  }
  
  if (multizone){
    SetHistoryOutputValue("BGS_DISP_X", log10(fea_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_DISP_Y", log10(fea_solver->GetRes_BGS(1)));
    if (nDim == 3) SetHistoryOutputValue("BGS_DISP_Z", log10(fea_solver->GetRes_BGS(2)));
  }
  
  SetHistoryOutputValue("VMS", fea_solver->GetTotal_CFEA());
  SetHistoryOutputValue("LOAD_INCREMENT", fea_solver->GetLoad_Increment());
  SetHistoryOutputValue("LOAD_RAMP", fea_solver->GetForceCoeff());
  
} 

void CElasticityOutput::SetHistoryOutputFields(CConfig *config){
  
  // Iteration numbers
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER"); 
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER");  
  AddHistoryOutput("INNER_ITER",   "Inner_Iter",  FORMAT_INTEGER, "ITER");

  // Misc.
  AddHistoryOutput("PHYS_TIME",   "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  // Residuals

  AddHistoryOutput("RMS_UTOL",   "rms[U]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("RMS_RTOL",   "rms[R]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("RMS_ETOL",   "rms[E]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);

  AddHistoryOutput("RMS_DISP_X", "rms[DispX]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("RMS_DISP_Y", "rms[DispY]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("RMS_DISP_Z", "rms[DispZ]", FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);

  AddHistoryOutput("BGS_DISP_X", "bgs[DispX]", FORMAT_FIXED,  "BGS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("BGS_DISP_Y", "bgs[DispY]", FORMAT_FIXED,  "BGS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("BGS_DISP_Z", "bgs[DispZ]", FORMAT_FIXED,  "BGS_RES", TYPE_RESIDUAL);

  AddHistoryOutput("VMS",            "VonMises", FORMAT_SCIENTIFIC, "VMS");
  AddHistoryOutput("LOAD_INCREMENT", "Load_Increment",  FORMAT_FIXED, "LOAD_INCREMENT");
  AddHistoryOutput("LOAD_RAMP",      "Load_Ramp",       FORMAT_FIXED, "LOAD_RAMP");
  
}

void CElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
 
  CVariable* Node_Struc = solver[FEA_SOL]->node[iPoint]; 
  CPoint*    Node_Geo  = geometry->node[iPoint];
          
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("DISPLACEMENT-X", iPoint, Node_Struc->GetSolution(0));
  SetVolumeOutputValue("DISPLACEMENT-Y", iPoint, Node_Struc->GetSolution(1));
  if (nDim == 3) SetVolumeOutputValue("DISPLACEMENT-Z", iPoint, Node_Struc->GetSolution(2));
  
  if(dynamic){
    SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Struc->GetSolution_Vel(0));
    SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Struc->GetSolution_Vel(1));
    if (nDim == 3) SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Struc->GetSolution_Vel(2));
  
    SetVolumeOutputValue("ACCELERATION-X", iPoint, Node_Struc->GetSolution_Accel(0));
    SetVolumeOutputValue("ACCELERATION-Y", iPoint, Node_Struc->GetSolution_Accel(1));
    if (nDim == 3) SetVolumeOutputValue("ACCELERATION-Z", iPoint, Node_Struc->GetSolution_Accel(2));
  }
  
  SetVolumeOutputValue("STRESS-XX", iPoint, Node_Struc->GetStress_FEM()[0]);
  SetVolumeOutputValue("STRESS-YY", iPoint, Node_Struc->GetStress_FEM()[1]);
  SetVolumeOutputValue("STRESS-XY", iPoint, Node_Struc->GetStress_FEM()[2]);
  if (nDim == 3){
    SetVolumeOutputValue("STRESS-ZZ", iPoint, Node_Struc->GetStress_FEM()[3]);
    SetVolumeOutputValue("STRESS-XZ", iPoint, Node_Struc->GetStress_FEM()[4]);
    SetVolumeOutputValue("STRESS-YZ", iPoint, Node_Struc->GetStress_FEM()[5]);
  }
  SetVolumeOutputValue("VON_MISES_STRESS", iPoint, Node_Struc->GetVonMises_Stress());
  
}

void CElasticityOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");

  AddVolumeOutput("DISPLACEMENT-X",    "Displacement_x", "SOLUTION");
  AddVolumeOutput("DISPLACEMENT-Y",    "Displacement_y", "SOLUTION");
  if (nDim == 3) AddVolumeOutput("DISPLACEMENT-Z", "Displacement_z", "SOLUTION");
  
  AddVolumeOutput("VELOCITY-X",    "Velocity_x", "VELOCITY");
  AddVolumeOutput("VELOCITY-Y",    "Velocity_y", "VELOCITY");
  if (nDim == 3) AddVolumeOutput("VELOCITY-Z", "Velocity_z", "VELOCITY");
  
  AddVolumeOutput("ACCELERATION-X",    "Acceleration_x", "ACCELERATION");
  AddVolumeOutput("ACCELERATION-Y",    "Acceleration_y", "ACCELERATION");
  if (nDim == 3) AddVolumeOutput("ACCELERATION-Z", "Acceleration_z", "ACCELERATION");
  
  AddVolumeOutput("STRESS-XX",    "Sxx", "STRESS");
  AddVolumeOutput("STRESS-YY",    "Syy", "STRESS");
  AddVolumeOutput("STRESS-XY",    "Sxy", "STRESS");
  
  if (nDim == 3) {
    AddVolumeOutput("STRESS-ZZ",    "Szz", "STRESS");
    AddVolumeOutput("STRESS-XZ",    "Sxz", "STRESS");
    AddVolumeOutput("STRESS-YZ",    "Syz", "STRESS");
  }
    
  AddVolumeOutput("VON_MISES_STRESS", "Von_Mises_Stress", "STRESS");
  
}

inline bool CElasticityOutput::WriteHistoryFile_Output(CConfig *config) { return true;}

inline bool CElasticityOutput::WriteScreen_Header(CConfig *config) {  
  
  bool write_header;
  if (nonlinear_analysis) write_header = (config->GetIntIter() == 0);
  else write_header = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));

  /*--- For multizone problems, print the header only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_header = (write_header && config->GetWrt_ZoneConv());

  return write_header;
  
}
inline bool CElasticityOutput::WriteScreen_Output(CConfig *config) {

  bool write_output = true;

  /*--- For multizone problems, print the header only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_output = (write_output && config->GetWrt_ZoneConv());

  return write_output;
}

