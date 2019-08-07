/*!
 * \file SU2_DOT.cpp
 * \brief Main file of the Gradient Projection Code (SU2_DOT).
 * \author F. Palacios, T. Economon
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

#include "../include/SU2_DOT.hpp"
using namespace std;
su2double *beta_fiml_grad;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  
  char config_file_name[MAX_STRING_SIZE], *cstr;
  ofstream Gradient_file;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  /*--- MPI initialization, and buffer setting ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Init(&argc,&argv);
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  MPI_Comm_rank(MPICommunicator,&rank);
  MPI_Comm_size(MPICommunicator,&size);
#else
  SU2_Comm MPICommunicator(0);
#endif
  
  /*--- Pointer to different structures that will be used throughout the entire code ---*/
  
  CConfig **config_container          = NULL;
  CGeometry **geometry_container      = NULL;
  CSurfaceMovement *surface_movement  = NULL;
  CVolumetricMovement *mesh_movement  = NULL;
  //CSolver **solver_container		  = NULL; //JRH 06202017
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }
  
  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  //solver_container = new CSolver*[nZone]; //JRH 06202017
  
  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    //solver_container[iZone]		  = NULL; //JRH 06202017
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_DOT, iZone, nZone, 0, VERB_HIGH);

    /*--- Set the MPI communicator ---*/
    config_container[iZone]->SetMPICommunicator(MPICommunicator);
        
    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
    
    CGeometry *geometry_aux = NULL;
    
    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
    
    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
    
    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
    
    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    
    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);
    
    /*--- Deallocate the memory of geometry_aux ---*/
    
    delete geometry_aux;
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetSendReceive(config_container[iZone]);
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetBoundaries(config_container[iZone]);
    
  }
  
  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/
  
#ifdef HAVE_MPI
  StartTime = MPI_Wtime();
#else
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
  geometry_container[ZONE_0]->SetPoint_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl;
  geometry_container[ZONE_0]->SetBoundVolume();
  geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
  geometry_container[ZONE_0]->Check_BoundElem_Orientation(config_container[ZONE_0]);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
  geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry_container[ZONE_0]->SetCoord_CG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);
  
  /*--- Load the surface sensitivities from file. This is done only
   once: if this is an unsteady problem, a time-average of the surface
   sensitivities at each node is taken within this routine. ---*/
  
  if (!config_container[ZONE_0]->GetDiscrete_Adjoint()){
    if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl;
    geometry_container[ZONE_0]->SetBoundSensitivity(config_container[ZONE_0]);
  }
  else {
    if (rank == MASTER_NODE) cout << "Reading volume sensitivities at each node from file." << endl;
    mesh_movement = new CVolumetricMovement(geometry_container[ZONE_0], config_container[ZONE_0]);
    geometry_container[ZONE_0]->SetSensitivity(config_container[ZONE_0]);
    
    if (rank == MASTER_NODE)
      cout << endl <<"---------------------- Mesh sensitivity computation ---------------------" << endl;
    mesh_movement->SetVolume_Deformation(geometry_container[ZONE_0], config_container[ZONE_0], false, true);
    
    COutput *output = new COutput();
    output->SetSensitivity_Files(geometry_container, config_container, nZone);
  }
  
  /*--- Definition of the Class for surface deformation ---*/
  
  surface_movement = new CSurfaceMovement();
  
  /*--- Copy coordinates to the surface structure ---*/
  
  surface_movement->CopyBoundary(geometry_container[ZONE_0], config_container[ZONE_0]);
  
  if (config_container[ZONE_0]->GetDesign_Variable(0) != NONE){
    if (rank == MASTER_NODE)
      cout << endl <<"---------- Start gradient evaluation using sensitivity information ----------" << endl;
    
    /*--- Write the gradient in a external file ---*/
    
    if (rank == MASTER_NODE) {
      cstr = new char [config_container[ZONE_0]->GetObjFunc_Grad_FileName().size()+1];
      strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str());
      Gradient_file.open(cstr, ios::out);
    }
    
    /*--- If AD mode is enabled we can use it to compute the projection,
     * otherwise we use finite differences. ---*/
    
    if (config_container[ZONE_0]->GetAD_Mode()){
      SetProjection_AD(geometry_container[ZONE_0], config_container[ZONE_0], surface_movement, Gradient_file);
    }else{
      SetProjection_FD(geometry_container[ZONE_0], config_container[ZONE_0], surface_movement, Gradient_file);
    }
    
    if (rank == MASTER_NODE)
      Gradient_file.close();
  }
  
  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifdef HAVE_MPI
  StopTime = MPI_Wtime();
#else
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_DOT) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}

void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, ofstream& Gradient_file){
  
  unsigned short iDV, nDV, iFFDBox, nDV_Value, iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double delta_eps, my_Gradient, **Gradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
  bool *UpdatePoint, MoveSurface, Local_MoveSurface;
  CFreeFormDefBox **FFDBox;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  cout << "JRH Debugging: In SetProjection_FD" << endl;
  nDV = config->GetnDV();
  
  /*--- Boolean controlling points to be updated ---*/
  
  UpdatePoint = new bool[geometry->GetnPoint()];
  
  /*--- Definition of the FFD deformation class ---*/
  
  unsigned short nFFDBox = MAX_NUMBER_FFD;
  FFDBox = new CFreeFormDefBox*[nFFDBox];
  
  /*--- Structure to store the gradient ---*/
  
  Gradient = new su2double*[nDV];
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (nDV_Value != 1 && config->GetDesign_Variable(0) != FIML){
      cout << "The projection using finite differences currently only supports a fixed direction of movement for FFD points." << endl;
      exit(EXIT_FAILURE);
    }
    Gradient[iDV] = new su2double[nDV_Value];
  }
  
  /*--- Continuous adjoint gradient computation ---*/
  
  if (rank == MASTER_NODE)
    cout << "Evaluate functional gradient using Finite Differences." << endl;
  
  for (iDV = 0; iDV < nDV; iDV++) {
    
    MoveSurface = true;
    Local_MoveSurface = true;
    
    /*--- Free Form deformation based ---*/
    
    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
        (config->GetDesign_Variable(iDV) == FFD_NACELLE) ||
        (config->GetDesign_Variable(iDV) == FFD_GULL) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST) ||
        (config->GetDesign_Variable(iDV) == FFD_ROTATION) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS) ||
        (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK)) {
      
      /*--- Read the FFD information in the first iteration ---*/
      
      if (iDV == 0) {
        
        if (rank == MASTER_NODE)
          cout << "Read the FFD information from mesh file." << endl;
        
        /*--- Read the FFD information from the grid file ---*/
        
        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());
        
        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition() && (rank == MASTER_NODE)) {
          cout << "The input grid doesn't have the entire FFD information!" << endl;
          cout << "Press any key to exit..." << endl;
          cin.get();
        }
        
        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
          
          if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
          if (rank == MASTER_NODE) cout << "Check the FFD box intersections with the solid surfaces." << endl;
          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
        }
        
        if (rank == MASTER_NODE)
          cout <<"-------------------------------------------------------------------------" << endl;
        
      }
      
      if (rank == MASTER_NODE) {
        cout << endl << "Design variable number "<< iDV <<"." << endl;
        cout << "Performing 3D deformation of the surface." << endl;
      }
      
      /*--- Apply the control point change ---*/
      
      MoveSurface = false;
      
      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
        
        /*--- Reset FFD box ---*/
        
        switch (config->GetDesign_Variable(iDV) ) {
          case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_NACELLE :          Local_MoveSurface = surface_movement->SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_GULL :             Local_MoveSurface = surface_movement->SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ANGLE_OF_ATTACK :  Gradient[iDV][0] = config->GetAoA_Sens(); break;
        }
        
        /*--- Recompute cartesian coordinates using the new control points position ---*/
        
        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }
        
      }
      
    }
    
    /*--- Hicks Henne design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }
    
    /*--- Surface bump design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SURFACE_BUMP) {
      surface_movement->SetSurface_Bump(geometry, config, iDV, true);
    }

    /*--- Kulfan (CST) design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == CST) {
      surface_movement->SetCST(geometry, config, iDV, true);
    }
    
    /*--- Displacement design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == TRANSLATION) {
      surface_movement->SetTranslation(geometry, config, iDV, true);
    }
    
    /*--- Angle of Attack design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) {
      Gradient[iDV][0] = config->GetAoA_Sens();
    }
    
    /*--- Scale design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == SCALE) {
      surface_movement->SetScale(geometry, config, iDV, true);
    }
    
    /*--- Rotation design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == ROTATION) {
      surface_movement->SetRotation(geometry, config, iDV, true);
    }
    
    /*--- NACA_4Digits design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == NACA_4DIGITS) {
      surface_movement->SetNACA_4Digits(geometry, config);
    }
    
    /*--- Parabolic design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == PARABOLIC) {
      surface_movement->SetParabolic(geometry, config);
    }
    
    /*--- Custom design variable ---*/

    else if (config->GetDesign_Variable(iDV) == CUSTOM){
      if (rank == MASTER_NODE)
        cout <<"Custom design variable will be used in external script" << endl;
    }
    
    /*--- Design variable not implement ---*/
    
    else {
      if (rank == MASTER_NODE)
        cout << "Design Variable not implement yet" << endl;
    }
    
    /*--- Load the delta change in the design variable (finite difference step). ---*/
    
    if ((config->GetDesign_Variable(iDV) != ANGLE_OF_ATTACK) &&
        (config->GetDesign_Variable(iDV) != FFD_ANGLE_OF_ATTACK)) {
      
      /*--- If the Angle of attack is not involved, reset the value of the gradient ---*/
      
      my_Gradient = 0.0; Gradient[iDV][0] = 0.0;
      
      if (MoveSurface) {
        
        delta_eps = config->GetDV_Value(iDV);
        
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry->GetnPointDomain()) && UpdatePoint[iPoint]) {
                
                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
                Sensitivity = geometry->vertex[iMarker][iVertex]->GetAuxVar();
                
                dS = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);
                
                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }
                
                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
        }
        
      }
      //else if (config->GetDesign_Variable(iDV)==FIML) my_Gradient =
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&my_Gradient, &Gradient[iDV][0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      Gradient[iDV][0] = my_Gradient;
#endif
      
    }
    
  }
  
  /*--- Print gradients to screen and file ---*/
  
  OutputGradient(Gradient, config, Gradient_file);
  
  for (iDV = 0; iDV  < nDV; iDV++){
    delete [] Gradient[iDV];
  }
  delete [] Gradient;
  delete [] UpdatePoint;
  
}


void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, ofstream& Gradient_file){
  
  su2double DV_Value, *VarCoord, Sensitivity, **Gradient, my_Gradient, *Normal, Area = 0.0;
  unsigned short iDV_Value = 0, iMarker, nMarker, iDim, nDim, nDV_Value;
  unsigned long iVertex, nVertex, iPoint,nDV, iDV;
  bool fiml = false; //JRH 04172017
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();
  
  VarCoord = NULL;
  
  /*--- Structure to store the gradient ---*/
  cout << "JRH Debugging: In SetProjection_AD" << endl;
  Gradient = new su2double*[nDV];
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    Gradient[iDV] = new su2double[nDV_Value];
    if (config->GetDesign_Variable(iDV) == FIML) fiml = true;
  }
  
  /*--- Discrete adjoint gradient computation ---*/
  
  if (rank == MASTER_NODE)
    cout << "Evaluate functional gradient using Algorithmic Differentiation." << endl;
  
  /*--- Start recording of operations ---*/
  cout << "JRH Debugging: Starting AD Recording" << endl;
  AD::StartRecording();

  
  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/
  
  if (config->GetKind_Turb_Model()==SA_FIML) {
	  cout << "JRH Debugging: FIML Case, setting fiml = true" << endl;
	  fiml = true;
  }
  cout << "JRH Debugging: Setting DV Values With config->SetDV_Value()" << endl;
  for (iDV = 0; iDV < nDV; iDV++){
    
    nDV_Value =  config->GetnDV_Value(iDV);
    
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
    	//JRH - Don't want to reset the value of the FIML design variable to zero because it is defined differently then deformation DVs
    	//JRH - But we do want to register the DVs as input, but keep the current values.
      /*--- Initilization with su2double resets the index ---*/
      
    	//if (!fiml)  {
    	//  DV_Value = 0.0; //JRH - 04172017
    	//}
      //else
      DV_Value = config->GetDV_Value(iDV,iDV_Value);
      
      AD::RegisterInput(DV_Value); //This was originally performed for all DVs - JRH 04192017
      config->SetDV_Value(iDV, iDV_Value, DV_Value);

    }
  }
  
  /*--- Call the surface deformation routine ---*/
  //if (!fiml) {
	  cout << "JRH Debugging: Calling SetSurface_Deformation() in SU2_GEO.cpp" << endl;
	  surface_movement->SetSurface_Deformation(geometry, config);
	  cout << "JRH Debugging: Back from SU2_GEO.cpp" << endl;
	  /*--- Stop the recording --- */
  //}
  AD::StopRecording();

  
	  /*--- Initialize the derivatives of the output of the surface deformation routine
	   * with the discrete adjoints from the CFD solution ---*/
	//if (!fiml) { // JRH 04282017
		  for (iMarker = 0; iMarker < nMarker; iMarker++) {
			if (config->GetMarker_All_DV(iMarker) == YES) {
			  nVertex = geometry->nVertex[iMarker];
			  for (iVertex = 0; iVertex <nVertex; iVertex++) {
				iPoint      = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				Normal      = geometry->vertex[iMarker][iVertex]->GetNormal();

				Area = 0.0;
				for (iDim = 0; iDim < nDim; iDim++){
				  Area += Normal[iDim]*Normal[iDim];
				}
				Area = sqrt(Area);

				for (iDim = 0; iDim < nDim; iDim++){
				  if (config->GetDiscrete_Adjoint()){
					Sensitivity = geometry->GetSensitivity(iPoint, iDim);
				  } else {
					Sensitivity = -Normal[iDim]*geometry->vertex[iMarker][iVertex]->GetAuxVar()/Area;
				  }
				  SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
				}
			  }
			}
		  }
	//}

  /*--- Compute derivatives and extract gradient ---*/
  cout << "JRH Debugging: Calling AD::ComputeAdjoint()" << endl;
  AD::ComputeAdjoint();
  if (fiml) {
	  GetBetaFimlGrad(geometry, config);
  }
  cout << "JRH Debugging: Starting to call GetDerivative()" << endl;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      DV_Value = config->GetDV_Value(iDV, iDV_Value);
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
#ifdef HAVE_MPI
      if (!fiml) SU2_MPI::Allreduce(&my_Gradient, &Gradient[iDV][iDV_Value], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      if (!fiml) Gradient[iDV][iDV_Value] = my_Gradient;
#endif
      
      /*--- Angle of Attack design variable (this is different,
       the value comes from the input file) ---*/
      
      if ((config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) ||
          (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK))  {
        Gradient[iDV][iDV_Value] = config->GetAoA_Sens();
      }
      if (fiml) {
    	  //NOTE THE NEGATIVE...TESTING JRH 09042017 - REMOVED NEGATIVE!! 09282017
    	  Gradient[iDV][iDV_Value] = beta_fiml_grad[iDV];
    	  //cout << "JRH DEBUGGING: Gradient[iDV] = " << Gradient[iDV][iDV_Value] << "beta_fiml_grad[iDV] = " << beta_fiml_grad[iDV] << endl;
      }
    }
  }
  
  /*--- Print gradients to screen and file ---*/
  cout << "JRH Debugging: Beginning to output gradient file" << endl;
  if (rank == MASTER_NODE) OutputGradient(Gradient, config, Gradient_file);
  
  for (iDV = 0; iDV  < nDV; iDV++){
    delete [] Gradient[iDV];
  }
  delete [] Gradient;
  delete [] beta_fiml_grad;
  cout << "JRH Debugging: Done in SU2_DOT::SetProjection_AD()" << endl;
}

void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file){
  
  unsigned short iDV_Value, nDV_Value;
  unsigned long nDV, iDV;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nDV = config->GetnDV();
  
  /*--- Loop through all design variables and their gradients ---*/
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (rank == MASTER_NODE){
      
      /*--- Print the kind of design variable on screen ---*/
      
      cout << endl << "Design variable (";
      for (std::map<string, ENUM_PARAM>::const_iterator it = Param_Map.begin(); it != Param_Map.end(); ++it ){
        if (it->second == config->GetDesign_Variable(iDV)){
          cout << it->first << ") number "<< iDV << "." << endl;
        }
      }
      
      /*--- Print the kind of objective function to screen ---*/
      
      for (std::map<string, ENUM_OBJECTIVE>::const_iterator it = Objective_Map.begin(); it != Objective_Map.end(); ++it ){
        if (it->second == config->GetKind_ObjFunc()){
          cout << it->first << " gradient : ";
          if (iDV == 0) Gradient_file << it->first << " gradient " << endl;
        }
      }
      
      /*--- Print the gradient to file and screen ---*/
      
      for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
        cout << Gradient[iDV][iDV_Value];
        if (iDV_Value != nDV_Value-1 ){
          cout << ", ";
        }
        Gradient_file << Gradient[iDV][iDV_Value] << endl;
      }
      cout << endl;
      cout <<"-------------------------------------------------------------------------" << endl;
    }
  }
}

void GetBetaFimlGrad(CGeometry *geometry, CConfig* config){
	string 	restart_filename = "beta_fiml_grad.dat";
	string   text_line;
	ifstream restart_file;

	unsigned long nDV = config->GetnDV();
	beta_fiml_grad = new su2double[nDV];

		  /*--- Open the restart file, and throw an error if this fails. ---*/

	restart_file.open(restart_filename.data(), ios::in);
	if (restart_file.fail()) {
//		if (rank == MASTER_NODE)
			cout << "There is no beta fiml grad restart file!! " << restart_filename.data() << "."<< endl;
			exit(EXIT_FAILURE);
	}


	for (unsigned long iDV = 0; iDV < nDV; iDV++) {
		getline (restart_file, text_line);
		istringstream point_line(text_line);
		point_line >> beta_fiml_grad[iDV];
	}
	restart_file.close();

//	  /*JRH 06232017 - Stole this code from solver_direct_mean.c CEulerSolver::LoadRestart(). Attempting to read in beta_fiml_grad
//	   * from restart file of discrete adjoint solution. This gradient *shouldn't* have to be modified in any way by SU2_DOT and should
//	   * be written directly to the gradient output file.
//	   */
//	  /*--- Restart the solution from file information ---*/
//	  unsigned short iDim, iVar, iMesh, iMeshFine;
//	  unsigned long iPoint, index, iChildren, Point_Fine;
//	  unsigned short turb_model = config->GetKind_Turb_Model();
//	  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine, dull_val;
//	  bool grid_movement  = config->GetGrid_Movement();
//	  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
//	                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
//	  bool steady_restart = config->GetSteadyRestart();
//	  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
//	  string UnstExt, text_line;
//	  ifstream restart_file;
//	  unsigned short nDim = geometry->GetnDim();
//	  unsigned short iZone = config->GetiZone();
//	  unsigned short nZone = geometry->GetnZone();
//	  string restart_filename = config->GetSolution_FlowFileName();
//	  su2double beta_fiml, sensitivity_x, sensitivity_y, Solution_0, Solution_1, Solution_2, Solution_3, Solution_4, Solution_5, Solution_6, sensitivity;
//	  su2double  *beta_fiml_grad_nodes;
//	  su2double *beta_fiml_grad_local;
//	  Coord = new su2double [nDim];
//	  for (iDim = 0; iDim < nDim; iDim++)
//	    Coord[iDim] = 0.0;
//
//	  int rank = MASTER_NODE;
//	  int size;
//	#ifdef HAVE_MPI
//	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	  MPI_Comm_size(MPI_COMM_WORLD, &size);
//	#endif
//
//	  if (turb_model != SA_FIML) cout << "JRH Error SU2_DOT::GetBetaFimlGrad(): Turbulence model not SA_FIML and this routine should only be called if FIML case" << endl;
//
//	  /*--- Multizone problems require the number of the zone to be appended. ---*/
//
//	  //if (nZone > 1)
//	  //  restart_filename = config->GetMultizone_FileName(restart_filename, iZone);
//	  //restart_filename = config->GetRestart_AdjFileName();
//	  //restart_filename = config->GetAdj_FileName();
//	  //restart_filename = config->GetObjFunc_Extension(restart_filename);
//	  restart_filename = "solution_adj_cd.dat";
//	  //restart_filename = "solution_" << config->GetAdj_FileName();
//	  //restart_filename = config->GetObjFunc_Extension(restart_filename);
//	  /*--- Modify file name for an unsteady restart ---*/
//
//	  if (dual_time || time_stepping) {
//		  cout << "JRH Error SU2_DOT::GetBetaFimlGrad(): Unsteady FIML not yet implemented!!" << endl;
//	      //restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
//	  }
//	  cout << "JRH Debugging: In SU2_DOT::GetBetaFimlGrad() Beginning to read in beta_fiml_grad from restart file " << restart_filename << endl;
//
//	  /*--- Open the restart file, and throw an error if this fails. ---*/
//
//	  restart_file.open(restart_filename.data(), ios::in);
//	  if (restart_file.fail()) {
//	    if (rank == MASTER_NODE)
//	      cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
//	    exit(EXIT_FAILURE);
//	  }
//
//	  /*--- In case this is a parallel simulation, we need to perform the
//	   Global2Local index transformation first. ---*/
//
//	  map<unsigned long,unsigned long> Global2Local;
//	  map<unsigned long,unsigned long>::const_iterator MI;
//
//	  /*--- Now fill array with the transform values only for local points ---*/
//	  long nPoint = geometry->GetnPointDomain();
//	  beta_fiml_grad = new su2double[geometry->GetGlobal_nPointDomain()];
//	  beta_fiml_grad_local = new su2double[geometry->GetnPoint()];
//	  for (iPoint = 0; iPoint < nPoint; iPoint++) {
//	    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
//	  }
//
//	  /*--- Read all lines in the restart file ---*/
//	  cout << "JRH Debugging: In SU2_DOT::GetBetaFimlGrad() Beginning to read in beta_fiml_grad from restart file " << restart_filename << endl;
//	  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
//
//	  /*--- The first line is the header ---*/
//
//	  getline (restart_file, text_line);
//	  beta_fiml_grad_nodes = new su2double[geometry->GetGlobal_nPointDomain()];
//	  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++) beta_fiml_grad_nodes[iPoint_Global] = 0.0;
//	  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
//
//	    getline (restart_file, text_line);
//
//	    istringstream point_line(text_line);
//
//	    /*--- Retrieve local index. If this node from the restart file lives
//	     on the current processor, we will load and instantiate the vars. ---*/
//
//	    MI = Global2Local.find(iPoint_Global);
//	    if (MI != Global2Local.end()) {
//
//	      iPoint_Local = Global2Local[iPoint_Global];
//
//	      if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution_0 >> Solution_1 >> Solution_2 >> Solution_3 >> Solution_4 >> sensitivity >> sensitivity_x >> sensitivity_y >> beta_fiml >> beta_fiml_grad_local[iPoint_Local];
//	      if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution_0 >> Solution_1 >> Solution_2 >> Solution_3 >> Solution_4 >> Solution_5 >> sensitivity >> sensitivity_x >> sensitivity_y >> beta_fiml >> beta_fiml_grad_local[iPoint_Local];
//	      if (rank == MASTER_NODE) cout << "JRH_Debugging: For iPoint_Global = " << index << " beta = " << beta_fiml << " beta_fiml_grad = " << beta_fiml_grad_local[iPoint_Local] << endl;
//	     //"PointID"	"x"	"y"	"Conservative_1"	"Conservative_2"	"Conservative_3"	"Conservative_4"	"Conservative_5"	"Surface_Sensitivity"	"Sensitivity_x"	"Sensitivity_y"	"Beta_fiml"	"Beta_fiml_grad"
//
//	      //node[iPoint_Local]->SetSolution(Solution);
//
//	      /*--- For dynamic meshes, read in and store the
//	       grid coordinates and grid velocities for each node. ---*/
//
//	      if (grid_movement) {
//	    	  cout << "JRH ERROR SU2_DOT::GetBetaFimlGrad(): Grid Movement Not Implemented" << endl;
//	      }
//
//	    }
//
//	  }
//	  unsigned long Fiml_Skip_Index = 0;
//#ifdef HAVE_MPI
//	  //if (config->GetKind_Turb_Model() == SA_FIML) {
//		  unsigned long *nPointFiml_Local;
//		  unsigned long *nPointFiml;
//		  nPointFiml_Local = new unsigned long[size];
//		  nPointFiml = new unsigned long[size];
//		  for (unsigned short iRank = 0; iRank < size ; iRank++) nPointFiml_Local[iRank] = 0;
//		  for (iPoint = 0; iPoint < nPoint; iPoint++) {
//			  if (geometry->node[iPoint]->GetDomain()) nPointFiml_Local[rank]++;
//		  }
//		  //SU2_MPI::Allreduce(&Local_Sens_Beta_Fiml[iDV],  &Total_Sens_Beta_Fiml[iDV],  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//		  for (unsigned short iRank = 0; iRank < size ; iRank++) {
//			  SU2_MPI::Allreduce(&nPointFiml_Local[iRank],  &nPointFiml[iRank],  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
//		  }
//		  for (unsigned short iRank = 0 ; iRank < rank ; iRank++) Fiml_Skip_Index += nPointFiml[iRank];
//		  unsigned long Total_Fiml_Nodes = 0;
//		  for (unsigned short iRank = 0 ; iRank < size ; iRank++) Total_Fiml_Nodes += nPointFiml[iRank];
//		  //Fiml_Skip_Index += rank;
//		  cout << "JRH Debugging: In SU2_DOT::GetBetaFimlGrad() " << rank << " nPoint = " << nPoint << " and begins with fiml variable " << Fiml_Skip_Index << endl;
//		  if (rank == MASTER_NODE) cout << "Total_Fiml_Nodes = " << Total_Fiml_Nodes << endl;
//	  //}
//#else
//#endif
//	  	  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++) {
//	  		  beta_fiml_grad[iPoint_Global] = 0.0;
//	  	  }
//
//		  unsigned long nDV_Local = 0;
//		  for (iPoint = 0; iPoint < nPoint; iPoint++ ) {
//		  /*--- Retrieve local index. If this node from the restart file lives
//		   on the current processor, we will load and instantiate the vars. ---*/
//
//		  //MI = Global2Local.find(iPoint_Global);
//		  //if (MI != Global2Local.end()) {
//		  if (geometry->node[iPoint]->GetDomain()) {
//				//iPoint_Local = Global2Local[iPoint_Global];
//					//node[iPoint_Global]->SetBetaFiml(config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0);
//					//beta_fiml_grad_nodes[nDV_Local+Fiml_Skip_Index] = beta_fiml_grad_local[iPoint_Global];
//					//if (rank == MASTER_NODE) cout << "JRH Debugging: In turb solver setting node " << iPoint_Global << " to iDV =  " << nDV_Local+Fiml_Skip_Index << " = " << config->GetDV_Value(nDV_Local+Fiml_Skip_Index,0)+1.0 << endl;
//
//			  beta_fiml_grad_nodes[nDV_Local+Fiml_Skip_Index] = beta_fiml_grad_local[iPoint];
//			  nDV_Local++;
//		  }
//
//
//	  }
//#ifdef HAVE_MPI
//		  for(unsigned long iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
//			  SU2_MPI::Allreduce(&beta_fiml_grad_nodes[iPoint],  &beta_fiml_grad[iPoint],  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
//			  if (rank == MASTER_NODE) cout << "JRH Debugging: In GetBetaFiml iPoint = " << iPoint << " and beta_fiml_grad = " << beta_fiml_grad[iPoint] << endl;
//		  }
//#else
//		  for(unsigned long iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
//			  beta_fiml_grad[iPoint] = beta_fiml_grad_nodes[iPoint];
//		  }
//#endif
//
//	  /*--- Close the restart file ---*/
//
//	  restart_file.close();
//
//	  delete [] Coord;
//	  delete [] beta_fiml_grad_nodes;
//	  delete [] beta_fiml_grad_local;
}
