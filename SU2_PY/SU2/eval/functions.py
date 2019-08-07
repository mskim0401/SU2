#!/usr/bin/env python

## \file functions.py
#  \brief python package for functions
#  \author T. Lukaczyk, F. Palacios
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from ..io import redirect_folder, redirect_output


# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def function( func_name, config, state=None ):
    """ val = SU2.eval.func(func_name,config,state=None)
    
        Evaluates the aerodynamics and geometry functions.
        
        Wraps:
            SU2.eval.aerodynamics()
            SU2.eval.geometry()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh need not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS is not empty.
        
        Executes in:
            ./DIRECT or ./GEOMETRY
        
        Inputs:
            func_name - SU2 objective function name or 'ALL'
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            If func_name is 'ALL', returns a Bunch() of 
            functions with keys of objective function names
            and values of objective function floats.
            Otherwise returns a float.
    """
    
    # initialize
    state = su2io.State(state)
    
    # check for multiple objectives
    multi_objective = (type(func_name)==list)
    # func_name_string is only used to check whether the function has already been evaluated. 
    func_name_string = func_name
    if multi_objective:   func_name_string = func_name[0]  

    # redundancy check
    if not state['FUNCTIONS'].has_key(func_name_string):

        # Aerodynamics
        if multi_objective or func_name == 'ALL' or func_name in su2io.optnames_aero + su2io.grad_names_directdiff:
            aerodynamics( config, state )
            
        # Stability
        elif func_name in su2io.optnames_stab:
            stability( config, state )
        
        # Geometry
        elif func_name in su2io.optnames_geo:
            geometry( func_name, config, state )
            
        else:
            raise Exception, 'unknown function name, %s' % func_name
        
    #: if not redundant

    # prepare output
    if func_name == 'ALL':
        func_out = state['FUNCTIONS']
    elif (multi_objective):
        # If combine_objective is true, use the 'combo' output.
        objectives=config.OPT_OBJECTIVE
        func_out = 0.0
        for func in func_name:
            sign = su2io.get_objectiveSign(func)
            func_out+=state['FUNCTIONS'][func]*objectives[func]['SCALE']*sign
        state['FUNCTIONS']['COMBO'] = func_out
        
    elif (config.MULTI_MESH == "YES") : #NEW for combination of diff obj funs 4/26/2019 JRH
        func_out = state['FUNCTIONS']['COMBO'] #Populated in aerodynamics() below
        sys.stdout.write('Multi Mesh Simulation, Composite OF = '+str(func_out)+'\n')
    else:
        func_out = state['FUNCTIONS'][func_name]
        
    
    return copy.deepcopy(func_out)

#: def function()


# ----------------------------------------------------------------------
#  Aerodynamic Functions
# ----------------------------------------------------------------------

def aerodynamics( config, state=None ):
    """ vals = SU2.eval.aerodynamics(config,state=None)
    
        Evaluates aerodynamics with the following:
	          SU2.run.deform()
            SU2.run.direct()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS is not empty.
            
        Executes in:
            ./DIRECT
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            Bunch() of functions with keys of objective function names
            and values of objective function floats.
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    if config.NUM_CASES == 0 :
    
        # console output
        if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
            log_direct = 'log_Direct.out'
        else:
            log_direct = None
        
        # ----------------------------------------------------    
        #  Update Mesh
        # ----------------------------------------------------
        
        # does decomposition and deformation
        info = update_mesh(config,state)
        
        # ----------------------------------------------------    
        #  Adaptation (not implemented)
        # ----------------------------------------------------
        
        #if not state.['ADAPTED_FUNC']:
        #    config = su2run.adaptation(config)
        #    state['ADAPTED_FUNC'] = True
        
        # ----------------------------------------------------    
        #  Direct Solution
        # ----------------------------------------------------    
        
        # redundancy check
        direct_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_aero[:9] ] )
        if direct_done:
            # return aerodynamic function values
            aero = su2util.ordered_bunch()
            for key in su2io.optnames_aero:
                if state.FUNCTIONS.has_key(key):
                    aero[key] = state.FUNCTIONS[key]
            return copy.deepcopy(aero)    
        #: if redundant
        
        # files to pull
        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files['MESH']
        name = su2io.expand_part(name,config)
        link.extend(name)
        
        # files: direct solution
        #JRH - solution_flow_#.dat I think
        if files.has_key('DIRECT'):
            name = files['DIRECT']
            name = su2io.expand_time(name,config) #JRH-Returns name if steady...
            link.extend( name )
            ##config['RESTART_SOL'] = 'YES' # don't override config file
        else:
            config['RESTART_SOL'] = 'NO'
            
        # files: target equivarea distribution
        if ( 'EQUIV_AREA' in special_cases and 
             'TARGET_EA' in files ) : 
            pull.append( files['TARGET_EA'] )
    
        # files: target pressure distribution
        if ( 'INV_DESIGN_CP' in special_cases and
             'TARGET_CP' in files ) :
            pull.append( files['TARGET_CP'] )
    
        # files: target heat flux distribution
        if ( 'INV_DESIGN_HEATFLUX' in special_cases and
             'TARGET_HEATFLUX' in files ) :
            pull.append( files['TARGET_HEATFLUX'] )
            
        # files: NN Weights File
        if 'TRAIN_NN' in special_cases:
            pull.append(files['TRAIN_NN'])   
            
        
    
        # output redirection
        with redirect_folder( 'DIRECT', pull, link ) as push:
            with redirect_output(log_direct):     
                
                # # RUN DIRECT SOLUTION # #
                info = su2run.direct(config)
                su2io.restart2solution(config,info)
                state.update(info)
                
                # direct files to push
                name = info.FILES['DIRECT']
                name = su2io.expand_time(name,config)
                push.extend(name)
                          
                
                # equivarea files to push
                if 'WEIGHT_NF' in info.FILES:
                    push.append(info.FILES['WEIGHT_NF'])
    
                # pressure files to push
                if 'TARGET_CP' in info.FILES:
                    push.append(info.FILES['TARGET_CP'])
    
                # heat flux files to push
                if 'TARGET_HEATFLUX' in info.FILES:
                    push.append(info.FILES['TARGET_HEATFLUX'])
                
                if 'TRAIN_NN' in info.FILES:
                    push.append(info.FILES['TRAIN_NN'])
                    
        #: with output redirection
        # return output 
        funcs = su2util.ordered_bunch()
        for key in su2io.optnames_aero + su2io.grad_names_directdiff:
            if state['FUNCTIONS'].has_key(key):
                funcs[key] = state['FUNCTIONS'][key]
                
        if 'OUTFLOW_GENERALIZED' in config.OBJECTIVE_FUNCTION:    
            import downstream_function
            state['FUNCTIONS']['OUTFLOW_GENERALIZED']=downstream_function.downstream_function(config,state)
    
    
    else : #MULTIPLE CONFIGS FOR WEIGHTS PROBLEM - JRH 10042018
        # console output
        for i in range(0,config.NUM_CASES) :
            #sys.stdout.write('JRH: In functions.py->aerodynamics running direct solution design # '+str(i)+'\n')
            if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
                log_direct = 'log_Direct_' + str(i)+ '.out'
            else:
                log_direct = None
            
            # ----------------------------------------------------    
            #  Update Mesh
            # ----------------------------------------------------
            
            # does decomposition and deformation
            if i == 0 : info = update_mesh(config,state)
            
            # ----------------------------------------------------    
            #  Adaptation (not implemented)
            # ----------------------------------------------------
            
            #if not state.['ADAPTED_FUNC']:
            #    config = su2run.adaptation(config)
            #    state['ADAPTED_FUNC'] = True
            
            # ----------------------------------------------------    
            #  Direct Solution
            # ----------------------------------------------------    
            
            # redundancy check
            direct_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_aero[:9] ] )
            if i == 0 :
                if direct_done:
                    # return aerodynamic function values
                    aero = su2util.ordered_bunch()
                    for key in su2io.optnames_aero:
                        if state.FUNCTIONS.has_key(key):
                            aero[key] = state.FUNCTIONS[key]
                    return copy.deepcopy(aero)    
                #: if redundant
            
            # files to pull
            files = state.FILES
            pull = []; link = []
            
            # files: mesh
            #name = files['MESH']
            if config.NUM_CASES > 1 and config.MULTI_MESH == "YES" :
                name = files['MESH'+'_'+str(i)]
            else :
                name = files['MESH']            
            name = su2io.expand_part(name,config)
            link.extend(name)
            
            if i == 0 : curr_config = config
            else : 
                curr_config = su2io.Config(files['CONFIG_'+str(i)])
                curr_config.GRADIENT_METHOD = config.GRADIENT_METHOD
                curr_config.NUMBER_PART = config.NUMBER_PART
                curr_config.DV_VALUE = config.DV_VALUE
                curr_config.DV_VALUE_OLD = config.DV_VALUE_OLD
                curr_config.DV_VALUE_NEW = config.DV_VALUE_NEW
                #curr_config.NPOIN = config.NPOIN
            #sys.stdout.write('Current design variables: '+str(curr_config.DV_VALUE)+'\n\n')
            curr_config.CONFIG_I= str(i)
            
            # files: direct solution
            #JRH - solution_flow_#.dat I think
            if files.has_key('DIRECT_'+str(i)):
                name = files['DIRECT_'+str(i)]
                name = su2io.expand_time(name,config) #JRH-Returns name if steady...
                link.extend( name )
                ##config['RESTART_SOL'] = 'YES' # don't override config file
            else:
                config['RESTART_SOL'] = 'NO'
                
            # files: target equivarea distribution
            if ( 'EQUIV_AREA' in special_cases and 
                 'TARGET_EA' in files ) : 
                pull.append( files['TARGET_EA_'+str(i)] )
        
            # files: target pressure distribution
            if ( 'INV_DESIGN_CP' in special_cases and
                 'TARGET_CP_'+str(i) in files ) :
                pull.append( files['TARGET_CP_'+str(i)] )
        
            # files: target heat flux distribution
            if ( 'INV_DESIGN_HEATFLUX' in special_cases and
                 'TARGET_HEATFLUX' in files ) :
                pull.append( files['TARGET_HEATFLUX'+str(i)] )
                
            # files: NN Weights File
            if 'TRAIN_NN' in special_cases:
                pull.append(files['TRAIN_NN'])       
            
            if i > 0 : pull.append('beta_fiml.dat')
            
            # output redirection
            with redirect_folder( 'DIRECT_'+str(i), pull, link ) as push:
                with redirect_output(log_direct):     
                    
                    # # RUN DIRECT SOLUTION # #
                    info = su2run.direct(curr_config)
                    su2io.restart2solution(curr_config,info)
                    state.update(info)
                    
                    # direct files to push
                    if i == 0 : name = info.FILES['DIRECT']
                    else : 
                        #name = info.FILES['DIRECT_'+str(i)]
                        name = info.FILES['DIRECT'] #Need to save info?
                    
                    name = su2io.expand_time(name,curr_config)
                    push.extend(name)
                              
                    
                    # equivarea files to push
                    if 'WEIGHT_NF_'+str(i) in info.FILES:
                        push.append(info.FILES['WEIGHT_NF_'+str(i)])
        
                    # pressure files to push
                    if 'TARGET_CP_'+str(i) in info.FILES:
                        push.append(info.FILES['TARGET_CP_'+str(i)])
        
                    # heat flux files to push
                    if 'TARGET_HEATFLUX_'+str(i) in info.FILES:
                        push.append(info.FILES['TARGET_HEATFLUX_'+str(i)])
                    
                    if 'TRAIN_NN' in info.FILES:
                        push.append(info.FILES['TRAIN_NN'])
                        
            #: with output redirection
            # return output 
            if i == 0 : funcs = su2util.ordered_bunch()
            for key in su2io.optnames_aero + su2io.grad_names_directdiff:
                if state['FUNCTIONS'].has_key(key):
                    if i == 0: funcs[key] = state['FUNCTIONS'][key]
                    else :
                        if key != "COMBO": 
                            funcs[key] += state['FUNCTIONS'][key]
                     #sys.stdout.write(str(funcs[key]))
                    state['FUNCTIONS'][key] = funcs[key]
                if config.MULTI_MESH == "YES" :
                    objective = curr_config.OPT_OBJECTIVE
                    if objective.has_key(key) :
                        #sys.stdout.write(str(objective))
                        sign = su2io.get_objectiveSign(key)
                        sys.stdout.write('Config '+str(i)+' '+key+' = '+str(state['FUNCTIONS'][key])+' * '+str(objective[key]['SCALE'])+' * '+str(sign)+'\n')
                        if i == 0: funcs['COMBO'] = state['FUNCTIONS'][key]*objective[key]['SCALE']
                        else : funcs['COMBO'] += state['FUNCTIONS'][key]*objective[key]['SCALE']
                        state['FUNCTIONS']['COMBO'] = funcs['COMBO']
                        
                     
            if 'OUTFLOW_GENERALIZED' in config.OBJECTIVE_FUNCTION:    
                import downstream_function
                if i == 0 : state['FUNCTIONS']['OUTFLOW_GENERALIZED']=downstream_function.downstream_function(config,state)
                else : state['FUNCTIONS']['OUTFLOW_GENERALIZED']=[sum(x) for x in zip(downstream_function.downstream_function(config,state),downstream_function.downstream_function(config,state))]
    #sys.stdout.write('JRH: In functions.py->aerodynamics Done With Direct Solutions # \n')
                
    return funcs

#: def aerodynamics()


# ----------------------------------------------------------------------
#  Stability Functions
# ----------------------------------------------------------------------

def stability( config, state=None, step=1e-2 ):
   
    
    folder = 'STABILITY' # os.path.join('STABILITY',func_name) #STABILITY/D_MOMENT_Y_D_ALPHA/
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_direct = 'log_Direct.out'
    else:
        log_direct = None
    
    # ----------------------------------------------------    
    #  Update Mesh
    # ----------------------------------------------------
  
    
    # does decomposition and deformation
    info = update_mesh(config,state) 
    
    # ----------------------------------------------------    
    #  CENTRAL POINT
    # ----------------------------------------------------    
    
    # will run in DIRECT/
    func_0 = aerodynamics(config,state)      
    
    
    # ----------------------------------------------------    
    #  Run Forward Point
    # ----------------------------------------------------   
    
    # files to pull
    files = state.FILES
    pull = []; link = []
    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)
    
    # files: direct solution
    if files.has_key('DIRECT'):
        name = files['DIRECT']
        name = su2io.expand_time(name,config)
        link.extend( name )
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config['RESTART_SOL'] = 'NO'
        
    # files: target equivarea distribution
    if ( 'EQUIV_AREA' in special_cases and 
         'TARGET_EA' in files ) : 
        pull.append( files['TARGET_EA'] )

    # files: target pressure distribution
    if ( 'INV_DESIGN_CP' in special_cases and
         'TARGET_CP' in files ) :
        pull.append( files['TARGET_CP'] )

    # files: target heat flux distribution
    if ( 'INV_DESIGN_HEATFLUX' in special_cases and
         'TARGET_HEATFLUX' in files ) :
        pull.append( files['TARGET_HEATFLUX'] )

    # pull needed files, start folder
    with redirect_folder( folder, pull, link ) as push:
        with redirect_output(log_direct):     
            
            konfig = copy.deepcopy(config)
            ztate  = copy.deepcopy(state)
            
            # TODO: GENERALIZE
            konfig.AoA = konfig.AoA + step
            ztate.FUNCTIONS.clear()
            
            func_1 = aerodynamics(konfig,ztate)
                        
            ## direct files to store
            #name = ztate.FILES['DIRECT']
            #if not state.FILES.has_key('STABILITY'):
                #state.FILES.STABILITY = su2io.ordered_bunch()
            #state.FILES.STABILITY['DIRECT'] = name
            
            ## equivarea files to store
            #if 'WEIGHT_NF' in ztate.FILES:
                #state.FILES.STABILITY['WEIGHT_NF'] = ztate.FILES['WEIGHT_NF']
    
    # ----------------------------------------------------    
    #  DIFFERENCING
    # ----------------------------------------------------
        
    for derv_name in su2io.optnames_stab:

        matches = [ k for k in su2io.optnames_aero if k in derv_name ]
        if not len(matches) == 1: continue
        func_name = matches[0]

        obj_func = ( func_1[func_name] - func_0[func_name] ) / step
        
        state.FUNCTIONS[derv_name] = obj_func
    

    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_stab:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]    
    
    return funcs
    
    
    
    
# ----------------------------------------------------------------------
#  Geometric Functions
# ----------------------------------------------------------------------

def geometry( func_name, config, state=None ):
    """ val = SU2.eval.geometry(config,state=None)
    
        Evaluates geometry with the following:
            SU2.run.deform()
            SU2.run.geometry()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS does not have func_name.
            
        Executes in:
            ./GEOMETRY
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            Bunch() of functions with keys of objective function names
            and values of objective function floats.
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_geom = 'log_Geometry.out'
    else:
        log_geom = None
    
    # ----------------------------------------------------
    #  Update Mesh (check with Trent)
    # ----------------------------------------------------
    
    # does decomposition and deformation
    #info = update_mesh(config,state)


    # ----------------------------------------------------    
    #  Geometry Solution
    # ----------------------------------------------------    
    
    # redundancy check
    geometry_done = state.FUNCTIONS.has_key(func_name)
    #geometry_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_geo ] )
    if not geometry_done:    
        
        # files to pull
        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files['MESH']
        name = su2io.expand_part(name,config)
        link.extend(name)
        
        # update function name
        ## TODO
        
        # output redirection
        with redirect_folder( 'GEOMETRY', pull, link ) as push:
            with redirect_output(log_geom):     
                
                # setup config
                config.GEO_PARAM = func_name
                config.GEO_MODE  = 'FUNCTION'
                
                # # RUN GEOMETRY SOLUTION # #
                info = su2run.geometry(config)
                state.update(info)
                
                # no files to push
                
        #: with output redirection
        
    #: if not redundant 
    
    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_geo:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]
    return funcs
    

#: def geometry()



def update_mesh(config,state=None):
    """ SU2.eval.update_mesh(config,state=None)
    
        updates mesh with the following:
	          SU2.run.deform()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            
        Executes in:
            ./DECOMP and ./DEFORM
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            nothing
            
        Modifies:
            config and state by reference
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_decomp = 'log_Decomp.out'
        log_deform = 'log_Deform.out'
    else:
        log_decomp = None
        log_deform = None
    
        
    # ----------------------------------------------------
    #  Deformation
    # ----------------------------------------------------
    
    # redundancy check
    deform_set  = config['DV_KIND'] == config['DEFINITION_DV']['KIND']
    deform_todo = not config['DV_VALUE_NEW'] == config['DV_VALUE_OLD']
    if deform_set and deform_todo:
    
        # files to pull
        pull = []
        link = config['MESH_FILENAME']
        link = su2io.expand_part(link,config)
        
        # output redirection
        with redirect_folder('DEFORM',pull,link) as push:
            with redirect_output(log_deform):
                
                # # RUN DEFORMATION # #
                info = su2run.deform(config)
                state.update(info)
                
                # data to push
                meshname = info.FILES.MESH
                names = su2io.expand_part( meshname , config )
                push.extend( names )
        
        #: with redirect output
        
    elif deform_set and not deform_todo:
        state.VARIABLES.DV_VALUE_NEW = config.DV_VALUE_NEW

    #: if not redundant

    return 

