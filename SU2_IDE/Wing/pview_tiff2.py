#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
for dsn in range(1 , 81) :
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    if dsn < 10:
        dsn_str = "00%s" % str(dsn)
    else :
        if dsn < 100 :
            dsn_str = "0%s" % str(dsn)
        else :
            dsn_str = "%s" % str(dsn)

    # create a new 'Legacy VTK Reader'
    adjointvtk = LegacyVTKReader(FileNames=["/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/DESIGNS/DSN_%s/ADJOINT_DRAG/adjoint.vtk" % dsn_str])
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1006, 563]
    
    # get color transfer function/color map for 'Conservative_1'
    conservative_1LUT = GetColorTransferFunction('Conservative_1')
    conservative_1LUT.RescaleOnVisibilityChange = 1
    conservative_1LUT.RGBPoints = [-0.9132158756256104, 0.0, 0.0, 1.0, 0.4524846076965332, 1.0, 0.0, 0.0]
    conservative_1LUT.ColorSpace = 'HSV'
    conservative_1LUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    conservative_1LUT.ScalarRangeInitialized = 1.0
    
    # show data in view
    adjointvtkDisplay = Show(adjointvtk, renderView1)
    # trace defaults for the display properties.
    adjointvtkDisplay.ColorArrayName = ['POINTS', 'Conservative_1']
    adjointvtkDisplay.LookupTable = conservative_1LUT
    adjointvtkDisplay.OSPRayScaleArray = 'Conservative_1'
    adjointvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    adjointvtkDisplay.SelectOrientationVectors = 'Conservative_1'
    adjointvtkDisplay.ScaleFactor = 20.0
    adjointvtkDisplay.SelectScaleArray = 'Conservative_1'
    adjointvtkDisplay.GlyphType = 'Arrow'
    adjointvtkDisplay.ScalarOpacityUnitDistance = 9.968587327669786
    adjointvtkDisplay.GaussianRadius = 10.0
    adjointvtkDisplay.SetScaleArray = ['POINTS', 'Conservative_1']
    adjointvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    adjointvtkDisplay.OpacityArray = ['POINTS', 'Conservative_1']
    adjointvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    adjointvtkDisplay.OSPRayScaleFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    adjointvtkDisplay.ScaleTransferFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    adjointvtkDisplay.OpacityTransferFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    #changing interaction mode based on data extents
    renderView1.CameraPosition = [0.0, 0.0, 10000.0]
    renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]
    
    # show color bar/color legend
    adjointvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get opacity transfer function/opacity map for 'Conservative_1'
    conservative_1PWF = GetOpacityTransferFunction('Conservative_1')
    conservative_1PWF.Points = [-0.9132158756256104, 0.0, 0.5, 0.0, 0.4524846076965332, 1.0, 0.5, 0.0]
    conservative_1PWF.ScalarRangeInitialized = 1
    
    # set scalar coloring
    ColorBy(adjointvtkDisplay, ('POINTS', 'beta_fiml'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(conservative_1LUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    adjointvtkDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    adjointvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'beta_fiml'
    beta_fimlLUT = GetColorTransferFunction('beta_fiml')
    beta_fimlLUT.LockDataRange = 1
    beta_fimlLUT.RescaleOnVisibilityChange = 1
    beta_fimlLUT.RGBPoints = [0.0, 0.0, 0.0, 1.0, 1.00457, 1.0, 0.0, 0.0]
    beta_fimlLUT.ColorSpace = 'HSV'
    beta_fimlLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    beta_fimlLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'beta_fiml'
    beta_fimlPWF = GetOpacityTransferFunction('beta_fiml')
    beta_fimlPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.00457, 1.0, 0.5, 0.0]
    beta_fimlPWF.ScalarRangeInitialized = 1
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5257879439125188, -0.03926151812399895, 10000.0]
    renderView1.CameraFocalPoint = [0.5257879439125188, -0.03926151812399895, 0.0]
    renderView1.CameraParallelScale = 0.3172364985013145
    
    # save screenshot
    SaveScreenshot('/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/Pics/Beta/beta.tif', magnification=1, quality=100, view=renderView1)
    
    # set scalar coloring
    ColorBy(adjointvtkDisplay, ('POINTS', 'beta_fiml_grad'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(beta_fimlLUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    adjointvtkDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    adjointvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'beta_fiml_grad'
    beta_fiml_gradLUT = GetColorTransferFunction('beta_fiml_grad')
    beta_fiml_gradLUT.LockDataRange = 1
    beta_fiml_gradLUT.RescaleOnVisibilityChange = 1
    beta_fiml_gradLUT.RGBPoints = [-7.04983e-08, 0.0, 0.0, 1.0, 1.42793e-05, 1.0, 0.0, 0.0]
    beta_fiml_gradLUT.ColorSpace = 'HSV'
    beta_fiml_gradLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    beta_fiml_gradLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'beta_fiml_grad'
    beta_fiml_gradPWF = GetOpacityTransferFunction('beta_fiml_grad')
    beta_fiml_gradPWF.Points = [-7.04983e-08, 0.0, 0.5, 0.0, 1.42793e-05, 1.0, 0.5, 0.0]
    beta_fiml_gradPWF.ScalarRangeInitialized = 1
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5257879439125188, -0.03926151812399895, 10000.0]
    renderView1.CameraFocalPoint = [0.5257879439125188, -0.03926151812399895, 0.0]
    renderView1.CameraParallelScale = 0.3172364985013145
    
    # save screenshot
    SaveScreenshot('/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/Pics/Beta_Grad/beta_grad.tif', magnification=1, quality=100, view=renderView1)
    
    # destroy adjointvtk
    Delete(adjointvtk)
    del adjointvtk
    
    # create a new 'Legacy VTK Reader'
    flowvtk = LegacyVTKReader(FileNames=["/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/DESIGNS/DSN_%s/DIRECT/flow.vtk" % dsn_str])
    
    # show data in view
    flowvtkDisplay = Show(flowvtk, renderView1)
    # trace defaults for the display properties.
    flowvtkDisplay.ColorArrayName = ['POINTS', 'Conservative_1']
    flowvtkDisplay.LookupTable = conservative_1LUT
    flowvtkDisplay.OSPRayScaleArray = 'Conservative_1'
    flowvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    flowvtkDisplay.SelectOrientationVectors = 'Conservative_1'
    flowvtkDisplay.ScaleFactor = 20.0
    flowvtkDisplay.SelectScaleArray = 'Conservative_1'
    flowvtkDisplay.GlyphType = 'Arrow'
    flowvtkDisplay.ScalarOpacityUnitDistance = 9.968587327669786
    flowvtkDisplay.GaussianRadius = 10.0
    flowvtkDisplay.SetScaleArray = ['POINTS', 'Conservative_1']
    flowvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    flowvtkDisplay.OpacityArray = ['POINTS', 'Conservative_1']
    flowvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    flowvtkDisplay.OSPRayScaleFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    flowvtkDisplay.ScaleTransferFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    flowvtkDisplay.OpacityTransferFunction.Points = [0.2691405117511749, 0.0, 0.5, 0.0, 0.6025006771087646, 1.0, 0.5, 0.0]
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    #changing interaction mode based on data extents
    renderView1.CameraPosition = [0.0, 0.0, 10000.0]
    renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]
    
    # show color bar/color legend
    flowvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # set scalar coloring
    ColorBy(flowvtkDisplay, ('POINTS', 'Mach'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(conservative_1LUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    flowvtkDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    flowvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'Mach'
    machLUT = GetColorTransferFunction('Mach')
    machLUT.RescaleOnVisibilityChange = 1
    machLUT.RGBPoints = [0.0, 0.0, 0.0, 1.0, 1.1917582750320435, 1.0, 0.0, 0.0]
    machLUT.ColorSpace = 'HSV'
    machLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    machLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'Mach'
    machPWF = GetOpacityTransferFunction('Mach')
    machPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1917582750320435, 1.0, 0.5, 0.0]
    machPWF.ScalarRangeInitialized = 1
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.671591305041831, 0.21068078892535724, 10000.0]
    renderView1.CameraFocalPoint = [0.671591305041831, 0.21068078892535724, 0.0]
    renderView1.CameraParallelScale = 0.5620038085214876
    
    # save screenshot
    SaveScreenshot('/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/Pics/Mach/Mach.tif', magnification=1, quality=100, view=renderView1)
    
    # set scalar coloring
    ColorBy(flowvtkDisplay, ('POINTS', 'Pressure'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(machLUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    flowvtkDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    flowvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'Pressure'
    pressureLUT = GetColorTransferFunction('Pressure')
    pressureLUT.RescaleOnVisibilityChange = 1
    pressureLUT.RGBPoints = [23411.984375, 0.0, 0.0, 1.0, 55217.24609375, 1.0, 0.0, 0.0]
    pressureLUT.ColorSpace = 'HSV'
    pressureLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    pressureLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'Pressure'
    pressurePWF = GetOpacityTransferFunction('Pressure')
    pressurePWF.Points = [23411.984375, 0.0, 0.5, 0.0, 55217.24609375, 1.0, 0.5, 0.0]
    pressurePWF.ScalarRangeInitialized = 1
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.671591305041831, 0.21068078892535724, 10000.0]
    renderView1.CameraFocalPoint = [0.671591305041831, 0.21068078892535724, 0.0]
    renderView1.CameraParallelScale = 0.5620038085214876
    
    # save screenshot
    SaveScreenshot('/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/Pics/Pressure/pressure.tif', magnification=1, quality=100, view=renderView1)
    
    # set scalar coloring
    ColorBy(flowvtkDisplay, ('POINTS', 'Eddy_Viscosity'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(pressureLUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    flowvtkDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    flowvtkDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'Eddy_Viscosity'
    eddy_ViscosityLUT = GetColorTransferFunction('Eddy_Viscosity')
    eddy_ViscosityLUT.RescaleOnVisibilityChange = 1
    eddy_ViscosityLUT.RGBPoints = [1.839728827751606e-31, 0.0, 0.0, 1.0, 0.0037628745194524527, 1.0, 0.0, 0.0]
    eddy_ViscosityLUT.ColorSpace = 'HSV'
    eddy_ViscosityLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    eddy_ViscosityLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'Eddy_Viscosity'
    eddy_ViscosityPWF = GetOpacityTransferFunction('Eddy_Viscosity')
    eddy_ViscosityPWF.Points = [1.839728827751606e-31, 0.0, 0.5, 0.0, 0.0037628745194524527, 1.0, 0.5, 0.0]
    eddy_ViscosityPWF.ScalarRangeInitialized = 1
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.8872091072454573, 0.029002825957487335, 10000.0]
    renderView1.CameraFocalPoint = [0.8872091072454573, 0.029002825957487335, 0.0]
    renderView1.CameraParallelScale = 0.5620038085214876
    
    # save screenshot
    SaveScreenshot('/home/jon/SU2/TestCases/optimization_rans/steady_rae2822/Pics/Eddy/eddy.tif', magnification=1, quality=100, view=renderView1)
    
    # destroy flowvtk
    Delete(flowvtk)
    del flowvtk
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.8872091072454573, 0.029002825957487335, 10000.0]
    renderView1.CameraFocalPoint = [0.8872091072454573, 0.029002825957487335, 0.0]
    renderView1.CameraParallelScale = 0.5620038085214876
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).