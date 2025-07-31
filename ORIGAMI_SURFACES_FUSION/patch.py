import adsk.core, adsk.fusion, traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        
        # Create a document.
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
 
        product = app.activeProduct
        design = adsk.fusion.Design.cast(product)

        # Get the root component of the active design.
        rootComp = design.rootComponent

        sketches = rootComp.sketches 
        xyPlane = rootComp.xYConstructionPlane 

        sketch1 = sketches.add(xyPlane)
        points1 = sketch1.sketchPoints
        lines1 = sketch1.sketchCurves.sketchLines 

        p1 = adsk.core.Point3D.create( 2, 0 , 0 )
        p2 = adsk.core.Point3D.create( 0, 2 , 0 )
        p3 = adsk.core.Point3D.create( 0, -2 , 0 )
        pt1 = points1.add( p1 )
        pt2 = points1.add( p2 )
        pt3 = points1.add( p3 )

        line1 = lines1.addByTwoPoints( pt1 ,pt2 )
        line2 = lines1.addByTwoPoints( pt2 ,pt3 )
        line3 = lines1.addByTwoPoints( pt3 ,pt1 )

        prof1 = sketch1.profiles.item(0)

        patches = rootComp.features.patchFeatures
        patchInput2 = patches.createInput( prof1 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        patch1 = patches.add(patchInput2)

        # replicate it with anotther patch

        sketch2 = sketches.add(xyPlane)
        points2 = sketch2.sketchPoints
        lines2 = sketch2.sketchCurves.sketchLines 

        p1 = adsk.core.Point3D.create( 2, 0 , 0 )
        p2 = adsk.core.Point3D.create( 0, -2 , 0 )
        p3 = adsk.core.Point3D.create( 2, -2 , 0 )

        pt1 = points2.add( p1 )
        pt2 = points2.add( p2 )
        pt3 = points2.add( p3 )

        line1 = lines2.addByTwoPoints( pt1 ,pt2 )
        line2 = lines2.addByTwoPoints( pt2 ,pt3 )
        line3 = lines2.addByTwoPoints( pt3 ,pt1 )

        prof2 = sketch2.profiles.item(0)

        patchInput3 = patches.createInput( prof2 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        patch2 = patches.add(patchInput3)

        # surface = patches.bodies.item(0)
        # surface2 = patches.bodies.item(1)
        surface1 = rootComp.bRepBodies.item(0)
        surface2 = rootComp.bRepBodies.item(1)
        surfaces = adsk.core.ObjectCollection.create()
        surfaces.add( surface1 )
        surfaces.add( surface2 )

        tolerance = adsk.core.ValueInput.createByReal(1.0)
        features = rootComp.features

        stitches = features.stitchFeatures
        stitchInput = stitches.createInput(surfaces, tolerance, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        
        # Create a stitch feature.
        stitch = stitches.add(stitchInput)

        # Create thiken feature
  
        thickenFeatures = rootComp.features.thickenFeatures
        inputSurfaces = adsk.core.ObjectCollection.create()
        bodies = rootComp.bRepBodies.item(0)
        inputSurfaces.add(bodies)
        thickness = adsk.core.ValueInput.createByReal(0.05)
        thickenInput = thickenFeatures.createInput(inputSurfaces, thickness, False,  adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        thickenFeatures.add(thickenInput)
        
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


# import adsk.core, adsk.fusion, traceback

# def run(context):
#     ui = None
#     try:
#         app = adsk.core.Application.get()
#         ui  = app.userInterface

#         # Create a document
#         doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        
#         design = app.activeProduct

#         # Get the root component of the active design.
#         rootComp = design.rootComponent
        
#         # Create sketch
#         sketches = rootComp.sketches
#         sketch = sketches.add(rootComp.xZConstructionPlane)
#         sketchCircles = sketch.sketchCurves.sketchCircles
#         centerPoint = adsk.core.Point3D.create(0, 0, 0)
#         sketchCircle = sketchCircles.addByCenterRadius(centerPoint, 3.0)
        
#         # Create surface
#         openProfile = rootComp.createOpenProfile(sketchCircle)
#         features = rootComp.features
#         extrudes = features.extrudeFeatures
#         extrudeInput = extrudes.createInput(openProfile, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
#         extrudeInput.isSolid = False
#         distance = adsk.core.ValueInput.createByReal(3.0)
#         extrudeInput.setDistanceExtent(False, distance)
#         extrude = extrudes.add(extrudeInput)
        
#         # Create thiken feature
#         thickenFeatures = features.thickenFeatures
#         inputSurfaces = adsk.core.ObjectCollection.create()
#         bodies = extrude.bodies
#         for body in bodies :
#             inputSurfaces.add(body)
#         thickness = adsk.core.ValueInput.createByReal(1.0)
#         thickenInput = thickenFeatures.createInput(inputSurfaces, thickness, False,  adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
#         thickenFeatures.add(thickenInput)
#     except:
#         if ui:
#             ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))