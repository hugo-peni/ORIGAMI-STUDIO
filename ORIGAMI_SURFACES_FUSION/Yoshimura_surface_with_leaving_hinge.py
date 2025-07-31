import adsk.core, adsk.fusion, traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        
        sides = 5

        radius = 2.5
        radius2 = 2.55
        num_repetition = 2
        x = 0.05
        shell_thickness  = 0.1 ; 
        y = shell_thickness + x ; 

    
        thickness = 0.3
        support_length = 1.8 

        # === Define User Parameters ===
        paramName = 'circleDiameter'
        paramValue = 2*radius
        existingParam = design.userParameters.itemByName(paramName)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName, adsk.core.ValueInput.createByReal(paramValue), 'cm', 'Diameter of polygon bounding circle')

        angleParamName = 'angleDeg'
        angleValue = 180
        existingParam = design.userParameters.itemByName(angleParamName)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(angleParamName, adsk.core.ValueInput.createByString(f'{angleValue} deg'), 'deg', 'Rotation angle for first polygon')

        percent = 50
        angleParamName2 = 'angleDeg2'
        angleValue2 = ( 180 - ( 360 / sides ) * (percent / 100) )
        existingParam = design.userParameters.itemByName(angleParamName2)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(angleParamName2, adsk.core.ValueInput.createByString(f'{angleValue2} deg'), 'deg', 'Rotation angle for second polygon')

                # === Define User Parameters ===
        paramName3 = 'dist1'
        paramValue3 = 0.1
        existingParam = design.userParameters.itemByName(paramName3)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName3, adsk.core.ValueInput.createByReal(paramValue3), 'cm', 'Diameter of polygon bounding circle')

        paramName4 = 'dist2'
        paramValue4 = 0.1
        existingParam = design.userParameters.itemByName(paramName4)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName4, adsk.core.ValueInput.createByReal(paramValue4), 'cm', 'Diameter of polygon bounding circle')


        paramName5 = 'dist2'
        paramValue5 = 0.1
        existingParam = design.userParameters.itemByName(paramName5)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName5, adsk.core.ValueInput.createByReal(paramValue5), 'cm', 'Diameter of polygon bounding circle')
  

        paramName6 = 'circleDiameter2'
        paramValue6 = 2*radius2 - 2*thickness ; 
        existingParam = design.userParameters.itemByName(paramName6)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName6, adsk.core.ValueInput.createByReal(paramValue6), 'cm', 'Diameter of polygon bounding circle')

        paramName7 = 'circleDiameter3'
        paramValue7 = 2*radius2 - shell_thickness*2 ; 
        existingParam = design.userParameters.itemByName(paramName7)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName7, adsk.core.ValueInput.createByReal(paramValue7), 'cm', 'Diameter of polygon bounding circle')


        # === First Sketch on XY Plane ===
        sketches = rootComp.sketches
        xyPlane = rootComp.xYConstructionPlane
        sketch1 = sketches.add(xyPlane)
        
        constraints1 = sketch1.geometricConstraints
        lines1 = sketch1.sketchCurves.sketchLines
        circles1 = sketch1.sketchCurves.sketchCircles
        dimensions1 = sketch1.sketchDimensions
        sketchPoints1 = sketch1.sketchPoints
        
        origin = adsk.core.Point3D.create(0, 0, 0)
        
        bounding_circle_1 = circles1.addByCenterRadius(origin, radius)
        bounding_circle_1.isConstruction = True
        centerPoint = bounding_circle_1.centerSketchPoint
        constraints1.addCoincident( centerPoint , sketch1.originPoint)
        
        # Create the polygon
        polygon = lines1.addScribedPolygon( centerPoint , sides , 0 , radius , True)
        
        for i in range( polygon.count ) :
            if( i > 0 ):
                polygon.item(i).isConstruction = False
                
        radius_line_1 = lines1.addByTwoPoints( polygon[0].startSketchPoint , sketch1.originPoint )
        radius_line_1.isConstruction = True
        
        constraints1.addCoincident( polygon[0].startSketchPoint , bounding_circle_1 )
        constraints1.addCoincident( polygon[1].startSketchPoint , bounding_circle_1 )
        constraints1.addCoincident( polygon[2].startSketchPoint , bounding_circle_1 )
        
        # add reference lines
        # add an horizontal reference line
        
        HorizRefLine1 = lines1.addByTwoPoints( sketch1.originPoint , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        constraints1.addHorizontal(HorizRefLine1)
        constraints1.addCoincident( HorizRefLine1.endSketchPoint , bounding_circle_1 )
        HorizRefLine1.isConstruction = True
        
        diameterDim = dimensions1.addDiameterDimension( bounding_circle_1 , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        diameterDim.parameter.expression = paramName
        
        angularDim1 = dimensions1.addAngularDimension( HorizRefLine1, radius_line_1, adsk.core.Point3D.create(3, 6, 0) , True )
        angularDim1.parameter.expression = angleParamName
        
        plane_offset = 1
        planes = rootComp.constructionPlanes
        offsetValue = adsk.core.ValueInput.createByReal(plane_offset)
        planeInput = planes.createInput()
        planeInput.setByOffset(xyPlane, offsetValue)
        offsetPlane = planes.add(planeInput)
        
        sketch2 = sketches.add(offsetPlane)
        sketchPoints2 = sketch2.sketchPoints
        lines2 = sketch2.sketchCurves.sketchLines
        circles2 = sketch2.sketchCurves.sketchCircles
        constraints2 = sketch2.geometricConstraints
        dimensions2 = sketch2.sketchDimensions
        
        radius2 = radius
        sides2 = sides
        
        bounding_circle_2 = circles2.addByCenterRadius( sketch2.originPoint , radius2 )
        bounding_circle_2.isConstruction = True
        
        polygon2 = lines2.addScribedPolygon( sketch2.originPoint , sides2 , 0 , radius2 , True )
        
        constraints2.addCoincident( polygon2[0].startSketchPoint , bounding_circle_2 )
        constraints2.addCoincident( polygon2[1].startSketchPoint , bounding_circle_2 )
        constraints2.addCoincident( polygon2[2].startSketchPoint , bounding_circle_2 )
        
        polygon2_radius_line = lines2.addByTwoPoints( sketch2.originPoint , polygon2[0].startSketchPoint )
        polygon2_radius_line.isConstruction = True
        
        HorizRefLine2 = lines2.addByTwoPoints( sketch2.originPoint , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        constraints2.addHorizontal( HorizRefLine2 )
        constraints2.addCoincident( HorizRefLine2.endSketchPoint , bounding_circle_2 )
        HorizRefLine2.isConstruction = True
        
        for i in range( polygon2.count ):
            if i >= 1:
                polygon2.item(i).isConstruction = True
                
        diameterDim2 = dimensions2.addDiameterDimension( bounding_circle_2 , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        diameterDim2.parameter.expression = paramName
        
        angularDim2 = dimensions2.addAngularDimension( HorizRefLine2 , polygon2_radius_line , adsk.core.Point3D.create(3, 6, 0) , True )
        angularDim2.parameter.expression = angleParamName2
        
        # create two planes for the kresling pattern
        
        planeInput.setByThreePoints( polygon[0].startSketchPoint , polygon[0].endSketchPoint , polygon2[0].startSketchPoint )
        kresling_plane_1 = planes.add( planeInput )
        
        planeInput.setByThreePoints( polygon2[0].startSketchPoint , polygon2[0].endSketchPoint , polygon[0].endSketchPoint )
        kresling_plane_2 = planes.add( planeInput )
        
        # Create construction points
        
        constructionPoints = rootComp.constructionPoints
        pointInput = constructionPoints.createInput()
                    
        construction_point_kresling_1 = []
        construction_point_kresling_2 = []
        
        p1 = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
        pt1 = sketchPoints1.add(p1)
        constraints1.addCoincident(polygon[0].startSketchPoint, pt1)
        pointInput.setByPoint(pt1)
        construction_point_kresling_1.append( constructionPoints.add( pointInput ) )
        
        p2 = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
        pt2 = sketchPoints1.add(p2)
        constraints1.addCoincident(polygon[0].endSketchPoint, pt2)
        pointInput.setByPoint(pt2)
        construction_point_kresling_1.append( constructionPoints.add( pointInput ) )
        
        p3 = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
        pt3 = sketchPoints2.add(p3)
        constraints2.addCoincident(polygon2[0].startSketchPoint, pt3)
        pointInput.setByPoint(pt3)
        construction_point_kresling_2.append( constructionPoints.add( pointInput ) )
        
        p4 = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
        pt4 = sketchPoints2.add(p4)
        constraints2.addCoincident(polygon2[0].endSketchPoint, pt4)
        pointInput.setByPoint(pt4)
        construction_point_kresling_2.append( constructionPoints.add( pointInput ) )
        
        # sketch on kresling 1
        
        sketch_kresling_1 = sketches.add( kresling_plane_1 )
        constraints_kresling_1 = sketch_kresling_1.geometricConstraints
        points_kresling_1 = sketch_kresling_1.sketchPoints
        lines_kresling_1 = sketch_kresling_1.sketchCurves.sketchLines
        
        p1 = adsk.core.Point3D.create(3,3,0)
        pt1 = points_kresling_1.add( p1 )
        p2 = adsk.core.Point3D.create(3,3,0)
        pt2 = points_kresling_1.add( p2 )
        p3 = adsk.core.Point3D.create(3,3,0)
        pt3 = points_kresling_1.add( p3 )
        
        projectedEntity1 = sketch_kresling_1.project( construction_point_kresling_1[0] )
        projectedEntity2 = sketch_kresling_1.project( construction_point_kresling_1[1] )
        projectedEntity3 = sketch_kresling_1.project( construction_point_kresling_2[0] )
            
        constraints_kresling_1.addCoincident( pt1 , projectedEntity1.item(0) )
        constraints_kresling_1.addCoincident( pt2 , projectedEntity2.item(0) )
        constraints_kresling_1.addCoincident( pt3 , projectedEntity3.item(0) )
        
        line1 = lines_kresling_1.addByTwoPoints( pt1 , pt2 )
        line2 = lines_kresling_1.addByTwoPoints( pt2 , pt3 )
        line3 = lines_kresling_1.addByTwoPoints( pt3, pt1 )
        
        prof_kresling_1 = sketch_kresling_1.profiles.item(0)
            
        patches = rootComp.features.patchFeatures
        patchInput1 = patches.createInput( prof_kresling_1 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        patches.add(patchInput1)
        
        # sketch kresling 2
        
        sketch_kresling_2 = sketches.add( kresling_plane_2 )
        constraints_kresling_2 = sketch_kresling_2.geometricConstraints
        points_kresling_2 = sketch_kresling_2.sketchPoints
        lines_kresling_2 = sketch_kresling_2.sketchCurves.sketchLines
        
        p4 = adsk.core.Point3D.create(3,3,0)
        pt4 = points_kresling_2.add( p4 )
        p5 = adsk.core.Point3D.create(3,3,0)
        pt5 = points_kresling_2.add( p5 )
        p6 = adsk.core.Point3D.create(3,3,0)
        pt6 = points_kresling_2.add( p6 )
        
        projectedEntity4 = sketch_kresling_2.project( construction_point_kresling_2[0] )
        projectedEntity5 = sketch_kresling_2.project( construction_point_kresling_2[1] )
        projectedEntity6 = sketch_kresling_2.project( construction_point_kresling_1[1] )
            
        constraints_kresling_2.addCoincident( pt4 , projectedEntity4.item(0) )
        constraints_kresling_2.addCoincident( pt5 , projectedEntity5.item(0) )
        constraints_kresling_2.addCoincident( pt6 , projectedEntity6.item(0) )
        
        line4 = lines_kresling_2.addByTwoPoints( pt4 , pt5 )
        line5 = lines_kresling_2.addByTwoPoints( pt5 , pt6 )
        line6 = lines_kresling_2.addByTwoPoints( pt6, pt4 )
        
        prof_kresling_2 = sketch_kresling_2.profiles.item(0)

        patches = rootComp.features.patchFeatures
        patchInput2 = patches.createInput( prof_kresling_2 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        patches.add(patchInput2)
            
         # Get the body created by extrusion
        body1 = rootComp.bRepBodies.item(0)
        body2 = rootComp.bRepBodies.item(1)
        
        # Create input entities for circular pattern
        inputEntites = adsk.core.ObjectCollection.create()
        zAxis = rootComp.zConstructionAxis
        circularFeats = rootComp.features.circularPatternFeatures
        
        inputEntites.add(body1)
        circularFeatInput = circularFeats.createInput(inputEntites, zAxis)
        circularFeatInput.quantity = adsk.core.ValueInput.createByReal(sides)
        circularFeatInput.totalAngle = adsk.core.ValueInput.createByString('360 deg')
        circularFeatInput.isSymmetric = False
        circularFeat = circularFeats.add(circularFeatInput)
        
        inputEntites.add(body2)
        circularFeatInput = circularFeats.createInput(inputEntites, zAxis)
        circularFeatInput.quantity = adsk.core.ValueInput.createByReal(sides)
        circularFeatInput.totalAngle = adsk.core.ValueInput.createByString('360 deg')
        circularFeatInput.isSymmetric = False
        circularFeat = circularFeats.add(circularFeatInput)

        # rectangular repetition : 

        mirrorFeats = rootComp.features.mirrorFeatures
        bodies = adsk.core.ObjectCollection.create()
        for i in range(len(rootComp.bRepBodies)):
            bodies.add(rootComp.bRepBodies.item(i))

        mirrorInput = mirrorFeats.createInput(bodies, offsetPlane )
        mirrorFeats.add(mirrorInput)

        inputEntites = adsk.core.ObjectCollection.create()

        for i in range( len( rootComp.bRepBodies ) ):
            inputEntites.add( rootComp.bRepBodies.item(i) )
        
        # Get x and y axes for rectangular pattern
        xAxis = rootComp.xConstructionAxis
        zAxis = rootComp.zConstructionAxis
        
        # Quantity and distance
        quantityOne = adsk.core.ValueInput.createByString('1')
        distanceOne = adsk.core.ValueInput.createByString('8 cm') 
        quantityTwo = adsk.core.ValueInput.createByString(f"{num_repetition}")
        distanceTwo = adsk.core.ValueInput.createByString(f"{plane_offset * 2 } cm")
        
        # Create the input for rectangular pattern
        rectangularPatterns = rootComp.features.rectangularPatternFeatures
        rectangularPatternInput = rectangularPatterns.createInput(inputEntites, xAxis, quantityOne, distanceOne, adsk.fusion.PatternDistanceType.SpacingPatternDistanceType)
        
        # Set the data for second direction
        rectangularPatternInput.setDirectionTwo(zAxis, quantityTwo, distanceTwo)
        
        # Create the rectangular pattern
        rectangularFeature = rectangularPatterns.add(rectangularPatternInput)

        ctorPlanes = rootComp.constructionPlanes
        ctorPlaneInput1 = ctorPlanes.createInput()
        offset = adsk.core.ValueInput.createByString("-0.3 cm")
        ctorPlaneInput1.setByOffset(xyPlane, offset)
        ctorPlane1 = ctorPlanes.add(ctorPlaneInput1)
        sketch3 = sketches.add(ctorPlane1)
        sketchCirclesObj3 = sketch3.sketchCurves.sketchCircles
        sketchCirclesObj3.addByCenterRadius( sketch3.originPoint, radius2)


        offset = adsk.core.ValueInput.createByString("-2 cm")
        ctorPlaneInput1.setByOffset(xyPlane, offset)
        ctorPlane2 = ctorPlanes.add(ctorPlaneInput1)
        sketch4 = sketches.add(ctorPlane2)
        sketchCirclesObj4 = sketch4.sketchCurves.sketchCircles
        sketchCirclesObj4.addByCenterRadius( sketch4.originPoint, radius2)

        profile2 = sketch1.profiles.item(0)
        profile3 = sketch3.profiles.item(0)

        loftFeats = rootComp.features.loftFeatures
        loftInput = loftFeats.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        loftSectionsObj = loftInput.loftSections
        loftSectionsObj.add(profile2)
        loftSectionsObj.add(profile3)
        loftInput.isSolid = False
        loftInput.isClosed = False
        loftInput.isTangentEdgesMerged = True
        loftInput.startLoftEdgeAlignment = adsk.fusion.LoftEdgeAlignments.FreeEdgesLoftEdgeAlignment;
        loftInput.endLoftEdgeAlignment = adsk.fusion.LoftEdgeAlignments.FreeEdgesLoftEdgeAlignment;
        loftFeats.add(loftInput)

        profile4 = sketch4.profiles.item(0)

        loftInput2 = loftFeats.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        loftSectionsObj2 = loftInput2.loftSections
        loftSectionsObj2.add(profile3)
        loftSectionsObj2.add(profile4)
        loftInput2.isSolid = False
        loftInput2.isClosed = False
        loftInput2.isTangentEdgesMerged = True
        loftInput2.startLoftEdgeAlignment = adsk.fusion.LoftEdgeAlignments.FreeEdgesLoftEdgeAlignment;
        loftInput2.endLoftEdgeAlignment = adsk.fusion.LoftEdgeAlignments.FreeEdgesLoftEdgeAlignment;
        loftFeats.add(loftInput2)

        midPlaneInput = planes.createInput()
        offset = adsk.core.ValueInput.createByString(f"{plane_offset*num_repetition} cm")
        midPlaneInput.setByOffset(xyPlane, offset)
        midPlane = planes.add(midPlaneInput)

        bodies = adsk.core.ObjectCollection.create()

        max = len( rootComp.bRepBodies )
        bodies.add(rootComp.bRepBodies.item(max-2))  # Add the extruded body to be mirrored
        bodies.add(rootComp.bRepBodies.item(max-1))  # Add the extruded body to be mirrored

        mirrorInput = mirrorFeats.createInput(bodies, midPlane )
        mirrorFeats.add(mirrorInput)

        # create one surface : 
        surfaces = adsk.core.ObjectCollection.create()

        for i in range(len(rootComp.bRepBodies)):
            surface = rootComp.bRepBodies.item(i)
            surfaces.add(surface)
        
        tolerance = adsk.core.ValueInput.createByReal(1.0)
        features = rootComp.features

        stitches = features.stitchFeatures
        stitchInput = stitches.createInput(surfaces, tolerance, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        
        # Create a stitch feature.
        stitch = stitches.add(stitchInput)

        thickenFeatures = rootComp.features.thickenFeatures
        inputSurfaces = adsk.core.ObjectCollection.create()
        bodies = rootComp.bRepBodies.item(0)
        inputSurfaces.add(bodies)
        thickness = adsk.core.ValueInput.createByReal(shell_thickness)
        thickenInput = thickenFeatures.createInput(inputSurfaces, thickness, False,  adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        thickenFeatures.add(thickenInput)

        # sketch on kresling 1
        
        sketch_kresling_3 = sketches.add( kresling_plane_1 )
        constraints_kresling_3 = sketch_kresling_3.geometricConstraints
        points_kresling_3 = sketch_kresling_3.sketchPoints
        lines_kresling_3 = sketch_kresling_3.sketchCurves.sketchLines
        dimensions_3 = sketch_kresling_3.sketchDimensions
        
        p1 = adsk.core.Point3D.create(3,3,0)
        pt1 = points_kresling_3.add( p1 )
        p2 = adsk.core.Point3D.create(3,3,0)
        pt2 = points_kresling_3.add( p2 )
        p3 = adsk.core.Point3D.create(3,3,0)
        pt3 = points_kresling_3.add( p3 )
        
        projectedEntity1 = sketch_kresling_3.project( construction_point_kresling_1[0] )
        projectedEntity2 = sketch_kresling_3.project( construction_point_kresling_1[1] )
        projectedEntity3 = sketch_kresling_3.project( construction_point_kresling_2[0] )
            
        constraints_kresling_3.addCoincident( pt1 , projectedEntity1.item(0) )
        constraints_kresling_3.addCoincident( pt2 , projectedEntity2.item(0) )
        constraints_kresling_3.addCoincident( pt3 , projectedEntity3.item(0) )
        
        line1 = lines_kresling_3.addByTwoPoints( pt1 , pt2 )
        # line1.isConstruction = True 
        line2 = lines_kresling_3.addByTwoPoints( pt2 , pt3 )
        # line2.isConstruction = True 
        line3 = lines_kresling_3.addByTwoPoints( pt3, pt1 )
        # line3.isConstruction = True 

        allConnected = adsk.core.ObjectCollection.create()
        connected1 = sketch_kresling_3.findConnectedCurves(line3)  # Use correct sketch
        for c in connected1:
            allConnected.add(c)

        offset_dim = -0.1

        if allConnected.count > 0:
            dirPoint = adsk.core.Point3D.create(-1, 0, 0)
            offsetCurves = sketch_kresling_3.offset(allConnected, dirPoint, offset_dim)  # Use correct sketch

        line1.isConstruction = True 
        line2.isConstruction = True 
        line3.isConstruction = True 

        prof_kresling_3 = sketch_kresling_3.profiles.item(0)

        #create an extrusion input

        extrudes = rootComp.features.extrudeFeatures
        extInput = extrudes.createInput( prof_kresling_3 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation )
        distance = adsk.core.ValueInput.createByReal(-x)
        extInput.setDistanceExtent( False , distance )
        extInput.isSolid = True
        extrude_kresling_3 = extrudes.add( extInput )

        extInput = extrudes.createInput( prof_kresling_3 , adsk.fusion.FeatureOperations.NewBodyFeatureOperation )
        distance = adsk.core.ValueInput.createByReal(y)
        extInput.setDistanceExtent( False , distance )
        extInput.isSolid = True
        extrude_kresling_5 = extrudes.add( extInput )

                # sketch on kresling 2 with extrusion
        sketch_kresling_4 = sketches.add(kresling_plane_2)
        constraints_kresling_4 = sketch_kresling_4.geometricConstraints
        points_kresling_4 = sketch_kresling_4.sketchPoints
        lines_kresling_4 = sketch_kresling_4.sketchCurves.sketchLines
        
        pt4 = points_kresling_4.add(adsk.core.Point3D.create(3, 3, 0))
        pt5 = points_kresling_4.add(adsk.core.Point3D.create(3, 3, 0))
        pt6 = points_kresling_4.add(adsk.core.Point3D.create(3, 3, 0))
        
        projectedEntity4 = sketch_kresling_4.project(construction_point_kresling_2[0])
        projectedEntity5 = sketch_kresling_4.project(construction_point_kresling_2[1])
        projectedEntity6 = sketch_kresling_4.project(construction_point_kresling_1[1])
        
        constraints_kresling_4.addCoincident(pt4, projectedEntity4.item(0))
        constraints_kresling_4.addCoincident(pt5, projectedEntity5.item(0))
        constraints_kresling_4.addCoincident(pt6, projectedEntity6.item(0))
        
        line4 = lines_kresling_4.addByTwoPoints(pt4, pt5)
        line5 = lines_kresling_4.addByTwoPoints(pt5, pt6)
        line6 = lines_kresling_4.addByTwoPoints(pt6, pt4)

        allConnected2 = adsk.core.ObjectCollection.create()
        for c in sketch_kresling_4.findConnectedCurves(line6):
            allConnected2.add(c)

        offset_dim_2 = -0.1
        if allConnected2.count > 0:
            dirPoint2 = adsk.core.Point3D.create(-1, 0, 0)
            offsetCurves2 = sketch_kresling_4.offset(allConnected2, dirPoint2, -offset_dim_2)

        # Mark original lines as construction
        line4.isConstruction = True
        line5.isConstruction = True
        line6.isConstruction = True

        prof_kresling_4 = sketch_kresling_4.profiles.item(0)

        extrudes = rootComp.features.extrudeFeatures
        extInput2 = extrudes.createInput(prof_kresling_4, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        distance2 = adsk.core.ValueInput.createByReal(x)
        extInput2.setDistanceExtent(False, distance2)
        extInput2.isSolid = True
        extrude_kresling_4 = extrudes.add(extInput2)

        extInput2 = extrudes.createInput(prof_kresling_4, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        distance2 = adsk.core.ValueInput.createByReal(-y)
        extInput2.setDistanceExtent(False, distance2)
        extInput2.isSolid = True
        extrude_kresling_6 = extrudes.add(extInput2)

        # sketch on kresling 2 with extrusion

        extruded_body_3 = extrude_kresling_3.bodies.item(0)
        extruded_body_4 = extrude_kresling_4.bodies.item(0)
        extruded_body_32 = extrude_kresling_5.bodies.item(0)
        extruded_body_42 = extrude_kresling_6.bodies.item(0)

        # === Add circular pattern for extruded kresling 3 ===
        inputEntities_3 = adsk.core.ObjectCollection.create()
        inputEntities_3.add( extruded_body_3 )
        inputEntities_3.add( extruded_body_32 )
        circularInput_3 = circularFeats.createInput( inputEntities_3 , zAxis )
        circularInput_3.quantity = adsk.core.ValueInput.createByReal( sides )
        circularInput_3.totalAngle = adsk.core.ValueInput.createByString('360 deg')
        circularInput_3.isSymmetric = False
        circularFeat_3 = circularFeats.add(circularInput_3)

        # === Add circular pattern for extruded kresling 4 ===
        inputEntities_4 = adsk.core.ObjectCollection.create()
        inputEntities_4.add( extruded_body_4 )
        inputEntities_4.add( extruded_body_42 )
        circularInput_4 = circularFeats.createInput(inputEntities_4, zAxis)
        circularInput_4.quantity = adsk.core.ValueInput.createByReal(sides)
        circularInput_4.totalAngle = adsk.core.ValueInput.createByString('360 deg')
        circularInput_4.isSymmetric = False
        circularFeat_4 = circularFeats.add(circularInput_4)

        bodies_to_mirror = adsk.core.ObjectCollection.create()
        for i in range(len(rootComp.bRepBodies)):
            if i >= 2:
                bodies_to_mirror.add(rootComp.bRepBodies.item(i))

        # Display size
        ui.messageBox(f"Number of bodies to mirror: {len(rootComp.bRepBodies)}")

        mirrorInput = mirrorFeats.createInput(bodies_to_mirror, offsetPlane )
        mirrorFeats.add(mirrorInput)

        inputEntites = adsk.core.ObjectCollection.create()

        for i in range(len(rootComp.bRepBodies)):
            if i >= 2:
                inputEntites.add(rootComp.bRepBodies.item(i))
        
        # Get x and y axes for rectangular pattern
        xAxis = rootComp.xConstructionAxis
        zAxis = rootComp.zConstructionAxis
        
        # Quantity and distance
        quantityOne = adsk.core.ValueInput.createByString('1')
        distanceOne = adsk.core.ValueInput.createByString('8 cm') 
        quantityTwo = adsk.core.ValueInput.createByString(f"{num_repetition}")
        distanceTwo = adsk.core.ValueInput.createByString(f"{plane_offset * 2 } cm")
        
        # Create the input for rectangular pattern
        rectangularPatterns = rootComp.features.rectangularPatternFeatures
        rectangularPatternInput = rectangularPatterns.createInput(inputEntites, xAxis, quantityOne, distanceOne, adsk.fusion.PatternDistanceType.SpacingPatternDistanceType)
        
        # Set the data for second direction
        rectangularPatternInput.setDirectionTwo(zAxis, quantityTwo, distanceTwo)
        
        # Create the rectangular pattern
        rectangularFeature = rectangularPatterns.add(rectangularPatternInput)

        ### SUPPORT 

        # === First Sketch on XY Plane ===

        sketch5 = sketches.add(ctorPlane2)
        
        constraints5 = sketch5.geometricConstraints
        lines5 = sketch5.sketchCurves.sketchLines
        circles5 = sketch5.sketchCurves.sketchCircles
        dimensions5 = sketch5.sketchDimensions
        sketchPoints5 = sketch5.sketchPoints
        
        origin = adsk.core.Point3D.create(0, 0, 0)
        
        circle_1 = circles5.addByCenterRadius(origin, radius)
        centerPoint = circle_1.centerSketchPoint
        constraints5.addCoincident( centerPoint , sketch5.originPoint)

        circle_2 = circles5.addByCenterRadius(origin, radius)
        centerPoint2 = circle_2.centerSketchPoint
        constraints5.addCoincident( centerPoint2 , sketch5.originPoint)
        
        diameterDim = dimensions5.addDiameterDimension( circle_1 , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        diameterDim.parameter.expression = paramName7

        diameterDim2 = dimensions5.addDiameterDimension( circle_2 , adsk.core.Point3D.create( 10 , 0 , 0 ) )
        diameterDim2.parameter.expression = paramName6

        prof_kresling_5 = sketch5.profiles.item(0)

        extrudes = rootComp.features.extrudeFeatures
        extInput3 = extrudes.createInput(prof_kresling_5, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        distance3 = adsk.core.ValueInput.createByReal(support_length)
        extInput3.setDistanceExtent(False, distance3)
        extInput3.isSolid = True
        extrude_kresling_3 = extrudes.add(extInput3)

        bodies_support = adsk.core.ObjectCollection.create()

        max = len( rootComp.bRepBodies )
        bodies_support.add(rootComp.bRepBodies.item(max-1))  # Add the extruded body to be mirrored

        mirrorInput = mirrorFeats.createInput(bodies_support, midPlane )
        mirrorFeats.add(mirrorInput)

    except: 
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))



