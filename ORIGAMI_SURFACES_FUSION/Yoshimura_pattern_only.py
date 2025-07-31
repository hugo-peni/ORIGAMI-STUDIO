import adsk.core, adsk.fusion, traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        
        sides = 6

        # === Define User Parameters ===
        paramName = 'circleDiameter'
        paramValue = 10
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
        
        # === First Sketch on XY Plane ===
        sketches = rootComp.sketches
        xyPlane = rootComp.xYConstructionPlane
        sketch1 = sketches.add(xyPlane)
        
        constraints1 = sketch1.geometricConstraints
        lines1 = sketch1.sketchCurves.sketchLines
        circles1 = sketch1.sketchCurves.sketchCircles
        dimensions1 = sketch1.sketchDimensions
        sketchPoints1 = sketch1.sketchPoints
        
        radius = 5
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
        
        plane_offset = 2.0
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
        num_repetition = 3 
        quantityTwo = adsk.core.ValueInput.createByString(f"{num_repetition}")
        distanceTwo = adsk.core.ValueInput.createByString(f"{plane_offset * 2 } cm")
        
        # Create the input for rectangular pattern
        rectangularPatterns = rootComp.features.rectangularPatternFeatures
        rectangularPatternInput = rectangularPatterns.createInput(inputEntites, xAxis, quantityOne, distanceOne, adsk.fusion.PatternDistanceType.SpacingPatternDistanceType)
        
        # Set the data for second direction
        rectangularPatternInput.setDirectionTwo(zAxis, quantityTwo, distanceTwo)
        
        # Create the rectangular pattern
        rectangularFeature = rectangularPatterns.add(rectangularPatternInput)

       
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

        
    except: 
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))



