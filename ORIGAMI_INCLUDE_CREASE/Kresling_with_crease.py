
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
        height = 2 

        # === Define User Parameters ===
        paramName = 'circleDiameter'
        paramValue = 9
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
        
        angularDim1 = dimensions1.addAngularDimension( HorizRefLine1, radius_line_1, adsk.core.Point3D.create(10, 0, 0), True)
        angularDim1.parameter.expression = angleParamName
        
        plane_offset = height
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
        
        angularDim2 = dimensions2.addAngularDimension( HorizRefLine2 , polygon2_radius_line , adsk.core.Point3D.create(10,0,0) , True )
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
        thickness = adsk.core.ValueInput.createByReal(0.05)
        thickenInput = thickenFeatures.createInput(inputSurfaces, thickness, False,  adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        thickenFeatures.add(thickenInput)



        highlight_body = rootComp.bRepBodies.item(0)
        num_edges = highlight_body.edges.count
        ui.activeSelections.clear()  # Clears previous selections
        edge = highlight_body.edges.item(1)  # Gets the second edge from the BRepBody
        ui.activeSelections.add(edge)       # Adds that edge to the active selection

        distance = adsk.core.ValueInput.createByReal(0.9)
        planeInput.setByDistanceOnPath(edge, distance)
        plane_living_hinge = planes.add(planeInput)

        distance = adsk.core.ValueInput.createByReal(0.1)
        planeInput.setByDistanceOnPath(edge, distance)
        plane_living_hinge_end = planes.add(planeInput)

        sketch_living_hinge = sketches.add(  plane_living_hinge )
        constraints_sketch_living_hinge = sketch_living_hinge.geometricConstraints
        points_sketch_living_hinge = sketch_living_hinge.sketchPoints
        lines_sketch_living_hinge = sketch_living_hinge.sketchCurves.sketchLines
        circles_sketch_living_hinge = sketch_living_hinge.sketchCurves.sketchCircles

        entities = []
        for i in range( highlight_body.faces.count ):
            entities.append(highlight_body.faces.item(i))
        
        sketchEntities = sketch_living_hinge.intersectWithSketchPlane(entities)


        # entity1 = sketchEntities[0]
        # entity1.isConstruction = True
        # entity2 = sketchEntities[1]
        # entity2.isConstruction = True

        # start1 = entity1.startSketchPoint.geometry
        # end1 = entity1.endSketchPoint.geometry

        # start2 = entity2.startSketchPoint.geometry
        # end2 = entity2.endSketchPoint.geometry

        # i = 0 
        # j = 0 

        # if end1.y >= 0 and end2.y <= 0:
        #     if start1.y >= 0 and start2.y <= 0:
        #         i = 0 
        #         j = 1
        #     elif start1.y <= 0 and start2.y >= 0:
        #         i = 1
        #         j = 0
        # elif end1.y <= 0 and end2.y >= 0:
        #     if start1.y >= 0 and start2.y <= 0:
        #         i = 0 
        #         j = 1
        #     elif start1.y <= 0 and start2.y >= 0:
        #         i = 1
        #         j = 0
        # else:
        #     i = 0
        #     j = 1

        # msg = f"""Entity 0:
        # Start: ({start1.x:.3f}, {start1.y:.3f}, {start1.z:.3f})
        # End:   ({end1.x:.3f}, {end1.y:.3f}, {end1.z:.3f})

        # Entity 1:
        # Start: ({start2.x:.3f}, {start2.y:.3f}, {start2.z:.3f})
        # End:   ({end2.x:.3f}, {end2.y:.3f}, {end2.z:.3f})
        # """
        # # ui.messageBox(msg)

        # # # Prepare the message
        # # msg = " "; 
        # # for i in range( len(sketchEntities )):
        # #     msg += f"Intersected {sketchEntities[i]} entities:\n"

        # # # Show in message box
        
        # HozRefLH = lines_sketch_living_hinge.addByTwoPoints(sketch_living_hinge.originPoint, adsk.core.Point3D.create(1, 0, 0))
        # constraints_sketch_living_hinge.addHorizontal(HozRefLH)
        # HozRefLH.isConstruction = True

        # VerRefLH = lines_sketch_living_hinge.addByTwoPoints(sketch_living_hinge.originPoint, adsk.core.Point3D.create(0, 1, 0))
        # constraints_sketch_living_hinge.addVertical(VerRefLH)
        # VerRefLH.isConstruction = True


        # circus = circles_sketch_living_hinge.addByCenterRadius( sketch_living_hinge.originPoint , 0.3 )
        # circus.isConstruction = True
        # circus2 = circles_sketch_living_hinge.addByCenterRadius( sketch_living_hinge.originPoint , 1 )
        # circus2.isConstruction = True 
        # circus3 = circles_sketch_living_hinge.addByCenterRadius( sketch_living_hinge.originPoint , 0.2 )
        # circus3.isConstruction = False
        # # circles_sketch_living_hinge.addByCenterRadius(sketch_living_hinge.originPoint, 0.3)

        # line1 = lines_sketch_living_hinge.addByTwoPoints(sketch_living_hinge.originPoint, adsk.core.Point3D.create(1,1,0))
        # line1.isConstruction = True 

        # # line2 = lines_sketch_living_hinge.addByTwoPoints(sketch_living_hinge.originPoint, adsk.core.Point3D.create(1,2,0))
        
        # # constraints_sketch_living_hinge.addCoincident(line1.endSketchPoint, sketchEntities[0].endSketchPoint )
   
        # line3 = lines_sketch_living_hinge.addByTwoPoints(adsk.core.Point3D.create(-0.5,-0.5,0), adsk.core.Point3D.create(-1,-1,0))
        # constraints_sketch_living_hinge.addCollinear(line3, sketchEntities[i] )
        # constraints_sketch_living_hinge.addCoincident( line3.startSketchPoint, circus )
        # constraints_sketch_living_hinge.addCoincident( line3.endSketchPoint, circus2 )
        
        # line4 = lines_sketch_living_hinge.addByTwoPoints(adsk.core.Point3D.create(0.5,0.5,0), adsk.core.Point3D.create(1,1,0))
        # constraints_sketch_living_hinge.addCollinear(line4, sketchEntities[j] )
        # constraints_sketch_living_hinge.addCoincident( line4.startSketchPoint, circus )
        # constraints_sketch_living_hinge.addCoincident( line4.endSketchPoint, circus2 )

        # # line4 = lines_sketch_living_hinge.addByTwoPoints(line3.startSketchPoint , adsk.core.Point3D.create(-0.1, 0.2, 0) )
        # # constraints_sketch_living_hinge.addPerpendicular(line4, line3)

        # # line4REF = lines_sketch_living_hinge.addByTwoPoints(line4.endSketchPoint , adsk.core.Point3D.create(-0.1, 0.2, 0) )
        # # constraints_sketch_living_hinge.addPerpendicular(line4REF, HozRefLH)
        # # constraints_sketch_living_hinge.addCoincident(line4REF.endSketchPoint, HozRefLH )

        # # Get the arcs collection from the sketch
        # arcs = sketch_living_hinge.sketchCurves.sketchArcs

        # # Add arc by center, start, and sweep angle
        # arc = arcs.addByCenterStartSweep(
        #     adsk.core.Point3D.create(0, 0, 0),
        #     adsk.core.Point3D.create(0.3, 0, 0),
        #     -90
        # )

        # # Add constraints: center and endpoints coincident
        # constraints_sketch_living_hinge.addCoincident(arc.centerSketchPoint, sketch_living_hinge.originPoint)
        # constraints_sketch_living_hinge.addCoincident(arc.startSketchPoint, line3.startSketchPoint)
        # constraints_sketch_living_hinge.addCoincident(arc.endSketchPoint, line4.startSketchPoint)

        # # Add fillets: use .geometry for SketchPoints
        # arc2 = sketch_living_hinge.sketchCurves.sketchArcs.addFillet(
        #     line3, line3.startSketchPoint.geometry,
        #     arc, arc.startSketchPoint.geometry,
        #     0.05
        # )

        # arc3 = sketch_living_hinge.sketchCurves.sketchArcs.addFillet(
        #     arc, arc.endSketchPoint.geometry,
        #     line4, line4.startSketchPoint.geometry,
        #     0.05
        # )

        # # Find connected curves to line3 (including arc + fillets if connected)
        # curves = sketch_living_hinge.findConnectedCurves(line3)

        # # # Offset the connected curves
        # if curves.count > 0:
        #     dirPoint = adsk.core.Point3D.create(-1, 0, 0)
        #     offsetCurves = sketch_living_hinge.offset(curves, dirPoint, 0.05)
        
        # # join_line_1 = lines_sketch_living_hinge.addByTwoPoints( offsetCurves.endSketchPoint , curves.endSketchPoint )

        # # Extract individual curves
        # offset_last = offsetCurves.item(offsetCurves.count - 1)
        # original_last = curves.item(curves.count - 1)

        # # Create the joining line
        # join_line_1 = lines_sketch_living_hinge.addByTwoPoints(
        #     offset_last.endSketchPoint,
        #     original_last.endSketchPoint
        # )

        # # Extract first curves
        # offset_first = offsetCurves.item(0)
        # original_first = curves.item(0)

        # # Create the joining line between the starting points
        # join_line_2 = lines_sketch_living_hinge.addByTwoPoints(
        #     offset_first.endSketchPoint,
        #     original_first.endSketchPoint
        # )

        # #constraints_sketch_living_hinge.addLineOnPlanarSurface(line1, highlight_body.faces.item(0))
        # #constraints_sketch_living_hinge.addLineOnPlanarSurface(line1, plane_living_hinge )
        # #constraints_sketch_living_hinge.addLineOnPlanarSurface(line2, highlight_body.faces.item(1))

        # # prof = sketch_living_hinge.profiles.item(1)
        # # prof2 = sketch_living_hinge.profiles.item(0)

        # # extrudes = rootComp.features.extrudeFeatures
        # # distance = adsk.core.ValueInput.createByReal(-1.5)
        # # extrude1 = extrudes.addSimple(prof, distance, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)

        # # # extrudeInput = extrudes.createInput(profVertical, adsk.fusion.FeatureOperations.CutFeatureOperation)

        # # extrude2 = extrudes.addSimple(prof2, distance, adsk.fusion.FeatureOperations.CutFeatureOperation)

        # # ui.messageBox(msg)

    except: 
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))








