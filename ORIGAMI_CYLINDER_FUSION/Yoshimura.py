import adsk.core, adsk.fusion, traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        
    
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        
        sides = 18

        # === Define User Parameters ===
        # Circle Diameter
        paramName = 'circleDiameter'
        paramValue = 10
        existingParam = design.userParameters.itemByName(paramName)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName, adsk.core.ValueInput.createByReal(paramValue), 'cm', 'Diameter of polygon bounding circle')

        # First polygon angle
        angleParamName = 'angleDeg'
        angleValue = 0
        existingParam = design.userParameters.itemByName(angleParamName)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(angleParamName, adsk.core.ValueInput.createByString(f'{angleValue} deg'), 'deg', 'Rotation angle for first polygon')

        # Second polygon angle
        angleParamName2 = 'angleDeg2'
        angleValue2 = 360 / ( 2 * sides )
        existingParam = design.userParameters.itemByName(angleParamName2)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(angleParamName2, adsk.core.ValueInput.createByString(f'{angleValue2} deg'), 'deg', 'Rotation angle for second polygon')
        

        # === First Sketch on XY Plane ===
        sketches = rootComp.sketches
        xyPlane = rootComp.xYConstructionPlane
        sketch = sketches.add(xyPlane)
        
        radius = 5
        origin = adsk.core.Point3D.create(0, 0, 0)
        
        circle = sketch.sketchCurves.sketchCircles.addByCenterRadius(origin, radius)
        circle.isConstruction = True
        
        constraints = sketch.geometricConstraints
        centerPoint = circle.centerSketchPoint
        constraints.addCoincident(centerPoint, sketch.originPoint)
        
        polygon = sketch.sketchCurves.sketchLines.addScribedPolygon(centerPoint, sides, 0, radius, True)
        
        constraints.addCoincident(polygon[0].startSketchPoint, circle)
        constraints.addCoincident(polygon[1].startSketchPoint, circle)
        constraints.addCoincident(polygon[2].startSketchPoint, circle)
        
        lines = sketch.sketchCurves.sketchLines
        PointOnHorizon = adsk.core.Point3D.create(-10, 0, 0)
        HorizontalLine = lines.addByTwoPoints(sketch.originPoint, PointOnHorizon)
        HorizontalLine.isConstruction = True
        constraints.addHorizontal(HorizontalLine)
        
        line1 = lines.addByTwoPoints(polygon[0].startSketchPoint, sketch.originPoint)
        line1.isConstruction = True
        
        dimensions = sketch.sketchDimensions
        diameterDim = dimensions.addDiameterDimension(circle, PointOnHorizon)
        diameterDim.parameter.expression = paramName
        
        angularDim1 = dimensions.addAngularDimension(HorizontalLine, line1, adsk.core.Point3D.create(3, 6, 0), True)
        angularDim1.parameter.expression = angleParamName

        sketchPoints = sketch.sketchPoints
        points1 = []
        
        constructionPoints = rootComp.constructionPoints
        pointInput = constructionPoints.createInput()
                    
        construction_point_bottom = []
        construction_point_up = []
        
        for i in range(sides):
            p = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
            pt = sketchPoints.add(p)
            constraints.addCoincident(polygon[i].startSketchPoint, pt)
            points1.append(pt)
            pointInput.setByPoint(points1[i])
            construction_point_bottom.append( constructionPoints.add( pointInput ) )

    
        # === Offset Plane and Second Sketch ===
        plane_offset = 5
        planes = rootComp.constructionPlanes
        offsetValue = adsk.core.ValueInput.createByReal(plane_offset)
        planeInput = planes.createInput()
        planeInput.setByOffset(xyPlane, offsetValue)
        offsetPlane = planes.add(planeInput)

        sketch2 = sketches.add(offsetPlane)
        constraints2 = sketch2.geometricConstraints
        circles2 = sketch2.sketchCurves.sketchCircles
        lines2 = sketch2.sketchCurves.sketchLines
        dimensions2 = sketch2.sketchDimensions
        
        circle2 = circles2.addByCenterRadius(origin, radius)
        circle2.isConstruction = True
        centerPoint2 = circle2.centerSketchPoint
        constraints2.addCoincident(centerPoint2, sketch2.originPoint)
        
        polygon2 = sketch2.sketchCurves.sketchLines.addScribedPolygon(centerPoint2, sides, 0, radius, True)
        constraints2.addCoincident(polygon2[0].startSketchPoint, circle2)
        constraints2.addCoincident(polygon2[1].startSketchPoint, circle2)
        constraints2.addCoincident(polygon2[2].startSketchPoint, circle2)
        
        line2 = lines2.addByTwoPoints(polygon2[0].startSketchPoint, sketch2.originPoint)
        line2.isConstruction = True
        
        diameterDim2 = dimensions2.addDiameterDimension(circle2, adsk.core.Point3D.create(radius, 0, 0))
        diameterDim2.parameter.expression = paramName
        
        PointOnHorizon2 = adsk.core.Point3D.create(-10, 0, 0)
        HorizontalLine2 = lines2.addByTwoPoints(sketch2.originPoint, PointOnHorizon2)
        HorizontalLine2.isConstruction = True
        constraints2.addHorizontal(HorizontalLine2)
        
        angularDim2 = dimensions2.addAngularDimension(HorizontalLine2, line2, PointOnHorizon2, True)
        angularDim2.parameter.expression = angleParamName2

        sketchPoints2 = sketch2.sketchPoints
        points2 = []
        for i in range(sides):
            p = adsk.core.Point3D.create(1, 2, 0)  # Example placeholder coordinates
            pt = sketchPoints2.add(p)
            constraints2.addCoincident(polygon2[i].startSketchPoint, pt)
            points2.append(pt)
            pointInput.setByPoint( points2[i] )
            construction_point_up.append( constructionPoints.add( pointInput ) )

        # === Construction Planes Between Points ===
        facets_bottom = []
        facets_top = []
        
        planeInput = planes.createInput()
        
        i = 0 ;
                    
        for i in range(sides - 1 ):
            
            planeInput.setByThreePoints(points1[i], points2[i], points2[i + 1])
            facets_bottom.append(planes.add(planeInput))

            planeInput.setByThreePoints(points1[i], points1[i + 1], points2[i+1])
            facets_top.append(planes.add(planeInput))

        
        # === Sketch on First Bottom Facet ===
        
        param_facet_2H_name = 'x_distance_point2'
        paramValue2H = -1.0
        existingParam = design.userParameters.itemByName(param_facet_2H_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_2H_name, adsk.core.ValueInput.createByReal(paramValue2H), 'cm', 'Distance in X of point1')
        
        param_facet_2V_name = 'y_distance_point2'
        paramValue2V = -1.0
        existingParam = design.userParameters.itemByName(param_facet_2V_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_2V_name, adsk.core.ValueInput.createByReal(paramValue2V), 'cm', 'Distance in X of point1')
        
        
        param_facet_3H_name = 'x_distance_point3'
        paramValue3H = 2.0
        existingParam = design.userParameters.itemByName(param_facet_3H_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_3H_name, adsk.core.ValueInput.createByReal(paramValue3H), 'cm', 'Distance in X of point2')
        
        param_facet_3V_name = 'y_distance_point3'
        paramValue3V = 2.0
        existingParam = design.userParameters.itemByName(param_facet_3V_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_3V_name, adsk.core.ValueInput.createByReal(paramValue3V), 'cm', 'Distance in X of point2')
        
        param_facet_4H_name = 'x_distance_point4'
        paramValue4H = -2.0
        existingParam = design.userParameters.itemByName(param_facet_4H_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_4H_name, adsk.core.ValueInput.createByReal(paramValue4H), 'cm', 'Distance in X of point3')
        
        param_facet_4V_name = 'y_distance_point4'
        paramValue4V = 2.0
        existingParam = design.userParameters.itemByName(param_facet_4V_name)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(param_facet_4V_name, adsk.core.ValueInput.createByReal(paramValue4V), 'cm', 'Distance in Y of point3')
        

        sketch_facet_bottom = []
        sketch_facet_bottom_constraints = []
        sketch_facet_bottom_dimensions = []
        sketch_facet_bottom_points = []
        sketch_facet_bottom_lines = []
        prof_bottom = []
        extrude_bottom = []

        sketch_facet_up = []
        sketch_facet_up_constraints = []
        sketch_facet_up_dimensions = []
        sketch_facet_up_points = []
        sketch_facet_up_lines = []
        prof_up = []
        extrude_up = []
        
        projectedEntities = []

        sides_test = sides-1
        
        extrudes = rootComp.features.extrudeFeatures
                    
        for i in range(sides_test):
            
            # === Sketch on bottom facet ===
            sketch_bottom = sketches.add(facets_bottom[i])
            sketch_facet_bottom.append(sketch_bottom)
            sketch_facet_bottom_constraints.append(sketch_bottom.geometricConstraints)
            sketch_facet_bottom_dimensions.append(sketch_bottom.sketchDimensions)
            sketch_facet_bottom_points.append(sketch_bottom.sketchPoints)
            sketch_facet_bottom_lines.append( sketch_bottom.sketchCurves.sketchLines )
            
            sketch_up = sketches.add(facets_top[i])
            sketch_facet_up.append(sketch_up)
            sketch_facet_up_constraints.append( sketch_up.geometricConstraints )
            sketch_facet_up_dimensions.append( sketch_up.sketchDimensions )
            sketch_facet_up_points.append( sketch_up.sketchPoints )
            sketch_facet_up_lines.append( sketch_up.sketchCurves.sketchLines )


            # Add points
            point1 = adsk.core.Point3D.create(1, 1, 0)
            point2 = adsk.core.Point3D.create(2, 1, 0)
            point3 = adsk.core.Point3D.create(3, 1, 0)
            point4 = adsk.core.Point3D.create(3,3,0)
            point5 = adsk.core.Point3D.create(3,3,0)
            point6 = adsk.core.Point3D.create(3,3,0)
            pt4 = sketch_facet_bottom_points[i].add( point4 )
            pt5 = sketch_facet_bottom_points[i].add( point5 )
            pt6 = sketch_facet_bottom_points[i].add( point6 )
            
            sketch_facet_bottom_constraints[i].addCoincident( sketch_bottom.originPoint, pt4 )
            
            pt1 = sketch_facet_bottom_points[i].add(point1)
            pt2 = sketch_facet_bottom_points[i].add(point2)
            pt3 = sketch_facet_bottom_points[i].add(point3)
            
            Hline = sketch_facet_bottom_lines[i].addByTwoPoints( pt4 , pt5 )
            sketch_facet_bottom_constraints[i].addHorizontal( Hline )
            Hline.isConstruction = True
            
            Vline = sketch_facet_bottom_lines[i].addByTwoPoints( pt4 , pt6 )
            sketch_facet_bottom_constraints[i].addVertical( Vline )
            Vline.isConstruction = True
            
                        
            # Add dimensions
        
            # dimension2H = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt1, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.HorizontalDimensionOrientation, adsk.core.Point3D.create(1, 1, 0) )
            
            # dimension2H.parameter.expression = param_facet_2H_name
            
            # dimension2V = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt1, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.VerticalDimensionOrientation, adsk.core.Point3D.create(1, 1, 0) )
            # dimension2V.parameter.expression = param_facet_2V_name

            # dimension3H = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt2, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.HorizontalDimensionOrientation, adsk.core.Point3D.create(2, 1, 0) )
            # dimension3H.parameter.expression = param_facet_3H_name
            
            # dimension3V = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt2, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.VerticalDimensionOrientation, adsk.core.Point3D.create(2, 1, 0) )
            # dimension3V.parameter.expression = param_facet_3V_name
            
            # dimension4H = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt3, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.HorizontalDimensionOrientation, adsk.core.Point3D.create(2, 1, 0) )
            # dimension4H.parameter.expression = param_facet_4H_name
            
            # dimension4V = sketch_facet_bottom_dimensions[i].addDistanceDimension( pt3, sketch_facet_bottom[i].originPoint, adsk.fusion.DimensionOrientations.VerticalDimensionOrientation, adsk.core.Point3D.create(2, 1, 0) )
            # dimension4V.parameter.expression = param_facet_4V_name
            
            
            #projectedEntities.append( sketch_bottom.project( construction_point_bottom[i] ).item(0) )
            #projectedEntities.append( sketch_bottom.project( construction_point_up[i] ).item(0) )
            #projectedEntities.append( sketch_bottom.project( construction_point_up[i+1] ).item(0) )
            
            #sketch_facet_bottom_constraints[i].addCoincident( pt1 , projectedEntities[i] )
            #sketch_facet_bottom_constraints[i].addCoincident( pt2 , projectedEntities[i+1] )
            #sketch_facet_bottom_constraints[i].addCoincident( pt3 , projectedEntities[i+2] )
            
            # planeInput.setByThreePoints(points1[i], points2[i], points2[i + 1])
            # facets_bottom.append(planes.add(planeInput))

            # planeInput.setByThreePoints(points1[i], points1[i + 1], points2[i])
            # facets_top.append(planes.add(planeInput))

            
            projectedEntity1 = sketch_bottom.project( construction_point_bottom[i] )
            projectedEntity2 = sketch_bottom.project( construction_point_up[i] )
            projectedEntity3 = sketch_bottom.project( construction_point_up[i+1] )
            
            sketch_facet_bottom_constraints[i].addCoincident( pt1 , projectedEntity1.item(0) )
            sketch_facet_bottom_constraints[i].addCoincident( pt2 , projectedEntity2.item(0) )
            sketch_facet_bottom_constraints[i].addCoincident( pt3 , projectedEntity3.item(0) )
            
            line1 = sketch_facet_bottom_lines[i].addByTwoPoints( pt1 , pt2 )
            line2 = sketch_facet_bottom_lines[i].addByTwoPoints( pt2 , pt3 )
            line3 = sketch_facet_bottom_lines[i].addByTwoPoints( pt3, pt1 )
            
            prof_bottom.append( sketch_facet_bottom[i].profiles.item(0) )
            
            # create an extrusion input
            extInput = extrudes.createInput( prof_bottom[i] ,
            adsk.fusion.FeatureOperations.NewBodyFeatureOperation )
            
            distance = adsk.core.ValueInput.createByReal(0.05)
            extInput.setDistanceExtent( True , distance )
            extInput.isSolid = True
            extrude_bottom.append( extrudes.add( extInput ) )
            
            # Do the same on the upper plane
            
            point7 = adsk.core.Point3D.create(1, 1, 0)
            point8 = adsk.core.Point3D.create(2, 1, 0)
            point9 = adsk.core.Point3D.create(3, 1, 0)
            point10 = adsk.core.Point3D.create(3,3,0)
            point11 = adsk.core.Point3D.create(3,3,0)
            point12 = adsk.core.Point3D.create(3,3,0)
            pt10 = sketch_facet_up_points[i].add( point4 )
            pt11 = sketch_facet_up_points[i].add( point5 )
            pt12 = sketch_facet_up_points[i].add( point6 )

            sketch_facet_up_constraints[i].addCoincident( sketch_up.originPoint, pt10 )

            pt7 = sketch_facet_up_points[i].add(point1)
            pt8 = sketch_facet_up_points[i].add(point2)
            pt9 = sketch_facet_up_points[i].add(point3)


            projectedEntity4 = sketch_up.project( construction_point_bottom[i] )
            projectedEntity5 = sketch_up.project( construction_point_bottom[i+1] )
            projectedEntity6 = sketch_up.project( construction_point_up[i+1] )

            sketch_facet_up_constraints[i].addCoincident( pt7 , projectedEntity4.item(0) )
            sketch_facet_up_constraints[i].addCoincident( pt8 , projectedEntity5.item(0) )
            sketch_facet_up_constraints[i].addCoincident( pt9 , projectedEntity6.item(0) )

            line4 = sketch_facet_up_lines[i].addByTwoPoints( pt7 , pt8 )
            line5 = sketch_facet_up_lines[i].addByTwoPoints( pt8 , pt9 )
            line6 = sketch_facet_up_lines[i].addByTwoPoints( pt9, pt7 )

            prof_up.append( sketch_facet_up[i].profiles.item(0) )

            # create an extrusion input
            extInput = extrudes.createInput( prof_up[i] , adsk.fusion.FeatureOperations.NewBodyFeatureOperation )

            distance = adsk.core.ValueInput.createByReal(0.1)
            extInput.setDistanceExtent( True , distance )
            extInput.isSolid = True
            extrude_up.append( extrudes.add( extInput ) )
      
            # === Sketch on top facet ===
 
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))



