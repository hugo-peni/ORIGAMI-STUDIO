import adsk.core, adsk.fusion, traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface

        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent

        offset_dim = 0.05

        paramName1 = 'circleDiameter1'
        paramValue1 = 0.35
        existingParam = design.userParameters.itemByName(paramName1)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName1, adsk.core.ValueInput.createByReal(paramValue1), 'cm', 'Diameter of polygon bounding circle')

        paramName2 = 'circleDiameter2'
        paramValue2 = paramValue1 - ( offset_dim * 2 )
        existingParam = design.userParameters.itemByName(paramName2)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName2, adsk.core.ValueInput.createByReal(paramValue2), 'cm', 'Diameter of polygon bounding circle')

        paramName3 = 'circleDiameter3'
        paramValue3 = 0.1
        existingParam = design.userParameters.itemByName(paramName3)
        if existingParam:
            existingParam.deleteMe()
        design.userParameters.add(paramName3, adsk.core.ValueInput.createByReal(paramValue3), 'cm', 'Diameter of polygon bounding circle')

        # === First Sketch on XY Plane ===
        sketches = rootComp.sketches
        xyPlane = rootComp.xYConstructionPlane
        sketch1 = sketches.add(xyPlane)

        constraints1 = sketch1.geometricConstraints
        dimensions1 = sketch1.sketchDimensions
        lines1 = sketch1.sketchCurves.sketchLines
        arcs1 = sketch1.sketchCurves.sketchArcs

        # Reference horizontal construction line from origin
        hozRef = lines1.addByTwoPoints(sketch1.originPoint, adsk.core.Point3D.create(1, 0, 0))
        constraints1.addHorizontal(hozRef)
        hozRef.isConstruction = True

        # Reference vertical construction line from origin
        verRef = lines1.addByTwoPoints(sketch1.originPoint, adsk.core.Point3D.create(0, 1, 0))
        constraints1.addVertical(verRef)
        verRef.isConstruction = True

        # Add symmetric angled lines
        line1 = lines1.addByTwoPoints(adsk.core.Point3D.create(0.5, 1, 0), adsk.core.Point3D.create(0.7, 0.3, 0))
        line2 = lines1.addByTwoPoints(line1.endSketchPoint, adsk.core.Point3D.create(1, 0.2, 0))
        line4 = lines1.addByTwoPoints(adsk.core.Point3D.create(0.5, -1, 0), adsk.core.Point3D.create(0.7, -0.3, 0))
        line3 = lines1.addByTwoPoints(line4.endSketchPoint, adsk.core.Point3D.create(1, -0.2, 0))

        constraints1.addSymmetry(line1, line4, hozRef)
        constraints1.addSymmetry(line2, line3, hozRef)

        # Add arc connecting line2 to line3
        arc = arcs1.addByCenterStartSweep(
            adsk.core.Point3D.create(1, 0, 0),
            line2.endSketchPoint,
            -90
        )
        constraints1.addTangent(arc, line2)
        constraints1.addTangent(arc, line3)
        constraints1.addCoincident( arc.startSketchPoint , line3.endSketchPoint )

        # Add fillets
        fillet1 = arcs1.addFillet(
            line1, line1.endSketchPoint.geometry,
            line2, line2.startSketchPoint.geometry,
            0.05
        )
        fillet2 = arcs1.addFillet(
            line4, line4.endSketchPoint.geometry,
            line3, line3.startSketchPoint.geometry,
            0.05
        )

        # Find all connected curves from line4 and line1 and merge them into one collection
        allConnected = adsk.core.ObjectCollection.create()
        connected1 = sketch1.findConnectedCurves(line4)
        for c in connected1:
            allConnected.add(c)
        connected2 = sketch1.findConnectedCurves(line1)
        for c in connected2:
            allConnected.add(c)

        if allConnected.count > 0:
            dirPoint = adsk.core.Point3D.create(-1, 0, 0)
            offsetCurves = sketch1.offset(allConnected, dirPoint, offset_dim )

            # Collect the 4 candidate points from the first and last offset curves
            offset_first = offsetCurves.item(0)
            offset_last = offsetCurves.item(offsetCurves.count - 1)

            candidates = [
                offset_first.startSketchPoint,
                offset_first.endSketchPoint,
                offset_last.startSketchPoint,
                offset_last.endSketchPoint
            ]

            # Helper to find the closest point from a list to a reference point
            def find_closest(refPoint, candidates):
                return min(candidates, key=lambda pt: pt.geometry.distanceTo(refPoint.geometry))

            # === Find and connect closest to line1.startSketchPoint ===
            closest_to_line1 = find_closest(line1.startSketchPoint, candidates)
            join_line_1 = lines1.addByTwoPoints(closest_to_line1, line1.startSketchPoint)

            # === Remove that point from candidates to avoid reusing it ===
            candidates.remove(closest_to_line1)

            # === Find and connect closest to line4.startSketchPoint ===
            closest_to_line4 = find_closest(line4.startSketchPoint, candidates)
            join_line_2 = lines1.addByTwoPoints(closest_to_line4, line4.startSketchPoint)

        # Reference horizontal construction line from origin
        hozRef_2 = lines1.addByTwoPoints( arc.centerSketchPoint , adsk.core.Point3D.create(2, 0, 0) )
        constraints1.addHorizontal(hozRef_2)
        hozRef_2.isConstruction = True

        # Reference vertical construction line from origin
        verRef_2 = lines1.addByTwoPoints( arc.centerSketchPoint, adsk.core.Point3D.create(1, 2, 0) )
        constraints1.addVertical(verRef_2)
        verRef_2.isConstruction = True

        # Reference vertical construction line from origin
        verRef_3 = lines1.addByTwoPoints( adsk.core.Point3D.create(0.5, 0, 0) , adsk.core.Point3D.create(0.5, 2, 0) )
        constraints1.addVertical(verRef_3)
        verRef_3.isConstruction = True

        arc2 = arcs1.addByCenterStartSweep(
            arc.centerSketchPoint,
            adsk.core.Point3D.create(1, 0.1, 0),
            -90
        )

        line5 = lines1.addByTwoPoints( arc2.endSketchPoint , adsk.core.Point3D.create(0, 0.2, 0) )
        line6 = lines1.addByTwoPoints( arc2.startSketchPoint , adsk.core.Point3D.create(0, -0.2, 0) )
        constraints1.addParallel( line5 , line2 )
        constraints1.addParallel( line6 , line3 )

        constraints1.addSymmetry( arc2.endSketchPoint , arc2.startSketchPoint , hozRef_2 )

        constraints1.addTangent( arc2 , line5 )
        constraints1.addTangent( arc2 , line6 )

        constraints1.addCoincident( line5.endSketchPoint , verRef_3 )
        constraints1.addCoincident( line6.endSketchPoint , verRef_3 )
    

        line7 = lines1.addByTwoPoints( line5.endSketchPoint , line6.endSketchPoint )

        diameterDim1 = dimensions1.addRadialDimension( arc , adsk.core.Point3D.create( 4 , 0 , 0 ) )
        diameterDim1.parameter.expression = paramName1

        diameterDim2 = dimensions1.addRadialDimension( arc2 , adsk.core.Point3D.create( 3 , 0 , 0 ) )
        diameterDim2.parameter.expression = paramName2

        diameterDim3 = dimensions1.addRadialDimension( fillet1 , adsk.core.Point3D.create( 3 , 0 , 0 ) )
        diameterDim3.parameter.expression = paramName3

        constraints1.addCoincident( arc2.centerSketchPoint , arc.centerSketchPoint )

        constraints1.addSymmetry( line1.startSketchPoint , line4.startSketchPoint , hozRef )
        constraints1.addSymmetry( fillet1 , fillet2 , hozRef )


    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))