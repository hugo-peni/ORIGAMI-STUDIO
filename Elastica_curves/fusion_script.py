
import adsk.core
import adsk.fusion
import traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        sketches = rootComp.sketches

        # START CODE GENERATION SECTION DO NOT REMOVE COMMENTS

        sketch0 = sketches.add(rootComp.xYConstructionPlane)
        sketch0.is3D = True
        pts0 = adsk.core.ObjectCollection.create()
        pts0.add(adsk.core.Point3D.create(1.9378248434212895, -0.4948079185090461, 0.0))
        pts0.add(adsk.core.Point3D.create(1.9436231366470833, -0.47151787101885445, 0.02400111309245556))
        pts0.add(adsk.core.Point3D.create(1.9576048924295575, -0.40961333613038114, 0.08746763134679503))
        pts0.add(adsk.core.Point3D.create(1.969615506024416, -0.3472963553338607, 0.15093414960113405))
        pts0.add(adsk.core.Point3D.create(1.9796428837618656, -0.2846296765465701, 0.21440066785547351))
        pts0.add(adsk.core.Point3D.create(1.9876769289225082, -0.22167639980202208, 0.27786718610981254))
        pts0.add(adsk.core.Point3D.create(1.9937095519038848, -0.15849991371357677, 0.341333704364152))
        pts0.add(adsk.core.Point3D.create(1.997734678366016, -0.09516383164748424, 0.40480022261849147))
        pts0.add(adsk.core.Point3D.create(1.9997482553477501, -0.03173192766961574, 0.4682667408728305))
        pts0.add(adsk.core.Point3D.create(2.0, 5.551115123125784e-17, 0.5000000000000002))
        spline0 = sketch0.sketchCurves.sketchFittedSplines.add(pts0)


        sketch1 = sketches.add(rootComp.xYConstructionPlane)
        sketch1.is3D = True
        pts1 = adsk.core.ObjectCollection.create()
        pts1.add(adsk.core.Point3D.create(2.0, 5.551115123125784e-17, -0.5000000000000002))
        pts1.add(adsk.core.Point3D.create(1.9997482553477501, 0.03173192766961618, -0.46826674087283005))
        pts1.add(adsk.core.Point3D.create(1.997734678366016, 0.09516383164748468, -0.404800222618491))
        pts1.add(adsk.core.Point3D.create(1.9937095519038848, 0.15849991371357722, -0.34133370436415156))
        pts1.add(adsk.core.Point3D.create(1.9876769289225082, 0.22167639980202208, -0.27786718610981254))
        pts1.add(adsk.core.Point3D.create(1.9796428837618654, 0.28462967654657056, -0.21440066785547307))
        pts1.add(adsk.core.Point3D.create(1.969615506024416, 0.3472963553338607, -0.15093414960113405))
        pts1.add(adsk.core.Point3D.create(1.9576048924295575, 0.4096133361303816, -0.08746763134679458))
        pts1.add(adsk.core.Point3D.create(1.9436231366470833, 0.47151787101885445, -0.02400111309245556))
        pts1.add(adsk.core.Point3D.create(1.9378248434212895, 0.4948079185090462, 0.0))
        spline1 = sketch1.sketchCurves.sketchFittedSplines.add(pts1)


        sketch2 = sketches.add(rootComp.xYConstructionPlane)
        sketch2.is3D = True
        pts2 = adsk.core.ObjectCollection.create()
        pts2.add(adsk.core.Point3D.create(1.9378248434212895, -0.4948079185090461, 0.0))
        pts2.add(adsk.core.Point3D.create(1.9436231366470833, -0.47151787101885445, -0.02400111309245556))
        pts2.add(adsk.core.Point3D.create(1.9576048924295575, -0.40961333613038114, -0.08746763134679503))
        pts2.add(adsk.core.Point3D.create(1.969615506024416, -0.3472963553338607, -0.15093414960113405))
        pts2.add(adsk.core.Point3D.create(1.9796428837618656, -0.2846296765465701, -0.21440066785547351))
        pts2.add(adsk.core.Point3D.create(1.9876769289225082, -0.22167639980202208, -0.27786718610981254))
        pts2.add(adsk.core.Point3D.create(1.9937095519038848, -0.15849991371357677, -0.341333704364152))
        pts2.add(adsk.core.Point3D.create(1.997734678366016, -0.09516383164748424, -0.40480022261849147))
        pts2.add(adsk.core.Point3D.create(1.9997482553477501, -0.03173192766961574, -0.4682667408728305))
        pts2.add(adsk.core.Point3D.create(2.0, 5.551115123125784e-17, -0.5000000000000002))
        spline2 = sketch2.sketchCurves.sketchFittedSplines.add(pts2)


        sketch3 = sketches.add(rootComp.xYConstructionPlane)
        sketch3.is3D = True
        pts3 = adsk.core.ObjectCollection.create()
        pts3.add(adsk.core.Point3D.create(2.0, 5.551115123125784e-17, 0.5000000000000002))
        pts3.add(adsk.core.Point3D.create(1.9997482553477501, 0.03173192766961618, 0.46826674087283005))
        pts3.add(adsk.core.Point3D.create(1.997734678366016, 0.09516383164748468, 0.404800222618491))
        pts3.add(adsk.core.Point3D.create(1.9937095519038848, 0.15849991371357722, 0.34133370436415156))
        pts3.add(adsk.core.Point3D.create(1.9876769289225082, 0.22167639980202208, 0.27786718610981254))
        pts3.add(adsk.core.Point3D.create(1.9796428837618654, 0.28462967654657056, 0.21440066785547307))
        pts3.add(adsk.core.Point3D.create(1.969615506024416, 0.3472963553338607, 0.15093414960113405))
        pts3.add(adsk.core.Point3D.create(1.9576048924295575, 0.4096133361303816, 0.08746763134679458))
        pts3.add(adsk.core.Point3D.create(1.9436231366470833, 0.47151787101885445, 0.02400111309245556))
        pts3.add(adsk.core.Point3D.create(1.9378248434212895, 0.4948079185090462, 0.0))
        spline3 = sketch3.sketchCurves.sketchFittedSplines.add(pts3)




        pathCurves3 = adsk.core.ObjectCollection.create()
        pathCurves3.add(spline0)
        pathCurves3.add(spline2)
        pathCurves3.add(spline1)
        pathCurves3.add(spline3)


        patches = rootComp.features.patchFeatures
        patchInput = patches.createInput(pathCurves3, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        patch3 = patches.add(patchInput)

        # END CODE GENERATION SECTION DO NOT REMOVE COMMENTS

        # create one surface : 
        # surfaces = adsk.core.ObjectCollection.create()

        # for i in range(len(rootComp.bRepBodies)):
        #     surface = rootComp.bRepBodies.item(i)
        #     surfaces.add(surface)
        
        # tolerance = adsk.core.ValueInput.createByReal(1.0)
        # features = rootComp.features

        # stitches = features.stitchFeatures
        # stitchInput = stitches.createInput(surfaces, tolerance, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        
        # # Create a stitch feature.
        # stitch = stitches.add(stitchInput)

        # thickenFeatures = rootComp.features.thickenFeatures
        # inputSurfaces = adsk.core.ObjectCollection.create()
        # bodies = rootComp.bRepBodies.item(0)
        # inputSurfaces.add(bodies)
        # thickness = adsk.core.ValueInput.createByReal(0.05)
        # thickenInput = thickenFeatures.createInput( inputSurfaces , thickness , False , adsk.fusion.FeatureOperations.NewBodyFeatureOperation )
        # thickenFeatures.add(thickenInput)


    except:
        if ui:
            ui.messageBox('Script failed:\n{}'.format(traceback.format_exc()))


