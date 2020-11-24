"""
    ###############################
    ## Subglacial Overdeepenings ##
    ###############################

    This is the main code file of the first work package of the project 'Future glacial lakes
    in High Mountain Asia - Modeling and Risk Analysis' (GLAMoR).
    The script was written by Wilhelm Furian.

    With this script it is possible to analyze the subglacial bedrock topography of
    all glaciers in a given RGI region and to find potential future glacial lakes.
    It also allows for the investigation of the morphology of surrounding mountain slopes
    in order to quantify the potential hazard of mass movement impacts into the future lakes.

    In order to run properly, this script requires different datasets to be downloaded:
    - glacier ice thickness data in raster format. For High Mountain Asia (HMA),
      we recommend the data provided by Farinotti et al. (2019)
      (https://doi.org/10.1038/s41561-019-0300-3)
    - glacier outlines in shapefile format, downloadable from the Randolph Glacier Inventory, v6
      (https://www.glims.org/RGI/)
    - a DEM raster for the whole study area. For HMA, we recommend using the ALOS World 3D - 30m (AW3D30),
      available at https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm.

    The model is written in Python 2.7 and has been tested with ArcGIS 10.7 and PyCharm CE 2019.3.1.
    It needs the following ArcGIS extensions to be enabled:
    3D Analyst, Spatial Analyst and Geostatistical Analyst.

    The code is available on github. https://github.com/cryotools/subglacial-overdeepenings
    For more information on this work package see the README.
    For more information on the project as a whole see https://hu-berlin.de/glamor.

    You are allowed to use and modify this code in a noncommercial manner and by
    appropriately citing the above mentioned developer. If you would like to share your own improvements,
    please fork the repository on GitHub, commit your changes, and file a merge request.

    Correspondence: furiawil@geo.hu-berlin.de
"""

# imports
import os
import os.path
import arcpy
from arcpy.sa import *
import shutil
import gc
import re

# checkout the necessary extensions for arcpy
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("GeoStats")
arcpy.env.overwriteOutput = True

# prepare input paths
# (at the moment, this has to be done manually for every RGI region)
input_path = "~"                            # path to folder with all the preprocessed glacier thickness data
input_path_raster = "~"                     # path to final glacier thickness rasters
RGI_path = "~/14_rgi60_SouthAsiaWest.shp"   # absolute path to the RGI shapefile of the current region
work_path = "~"                             # where should the files be saved

# First Step: Find overdeepenings
for raster in os.listdir(input_path_raster):
    if raster.startswith("RGI60") and raster.endswith(".tif"):
        glacierID = raster[:-14]
        glacierLocation = input_path_raster + "/" + raster
        print glacierID

        # define several work folder connections
        current_path = work_path + "/" + glacierID
        subfolder_SHP = current_path + "/" + "single_sink_shapes"
        subfolder_RGI = current_path + "/single_RGI_shapes"
        subfolder_TIN = current_path + "/" + "single_TINs"
        subfolder_TIF = current_path + "/intersecting_TIFs"
        subfolder_TIF_projected = work_path + "/intersecting_TIFs_projected"
        subfolder_SINKS = current_path + "/single_sinks"

        print "\n"
        print "=============================================================================="
        print "================== Working on raster " + glacierID + " =========================="
        print "=============================================================================="
        print "\n"

        # create folder for current shapefile
        try:
            os.mkdir(current_path)
        except OSError:
            print("Creation of the directory %s failed" % current_path)

        try:
            os.mkdir(subfolder_SHP)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_SHP)

        print "extracting single glacier-shp from RGI-shapes"
        where_clause = str(glacierID[0:8] + "." + glacierID[9:])
        arcpy.Select_analysis(RGI_path, current_path + "/glacier.shp",
                              where_clause=' "RGIId" = \'{}\''.format(where_clause))

        # exclude glaciers with less than 1km2 area
        with arcpy.da.SearchCursor(current_path + "/glacier.shp", "Area") as cursor:
            for row in cursor:
                area = row[0]
        del cursor
        if area < 1.0:
            print "Glacier " + glacierID + " is to small (" + str(row[0]) + " km2). Skipping iteration!"
            shutil.rmtree(current_path)
            continue

        print "buffer the single shapefile for intersect program"
        arcpy.Buffer_analysis(in_features=current_path + "/glacier.shp",
                              out_feature_class=current_path + "/glacier_buffer.shp",
                              buffer_distance_or_field="100 Meters")

        print "find intersecting glaciers"
        arcpy.Intersect_analysis(in_features=[RGI_path, current_path + "/glacier_buffer.shp"],
                                 out_feature_class=current_path + "/intersect",
                                 join_attributes="ALL")

        print "select corresponding shapefiles and save them to single_RGI_shapes"
        try:
            os.mkdir(subfolder_RGI)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_RGI)

        with arcpy.da.SearchCursor(current_path + "/intersect.shp", "RGIId") as cursor:
            for row in cursor:
                print row
                output_name = str(str(row)[9:11] + "_" + str(row)[12:17])
                output = subfolder_RGI + "/{}.shp".format(output_name)
                where_clause = str(row)[2:17]
                arcpy.Select_analysis(RGI_path, output,
                                      where_clause=' "RGIId" = {}\''.format(where_clause))

        print "merging the intersecting shapefiles"
        arcpy.env.workspace = subfolder_RGI
        try:
            for i in range(len(os.listdir(subfolder_RGI))):
                shp_paths = []
                list_shps = arcpy.ListFeatureClasses("*{}*".format(".shp"))
                if list_shps:
                    for shp in list_shps:
                        desc = arcpy.Describe(shp)
                        shppath = desc.catalogPath
                        shp_paths.append(shppath)
                    arcpy.Merge_management(shp_paths, "merged.shp")
        except arcpy.ExecuteError:
            print("Merged shapefile should have been created at " + subfolder_RGI)
        arcpy.ResetEnvironments()

        # buffer the multiple shapefiles for ALOS extent
        print "Buffering"
        arcpy.Buffer_analysis(in_features=subfolder_RGI + "/merged.shp",
                              out_feature_class=subfolder_RGI + "/merged_buffer.shp",
                              buffer_distance_or_field="2200 Meters")

        # dissolve the buffers for ALOS extent
        print "Dissolving"
        arcpy.Dissolve_management(in_features=subfolder_RGI + "/merged_buffer.shp",
                                  out_feature_class=subfolder_RGI + "/dissolved_buffer.shp")

        # loading ALOS clipped to buffer extent
        print "Clipping"
        arcpy.Clip_management(in_raster="D:/03_HMA_DEM/HMA_alos-jaxa.tif",
                              out_raster=current_path + "/DEM_area.tif",
                              in_template_dataset=subfolder_RGI + "/merged_buffer.shp",
                              clipping_geometry="ClippingGeometry")

        #####################################################################################

        # create subfolder_TIF in current_path
        try:
            os.mkdir(subfolder_TIF)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_TIF)

        print "copy intersecting tifs and glacier bedrock tifs into subfolder_TIF"
        for filename in os.listdir(subfolder_RGI):
            if filename.startswith("1") and filename.endswith(".shp"):
                RGI_nr = str(filename)[0:2] + "." + str(filename)[3:8]
                for tif in os.listdir(input_path + "/03_substracted"):
                    if str(tif)[6:14] == RGI_nr:
                        arcpy.CopyRaster_management(in_raster=input_path + "/03_substracted/" + tif,
                                                    out_rasterdataset=subfolder_TIF + "/" + str(tif)[6:8] + "_" + str(
                                                        tif)[9:14] + str(tif)[24:28])

        print "project main area DEM raster"
        arcpy.ProjectRaster_management(in_raster=current_path + "/DEM_area.tif",
                                       out_raster=current_path + "/DEM_area_proj.tif",
                                       out_coor_system=glacierLocation,
                                       resampling_type="BILINEAR",
                                       cell_size=32)

        # create subfolder_TIF_projected in work_path
        # (so that every glacier can check there if its needed file is already projected)
        try:
            os.mkdir(subfolder_TIF_projected)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_TIF_projected)

        print "project all intersecting rasters and main glacier raster"
        for filename in os.listdir(subfolder_TIF):
            if filename.endswith(".tif"):
                if os.path.isfile(subfolder_TIF_projected + "/" + filename):
                    print("File exists")
                else:
                    print "Projecting " + filename
                    arcpy.ProjectRaster_management(in_raster=subfolder_TIF + "/" + filename,
                                                   out_raster=subfolder_TIF_projected + "/" + filename,
                                                   out_coor_system=glacierLocation,
                                                   resampling_type="BILINEAR",
                                                   cell_size=32)

        print "clipmerge the projected bedrock rasters into the projected DEM"
        raster_list = []
        for filename in os.listdir(subfolder_TIF):  # get list of rasters to merge
            if filename.endswith(".tif"):
                raster_list.append(subfolder_TIF_projected + "/" + filename)
        raster_list.append(current_path + "/DEM_area_proj.tif")  # append background DEM

        print "get coordinate system from the original projected glacier bedrock file"
        out_coor_system = arcpy.Describe(glacierLocation).spatialReference
        arcpy.MosaicToNewRaster_management(input_rasters=raster_list,
                                           output_location=current_path,
                                           raster_dataset_name_with_extension="area_bedrock.tif",
                                           coordinate_system_for_the_raster=out_coor_system,
                                           pixel_type="32_BIT_FLOAT",
                                           cellsize="32",
                                           number_of_bands="1",
                                           mosaic_method="MINIMUM")

        ####################################################################################################

        print "identify only those slopes that would directly hit the lake"
        outFill_area = Fill(in_surface_raster=current_path + "/area_bedrock.tif")
        outFill_area.save(current_path + "/area_bedrock_filled.tif")
        del outFill_area

        print "create area slope raster"
        outSlope_area = Slope(in_raster=current_path + "/area_bedrock.tif", output_measurement="DEGREE")
        outSlope_area.save(current_path + "/area_slope.tif")
        del outSlope_area

        print "create weight raster for flow accumulation"
        outWeight_area = Con(current_path + "/area_slope.tif", 1, 0, "VALUE < 20")
        outWeight_area.save(current_path + "/weight_raster.tif")
        del outWeight_area

        print "calculate area flowDir raster"
        outFlowDir_area = FlowDirection(current_path + "/area_bedrock_filled.tif")
        outFlowDir_area.save(current_path + "/area_flowDir.tif")
        del outFlowDir_area

        print "calculate flow accumulation with weight raster"
        outFlowAcc_area = FlowAccumulation(in_flow_direction_raster=current_path + "/area_flowDir.tif",
                                           in_weight_raster=current_path + "/weight_raster.tif")
        outFlowAcc_area.save(current_path + "/area_flowAcc.tif")
        del outFlowAcc_area

        print "calculate binary stream raster"
        outStreams_area = Con(current_path + "/area_flowAcc.tif", 1, 0, "VALUE > 100")
        outStreams_area.save(current_path + "/area_streams.tif")
        del outStreams_area

        ################################

        print "find sinks with newly calculated hydrological maps"
        bedrock = Raster(current_path + "/area_bedrock.tif")
        filled = Raster(current_path + "/area_bedrock_filled.tif")
        outSinks_area = Con(filled > bedrock, 1)
        arcpy.Clip_management(in_raster=outSinks_area,
                              out_raster=current_path + "/area_sinks.tif",
                              in_template_dataset=current_path + "/glacier.shp",
                              clipping_geometry="ClippingGeometry")
        del outSinks_area

        n = gc.collect()
        print("Number of unreachable objects collected by GC: ", n)

"""
Next, we need to region group all the overdeepenings in the bed of each glacier.
Normally, this should be possible with arcpy, but there are still unresolved bugs in the region group routine.
Therefore, we recommend to stop this script temporarily and switch to R for the next step.
There, using the clump package will accomplish the same.
"""

# Second Step: Calculate overdeepening properties
glaciersToDo = [name for name in os.listdir(work_path)]
for i in glaciersToDo:
    if i.startswith("RGI"):
        glacierID = i
        glacierLocation = input_path + "/01_original/" + glacierID[:8] + "." + glacierID[9:] + "_thickness.tif"

        # define several work folder connections
        current_path = work_path + "/" + glacierID
        subfolder_SHP = current_path + "/" + "single_sink_shapes"
        subfolder_RGI = current_path + "/single_RGI_shapes"
        subfolder_TIN = current_path + "/" + "single_TINs"
        subfolder_TIF = current_path + "/intersecting_TIFs"
        subfolder_TIF_projected = work_path + "/intersecting_TIFs_projected"
        subfolder_SINKS = current_path + "/single_sinks"

        print "convert to Int-raster and calculate attribute table for R-import"
        outInt = Int(current_path + "/area_sinksRegion.tif")  # this is the file produced by R
        outInt.save(current_path + "/area_sinksRegionInt.tif")
        del outInt
        arcpy.BuildRasterAttributeTable_management(in_raster=current_path + "/area_sinksRegionInt.tif",
                                                   overwrite="OVERWRITE")

        print "== Creating and saving region grouped polygon file =="
        polygonGrp = current_path + "/polygonGrp.shp"
        arcpy.RasterToPolygon_conversion(in_raster=current_path + "/area_sinksRegionInt.tif",
                                         out_polygon_features=polygonGrp,
                                         simplify="NO_SIMPLIFY",
                                         raster_field="VALUE")

        print "== Checking for sinks "
        rows = [row for row in arcpy.da.SearchCursor(polygonGrp, field_names="*")]
        if len(rows) == 0:
            print "======================"
            print "== No  sinks found! =="
            print "== Skipping glacier =="
            print "== deleting folder ==="
            print "======================"
            del polygonGrp
            shutil.rmtree(current_path)
            output = "Glacier folder " + glacierID + " was deleted."
            text = open(work_path + "/glacier_" + glacierID + "_was_deleted.txt", "a")
            text.write(output)
            text.close()
            continue

        print "== Creating and saving zonal statistics =="
        bedrockRaster = current_path + "/area_bedrock.tif"
        zonStatOutput = current_path + "/maxHeight.dbf"
        ZonalStatisticsAsTable(in_zone_data=polygonGrp,
                               zone_field="Id",
                               in_value_raster=bedrockRaster,
                               out_table=zonStatOutput,
                               statistics_type="ALL",
                               ignore_nodata=True)

        print "== Joining shapefile and zonal statistics =="
        arcpy.JoinField_management(in_data=polygonGrp,
                                   in_field="Id",
                                   join_table=zonStatOutput,
                                   join_field="Id",
                                   fields=["MIN", "MAX", "RANGE", "MEAN", "AREA", "SUM"])

        print "== Rounding MAX height =="
        arcpy.AddField_management(polygonGrp, field_name="MAX_round")
        with arcpy.da.UpdateCursor(polygonGrp, ["MAX", "MAX_round"]) as cursor:
            for row in cursor:
                row[1] = round(row[0], 0)
                cursor.updateRow(row)

        print "== Applying thresholds (>5m depth, >10.000 qm)"
        with arcpy.da.UpdateCursor(polygonGrp, ["RANGE"]) as cursor:
            for row in cursor:
                if row[0] < 5:
                    cursor.deleteRow()
        with arcpy.da.UpdateCursor(polygonGrp, ["AREA"]) as cursor:
            for row in cursor:
                if row[0] < 10000:
                    cursor.deleteRow()

        print "== Checking for remaining sinks "
        rows = [row for row in arcpy.da.SearchCursor(polygonGrp, field_names="*")]
        if len(rows) == 0:
            print "==============================="
            print "== No remaining sinks found! =="
            print "== Skipping glacier " + glacierID + " and starting next iteration"
            print "== deleting glacier folder =="
            print "==============================="
            del polygonGrp
            del bedrockRaster
            del zonStatOutput
            try:
                shutil.rmtree(current_path)
            except WindowsError:
                print "deleted everything save the LOCK-file"
                print "manual check required at " + current_path
            output = "Glacier folder " + glacierID + " was deleted."
            text = open(work_path + "/glacier_" + glacierID + "_was_deleted.txt", "a")
            text.write(output)
            text.close()
            continue
        elif len(rows) >= 1:
            print "== Found " + str(rows) + " remaining sinks =="
            print "== Clipping bedrock raster to sinks =="
            bedrockSinks = current_path + "/bedrock_sinks.tif"
            arcpy.Clip_management(bedrockRaster, out_raster=bedrockSinks,
                                  in_template_dataset=polygonGrp,
                                  clipping_geometry="ClippingGeometry",
                                  maintain_clipping_extent="NO_MAINTAIN_EXTENT")

            print "== Creating contours =="
            contours = current_path + "/contours.shp"
            Contour(in_raster=bedrockSinks,
                    out_polyline_features=contours,
                    contour_interval=5)

            print "== Creating TINs =="
            TINs = current_path + "/TINs"
            arcpy.CreateTin_3d(out_tin=TINs,
                               spatial_reference=arcpy.Describe(contours).spatialReference,
                               in_features=contours)

            print "== Clipping TINs to sink extents =="
            arcpy.EditTin_3d(in_tin=TINs,
                             in_features=polygonGrp,
                             constrained_delaunay="DELAUNAY")

            print "== Calculating separate TIN volume and area =="
            arcpy.PolygonVolume_3d(in_surface=TINs,
                                   in_feature_class=polygonGrp,
                                   in_height_field="MAX",
                                   out_volume_field="Volume",
                                   surface_area_field="SArea",
                                   pyramid_level_resolution="0")
            print "== Completed =="

        ################################################################

        # create folder for current TIN and copy TINS
        try:
            os.mkdir(subfolder_TIN)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_TIN)

        print "editing single TINs with each row of the shp"
        print "and copying the respective shp-files into subfolder_SHP"
        rows = [row for row in arcpy.da.SearchCursor(current_path + "/polygonGrp.shp", field_names="*")]
        for j in range(len(rows)):
            subsubfolder_TIN = subfolder_TIN + "/" + "TIN_" + str(j)
            try:
                os.mkdir(subsubfolder_TIN)
            except OSError:
                print("Creation of the directory %s failed" % subsubfolder_TIN)
            # copy the whole TIN-set in the subfolder:
            arcpy.CopyTin_3d(in_tin=current_path + "/TINs", out_tin=subsubfolder_TIN)
            # create copy of "polygonGrp" to edit the TIN-set into the single TINs:
            arcpy.CopyFeatures_management(current_path + "/polygonGrp.shp", current_path + "/sinks.shp")
            # delete everything save the current row
            with arcpy.da.UpdateCursor(current_path + "/sinks.shp", "FID") as cursor:
                for row in cursor:
                    if row[0] > j or row[0] < j:
                        cursor.deleteRow()
            # copy the remaining shp into the subfolder for the sink-shapes
            arcpy.CopyFeatures_management(current_path + "/sinks.shp", subfolder_SHP + "/sinks" + str(i) + ".shp")
            # clip the TIN-set with the one-shape-clipped result of polygonGrp
            arcpy.EditTin_3d(in_tin=subfolder_TIN + "/" + "TIN_" + str(i),
                             in_features=current_path + "/sinks.shp",
                             constrained_delaunay="DELAUNAY")

        # extract only relevant sinks
        try:
            os.mkdir(subfolder_SINKS)
        except OSError:
            print("Creation of the directory %s failed" % subfolder_SINKS)

        print "copy sink shapefile into new sink folder and create sink-raster"
        for sinkshape in os.listdir(subfolder_SHP):
            if sinkshape.endswith(".shp"):
                sinkNr = re.findall(r'\d+', sinkshape)[0]
                subsubfolder_SINKS = subfolder_SINKS + "/sink" + str(sinkNr)
                try:
                    os.mkdir(subsubfolder_SINKS)
                except OSError:
                    print("Creation of the directory %s failed" % subsubfolder_SINKS)
                arcpy.CopyFeatures_management(subfolder_SHP + "/" + sinkshape, subsubfolder_SINKS + "/sink.shp")
                with arcpy.da.SearchCursor(subsubfolder_SINKS + "/sink.shp", "gridcode") as cursor:
                    for row in cursor:
                        gridcode = row[0]
                extract = ExtractByAttributes(current_path + "/area_sinksRegionInt.tif",
                                              ' "Value" = {}'.format(gridcode))
                extract.save(subsubfolder_SINKS + "/sink.tif")
                del extract

        ############################################################################################################

        # Step Three: Find relevant slopes for each overdeepening

        print "calculate relevant slopes for every sink"
        for j in range(len(os.listdir(subfolder_SINKS))):
            sinkNr = j
            sinkfolder = subfolder_SINKS + "/sink" + str(sinkNr)
            for lakeraster in os.listdir(sinkfolder):
                if lakeraster.startswith("sink") and lakeraster.endswith("tif"):
                    lakeraster = Raster(sinkfolder + "/" + lakeraster)
                    # expand the lake raster:
                    lakeValue = lakeraster.maximum
                    outExpand = Expand(in_raster=lakeraster,
                                       number_cells=1,
                                       zone_values=lakeValue)
                    outExpand.save(sinkfolder + "/expandLake.tif")
                    del outExpand
                    del lakeraster
                    del lakeValue

                    # find inlets and outlets
                    lake = Raster(sinkfolder + "/sink.tif")
                    expand = Raster(sinkfolder + "/expandLake.tif")
                    streams = Raster(current_path + "/area_streams.tif")
                    flowACC = Raster(current_path + "/area_flowAcc.tif")
                    # calculate raster
                    IOS = Con(IsNull(lake), streams * flowACC * expand / expand)
                    outletValue = IOS.maximum
                    IOS_outlet = Con(IOS, 1, None, "VALUE = {}".format(outletValue))
                    arcpy.RasterToPoint_conversion(in_raster=IOS_outlet,
                                                   out_point_features=sinkfolder + "/outlet.shp")
                    # exclude small tributary flows
                    if outletValue >= 200:
                        IOS_inlets = Con(IOS, 1, None, "VALUE > 200 AND VALUE < {}".format(outletValue))
                        arcpy.RasterToPoint_conversion(in_raster=IOS_inlets,
                                                       out_point_features=sinkfolder + "/inlets.shp")
                    elif outletValue < 200:
                        IOS_inlets = Con(IOS, 1, None, "VALUE < {}".format(outletValue))
                        arcpy.RasterToPoint_conversion(in_raster=IOS_inlets,
                                                       out_point_features=sinkfolder + "/inlets.shp")
                        with arcpy.da.UpdateCursor(sinkfolder + "/inlets.shp", "FID") as cursor:
                            for row in cursor:
                                cursor.deleteRow()

                    del IOS_inlets
                    del IOS_outlet

                    # check if inlets exist and calculate watershed for whole area around lake
                    if int(arcpy.GetCount_management(sinkfolder + "/inlets.shp").getOutput(0)) == 0:
                        print "Sink " + str(sinkNr) + " has no inlets."
                        # arcpy.RasterToPolygon_conversion(in_raster=lake,
                        #                                  out_polygon_features=sinkfolder + "/sink.shp",
                        #                                  simplify="NO_SIMPLIFY",
                        #                                  create_multipart_features="MULTIPLE_OUTER_PART")
                        arcpy.Buffer_analysis(in_features=sinkfolder + "/sink.shp",
                                              out_feature_class=sinkfolder + "/sink_buffer.shp",
                                              buffer_distance_or_field="2 Kilometers")
                        arcpy.Clip_management(in_raster=current_path + "/area_flowDir.tif",
                                              out_raster=sinkfolder + "/flowDir_clip.tif",
                                              in_template_dataset=sinkfolder + "/sink_buffer.shp",
                                              clipping_geometry="ClippingGeometry")
                        ExtractValuesToPoints(in_point_features=sinkfolder + "/outlet.shp",
                                              in_raster=current_path + "/area_flowAcc.tif",
                                              out_point_features=sinkfolder + "/outletPointVAL.shp")
                        watershedFinal = Watershed(in_flow_direction_raster=sinkfolder + "/flowDir_clip.tif",
                                                   in_pour_point_data=sinkfolder + "/outletPointVAL.shp")
                        watershedFinal.save(sinkfolder + "/watershed_final.tif")
                        del watershedFinal
                        arcpy.RasterToPolygon_conversion(in_raster=sinkfolder + "/watershed_final.tif",
                                                         out_polygon_features=sinkfolder + "/watershed_final.shp",
                                                         simplify="NO_SIMPLIFY",
                                                         create_multipart_features="MULTIPLE_OUTER_PART")

                    # if inlets exist, calculate and disregard their watershed areas
                    elif int(arcpy.GetCount_management(sinkfolder + "/inlets.shp").getOutput(0)) > 0:
                        print "Sink " + str(sinkNr) + " has " + str(
                            int(arcpy.GetCount_management(sinkfolder + "/inlets.shp").getOutput(0))) + " inlets."

                        with arcpy.da.SearchCursor(sinkfolder + "/inlets.shp", "FID") as cursor:
                            for row in cursor:
                                FID = row[0]
                                arcpy.Select_analysis(in_features=sinkfolder + "/inlets.shp",
                                                      out_feature_class=sinkfolder + "/inletPoint{}.shp".format(
                                                          str(FID)),
                                                      where_clause=' "FID" = {}'.format(str(FID)))
                                ExtractValuesToPoints(in_point_features=sinkfolder + "/inletPoint{}.shp".
                                                      format(str(FID)),
                                                      in_raster=current_path + "/area_flowAcc.tif",
                                                      out_point_features=sinkfolder + "/inletPoint{}_VAL.shp".
                                                      format(str(FID)))

                                # calculate watershed
                                outWatershed = Watershed(in_flow_direction_raster=current_path + "/area_flowDir.tif",
                                                         in_pour_point_data=sinkfolder + "/inletPoint{}_VAL.shp".format(
                                                             str(FID)))
                                outWatershed.save(sinkfolder + "/watershed{}.tif".format(str(FID)))
                                del outWatershed

                                # create polygon file
                                arcpy.RasterToPolygon_conversion(
                                    in_raster=sinkfolder + "/watershed{}.tif".format(str(FID)),
                                    out_polygon_features=sinkfolder + "/watershed{}.shp".format(str(FID)),
                                    simplify="NO_SIMPLIFY",
                                    create_multipart_features="MULTIPLE_OUTER_PART")

                                # buffer sinks to use for area clipping
                                # arcpy.RasterToPolygon_conversion(in_raster=sinkfolder + "/sink.tif",
                                #                                  out_polygon_features=sinkfolder + "/sink.shp",
                                #                                  simplify="NO_SIMPLIFY",
                                #                                  create_multipart_features="MULTIPLE_OUTER_PART")
                                arcpy.Buffer_analysis(in_features=sinkfolder + "/sink.shp",
                                                      out_feature_class=sinkfolder + "/sink_buffer.shp",
                                                      buffer_distance_or_field="2 Kilometers")

                        # find relevant slope areas
                        watershed_list = []
                        for watershed in os.listdir(sinkfolder):
                            if watershed.startswith("watershed") and watershed.endswith("shp"):
                                watershed_list.append(sinkfolder + "/" + watershed)
                                # print watershed_list
                        if len(watershed_list) == 1:
                            print "Only 1 inlet and watershed"
                            arcpy.Clip_analysis(in_features=sinkfolder + "/watershed0.shp",
                                                clip_features=sinkfolder + "/sink_buffer.shp",
                                                out_feature_class=sinkfolder + "/watershed_clip.shp")

                        elif len(watershed_list) > 1:
                            print "{} inlets and subsequent watersheds".format(str(len(watershed_list)))
                            arcpy.Merge_management(watershed_list, sinkfolder + "/watersh_merge.shp")
                            arcpy.Dissolve_management(in_features=sinkfolder + "/watersh_merge.shp",
                                                      out_feature_class=sinkfolder + "/watersh_diss.shp")
                            arcpy.Clip_analysis(in_features=sinkfolder + "/watersh_diss.shp",
                                                clip_features=sinkfolder + "/sink_buffer.shp",
                                                out_feature_class=sinkfolder + "/watershed_clip.shp")

                        arcpy.Erase_analysis(in_features=sinkfolder + "/sink_buffer.shp",
                                             erase_features=sinkfolder + "/watershed_clip.shp",
                                             out_feature_class=sinkfolder + "/remainingArea.shp")

                        # calculate watershed for outlet point within relevant regions
                        arcpy.Clip_management(in_raster=current_path + "/area_flowDir.tif",
                                              out_raster=sinkfolder + "/flowDir_clip.tif",
                                              in_template_dataset=sinkfolder + "/remainingArea.shp",
                                              clipping_geometry="ClippingGeometry")
                        ExtractValuesToPoints(in_point_features=sinkfolder + "/outlet.shp",
                                              in_raster=current_path + "/area_flowAcc.tif",
                                              out_point_features=sinkfolder + "/outletPointVAL.shp")

                        watershedFinal = Watershed(in_flow_direction_raster=sinkfolder + "/flowDir_clip.tif",
                                                   in_pour_point_data=sinkfolder + "/outletPointVAL.shp")
                        watershedFinal.save(sinkfolder + "/watershed_final.tif")
                        arcpy.RasterToPolygon_conversion(in_raster=sinkfolder + "/watershed_final.tif",
                                                         out_polygon_features=sinkfolder + "/watershed_final.shp",
                                                         simplify="NO_SIMPLIFY",
                                                         create_multipart_features="MULTIPLE_OUTER_PART")

        #############################################################################################

        print "calculate slope raster for watershed area"
        for j in range(len(os.listdir(subfolder_SINKS))):
            sinkNr = j
            sink_path = subfolder_SINKS + "/sink" + str(sinkNr)
            arcpy.Clip_management(in_raster=current_path + "/area_slope.tif",
                                  out_raster=sink_path + "/slope.tif",
                                  in_template_dataset=sink_path + "/watershed_final.shp",
                                  clipping_geometry="ClippingGeometry")

        # define function for slope detection
        def single_slope_calc(sinknr, threshold_low, threshold_high):
            sink_folder = subfolder_SINKS + "/sink" + str(sinknr)
            # check for areas within threshold
            slopes = Raster(sink_folder + "/slope.tif")
            if slopes.maximum < threshold_low:
                print "###############################################"
                print "Sink " + str(sinknr) + " has no slopes between " + str(threshold_low) + " and " + str(
                    threshold_high) + " degrees!"
                print "###############################################"
                return
            elif slopes.maximum >= threshold_low:
                print "Sink " + str(sinknr) + " has some contributing slope areas between " + str(
                    threshold_low) + " and " + str(threshold_high)

            # create new subsubfolder for the slope class
            subfolder_slopes = sink_folder + "/slopes" + str(threshold_low)
            try:
                os.mkdir(subfolder_slopes)
            except OSError:
                print("Creation of the directory %s failed" % subfolder_slopes)

            # generate 0-1-raster for all critical slopes
            slopes = Con(sink_folder + "/slope.tif", 1, 0,
                         "VALUE >= {} AND VALUE < {}".format(threshold_low, threshold_high))
            slopes.save(subfolder_slopes + "/slopes_prelim.tif")

            # exclude the lake area from the binary slope raster
            slopes_binary = Con(IsNull(sink_folder + "/sink.tif"), subfolder_slopes + "/slopes_prelim.tif")
            slopes_binary.save(subfolder_slopes + "/slopes_binary.tif")
            del slopes
            del slopes_binary

            # generate polygon shapefile for slopes
            arcpy.RasterToPolygon_conversion(in_raster=subfolder_slopes + "/slopes_binary.tif",
                                             out_polygon_features=subfolder_slopes + "/slopes.shp",
                                             simplify="NO_SIMPLIFY",
                                             create_multipart_features="MULTIPLE_OUTER_PART")
            with arcpy.da.UpdateCursor(subfolder_slopes + "/slopes.shp", ["gridcode"]) as crsr:
                for rw in crsr:
                    if rw[0] == 0:
                        crsr.deleteRow()

            # new check to see if deleting the lake area also deleted all the relevant slopes
            rws = [rw for rw in arcpy.da.SearchCursor(subfolder_slopes + "/slopes.shp", field_names="*")]
            if len(rws) == 0:
                print "###############################################"
                print "Sink " + str(sinknr) + " has no remaining slopes between " + str(threshold_low) + " and " + str(
                    threshold_high) + " degrees!"
                print "###############################################"
                return

            # get slope and (rounded) elevation values for critical slopes
            arcpy.Clip_management(in_raster=current_path + "/area_slope.tif",
                                  out_raster=subfolder_slopes + "/slope_clip.tif",
                                  in_template_dataset=subfolder_slopes + "/slopes.shp",
                                  clipping_geometry="ClippingGeometry")
            arcpy.Clip_management(in_raster=current_path + "/area_bedrock.tif",
                                  out_raster=subfolder_slopes + "/elevation_clipFL.tif",
                                  in_template_dataset=subfolder_slopes + "/slopes.shp",
                                  clipping_geometry="ClippingGeometry")
            out_int = Int(subfolder_slopes + "/elevation_clipFL.tif")
            out_int.save(subfolder_slopes + "/elevation_clipINT.tif")
            del out_int

            # find single slope areas
            # create single-pixel polygon from integer elevation raster
            arcpy.RasterToPolygon_conversion(in_raster=subfolder_slopes + "/elevation_clipINT.tif",
                                             out_polygon_features=subfolder_slopes + "/slope_poly.shp",
                                             simplify="NO_SIMPLIFY",
                                             create_multipart_features="SINGLE_OUTER_PART")
            # dissolve into single slope areas
            arcpy.Dissolve_management(in_features=subfolder_slopes + "/slope_poly.shp",
                                      out_feature_class=subfolder_slopes + "/slope_dissolved.shp",
                                      multi_part="SINGLE_PART")
            # delete small slope areas
            arcpy.AddGeometryAttributes_management(Input_Features=subfolder_slopes + "/slope_dissolved.shp",
                                                   Geometry_Properties="AREA",
                                                   Area_Unit="SQUARE_METERS")
            row_list = []
            with arcpy.da.UpdateCursor(subfolder_slopes + "/slope_dissolved.shp", "POLY_AREA") as crsr:
                for rw in crsr:
                    if rw[0] < 5000:
                        crsr.deleteRow()
                    elif rw[0] >= 5000:
                        row_list.append(rw)

            if len(row_list) == 0:
                print "###############################################"
                print "Sink " + str(sinknr) + " has only very small slopes between " + str(
                    threshold_low) + " and " + str(threshold_high) + " degrees!"
                print "###############################################"
                # return
            elif len(row_list) > 0:
                print "Sink " + str(sinknr) + " has " + str(len(row_list)) + " contributing slope areas."

            # split slopes shapefile into a separate file for every slope
            single_slopes = subfolder_slopes + "/single_slopes"
            try:
                os.mkdir(single_slopes)
            except OSError:
                print("Creation of the directory %s failed" % single_slopes)

            with arcpy.da.SearchCursor(subfolder_slopes + "/slope_dissolved.shp", ["FID"]) as crsr:
                for rw in crsr:
                    try:
                        slope_nr = single_slopes + "/" + "slope" + str(rw[0])
                        os.mkdir(slope_nr)
                    except OSError:
                        print("Creation of the directory %s failed" % slope_nr)
                    arcpy.Select_analysis(in_features=subfolder_slopes + "/slope_dissolved.shp",
                                          out_feature_class=slope_nr + "/slope" + str(rw[0]),
                                          where_clause="""{0} = {1}""".format("FID", int(rw[0])))

            # iterate through all the slopes for each sink and calculate everything:
            # the rounded distance to sink, mean elevation and mean slope
            # and merge all the shapefiles
            for k in range(len(os.listdir(single_slopes))):
                slopenr = k
                slopefolder = single_slopes + "/slope" + str(slopenr)
                for subdir, dirs, files in os.walk(slopefolder):
                    for file_name in files:
                        if file_name.endswith(".shp"):

                            # calculate and round distance
                            arcpy.Near_analysis(in_features=os.path.join(subdir, file_name),
                                                near_features=sink_folder + "/sink.shp")
                            with arcpy.da.UpdateCursor(os.path.join(subdir, file_name), ["NEAR_DIST"]) as crsr:
                                for rw in crsr:
                                    rw[0] = round(rw[0], 0)
                                    cursor.updateRow(rw)

                            # clip bedrock raster to slope extent
                            arcpy.Clip_management(in_raster=current_path + "/area_bedrock.tif",
                                                  out_raster=subdir + "/DEM_slope.tif",
                                                  in_template_dataset=os.path.join(subdir, file_name),
                                                  clipping_geometry="ClippingGeometry")

                            # calculate MEAN elevation for slope
                            zonstat_output = subdir + "/slopeStats.dbf"
                            ZonalStatisticsAsTable(in_zone_data=os.path.join(subdir, file_name),
                                                   zone_field="NEAR_FID",
                                                   in_value_raster=subdir + "/DEM_slope.tif",
                                                   out_table=zonstat_output,
                                                   statistics_type="ALL",
                                                   ignore_nodata=True)

                            # join results with slope polygon
                            arcpy.JoinField_management(in_data=os.path.join(subdir, file_name),
                                                       in_field="NEAR_FID",
                                                       join_table=zonstat_output,
                                                       join_field="NEAR_FID",
                                                       fields=["MIN", "MAX", "MEAN"])

                            # clip slope raster to slope extent
                            arcpy.Clip_management(in_raster=current_path + "/area_slope.tif",
                                                  out_raster=subdir + "/slope_Abhang.tif",
                                                  in_template_dataset=os.path.join(subdir, file_name),
                                                  clipping_geometry="ClippingGeometry")

                            # calculate MEAN elevation for slope
                            zonstat_output = subdir + "/slopeStats_Abhang.dbf"
                            ZonalStatisticsAsTable(in_zone_data=os.path.join(subdir, file_name),
                                                   zone_field="NEAR_FID",
                                                   in_value_raster=subdir + "/slope_Abhang.tif",
                                                   out_table=zonstat_output,
                                                   statistics_type="ALL",
                                                   ignore_nodata=True)

                            # join results with slope polygon
                            arcpy.JoinField_management(in_data=os.path.join(subdir, file_name),
                                                       in_field="NEAR_FID",
                                                       join_table=zonstat_output,
                                                       join_field="NEAR_FID",
                                                       fields=["MIN", "MAX", "MEAN"])

                            # calculate area estimate regarding different possibilities of sloping terrain
                            arcpy.AddField_management(os.path.join(subdir, file_name), field_name="DEM_AREA")
                            expression = "((32*(math.sqrt(((math.tan(math.radians(!MEAN_1!))*32)**2)+32**2)) + " \
                                         "((math.tan(math.radians(!MEAN_1!))*32)**2)+32**2) / 2) * (!POLY_AREA!/1024)"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="DEM_AREA",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")

                            # add height and area fields to each slope
                            lake_table = sink_folder + "/sink.shp"
                            arcpy.JoinField_management(in_data=os.path.join(subdir, file_name),
                                                       in_field="FID",
                                                       join_table=lake_table,
                                                       join_field="FID",
                                                       fields=["MAX_round", "AREA"])

                            # Calculate the hazard scores:
                            # height difference
                            arcpy.AddField_management(os.path.join(subdir, file_name),
                                                      field_name="HZRD_HGHT", field_type="FLOAT")
                            expression = "(!MAX! - !MAX_round!) / 1000"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="HZRD_HGHT",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")
                            with arcpy.da.UpdateCursor(os.path.join(subdir, file_name), ["HZRD_HGHT"]) as crsr:
                                for rw in crsr:
                                    if rw[0] > 1:
                                        rw[0] = 1
                                        crsr.updateRow(rw)

                            # slope
                            arcpy.AddField_management(os.path.join(subdir, file_name),
                                                      field_name="HZRD_SLOPE", field_type="FLOAT")
                            expression = "(!MEAN_1! - 20) / (70 - 20)"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="HZRD_SLOPE",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")
                            with arcpy.da.UpdateCursor(os.path.join(subdir, file_name), ["HZRD_SLOPE"]) as crsr:
                                for rw in crsr:
                                    if rw[0] > 1:
                                        rw[0] = 1
                                        crsr.updateRow(row)

                            # distance
                            arcpy.AddField_management(os.path.join(subdir, file_name),
                                                      field_name="HZRD_DIST", field_type="FLOAT")
                            expression = "1 - !NEAR_DIST! / 2000"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="HZRD_DIST",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")

                            # area
                            arcpy.AddField_management(os.path.join(subdir, file_name),
                                                      field_name="HZRD_AREA", field_type="FLOAT")
                            expression = "!DEM_AREA! / !AREA!"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="HZRD_AREA",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")
                            with arcpy.da.UpdateCursor(os.path.join(subdir, file_name), ["HZRD_AREA"]) as crsr:
                                for rw in cursor:
                                    if rw[0] > 1:
                                        rw[0] = 1
                                        crsr.updateRow(rw)

                            # add field for final hazard classification of each slope
                            arcpy.AddField_management(os.path.join(subdir, file_name),
                                                      field_name="HAZARD_NEW", field_type="FLOAT")
                            expression = "((!HZRD_AREA!*2) + !HZRD_DIST! + (!HZRD_SLOPE!*2) + (!HZRD_HGHT!/2)) / 5.5"
                            arcpy.CalculateField_management(in_table=os.path.join(subdir, file_name),
                                                            field="HAZARD_NEW",
                                                            expression=expression,
                                                            expression_type="PYTHON_9.3")
                            with arcpy.da.UpdateCursor(os.path.join(subdir, file_name), ["HAZARD_NEW"]) as crsr:
                                for rw in cursor:
                                    if rw[0] > 1:
                                        rw[0] = 1
                                        crsr.updateRow(rw)

            # merge all slopes of the current class
            matches_list = []
            for subdir, dirs, files in os.walk(single_slopes):
                for file_name in files:
                    if file_name.endswith(".shp"):
                        mtch = (os.path.join(subdir, file_name))
                        matches_list.append(mtch)
            if matches_list:
                slopes_merged = sink_folder + "/merged_slopes"
                try:
                    os.mkdir(slopes_merged)
                except OSError:
                    print("Creation of the directory %s failed" % slopes_merged)
                arcpy.Merge_management(matches_list,
                                       slopes_merged + "/sink{}_merged_{}_to_{}.shp".format(sinknr, threshold_low,
                                                                                            threshold_high))


        # Step Four: Calculate hazards for all slopes
        # detect all slopes
        for j in range(len(os.listdir(subfolder_SINKS))):
            single_slope_calc(sinknr=j, threshold_low=60, threshold_high=90)
            single_slope_calc(sinknr=j, threshold_low=40, threshold_high=60)
            single_slope_calc(sinknr=j, threshold_low=30, threshold_high=40)

        print "Merge slopes of all classes for each sink"
        for j in range(len(os.listdir(subfolder_SINKS))):
            sinkfolder = subfolder_SINKS + "/sink" + str(j)

            # see if the sink has critical slopes in the respective threshold category
            if os.path.exists(sinkfolder + "/merged_slopes"):
                matches = []
                for filename in os.listdir(sinkfolder + "/merged_slopes"):
                    if filename.startswith("sink") and filename.endswith("shp"):
                        match = (os.path.join(sinkfolder + "/merged_slopes/" + filename))
                        matches.append(match)
            else:
                print "Sink Nr. " + str(j) + " has no critical slopes."
                output = "Sink Nr. " + str(j) + " has no critical slopes.\n"
                text = open(current_path + "/error_messages.txt", "w")
                text.write(output)
                text.close()
                continue

            # merge the different slope categories
            arcpy.Merge_management(matches, sinkfolder + "/merged_slopes/all_merged.shp")

        ############################################################################################

        # calculate hazard approximation mean
        for j in range(len(os.listdir(subfolder_SINKS))):
            print "working on sink " + str(j)
            sinkfolder = subfolder_SINKS + "/sink" + str(j)
            arcpy.AddField_management(sinkfolder + "/sink.shp", field_name="HAZARD_NEW", field_type="FLOAT")
            arcpy.AddField_management(sinkfolder + "/sink.shp", field_name="HAZARD_MAX", field_type="FLOAT")
            arcpy.AddField_management(sinkfolder + "/sink.shp", field_name="IMPCT_AREA", field_type="FLOAT")
            hazard_sum = 0
            count = 0
            impact_area = 0

            # calculate normal hazard classification
            try:
                with arcpy.da.SearchCursor(sinkfolder + "/merged_slopes/all_merged.shp", "HAZARD_NEW") as cursor:
                    for row in cursor:
                        hazard_sum = hazard_sum + row[0]
                        count = count + 1
                with arcpy.da.UpdateCursor(sinkfolder + "/sink.shp", "HAZARD_NEW") as cursor:
                    for row in cursor:
                        row[0] = hazard_sum / count
                        cursor.updateRow(row)
                    for row in cursor:
                        if row[0] > 1:
                            row[0] = 1
                            cursor.updateRow(row)

                # calculate maximum hazard
                hazards = [j[0] for j in
                           arcpy.da.SearchCursor(sinkfolder + "/merged_slopes/all_merged.shp", "HAZARD_NEW")]
                max_hazard = max(hazards)
                with arcpy.da.UpdateCursor(sinkfolder + "/sink.shp", "HAZARD_MAX") as cursor:
                    for row in cursor:
                        row[0] = max_hazard
                        cursor.updateRow(row)

                # calculate lake impact predisposition area
                with arcpy.da.SearchCursor(sinkfolder + "/merged_slopes/all_merged.shp", "DEM_AREA") as cursor:
                    for row in cursor:
                        impact_area = impact_area + row[0]
                with arcpy.da.UpdateCursor(sinkfolder + "/sink.shp", "IMPCT_AREA") as cursor:
                    for row in cursor:
                        row[0] = impact_area
                        cursor.updateRow(row)
            except RuntimeError:
                print "Setting hazard level for sink Nr. " + str(j) + " to 0."
                output = "Setting hazard level for sink Nr. " + str(j) + " to 0."
                text = open(current_path + "/error_messages.txt", "a")
                text.write(output)
                text.close()
                with arcpy.da.UpdateCursor(sinkfolder + "/sink.shp", "HAZARD_NEW") as cursor:
                    for row in cursor:
                        row[0] = 0
                        cursor.updateRow(row)
