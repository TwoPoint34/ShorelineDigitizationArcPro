# Import required libraries and modules
import arcpy
from arcpy import env
from arcpy.sa import *
import os
import datetime

# This version of the tool contains little error handling (e.g. try/except statements).  The tool usage instructions
# contain hints on common errors the user might encounter and suggests possible work-arounds for such errors.

arcpy.env.overwriteOutput = True   # Output from running the tool multiple times replaces the output from previous runs

aprx = arcpy.mp.ArcGISProject("CURRENT")   # Saves the current ArcGIS project file path to a variable so that it can later
# be used to automatically add the script tool's output layers to the current project's table of contents

# Adding a timestamp to the output file names is useful for output versioning and code debugging purposes
now = datetime.datetime.now()

timestamp = str(now.strftime("%m%d%H%M"))

image_raster_path = arcpy.env.workspace = arcpy.GetParameterAsText(0)  # Set workspace, which is the folder containing the input
# rasters/images. Code saves the initial user input workspace setting to a variable so that its file path can be referenced for raster
# access even after the workspace has been changed to a new (output) location below.  The workspace is changed so that the outputs of
# the tool do not go into the same folder as the raster inputs.

raster_list = arcpy.ListRasters("*")   # List all the rasters in the workspace folder
imagerasterlist = []   # Create a new list so that the rasters in the original workspace folder can be accessed even after the
# workspace is changed.  This is possibly redundant, as I have not yet checked as to whether raster_list above has already done so.
for raster in raster_list:
    imagerasterlist.append(raster)   # Loop/iterate through each raster in the list/folder and add each raster to imagerasterlist

DEM_raster_path = arcpy.env.workspace = arcpy.GetParameterAsText(1)  # Reset workspace to access DEM rasters (if input by user)

if len(DEM_raster_path) != 0:
    DEM_raster_list = arcpy.ListRasters("*")   # List all the rasters in the workspace folder
    DEMrasterlist = []
    for demraster in DEM_raster_list:
        DEMrasterlist.append(demraster)

BlueBandNumber = arcpy.GetParameterAsText(2)   # Ask tool user to define sequential number of blue band within multi-band image
if len(BlueBandNumber) != 0:
    BlueBandNumber = int(BlueBandNumber)   # If blue band parameter is not blank, then convert the input to an integer in order to
    # make a blue band raster layer later in the code

GreenBandNumber = arcpy.GetParameterAsText(3)
if len(GreenBandNumber) != 0:
    GreenBandNumber = int(GreenBandNumber)

RedBandNumber = arcpy.GetParameterAsText(4)
if len(RedBandNumber) != 0:
    RedBandNumber = int(RedBandNumber)

NIRBandNumber = arcpy.GetParameterAsText(5)
if len(NIRBandNumber) != 0:
    NIRBandNumber = int(NIRBandNumber)

SWIRBandNumber = arcpy.GetParameterAsText(6)
if len(SWIRBandNumber) != 0:
    SWIRBandNumber = int(SWIRBandNumber)

AdditionToOutputFileName = arcpy.GetParameterAsText(7)   # Ask script tool user to name script tool output.  This parameter
# is optional but can help the user keep track of multiple tool runs

arcpy.env.workspace = arcpy.GetParameterAsText(8)   # Change workspace so that output files are saved to a different
# location than the input raster folder.  Saving as a file geodatabase (versus shapefile) allows for easy water surface area
# measurement through the shape_area field.  Also, for some reason, the tool does not work properly if the workspace is set as
# a folder rather than a geodatabase.  It's worth noting too that the direction of the script tool parameter works only as input
# and not output.

ImageDEMWLTable = arcpy.GetParameterAsText(9)
# Inputting a DEM allows for the creation of contour lines at ground reference water levels.  The user can compare these contours
# to the shoreline features created from multi-band imagery by the script tool.  The user may identify an existing DEM that
# approximately matches the time period and tide conditions of the aerial images being analyzed.  Alternatively, users may make
# their own DEMs using a variety of approaches.  For example, custom DEMs might be made through UAV-mounted LiDAR sensors or
# UAV-based stereo imagery.  A basic DEM could alternatively be created by interpolating a surface between measured beach elevation
# points.  Possible interpolation techniques for points include spline, kriging, and TIN creation. If the user has contour lines
# rather than elevation points, the ArcGIS Topo to Raster tool can interpolate a crude DEM.
# With this parameter, the user can input an Excel file that contains the relationship between the aerial image, the DEM, and
# water levels at the time of aerial image capture.  In this way, the user can input water levels under various tide (e.g. MHHW)
# or sea/water-level rise conditions.  The shoreline script tool will use the water level input value to create a contour line
# (shoreline) on the input DEM.  Assuming the water level at time of aerial image collection is above the DEM water level, the
# script tool derived aerial image shoreline can be compared to the DEM-based contour and thereby aid in the identification of
# shoreline change.  NOAA water level data can be found at:  https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels

ThresholdInputString = arcpy.GetParameterAsText(10)   # Allows user to specify the NDWI land/water threshold values, given that
# different values delineate shoreline more effectively for different aerial images.  Input format must be as follows (though
# values themselves can be changed):  0 0.1 0.2 0.3

MinWaterSurfaceArea = arcpy.GetParameterAsText(11)   # The act of converting to an integer automatically rounds positive numbers down to the
# nearest whole number such that the user input value will be included when values greater than the integer of that value are selected.
if len(MinWaterSurfaceArea) != 0:
    MinWaterSurfaceArea = int(MinWaterSurfaceArea)

# Write a function to avoid code repetition when adding outputs to ArcPro table of contents
def AddOutput2TOC():
    desc = arcpy.Describe(AddToTOC)
    mappyfilepath = desc.catalogPath   # Get the file path of the features to be added to table of contents
    aprxMap = aprx.listMaps("Map")   # Access the current ArcPro map
    for mappy in aprxMap:
        mappy.addDataFromPath(mappyfilepath)   # Add the output feature for each table row to the map table of contents

if len(ImageDEMWLTable) != 0:   # Assuming the user inputs an Excel file, execute the code below
    descextension = arcpy.Describe(ImageDEMWLTable)  # Create a describe object of the Excel file to extract properties below
    WLVTextension = descextension.extension  # Extract the file extension of the input Excel file
    if WLVTextension == '.xls' or '.xlsx':  # Make sure that input spreadsheet has proper file extension
        Input_Excel_File = ImageDEMWLTable
        Output_Table = 'FGDBTable' + '_' + AdditionToOutputFileName + '_' + timestamp
        SpreadsheetToFileGDBTable = arcpy.conversion.ExcelToTable(Input_Excel_File, Output_Table)  # Convert Excel file to FGDB table
        FGDBsearchcursor = arcpy.SearchCursor(SpreadsheetToFileGDBTable)  # Establish search cursor for analysis of table

        rowcount = 0
        tableimagerasterlist = []
        for row in FGDBsearchcursor:  # Get values below from FGDB table
            rowcount = rowcount + 1
            imagebasenameobject = row.getValue('ImageBaseName')
            DEMbasenameobject = row.getValue('DEMBaseName')
            WaterLevel = row.getValue('WaterLevel')

            for DEMraster in DEMrasterlist:  # Access folder containing DEM rasters
                DEMraster = arcpy.sa.Raster(os.path.join(DEM_raster_path, DEMraster))  # Create raster object from raster path
                descfolderDEM = arcpy.Describe(DEMraster)  # Describe raster object to allow extraction of baseName and catalogPath
                folderDEMbasename = descfolderDEM.baseName
                if folderDEMbasename == DEMbasenameobject:  # Match folder DEM file name with table DEM file name so that the
                # DEM file specified in the Excel/FGDB table is further processed by the subsequent code
                    DEMrasterpath = descfolderDEM.catalogPath  # Save the folder DEM file path to a variable
                    break  # Only a single matching DEM value is desired, so break DEMrasterlist

            for imageraster in imagerasterlist:  # Repeat process described for DEM above but this time for aerial image
                imageraster = arcpy.sa.Raster(os.path.join(image_raster_path, imageraster))
                descimageraster = arcpy.Describe(imageraster)
                imagerasterbasename = descimageraster.baseName
                if imagerasterbasename == imagebasenameobject:
                    imagerasterpath = descimageraster.catalogPath
                    break

            tableimagerasterlist.append(imagerasterpath)

            # Use the water level value from the input Excel/FGDB table to draw a contour line on the corresponding DEM
            out_polyline_features = imagerasterbasename + '_' + "GroundReferenceWaterLevel" + '_' + AdditionToOutputFileName + '_' + timestamp
            contour_value = [WaterLevel]  # Input type is a list (even if only contains a single item)
            GroundReferenceWaterLevel = AddToTOC = ContourList(DEMrasterpath, out_polyline_features, contour_value)  # Draw a contour line at
            # the user input water level

            if rowcount < 6:  # Call function to add the first five water level contours to the table of contents
                AddOutput2TOC()

        imagerasterlist = tableimagerasterlist  # If the user does not input an Excel file, the script tool simply loops through the folder with
        # the image rasters.  With an Excel file input, however, the script tool loops through the table and accesses the image rasters based on
        # their order in the table

        del FGDBsearchcursor  # Delete the search cursor just as a matter of best practice
        arcpy.management.Delete(SpreadsheetToFileGDBTable)  # Delete intermediate FGDB table

imagecount = 0
for imageraster in imagerasterlist:  # Loop through list of image rasters to digitize shorelines
    binaryrasterpathlist = []  # Create list of binary raster paths for subsequent looping and shoreline analysis
    binrasterbnlist = []  # Create list of binary raster base names, which will be used to match table inputs to files in folder
    mergepolyinputlist = []  # Create list to be populated with output polygon features and later merged
    mergeinputlist = []   # Create list to be populated with output line features and later merged
    imagecount = imagecount + 1    # Add "1" to the count value for each iteration of the loop (i.e. each raster)
    imageraster = arcpy.sa.Raster(os.path.join(image_raster_path, imageraster))
    descimageraster = arcpy.Describe(imageraster)
    rasterbasename = descimageraster.baseName

    # Unsupervised classification (machine learing) output
    Input_raster_bands = imageraster
    unsupclassification = arcpy.sa.IsoClusterUnsupervisedClassification(Input_raster_bands, 2, "", 1)
    out_unsup = 'UnsupClass'
    unsupclasscopy = arcpy.management.CopyRaster(unsupclassification, out_unsup)
    descunsupraster = arcpy.Describe(unsupclasscopy)
    unsuprasterpath = descunsupraster.catalogPath
    unsuprastername = descunsupraster.baseName
    binaryrasterpathlist.append(unsuprasterpath)  # Add catalog path of iso cluster binary raster to list
    binrasterbnlist.append(unsuprastername)

    # The script tool can digitize shoreline without the Image Analyst extension, but in such a case the only
    # algorithm used is unsupervised classification.  Other land/water interface detection algorithms require
    # image analyst.
    if arcpy.CheckExtension("ImageAnalyst") == "Available":
        from arcpy.ia import *
        # As an essential component of the (M)NDWI calculation, the green band is the only mandatory band
        Green_rasterlayer = 'GreenBand'
        GreenRasterLayer = arcpy.management.MakeRasterLayer(imageraster, Green_rasterlayer, "", "", GreenBandNumber)  # Use GreenBandNumber
        # parameter to make a green band layer from the multi-band image
        out_greenraster = 'greenraster'
        GreenBand = arcpy.management.CopyRaster(GreenRasterLayer, out_greenraster)   # CopyRaster creates a permanent raster layer
        # from the temporary layer created by the MakeRasterLayer tool

        NDWIbandlist = []  # List for NIR and/or SWIR layers.  Looping through this list in subsequent steps allows for both NDWI
        # and MNDWI outputs to be produced if both NIR and SWIR bands are input to the tool
        if len(str(NIRBandNumber)) != 0:  # If the user inputs a NIR band number
            NIR_rasterlayer = 'NIRBand'
            NIRRasterLayer = arcpy.management.MakeRasterLayer(imageraster, NIR_rasterlayer, "", "", NIRBandNumber)
            out_nirraster = 'NDWI'
            NIRBand = arcpy.management.CopyRaster(NIRRasterLayer, out_nirraster)
            NIRRasterPath = os.path.join(arcpy.env.workspace, 'NDWI')  # For later usage in the RasterCalculator tool, it works
            # better to save the raster path rather than trying to use the raster object (NIRBand) itself
            NDWIbandlist.append(NIRBand)

        if len(str(SWIRBandNumber)) != 0:
            SWIR_rasterlayer = 'SWIRBand'
            SWIRRasterLayer = arcpy.management.MakeRasterLayer(imageraster, SWIR_rasterlayer, "", "", SWIRBandNumber)
            out_swirraster = 'MNDWI'
            SWIRBand = arcpy.management.CopyRaster(SWIRRasterLayer, out_swirraster)
            SWIRRasterPath = os.path.join(arcpy.env.workspace, 'MNDWI')
            NDWIbandlist.append(SWIRBand)

        if len(str(NIRBandNumber)) !=0 or len(str(SWIRBandNumber)) !=0:
            Bandcount = 0
            for Band in NDWIbandlist:
                Bandcount = Bandcount + 1
                IRBand = Band
                desc = arcpy.Describe(Band)
                basename = desc.baseName  # Getting the basename serves to retain the label of whether the band is NIR or SWIR
                # Formula for Normalized Difference Water Index is:  NDWI = (GreenBand - NIRBand) / (GreenBand + NIRBand)
                num = arcpy.sa.Float(GreenBand) - arcpy.sa.Float(IRBand)   # Define NDWI numerator. Convert to Float to enable calculation.
                denom = arcpy.sa.Float(GreenBand) + arcpy.sa.Float(IRBand)   # Define NDWI denominator. Convert to Float to enable calculation.
                NDWI = arcpy.sa.Divide(num, denom)  # Calculate NDWI by dividing numerator by denominator.

                Otsu_raster = arcpy.ia.Threshold(NDWI)  # I got the idea to use Otsu thresholding from the CoastSat tool
                out_otsu = basename + '_' + 'Otsu'
                OtsuRasterCopy = arcpy.management.CopyRaster(Otsu_raster, out_otsu)
                descotsuraster = arcpy.Describe(OtsuRasterCopy)
                otsurasterpath = descotsuraster.catalogPath
                otsurastername = descotsuraster.baseName
                binaryrasterpathlist.append(otsurasterpath)
                binrasterbnlist.append(otsurastername)

                count1 = 0    # Define a counter variable for the nested loop below.  This counter will give output polygons unique names
                # based on the cutoff value used in their production

                if len(ThresholdInputString) != 0:
                    thresholdlist = ThresholdInputString.split()   # Turns cutoffstring variable into a list that allows for iteration (below) and
                    # creation of output files for each of the cutoff values
                    for i in thresholdlist:  # Loops through all the items in cutofflist so that each raster in the raster list has various land/water
                    # cutoff values calculated such that the user can visually identify the polygon that best represents the shoreline
                        count1 = count1 + 1
                        threshold = float(i)  # Define a variable named parameter that converts each value in cutoff list to a float so that it
                        # can be used in the RemapRange tool below
                        remap = RemapRange([[-1, threshold, 0], [threshold, 1, 1]])  # Define the cutoff ranges for the binary (two value) raster
                        # to be created by the Reclassify tool below
                        reclass_field = "Value"  # Define the raster pixel value as the field to base reclassification on
                        NDWIbinary = Reclassify(NDWI, reclass_field, remap)  # Create a binary raster from each NDWI raster produced in outer for loop
                        out_binraster = basename + 'threshold' + str(count1)
                        NDWIbinaryCopy = arcpy.management.CopyRaster(NDWIbinary, out_binraster)
                        descbinraster = arcpy.Describe(NDWIbinaryCopy)
                        binrasterpath = descbinraster.catalogPath
                        brbasename = descbinraster.baseName
                        binaryrasterpathlist.append(binrasterpath)
                        binrasterbnlist.append(brbasename)

                if len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
                    Blue_rasterlayer = 'BlueBand'   # Define a unique layer output name
                    BlueRasterLayer = arcpy.management.MakeRasterLayer(imageraster, Blue_rasterlayer, "", "", BlueBandNumber)  # Provide parameters
                    # for MakeRasterLayer tool in order to separate the green band from the four band composite image
                    out_blueraster = 'blueraster'   # Name CopyRaster output.  By default, ArcPro saves the output raster to the project's geodatabase
                    BlueBand = arcpy.management.CopyRaster(BlueRasterLayer, out_blueraster)   # CopyRaster creates a permanent raster layer
                    # from the temporary layer created by the MakeRasterLayer tool above
                    BlueRasterPath = os.path.join(arcpy.env.workspace, 'blueraster')

                    Red_rasterlayer = 'RedBand'   # Define a unique layer output name
                    RedRasterLayer = arcpy.management.MakeRasterLayer(imageraster, Red_rasterlayer, "", "", RedBandNumber)  # Provide parameters
                    # for MakeRasterLayer tool in order to separate the NIR band from the four band composite image
                    out_redraster = 'redraster'   # Name CopyRaster output.  By default, ArcPro saves the output raster to the project's geodatabase
                    RedBand = arcpy.management.CopyRaster(RedRasterLayer, out_redraster)   # CopyRaster creates a permanent raster layer
                    # from the temporary layer created by the MakeRasterLayer tool above
                    RedRasterPath = os.path.join(arcpy.env.workspace, 'redraster')

                if len(str(NIRBandNumber)) != 0 and len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
                    rasters = [RedRasterPath, BlueRasterPath, NIRRasterPath]
                    input_names = ['x', 'y', 'z']
                    expression = "(y > z) & (y > x)"  # Band ratio formula
                    BandRatioNIR = arcpy.sa.RasterCalculator(rasters, input_names, expression)
                    output_raster = 'BRNIR'
                    BandRatioRasterNIR = arcpy.management.CopyRaster(BandRatioNIR, output_raster)
                    descBRRNIR = arcpy.Describe(BandRatioRasterNIR)
                    BRRNIRpath = descBRRNIR.catalogPath
                    BRRNIRbasename = descBRRNIR.baseName
                    binaryrasterpathlist.append(BRRNIRpath)
                    binrasterbnlist.append(BRRNIRbasename)

                    # Adding the band ratio binary raster created in the previous step to the otsu binary raster
                    # sometimes results in superior differentiation between water and shoreline impervious surfaces
                    rasters = [BRRNIRpath, otsurasterpath]
                    input_names = ['x', 'y']
                    expression = "x + y"
                    BRRNIRPlusOtsuNDWI = arcpy.sa.RasterCalculator(rasters, input_names, expression)
                    output_raster = basename + '_' + 'BRNIROtsuNDWI'
                    BRNIROtsuNDWIRaster = arcpy.management.CopyRaster(BRRNIRPlusOtsuNDWI, output_raster)
                    descBRNIRONDWIR = arcpy.Describe(BRNIROtsuNDWIRaster)
                    BRNIRONDWIRpath = descBRNIRONDWIR.catalogPath
                    BRNIRONDWIRbasename = descBRNIRONDWIR.baseName
                    binaryrasterpathlist.append(BRNIRONDWIRpath)
                    binrasterbnlist.append(BRNIRONDWIRbasename)

                # The SWIR band does not appear to perform well in the band ratio calculation by itself, but it sometimes yields
                # good results when added to the Otsu NDWI raster.  The code below contains some repetition of the NIR band ratio
                # calculation.  Some of the code in the following block could be removed for the sake of conciseness.
                if len(str(SWIRBandNumber)) != 0 and len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
                    rasters = [RedRasterPath, BlueRasterPath, SWIRRasterPath]
                    input_names = ['x', 'y', 'z']
                    expression = "(y > z) & (y > x)"
                    BandRatioSWIR = arcpy.sa.RasterCalculator(rasters, input_names, expression)
                    output_raster = 'BRSWIR'
                    BandRatioRasterSWIR = arcpy.management.CopyRaster(BandRatioSWIR, output_raster)
                    descBRRSWIR = arcpy.Describe(BandRatioRasterSWIR)
                    BRRSWIRpath = descBRRSWIR.catalogPath

                    rasters = [BRRSWIRpath, otsurasterpath]
                    input_names = ['x', 'y']
                    expression = "x + y"
                    BRRSWIRPlusOtsuNDWI = arcpy.sa.RasterCalculator(rasters, input_names, expression)
                    output_raster = basename + '_' + 'BRSWIROtsuNDWI'
                    BRSWIROtsuNDWIRaster = arcpy.management.CopyRaster(BRRSWIRPlusOtsuNDWI, output_raster)
                    arcpy.management.Delete(BandRatioRasterSWIR)
                    descBRSWIRONDWIR = arcpy.Describe(BRSWIROtsuNDWIRaster)
                    BRSWIRONDWIRpath = descBRSWIRONDWIR.catalogPath
                    BRSWIRONDWIRbasename = descBRSWIRONDWIR.baseName
                    binaryrasterpathlist.append(BRSWIRONDWIRpath)
                    binrasterbnlist.append(BRSWIRONDWIRbasename)

    binarycount = 0
    for binaryraster in binaryrasterpathlist:  # Loop through the binary rasters created above in order to
    # produce shoreline features and water body polygons in the steps below
        binarycount = binarycount + 1
        binrasterbncount = 0
        for binrasterbn in binrasterbnlist:  # For each binary raster produced earlier in the code, a catalog
        # path and a base name were added to their respective lists.  Therefore, place of a given raster in
        # both lists is precisely the same.  For example, the fourth binary raster in the binaryrasterpathlist
        # will correspond to the fourth binary raster base name in the binrasterbnlist.  Therefore, when the
        # count during the binaryrasterpathlist loop matches the count during the binrasterbnlist loop, the
        # basename and the catalog path both correspond to the same input raster.
            binrasterbncount = binrasterbncount + 1
            if binarycount == binrasterbncount:
                binrasterbasename = binrasterbn
                break
        # Conversion of the binary rasters to polygons is a step toward creating the water body polygons and
        # toward creating the shoreline features
        TempBinPoly = 'BinaryPolygon' + '_' + timestamp + '_' + str(imagecount) + '_' + str(binarycount)
        temp_binary_polygons = arcpy.conversion.RasterToPolygon(binaryraster, TempBinPoly, "NO_SIMPLIFY")

        where_clause = 'gridcode = 1'  # Upon converting a raster to a polygon, ArcGIS assigns a gridcode
        # to the resulting polygons.  Based on the binary rasters created by the shoreline digitization tool,
        # a gridcode of 1 corresponds to water.
        WaterPolygons = arcpy.management.SelectLayerByAttribute(temp_binary_polygons, "", where_clause)

        WPcursor = arcpy.SearchCursor(WaterPolygons)  # Establishing a search cursor on the WaterPolygons
        # attribute table allows for subsequent access of the Shape_Area field.  The next block of code
        # below serve to select the water polygon with the largest Shape_Area
        value = 0
        for row in WPcursor:
            rowvalue = row.getValue('Shape_Area')
            if rowvalue > value:
                value = rowvalue
        del WPcursor
        maxvalue = value

        if len(str(MinWaterSurfaceArea)) == 0:  # If no MinWaterSurfaceArea is input by the user, then only the largest water
        # polygon is used for shoreline digitization and for water body surface area estimation
            MaxWaterInt = int(maxvalue)   # The act of converting to an integer automatically rounds positive numbers down to the
            # nearest whole number.  The RasterToPolygon tool appears to automatically round its Shape_Area values to six decimal
            # digits, but the full value appears to contain up to nine digits.  Therefore, the attempt to SelectLayerByAttribute using
            # the exact WaterValue fails due to the value being rounded/truncated to six decimal digits.  A work-around involves the
            # MaxWater value down to the nearest integer then writing an SQL expression to select the value greater than this integer.
            MaxWaterIntString = str(MaxWaterInt)
            where_clause = 'Shape_Area > ' + MaxWaterIntString + ' AND gridcode = 1'   # Identify largest water polygon as the only feature
            # larger than int(MaxWater).  For some reason, the selection of gridcode =1 is lost prior to this step, so I have added
            # gridcode = 1 to the SQL where_clause once more.
            WaterBodySelection = arcpy.management.SelectLayerByAttribute(WaterPolygons, "", where_clause)
        else:  # If the user inputs a MinWaterSurfaceArea value, then all water polygons that are at least as large as MinWaterSurfaceArea
        # will be used for subsequent shoreline digitization and water body surface area estimation
            if maxvalue < MinWaterSurfaceArea:  # Some of the shoreline digitization tool's water/land distinction algorithms will result in
            # better land/water distinction for some input images than will other algorithms.  In cases where the algorithm poorly distinguishes
            # the land/water interface, the resultant polygons may be significantly smaller than the water bodies they represent.  In some cases,
            # no water polygons will be created at all.  In these cases, all water polygons that are smaller than MinWaterSurfaceArea will be
            # deleted and not used for shoreline digitization or water surface area estimation.
                arcpy.management.Delete(temp_binary_polygons)
                arcpy.management.Delete(WaterPolygons)
            if maxvalue >= MinWaterSurfaceArea:  # All water polygons larger than MinWaterSurfaceArea are selected
                MinWaterSurfaceAreaStr = str(MinWaterSurfaceArea)
                where_clause = 'Shape_Area > ' + MinWaterSurfaceAreaStr + ' AND gridcode = 1'
                WaterBodySelection = arcpy.management.SelectLayerByAttribute(WaterPolygons, "", where_clause)

        try:  # It is possible that no appropriate water polygons will be created for a given binary raster.  In such a case,
        # no features will be selected and saved to the WaterBodySelection variable.  If there are no features in the
        # WaterBodySelection variable, the CopyFeatures step below will fail.
            out_feature_class = 'WaterPolygonsofInterest' + '_' + timestamp + '_' + str(imagecount) + '_' + str(binarycount)
            WaterBodySelectionLayer = AddToTOC = arcpy.management.CopyFeatures(WaterBodySelection, out_feature_class)
            arcpy.management.Delete(temp_binary_polygons)
        except:
            continue

        # Convert the water polygons to a shoreline features
        out_line_feature = rasterbasename + '_' + binrasterbasename + '_' + 'shoreline' + '_' + AdditionToOutputFileName + '_' + timestamp
        Shoreline = arcpy.management.PolygonToLine(WaterBodySelectionLayer, out_line_feature)

        # The script tool's initial shoreline and water polygon outputs may be fragmented into more pieces than the number of water bodies they
        # actually represent.  Dissolving the multiple output features into a single feature yields a much cleaner output.  If the single dissolved
        # feature actually represents different water bodies, the user can simply use the clip tool to select the water body of interest.
        try:  # Attempting to dissolve when no feature is created by a given water/land detection algorithm will result in an
        # error.
            out_feature_class = rasterbasename + '_' + binrasterbasename + '_' + timestamp
            WaterPolygonsDissolved = arcpy.management.Dissolve(WaterBodySelectionLayer, out_feature_class)
            arcpy.management.Delete(WaterBodySelectionLayer)
            out_feature_class = rasterbasename + '_' + binrasterbasename + '_' + 'DissolveShoreline' + '_' + timestamp
            ShorelineDissolved = arcpy.management.Dissolve(Shoreline, out_feature_class)
            mergepolyinputlist.append(WaterPolygonsDissolved)
        except:
            continue

        # The code block below gets the DEMraster cell width and saves it to a variable to be used as the tolerance for the
        # simplifyLine tool. As I understand it, a tolerance set to cell width would mean that the simplifyline algorithm used
        # on an input raster with unique pixel values (i.e. pixelated) would yield a simplified line that is at least as accurate
        # as the actual physical values (that the raster represents) as the pixalated input raster.
        CellSizeX = arcpy.management.GetRasterProperties(binaryraster, "CELLSIZEX")
        tolerance = CellSizeX
        simplified_line = rasterbasename + '_' + binrasterbasename + '_' + 'Line'
        SimplifiedShoreline = arcpy.cartography.SimplifyLine(ShorelineDissolved, simplified_line, "", tolerance, "", "NO_KEEP")
        mergeinputlist.append(SimplifiedShoreline)  # Add the simplified shoreline for each image raster in the folder to a
        # list for later merging.  This merge will help the user to organize and visualize the final output as a single feature class.

        # Delete intermediate products
        arcpy.management.Delete(Shoreline)
        arcpy.management.Delete(ShorelineDissolved)

    # Delete all intermediate objects in binaryrasterpathlist
    for item in binaryrasterpathlist:
        arcpy.management.Delete(item)

    # Merge the features in mergeinputlist
    linefeatures = mergeinputlist
    output = rasterbasename + '_' + 'FinalShoreline' + '_' + AdditionToOutputFileName + '_' + timestamp
    AddToTOC = arcpy.management.Merge(linefeatures, output, "", "ADD_SOURCE_INFO")  # "ADD_SOURCE_INFO" retains the
    # identity of the utilized algorithm for each line feature in the merged feature class
    for item in mergeinputlist:
        arcpy.management.Delete(item)

    if imagecount < 6:
        AddOutput2TOC()

    polygonfeatures = mergepolyinputlist  # Merges polygon features.  Retaining the polygon features as output
    # allows the user to easily calculate water surface area through the "Shape_Area" field
    output = rasterbasename + '_' + 'FinalPolygon' + '_' + AdditionToOutputFileName + '_' + timestamp
    AddToTOC = arcpy.management.Merge(polygonfeatures, output, "", "ADD_SOURCE_INFO")
    for item in mergepolyinputlist:
        arcpy.management.Delete(item)

    if imagecount < 6:
        AddOutput2TOC()

    # Delete intermediate products from the output folder
    if arcpy.CheckExtension("ImageAnalyst") == "Available":
        arcpy.management.Delete(binaryrasterpathlist)
        arcpy.management.Delete(GreenBand)
        if len(str(SWIRBandNumber)) != 0 or len(str(NIRBandNumber)) != 0:
            arcpy.management.Delete(IRBand)
            if len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
                arcpy.management.Delete(BlueBand)
                arcpy.management.Delete(RedBand)
        if len(str(NIRBandNumber)) != 0:
            arcpy.management.Delete(NIRBand)
        if len(str(SWIRBandNumber)) != 0:
            arcpy.management.Delete(SWIRBand)
