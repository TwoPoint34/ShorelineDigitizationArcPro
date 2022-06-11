#Github test
# Import required libraries and modules

import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy.ia import *
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

DEM_raster_path = arcpy.env.workspace = arcpy.GetParameterAsText(1)  # Re set workspace to access DEM rasters (if input by user)

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

outfilepath = arcpy.env.workspace = arcpy.GetParameterAsText(8)   # Change workspace so that output files are saved to a different
# location than the input raster folder.  Saving as a file geodatabase (versus shapefile) allows for easy water surface area
# measurement through the shape_area field.  Also, for some reason, the tool does not work properly if the workspace is set as
# a folder rather than a geodatabase.  It's worth noting too that the direction of the script tool parameter works only as input
# and not output.

ImageDEMWLTable = arcpy.GetParameterAsText(9)
# A DEM is the shoreline script tool's only acceptable input type for
# tide correction of the shoreline output.  The user may identify an existing DEM that approximately matches the time period
# and tide conditions of the aerial images being analyzed.  Alternatively, users may make their own DEMs using a variety of
# approaches.  For example, custom DEMs might be made through UAV-mounted LiDAR sensors or UAV-based stereo imagery.  A basic DEM
# could alternatively be created by interpolating a surface between measured beach elevation points.  Possible interpolation
# techniques for points include spline, kriging, and TIN creation. If the user has contour lines rather than elevation points,
# the ArcGIS Topo to Raster tool can interpolate a crude DEM.
# With this parameter, the user can input an Excel file that contains the relationship between the aerial image, the DEM, and
# water levels at the time of aerial image capture.  In this way, the user can input water levels under various tide (e.g. MHHW)
# or sea/water-level rise conditions.  The shoreline script tool will use the water level input value to create a contour line
# (shoreline) on the input DEM.  Assuming the water level at time of aerial image collection is above the DEM water level, the
# script tool derived aerial image shoreline can be compared to the DEM-based contour and thereby aid in the identification of
# shoreline change.  NOAA water level data can be found at:  https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels

ThresholdInputString = arcpy.GetParameterAsText(10)   # Allows user to specify the NDWI land/water threshold values, given that
# different values delineate shoreline more effectively for different aerial images.  Input format must be as follows (though
# values themselves can be changed):  0 0.1 0.2 0.3

if len(ImageDEMWLTable) != 0:   # Assuming the user inputs an Excel file, execute the code below
    descextension = arcpy.Describe(ImageDEMWLTable)  # Create a describe object of the Excel file to extract properties below
    WLVTextension = descextension.extension  # Extract the file extension of the input Excel file
    if WLVTextension == '.xls' or '.xlsx':  # Make sure that input spreadsheet has proper file extension
        Input_Excel_File = ImageDEMWLTable
        Output_Table = 'FGDBTable' + '_' + AdditionToOutputFileName + '_' + timestamp
        SpreadsheetToFileGDBTable = arcpy.conversion.ExcelToTable(Input_Excel_File, Output_Table)  # Convert Excel file to FGDB table
        FGDBsearchcursor = arcpy.SearchCursor(SpreadsheetToFileGDBTable)  # Establish search cursor for analysis of table
        scaledimagerasterlist = []  # Create list for analyzing only shared aerial image/DEM extent when a DEM is provided by user
        rowcount = 0
        for row in FGDBsearchcursor:  # Get values below from FGDB table
            rowcount = rowcount + 1
            imagebasenameobject = row.getValue('ImageBaseName')
            DEMbasenameobject = row.getValue('DEMBaseName')
            WaterLevel = row.getValue('WaterLevel')

            for DEMraster in DEMrasterlist:  # Access folder containing DEM rasters
                DEMraster = arcpy.sa.Raster(os.path.join(DEM_raster_path, DEMraster))
                descfolderDEM = arcpy.Describe(DEMraster)
                folderDEMbasename = descfolderDEM.baseName
                if folderDEMbasename == DEMbasenameobject:  # Match folder DEM file name with table DEM file name so that the
                # DEM file specified in the Excel/FGDB table is further processed by the subsequent code
                    DEMrasterpath = descfolderDEM.catalogPath  # Save the folder DEM file path to a variable
                    break  # Only a single DEM value is desired, so break DEMrasterlist

            for imageraster in imagerasterlist:  # Repeat process described for DEM above but this time for aerial image
                imageraster = arcpy.sa.Raster(os.path.join(image_raster_path, imageraster))
                descimageraster = arcpy.Describe(imageraster)
                imagerasterbasename = descimageraster.baseName
                if imagerasterbasename == imagebasenameobject:
                    imagerasterpath = descimageraster.catalogPath
                    break

            rasters = [imagerasterpath]   # Raster Calculator input is a list (even if it contains only a single raster)
            input_names = ['X']   # Raster Calculator input is a list of developer-defined variable names
            expression = "X / 10000000"  # Dividing the image raster by a very large number effectively guarantees that the
            # pixel values will all round down to zero, which makes for a nice integer raster and single polygon in the steps below
            Imageraster10millionth = arcpy.sa.RasterCalculator (rasters, input_names, expression)
            ImageExtentRasterInteger = Int(Imageraster10millionth)   # RasterToPolygon requires an integer raster input (e.g. "0")
            out_polygon_features = imagerasterbasename + '_' + 'TempImagePolygon_10mill' + '_' + timestamp
            ImageExtentPolygon = arcpy.conversion.RasterToPolygon(ImageExtentRasterInteger, out_polygon_features, "NO_SIMPLIFY")
            # Convert the integer raster to a single polygon that matches the extent of the input aerial image. Selecting
            # the "NO_SIMPLIFY" option appears to reduce information loss from the original aerial image

            # Create a DEM raster extent polygon in the same way the image extent polygon was created above
            demrasters = [DEMrasterpath]
            input_names = ['X']
            expression = "X / 10000000"
            DEMraster10millionth = arcpy.sa.RasterCalculator (demrasters, input_names, expression)
            DEMraster10millionthInt = Int(DEMraster10millionth)
            out_polygon_features = imagerasterbasename + '_' + 'TempDEMPolygon_10mill' + '_' + timestamp
            DEMExtentPolygon = arcpy.conversion.RasterToPolygon(DEMraster10millionthInt, out_polygon_features, "NO_SIMPLIFY")

            # Find the intersect of the aerial image extent and DEM extent such that all subsequent processing will only be
            # carried out on the intersect.  This step can significantly reduce tool processing time.  When, no Excel file in
            # provided, only the aerial image will be analyzed, and this analysis will include the entire image extent
            in_features = [DEMExtentPolygon, ImageExtentPolygon]
            out_feature_class = imagerasterbasename + '_' + "IntersectImageDEM" + '_' + timestamp
            ImageDEMIntersectPolygon = arcpy.analysis.PairwiseIntersect(in_features, out_feature_class)   # The polygon output
            # will be used as a mask for the ExtractByMask tool to extract shared image/DEM pixels

            imagerasterobject = arcpy.sa.Raster(imagerasterpath)  # Create raster object from raster path accessed from table row
            ExtractedImageraster = arcpy.sa.ExtractByMask(imagerasterobject, ImageDEMIntersectPolygon)  # Extract the image raster
            # pixels that correspond to the shared image/DEM extent
            out_rasterdata = imagerasterbasename
            SavedImageraster = arcpy.management.CopyRaster(ExtractedImageraster, out_rasterdata, "", "", "", "", "", "", "", "", "GRID")
            # Not sure if the above CopyRaster is necessary or not
            descSIr= arcpy.Describe(SavedImageraster)
            SIrCatPath = descSIr.catalogPath
            scaledimagerasterlist.append(SIrCatPath)  # Add the file path of the shared extent image raster to a list

            ExtractedDEMraster = arcpy.sa.ExtractByMask(DEMrasterpath, ImageDEMIntersectPolygon)  # Extract the DEM raster pixels
            # that correspond to the shared image/DEM extent

            # Delete intermediate products
            arcpy.management.Delete(ImageExtentPolygon)
            arcpy.management.Delete(DEMExtentPolygon)
            arcpy.management.Delete(ImageDEMIntersectPolygon)

            # Use the water level value from the input Excel/FGDB table to draw a contour line on the corresponding DEM
            out_polyline_features = imagerasterbasename + '_' + "GroundReferenceWaterLevel" + '_' + AdditionToOutputFileName + '_' + timestamp
            contour_value = [WaterLevel]  # Input value is a list (even if of a single item)
            GroundReferenceWaterLevel = ContourList(ExtractedDEMraster, out_polyline_features, contour_value)

            # In order to keep the ArcPro table of contents from getting too messy, only the results from the first five input
            # images will be automatically added.  The remaining results will be available in the output geodatabase.
            if rowcount < 6:
                desc = arcpy.Describe(GroundReferenceWaterLevel)
                contourfilepath = desc.catalogPath   # Get the file path of the contour line feature
                aprxMap = aprx.listMaps("Map")   # Access the current ArcPro map
                for mappy in aprxMap:
                    mappy.addDataFromPath(contourfilepath)   # Add the contour for each table row to the map table of contents

        del FGDBsearchcursor  # Delete the cursor just as a matter of best practice
        arcpy.management.Delete(SpreadsheetToFileGDBTable)
        imagerasterlist = scaledimagerasterlist  # The imagerasterlist below loops through image rasters whether or not the user
        # has input a DEM.  With a DEM, imagerasterlist is scaled (code above) to match the shared image/DEM extent.  Without the
        # DEM, the imagerasterlist was created simply by looping through the input image folder.

imagecount = 0
for raster in imagerasterlist:
    binaryrasterpathlist = []  # Create list of binary raster paths for subsequent looping and shoreline analysis
    binrasterbnlist = []  # Create list of binary raster base names, which will be used to match table inputs to files in folder
    mergepolyinputlist = []  # Create list to be populated with output polygon features and later merged
    mergeinputlist = []   # Create list to be populated with output line features and later merged
    imagecount = imagecount + 1    # Add "1" to the count value for each iteration of the loop (i.e. each raster)
    if len(ImageDEMWLTable) != 0:
        imageraster = arcpy.sa.Raster(raster)  # If the user inputs an Excel file, access the raster and save to variable
    else:
        imageraster = arcpy.sa.Raster(os.path.join(image_raster_path, raster))  # If no Excel, image rasters accessed from folder

    descimageraster = arcpy.Describe(imageraster)
    rasterbasename = descimageraster.baseName

    # As an essential component of the (M)NDWI calculation, the green band is the only mandatory band
    Green_rasterlayer = 'GreenBand'
    GreenRasterLayer = arcpy.management.MakeRasterLayer(imageraster, Green_rasterlayer, "", "", GreenBandNumber)  # Use GreenBandNumber
    # parameter to make a green band layer from the multi-band image
    out_greenraster = 'greenraster'   # By default, ArcPro saves the output raster to the project's geodatabase
    GreenBand = arcpy.management.CopyRaster(GreenRasterLayer, out_greenraster)   # CopyRaster creates a permanent raster layer
    # from the temporary layer created by the MakeRasterLayer tool

    if len(str(NIRBandNumber)) != 0:  # If the user inputs a NIR band number
        NIR_rasterlayer = 'NIRBand'
        NIRRasterLayer = arcpy.management.MakeRasterLayer(imageraster, NIR_rasterlayer, "", "", NIRBandNumber)
        out_nirraster = 'nirraster'
        NIRBand = arcpy.management.CopyRaster(NIRRasterLayer, out_nirraster)
        NIRRasterPath = os.path.join(arcpy.env.workspace, 'NIRraster')  # For later usage in the RasterCalculator tool, it works
        # better to save the raster path than trying to use the raster object (NIRBand) itself

    if len(str(SWIRBandNumber)) != 0:
        SWIR_rasterlayer = 'SWIRBand'
        SWIRRasterLayer = arcpy.management.MakeRasterLayer(imageraster, SWIR_rasterlayer, "", "", SWIRBandNumber)
        out_swirraster = 'SWIRraster'
        SWIRBand = arcpy.management.CopyRaster(SWIRRasterLayer, out_swirraster)
        SWIRRasterPath = os.path.join(arcpy.env.workspace, 'SWIRraster')
        IRBand = SWIRBand  # If a SWIR band is input, then the MNDWI will be calculated
    else:
        IRBand = NIRBand  # If no SWIR band is input, NDWI will be calculated

    # Formula for Normalized Difference Water Index is:  NDWI = (GreenBand - NIRBand) / (GreenBand + NIRBand)
    num = arcpy.sa.Float(GreenBand) - arcpy.sa.Float(IRBand)   # Define NDWI numerator. Convert to Float to enable calculation.
    denom = arcpy.sa.Float(GreenBand) + arcpy.sa.Float(IRBand)   # Define NDWI denominator. Convert to Float to enable calculation.
    NDWI = arcpy.sa.Divide(num, denom)  # Calculate NDWI by dividing numerator by denominator.

    out_NDWI_copy = 'out_NDWI_copy' + '_' + timestamp
    NDWIcopy = arcpy.management.CopyRaster(NDWI, out_NDWI_copy)
    descNDWI = arcpy.Describe(NDWIcopy)
    NDWIpath = descNDWI.catalogPath   # For some reason, the NDWI raster file path was getting lost when the script attempted
    # to create shorelines from both the Otsu threshold and the user input thresholds.  The script tool returned an error
    # (ERROR 160333: The table was not found).  Writing the catalogPath to a variable somehow bypasses this error.  I'm not
    # sure if this step is still necessary.

    Otsu_raster = arcpy.ia.Threshold(NDWI)  # I got the idea to use Otsu thresholding from the CoastSat tool
    if len(str(SWIRBandNumber)) != 0:
        out_otsu = 'OtsuMNDWI'
    else:
        out_otsu = 'OtsuNDWI'
    OtsuRasterCopy = arcpy.management.CopyRaster(Otsu_raster, out_otsu)
    descotsuraster = arcpy.Describe(OtsuRasterCopy)
    otsurasterpath = descotsuraster.catalogPath
    otsurastername = descotsuraster.baseName
    binaryrasterpathlist.append(otsurasterpath)
    binrasterbnlist.append(otsurastername)

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
        expression = "(y > z) & (y > x)"
        BandRatioNIR = arcpy.sa.RasterCalculator(rasters, input_names, expression)
        output_raster = 'BRNIR' + '_' + str(imagecount)
        BandRatioRasterNIR = arcpy.management.CopyRaster(BandRatioNIR, output_raster)
        descBRRNIR = arcpy.Describe(BandRatioRasterNIR)
        BRRNIRpath = descBRRNIR.catalogPath

        rasters = [BRRNIRpath, otsurasterpath]
        input_names = ['x', 'y']
        expression = "x + y"
        BRRNIRPlusOtsuNDWI = arcpy.sa.RasterCalculator(rasters, input_names, expression)
        output_raster = 'BRNIROtsuNDWI' + '_' + str(imagecount)
        BRNIROtsuNDWIRaster = arcpy.management.CopyRaster(BRRNIRPlusOtsuNDWI, output_raster)
        descBRNIRONDWIR = arcpy.Describe(BRNIROtsuNDWIRaster)
        BRNIRONDWIRpath = descBRNIRONDWIR.catalogPath
        BRNIRONDWIRbasename = descBRNIRONDWIR.baseName
        binaryrasterpathlist.append(BRNIRONDWIRpath)
        binrasterbnlist.append(BRNIRONDWIRbasename)

    # The SWIR band does not appear to perform well in the band ratio calculation by itself, but it sometimes yields
    # good results when added to the Otsu NDWI raster
    if len(str(SWIRBandNumber)) != 0 and len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
        rasters = [RedRasterPath, BlueRasterPath, SWIRRasterPath]
        input_names = ['x', 'y', 'z']
        expression = "(y > z) & (y > x)"
        BandRatioSWIR = arcpy.sa.RasterCalculator(rasters, input_names, expression)
        output_raster = 'BRSWIR' + '_' + str(imagecount)
        BandRatioRasterSWIR = arcpy.management.CopyRaster(BandRatioSWIR, output_raster)
        descBRRSWIR = arcpy.Describe(BandRatioRasterSWIR)
        BRRSWIRpath = descBRRSWIR.catalogPath

        rasters = [BRRSWIRpath, otsurasterpath]
        input_names = ['x', 'y']
        expression = "x + y"
        BRRSWIRPlusOtsuNDWI = arcpy.sa.RasterCalculator(rasters, input_names, expression)
        output_raster = 'BRSWIROtsuNDWI' + '_' + str(imagecount)
        BRSWIROtsuNDWIRaster = arcpy.management.CopyRaster(BRRSWIRPlusOtsuNDWI, output_raster)
        descBRSWIRONDWIR = arcpy.Describe(BRSWIROtsuNDWIRaster)
        BRSWIRONDWIRpath = descBRSWIRONDWIR.catalogPath
        BRSWIRONDWIRbasename = descBRSWIRONDWIR.baseName
        binaryrasterpathlist.append(BRSWIRONDWIRpath)
        binrasterbnlist.append(BRSWIRONDWIRbasename)

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
            out_binraster = 'NDWIThreshold' + '_' + str(count1)
            NDWIbinaryCopy = arcpy.management.CopyRaster(NDWIbinary, out_binraster)
            descbinraster = arcpy.Describe(NDWIbinaryCopy)
            binrasterpath = descbinraster.catalogPath
            brbasename = descbinraster.baseName
            binaryrasterpathlist.append(binrasterpath)
            binrasterbnlist.append(brbasename)

# GET BASENAME FOR EACH BINARY RASTER INPUT AND ADD TO OUTPUT FILE NAMES
    binarycount = 0
    for binaryraster in binaryrasterpathlist:
        binarycount = binarycount + 1
        binrasterbncount = 0
        for binrasterbn in binrasterbnlist:
            binrasterbncount = binrasterbncount + 1
            if binarycount == binrasterbncount:
                binrasterbasename = binrasterbn
                break
        TempBinPoly = 'BinaryPolygon' + '_' + timestamp + '_' + str(imagecount) + '_' + str(binarycount)
        temp_binary_polygons = arcpy.conversion.RasterToPolygon(binaryraster, TempBinPoly, "NO_SIMPLIFY")

        # The script tool assumes that the largest water body is the water body the user wants the shoreline
        # to be delineated for
        MaxCursor = arcpy.SearchCursor(temp_binary_polygons)
        value = 0
        for row in MaxCursor:
            rowvalue = row.getValue('Shape_Area')
            if rowvalue > value:
                value = rowvalue
        del MaxCursor
        MaxWater = value
        MaxWaterInt = int(MaxWater)   # The act of converting to an integer automatically rounds positive numbers down to the
            # nearest whole number.  The RasterToPolygon tool appears to automatically round its Shape_Area values to six decimal
            # digits, but the full value appears to contain up to nine digits.  Therefore, the attempt to SelectLayerByAttribute using
            # the exact WaterValue fails due to the value being rounded/truncated to six decimal digits.  A work-around involves the
            # MaxWater value down to the nearest integer then writing an SQL expression to select the value greater than this integer.

        where_clause = 'gridcode = 1'   # A gridcode field is added when a raster is converted to a polygon.  The (M)NDWI-derived
        # binary rasters have all been reclassified such that values below the threshold (i.e. land) become "0" and values above the
        # threshold (i.e. water) result in a gridcode of "1" corresponds to water in the binary raster
        Water_Selection = arcpy.management.SelectLayerByAttribute(temp_binary_polygons, "", where_clause)  # Select all water polygons

        MaxWaterIntString = str(MaxWaterInt)
        where_clause = 'Shape_Area > ' + MaxWaterIntString   # Identify largest water polygon as only feature larger than int(MaxWater)
        Largest_Water_Selection = arcpy.management.SelectLayerByAttribute(Water_Selection, "", where_clause)

            # Eliminate any polygon part that is smaller than 99% the area of the water polygon with largest area.  An alternative to this
            # approach would be to create a layer from the selected feature.
        out_feature_class = 'ElimPolyPart' '_' + timestamp + '_' + str(imagecount)
        WaterPolygonSelection = arcpy.management.EliminatePolygonPart(Largest_Water_Selection, out_feature_class, "PERCENT", "", 99, "ANY")

        out_path = arcpy.env.workspace
        out_name = rasterbasename + '_' + binrasterbasename + '_' + 'Polygon'
        SavedWaterPolygonSelection = arcpy.conversion.FeatureClassToFeatureClass(WaterPolygonSelection, out_path, out_name)
        mergepolyinputlist.append(SavedWaterPolygonSelection)

        # Convert the largest water polygon to a shoreline feature
        out_line_feature = rasterbasename + '_' + binrasterbasename + '_' + 'shoreline' + '_' + AdditionToOutputFileName + '_' + timestamp
        Shoreline = arcpy.management.PolygonToLine(SavedWaterPolygonSelection, out_line_feature)

        # Since we're only interested in a single shoreline feature, we can dissolve all Shoreline features into a single feature
        out_feature_class = rasterbasename + '_' + binrasterbasename + '_' + 'DissolveShorelineMask' + '_' + timestamp
        ShorelineDissolved = arcpy.management.Dissolve(Shoreline, out_feature_class)

        # Gets the DEMraster cell width and saves it to a variable to be used as the tolerance for the simplifyLine tool.
        # As I understand it, a tolerance of cell width would mean that the SimplifyLine algorithm used on an input raster
        # with unique pixel values (i.e. pixelated) would yield a simplified line that is at least as accurate as the
        # actual physical values (that the raster represents) as the pixalated input raster.
        CellSizeX = arcpy.management.GetRasterProperties(binaryraster, "CELLSIZEX")
        tolerance = CellSizeX
        simplified_line = rasterbasename + '_' + binrasterbasename + '_' + 'Line'
        SimplifiedShoreline = arcpy.cartography.SimplifyLine(ShorelineDissolved, simplified_line, "", tolerance, "", "NO_KEEP")
        mergeinputlist.append(SimplifiedShoreline)  # Add the simplified shoreline for each image raster in the folder to a
        # list for later merging.  This merge will help the user to organize and visualize the final output.

        # Delete intermediate products
        arcpy.management.Delete(temp_binary_polygons)
        arcpy.management.Delete(Shoreline)
        arcpy.management.Delete(ShorelineDissolved)

    linefeatures = mergeinputlist  # Merge the features in the merge list created previously in the binaryraster loop
    output = rasterbasename + '_' + 'FinalShoreline' + '_' + AdditionToOutputFileName + '_' + timestamp
    MergeSimpShore = arcpy.management.Merge(linefeatures, output, "", "ADD_SOURCE_INFO")  # Adding Source Info retains the
    # file names of the shoreline features being merged.  In this way, the tool user will be able to identify the binary
    # raster (e.g. OtsuNDWI, Threshold, etc.) used in the creation of the shoreline

    polygonfeatures = mergepolyinputlist  # Merges polygon features like line features are merged above.  Retaining the
    # polygon features as output allows the user to easily calculate water surface area through the "Shape_Area" field
    output = rasterbasename + '_' + 'FinalPolygon' + '_' + AdditionToOutputFileName + '_' + timestamp
    MergeSimpPoly = arcpy.management.Merge(polygonfeatures, output, "", "ADD_SOURCE_INFO")

    # Get the file path of the shoreline features
    if imagecount < 6:
        desc = arcpy.Describe(MergeSimpShore)
        simplinefilepath = desc.catalogPath   # Get the file path of simplified shoreline features
        aprxMap = aprx.listMaps("Map")
        for mappy in aprxMap:
            mappy.addDataFromPath(simplinefilepath)  # Add the shoreline for each row to the map table of contents

    # Get the file path of the water polygon features
    if imagecount < 6:
        desc = arcpy.Describe(MergeSimpPoly)
        contourfilepath = desc.catalogPath
        aprxMap = aprx.listMaps("Map")
        for mappy in aprxMap:
            mappy.addDataFromPath(contourfilepath)

# Delete intermediate products from the output folder
    arcpy.management.Delete(WaterPolygonSelection)
    arcpy.management.Delete(SavedWaterPolygonSelection)
    arcpy.management.Delete(binaryrasterpathlist)
    arcpy.management.Delete(GreenBand)
    arcpy.management.Delete(NDWIcopy)
    arcpy.management.Delete(OtsuRasterCopy)
    arcpy.management.Delete(mergeinputlist)
    arcpy.management.Delete(mergepolyinputlist)
    if len(ThresholdInputString) != 0:
        arcpy.management.Delete(NDWIbinaryCopy)
    if len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
        arcpy.management.Delete(IRBand)
        arcpy.management.Delete(BlueBand)
        arcpy.management.Delete(RedBand)
    if len(str(NIRBandNumber)) != 0:
        arcpy.management.Delete(NIRBand)
    if len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0 and len(str(NIRBandNumber)) != 0:
        arcpy.management.Delete(BandRatioRasterNIR)
        arcpy.management.Delete(BRNIROtsuNDWIRaster)
    if len(str(SWIRBandNumber)) != 0:
        arcpy.management.Delete(SWIRBand)
    if len(str(SWIRBandNumber)) != 0 and len(str(BlueBandNumber)) != 0 and len(str(RedBandNumber)) !=0:
        arcpy.management.Delete(BandRatioRasterSWIR)
        arcpy.management.Delete(BRSWIROtsuNDWIRaster)

if len(ImageDEMWLTable) != 0:
    arcpy.management.Delete(SavedImageraster)
arcpy.management.Delete(imagerasterlist)
arcpy.management.Delete(binrasterbnlist)
