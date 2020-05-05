#!/usr/bin/python

"""
extractViewAngle.py
Scope: export points or raster of viewing incidences angles from a Theia L2A product (rasters are scaled by 100 as UInt16) 
Author: simon.gascoin@cesbio.cnes.fr
"""

import csv
import gdal
import numpy as np
import ogr
import os
import osr
import sys
import xml.etree.ElementTree as ET


# function to read points file as lon lat values delimited by tab without header line
def readPoints(f):
    with open(f,'r') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        data = [r for r in reader]
        return data


# function to write points values as csv
def writePoints(newPointsFn,outDictList):
    with open(newPointsFn, 'w') as csvfile:
        fieldnames = list(outDictList[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for ouDict in outDictList:
            writer.writerow(ouDict)


# function to write an array to a (multiband) geotiff 
def array2geotiff(newRasterFn,geoTransform,array,noData,outSpatialRef,dataType=gdal.GDT_Float64):
    cols = array.shape[1]
    rows = array.shape[0]
    bands = array.shape[2]
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterFn, cols, rows, bands, dataType, options=['COMPRESS=DEFLATE'])
    outRaster.SetGeoTransform(geoTransform)
    # write bands
    for i in range(bands):
        outband = outRaster.GetRasterBand(i+1) # 1-based index
        outband.WriteArray(array[:,:,i])
        outband.SetNoDataValue(noData)
    outRaster.SetProjection(outSpatialRef.ExportToWkt())
    outRaster.FlushCache()


# function to get mask file name and bit number to test which detector was used
def getDetector(productFolder,root,bandId,detectorId):
    # find node containing detector metadata based on the presence of attribute "detector_id" in subnodes
    n = root.find(".//Product_Organisation//*[@detector_id]/..")
    if n is None:
        print('this product version does not provide detector mask')
        maskFn = bitNumber = None
    else:
        # get MASK_FILE element for target band and detector
        s = "./MASK_FILE/[@band_id='{}'][@detector_id='{}']".format(bandId,detectorId)
        element = n.find(s)
        # get detector mask file from element value 
        maskFn = os.path.join(productFolder,element.text)
        # get detector bit number from element attribute
        bitNumber = int(element.attrib['bit_number'])
    return maskFn, bitNumber


# function to test if detector was used at this point
def testDetector(point,maskFn,bitNumber):
    # open the raster file
    ds = gdal.Open(maskFn,gdal.GA_ReadOnly)
    if ds is None:
        print('Could not open the mask file')
        sys.exit(1)
    band = ds.GetRasterBand(1) # 1-based index
    data = band.ReadAsArray() # we could save memory and time by reading only the pixel using ReadRaster?
    geoTransform = ds.GetGeoTransform()
    # get position in array
    col,row = pix2map(point.GetX(),point.GetY(),geoTransform)
    # check if point is outside the mask
    if (col < 0 or row < 0 or col > band.XSize or row > band.YSize):
        print('Point is outside the product mask extent')
        test = False
    else:
        value = data[int(col)][int(row)]
        test = testBit(value, bitNumber)
    return test


# function which returns True if the bit number n is 1 in an integer value of base 10.
def testBit(value, n):
    mask = 1 << (n - 1) # bitNumber is 1-based index
    return(value & mask > 0)


# find position of x,y coordinates in georeferenced array with the same projection system
def pix2map(x,y,geoTransform):
    col = np.floor((x - geoTransform[0]) / geoTransform[1]) #x pixel
    row = np.floor((y - geoTransform[3]) / geoTransform[5]) #y pixel
    return col,row


# main function
def main(productFolder,outputFolder,points=None):
    # scale factor to export angles
    scale = 100
    # set no data value for UInt16 export
    noDataRaster = np.iinfo(np.uint16).max
    # set no data value for csv export
    noDataCsv = -10000    

    # MTD angle grid always have a 5 km resolution
    colstep = 5000
    rowstep = -5000
    # MTD angle grid always have an size of 23x23
    nx = ny = 23

    # open metadata file
    MTDFile = os.path.join(productFolder,os.path.basename(os.path.abspath(productFolder)+'_MTD_ALL.xml'))
    tree = ET.parse(MTDFile)
    root = tree.getroot()
    
    # get product id
    productId = root.find(".//PRODUCT_ID").text
    # get EPSG code
    epsg = root.find(".//HORIZONTAL_CS_CODE").text
    # get grid corners coordinates (warning in array geometry the lower left corner is the upper left in raster geometry)
    ulx = float(root.find(".//*[@name='upperLeft']/X").text)
    uly = float(root.find(".//*[@name='upperLeft']/Y").text)
    lrx = float(root.find(".//*[@name='lowerRight']/X").text)
    lry = float(root.find(".//*[@name='lowerRight']/Y").text)    

    # We assume that the above coordinates correspond to the *centers* of corner pixels
    # otherwise the 23x23 grid would have an extra row and column somewhere
    ulxMTD = ulx - colstep/2
    ulyMTD = uly - rowstep/2

    # define the affine transformation coefficients
    geoTransform = (ulxMTD, colstep, 0, ulyMTD, 0, rowstep)

    # create output spatial reference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(int(epsg))

    if points is not None:
        # create coordinate transformation
        inSpatialRef = osr.SpatialReference()
        inSpatialRef.ImportFromEPSG(4326)
        coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        #bandList = [element.text for element in root.find(".//*Band_Global_List").iter('BAND_ID')]

    # loop through angle definition
    for angle in ('Azimuth','Zenith'):

        # initialize output list of dictionnaries for points
        if points is not None:
            outDictList = list()
            [outDictList.append(dict()) for i in points]

        # loop through bands
        for band in root.iter('Band_Viewing_Incidence_Angles_Grids_List'):
            # init stack of grids
            Zd = np.array([], dtype=float).reshape(nx,ny,0)
            # loop through detectors
            for detector in band.iter('Viewing_Incidence_Angles_Grids'):
                rows = detector.find(angle).findall('.//VALUES')
                grid = ''
                # loop through grid rows to read grid values as a string
                for row in iter(rows):
                    grid = grid + row.text + '\n'
                # array with grid values 
                Z = np.fromstring(grid, dtype=float, sep=' ')
                # reshape to 2D array
                Z = Z.reshape((len(rows),-1))
                # add to the stack of detector grids
                Zd = np.dstack((Zd,Z))

            # display mean value for this angle and band
            bandId = band.attrib.get('band_id')
            print('{:s} {:s} mean value: {:g}'.format(bandId,angle,np.nanmean(Zd)))

            # export as multiband geotiff (we don't flatten the stack since the detector arrays overlap)
            if points is None:
                newRasterFn = os.path.join(\
                    outputFolder,'{:s}_{:s}_{:s}{:d}.tif'.format(productId,bandId,angle,scale))
                # scale 
                Zd = scale * Zd
                # set no data
                Zd[np.isnan(Zd)] = noDataRaster
                # write to disk
                array2geotiff(newRasterFn,geoTransform,Zd,noDataRaster,outSpatialRef,gdal.GDT_UInt16)

            # find values at points
            else:
                for ipoint,pointCoord in enumerate(points):
                    lon,lat = float(pointCoord[0]),float(pointCoord[1])
                    # create a geometry from coordinates
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(lat, lon)
                    # transform point
                    point.Transform(coordTransform)
                    # find position in array
                    col,row = pix2map(point.GetX(),point.GetY(),geoTransform)
                    # check if point is out of the grid
                    if (col < 0 or row < 0 or col > nx or row > ny):
                        v = noDataCsv

                    # otherwise retrieve the values in all bands
                    else:
                        vd = Zd[int(row),int(col),:]
                        # select the non-NaN value(s)
                        v = vd[np.isfinite(vd)]

                        # check if point is in no data area
                        if len(v) == 0:
                            v = noDataCsv

                        # check if more than one value is found in the stack 
                        # this can occur because angle grids overlap due to their coarse resolution
                        elif len(v) > 1:
                            print('solving an ambiguity for band = ' + bandId + ' at point ' + str(pointCoord))
                            detectorList = [d.attrib for d in band.iter('Viewing_Incidence_Angles_Grids')]
                            # indices where are the finite values
                            indexList = np.argwhere(np.isfinite(vd))
                            # look into the detector mask files to find which detector has measured this point 
                            test = False
                            for ix in indexList :
                                detectorId = detectorList[int(ix)]['detector_id']
                                print('testing detector = ' + detectorId)
                                maskFn,bitNumber = getDetector(productFolder,root,bandId,detectorId)
                                # if the detector mask file is provided then we assign the first value
                                if maskFn is None :
                                    print('takes first detector value by default')
                                    test = True
                                test = testDetector(point,maskFn,bitNumber)
                                if test: 
                                    print('found it!')
                                    v = vd[ix]
                                    break
                            # if test always false (point outside the mask) returns no data
                            if test is False: 
                                v = noDataCsv

                    # add this value to the output dictionnary 
                    if bandId in outDictList[ipoint]:
                        outDictList[ipoint][bandId].append(float(v))
                    else:
                        outDictList[ipoint][bandId] = float(v)

        # dump data to text file for this angle and band
        if points is not None:
            newPointsFn = os.path.join(\
                    outputFolder,'{:s}_{:s}.csv'.format(productId,angle))
            writePoints(newPointsFn,outDictList)

if __name__ == "__main__":

    # check arguments
    if len(sys.argv) == 4:
        print("Point mode")
        pointFile = sys.argv[3]
        # check if input file exists
        if not(os.path.exists(pointFile)):
            print("Error: input point file does not exists")
            sys.exit(1)
        points = readPoints(pointFile)    

    elif len(sys.argv) == 3:
        print("Raster mode")
        points = None

    else:
        print("Error: missing arguments\n")
        print("usage in raster mode: extractViewAngle.py productFolder outputFolder\n")
        print("usage in point mode: extractViewAngle.py productFolder outputFolder point_table_as_lon_lat.csv\n")
        print("example: python extractViewAngle.py SENTINEL2A_20180224-103018-463_L2A_T31TGK_C_V2-2 angles\n")
        print("example: python extractViewAngle.py SENTINEL2A_20180224-103018-463_L2A_T31TGK_C_V2-2 angles points.csv\n")
        sys.exit(1)

    # check if input file exists
    productFolder = sys.argv[1]
    if not(os.path.exists(productFolder)):
        print ("Error: input folder does not exists")
        sys.exit(1)

    # check if folder can be created
    outputFolder = sys.argv[2]
    try:
        os.makedirs(outputFolder,exist_ok=True)
    except OSError:
        print ("Error: cannot create output folder")
        sys.exit(1)
    else:
        main(productFolder,outputFolder,points)
 
