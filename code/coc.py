#!/usr/bin/python
'''Data Source:
http://guides.lib.purdue.edu/c.php?g=353290&p=2378621#Mississippi
'''

import os
from osgeo import ogr, osr

# Import the necessary modules
epsgdic = {'nad83':4269,'wgs84':4326,'pseudoutm':3857,'worldmercater':3395,'wgs84':42310}
#home =r'/home/naresh-1/workspace/machinelrn/data/gis/' in linux
home = "D:/workspace/sqdm-repo/sqdm/"
#srcshpfile = home + r'uscounties/cb_2016_us_county_500k.shp'
#outputShapefile = home+ r'out/uscounties.shp'

home ="D:/workspace/sqdm-repo/sqdm/out/tmp/redist/census_blocks_by_states/tabblock2010_30_pophu/"
srcshpfile = home + r'tabblock2010_30_pophu.shp'
outputShapefile = home+ r'prjtabblock2010_30_pophu.shp'


driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(srcshpfile)

# Get Projection from layer
layer = shp.GetLayer()
sourceSpatialRef = layer.GetSpatialRef()
outgcsref= osr.SpatialReference()
outgcsref2= osr.SpatialReference()
outprjref = osr.SpatialReference()
outprjrefwm = osr.SpatialReference()
if sourceSpatialRef.IsProjected:
    datum=sourceSpatialRef.GetAttrValue('GEOGCS|DATUM')
    if 'North_American_Datum_1983' in datum:
        print("Geogcs/datum")
        #outgcsref.ImportFromEPSG(epsgdic['nad83']) #NAD83 --> WGS84/4326
        #outgcsref2.ImportFromEPSG(epsgdic['wgs84']) #NAD83 --> WGS84/4326
        #outprjref.ImportFromEPSG(epsgdic['pseudoutm']) #sthGcs -->Proj/3857,also see 3395(v. imp)
    outprjrefwm.ImportFromEPSG(epsgdic['wgs84'])
transform= osr.CoordinateTransformation(sourceSpatialRef,outprjrefwm)
print(transform.TransformPoint(-92.355389,41.509646))
print(outprjrefwm)
#saving layer
if os.path.exists(outputShapefile):
    driver.DeleteDataSource(outputShapefile)
outDataSet = driver.CreateDataSource(outputShapefile)

outLayer = outDataSet.CreateLayer("uscounty", geom_type=ogr.wkbMultiPolygon)
# add fields
inLayerDefn = layer.GetLayerDefn()
for i in range(0, inLayerDefn.GetFieldCount()):
    fieldDefn = inLayerDefn.GetFieldDefn(i)
    outLayer.CreateField(fieldDefn)

outLayerDefn = outLayer.GetLayerDefn()
inFeature = layer.GetNextFeature()

inputGeom = inFeature.GetGeometryRef()
inputGeom.TransformTo(outprjrefwm)
while inFeature:
    inputGeom = inFeature.GetGeometryRef()
    # get the cover attribute for the input feature
    #inputGeom.TransformTo(outgcsref)  #to nad83
    #inputGeom.TransformTo(outgcsref2) #to WGS84
    inputGeom.TransformTo(outprjrefwm) #to WGS84/Pseudo Mercator
    
    #if inFeature.GetField('STATEFP') == '28' and cnt==0:
    #inFeature=layer.GetNextFeature()
    
    outFeature = ogr.Feature(outLayerDefn)
    outFeature.SetGeometry(inputGeom)
    for i in range(0, outLayerDefn.GetFieldCount()):
        outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
    outLayer.CreateFeature(outFeature)
    inFeature = None
    inFeature = layer.GetNextFeature()

prjfile = open(os.path.splitext(outputShapefile)[0]+'.prj','w')
prjfile.write(outprjrefwm.ExportToPrettyWkt())
prjfile.close()
shp = None
outDataSet = None
print("completed cordinate transform.")
