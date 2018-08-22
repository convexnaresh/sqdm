
# coding: utf-8

# In[121]:

import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
import math
import time
import pickle
import numpy as np
import sys
import time
from osgeo import osr

def test(geom_lists):
    count = 0
    for p in geom_lists:
        print p.GetGeometryName(), p.GetBoundary().GetPoints()[0:1]
        #print count
        count +=1
    print count, len(geom_lists)

def getPolyObjForLinearRing(wkt0):
    polyref =ogr.CreateGeometryFromWkt('POLYGON ('+wkt0.strip('LINEARRING')+')')
    return polyref

def getPolygons(datatype='datatype',dofilter=False):
    all_polygons = [] #list of all the polygons
    feature_names = []
    fieldname = ''
    if datatype != 'sample':
        
	driver = ogr.GetDriverByName('ESRI Shapefile')
    	if datatype=='mscounty': #ok
            shp = driver.Open(r'../data/gis/ms_county/stco.shp')
	    layer = shp.GetLayer()
	    howmany = layer.GetFeatureCount()
    	    fieldname = 'CONAME'
	
        if datatype=='mszip_ge':
            shp = driver.Open(r'../data/gis/mszip/mszip10.shp')
            layer = shp.GetLayer()
            howmany= layer.GetFeatureCount()
            removeid = []
            idlist = [fid for fid in range(0,howmany) if fid not in removeid]
	    fieldname = 'CONAME'
                    
        if datatype=='mszip_pr':
            shp = driver.Open(r'../data/gis/mszip/mstm/mszip10.shp')    
            layer = shp.GetLayer()
            howmany= layer.GetFeatureCount()
            removeid = []
	    fieldname = 'CONAME'
            idlist = [fid for fid in range(0,howmany) if fid not in removeid]
                    
        if datatype == 'cong':#ok
            shp = driver.Open(r'../data/gis/congress/cb_2015_us_cd114_500k.shp')
            layer = shp.GetLayer()
            howmany= layer.GetFeatureCount()
            removeid = []
            idlist = [fid for fid in range(0,howmany) if fid not in removeid]
            fieldname = 'GEOID'

        if datatype == 'uscounty':#ok
            shp = driver.Open(r'../data/gis/uscounty/tl_2016_us_county.shp')
            layer = shp.GetLayer()
            howmany= layer.GetFeatureCount()
            removeid = []
            idlist = [fid for fid in range(0,howmany) if fid not in removeid]
            fieldname = 'NAME'
        
        if datatype == 'usstate': #ok
            shp = driver.Open(r'../data/gis/usstate/tl_2016_us_state.shp')
            layer = shp.GetLayer()
            howmany= layer.GetFeatureCount()
            removeid = [31,35,40]
            idlist = [fid for fid in range(0,howmany) if fid not in removeid]
            fieldname = 'NAME'
            
        #Initialize cordinates systems
        if layer.GetSpatialRef().ExportToWkt()[0:10].startswith("GEOGCS"):
            print("Detected LAT/LON cordinate system. Projecting to Northing/Easting.")
        
        inputSptialRef = layer.GetSpatialRef()
        projSptialRef= osr.SpatialReference()
        projSptialRef.ImportFromEPSG(3857) #3814  3857#mercator #GCS GWS84 4326
	
	required = [i for i in range(0,howmany)]
	if dofilter and datatype == 'usstate':
	    required=filtering(howmany,datatype)

        print("#feats"),howmany
        for fid in required: #[i for i in range(0,howmany)]:
            inFeature = layer.GetFeature(fid)
            geom = inFeature.GetGeometryRef() #it is polygon now
	    feature_names +=[inFeature.GetField(fieldname)]
            if inputSptialRef.ExportToWkt()[0:10].startswith("GEOGCS"):
                #convert to PROJCS
                coordTransform = osr.CoordinateTransformation(inputSptialRef,projSptialRef)
                geom.Transform(coordTransform)
                pass
            if geom.GetBoundary().GetPoints() is None:
		if dofilter and datatype in ['cong','uscounty'] :
			continue
                for inner_geom in geom:
                    #geom is MULTIPOLYGON
                    if inner_geom.GetGeometryName() == 'POLYGON':
                        #POLYGON has rings.
                        if inner_geom.GetBoundary().GetPoints() is None:
                            for linearring in polyref:
                                wkt0 = linearring.ExportToWkt()
                                polyref = getPolyObjForLinearRing(wkt0)
                                all_polygons += [polyref]
                        else:
                            wkt0 = inner_geom.ExportToWkt()
                            polyref = ogr.CreateGeometryFromWkt(wkt0)
                            all_polygons += [polyref]
                    #geom is POLYGON
                    elif inner_geom.GetGeometryName() == 'LINEARRING':
                        wkt0 = inner_geom.ExportToWkt()
                        polyref = getPolyObjForLinearRing(wkt0)
                        all_polygons += [polyref]
            else:
                wkt0 = geom.ExportToWkt()
                polyref = ogr.CreateGeometryFromWkt(wkt0)
                all_polygons +=[polyref]

        print "Data Source returned list of references to geometric objects.End."
	if datatype == 'cong':
	    feature_names =['CD-GEOID-'+str(name) for name in feature_names]
 
	return all_polygons,inputSptialRef,feature_names
    #sample polygon           
    else:
        #Manual Polygons
        wkt0 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        wkt1 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 4 1,2 0, 1 1,1 3,2 5 ))"
        wkt2 = "POLYGON ((5 3, 7 1, 10 1, 11 2,11 4,7 3,5 3 ))"
        poly0 = ogr.CreateGeometryFromWkt(wkt0)
        poly1 = ogr.CreateGeometryFromWkt(wkt1)
        poly2 = ogr.CreateGeometryFromWkt(wkt2)
        polyref = ogr.CreateGeometryFromWkt(wkt0)
                    
        all_polygons = [poly0,poly1,poly2] #list of all the polygonsA
        
        return all_polygons,None,['p0','p1','p2']


def filtering(howmany,datatype):
    reqdgeom = []
    removable = []
    if datatype== 'usstate':
    	removable1 = [1, 2, 6, 7, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 78, 79, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 142, 143, 144]
	removable = [1,5,13,14,31,33,34,35,38,40,41,49]
    for i in range(0,howmany):
    	if i not in removable:
            reqdgeom +=[i]
    return reqdgeom


if __name__ == "__main__":
    #Load Data.
    geom_lists,geoinputSptialRef,feature_names= getPolygons(datatype='cong',dofilter=True)
    print("hello"), len(geom_lists)
    for i in range(len(geom_lists)):#geom in geom_lists:
	print feature_names[i],geom_lists[i].GetGeometryName(),len(geom_lists[i].GetBoundary().GetPoints())


