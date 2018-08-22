
'''
Algorithm:
M = [(x1,y1) for all the polygons in map]
X = sorted([x1,x32,x3,...,xn]) #list of unique x-points of all the vertices of points in polygon-map M
P1 = [(p1,q2),(p2,q2)...] # all the points in polygon P1 
Bp1p2 = [xa,xb,xc] #between p1 and p2, there lies other x-ordinates in X.
Find points (xa,ya),(xb,yb),(xc,yc) on line (p1,q2) and (p2,q2) using y-mx formula.
for a tuple: (p1,xa,T1), find T1=[ya,yb]
             (xa,xb,T2), find T2 = [ya,ybS]
'''

import MapCommonSides101 as Mcs
import numpy as np
from  osgeo import ogr, osr

driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')
layer = shp.GetLayer()
neighbormap = Mcs.getNeighborsT(layer)

'''Test Data'''
table0= [(2,5),(3,4),(4,4),(5,3),(7,3),(11,4),(11,2),(13,2),(13,7),(10,10),(6,10)]
map0 = {}
map0[0] = [(0,10)]#polygon_A
map0['0'] = [(0,10)]#outerpolygon


'''Test Data'''
table1= [(2,5),(3,4),(4,4),(5,3),(7,3),(11,4),(11,2),(13,2),(13,7),(10,10),(6,10),
		(1,3),(1,1),(2,0),(4,1)]
map1 = {}
map1[0] = [(0,3),(4,10)]#polygon_A
map1[1]=[(0,3),(11,15)] #polygon_B
map1['0'] = [(4,15)]#outerpolygon




pointstable = table1
polygon_to_ptsindexmapping = map1
#construct a table consisting of polygon points, and map consisting of polygion and indices of the points in the table.
pointstable, polygon_to_ptsindexmapping = Mcs.constructSharedPointsDatabase(layer,neighbormap)
#save

tablefile = 'GisSharedPoints2.txt'
mapfile = 'GisPolygonPtsIndexMap2.txt'
import pickle
with open(tablefile, 'wb') as fp:
    pickle.dump(pointstable, fp)
with open(mapfile, 'wb') as fp:
    pickle.dump(polygon_to_ptsindexmapping, fp)
