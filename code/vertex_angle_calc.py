
# coding: utf-8

# In[17]:

import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
from descartes import PolygonPatch
import treeformation101 as Tf
import MapCommonSides101 as Mcs
import math
import pickle

get_ipython().magic(u'matplotlib inline')


# In[64]:

get_ipython().magic(u'matplotlib inline')
def getPolygonPoints(datatype='sample', pid=0, debug = False):    
    saveFiles = True
    all_polygons = [] #list of all the polygons
    polygon = []
    #======================================================
    if datatype =="sample":
        #Manual Polygons
        wkt0 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        wkt1 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 4 1,2 0, 1 1,1 3,2 5 ))"
        wkt2 = "POLYGON ((5 3, 7 1, 10 1, 11 2,11 4,7 3,5 3 ))"
        poly0 = ogr.CreateGeometryFromWkt(wkt0)
        poly1 = ogr.CreateGeometryFromWkt(wkt1)
        poly2 = ogr.CreateGeometryFromWkt(wkt2)
        all_polygons = [poly0,poly1,poly2] #list of all the polygons
        polygon = all_polygons[pid]
    #======================================================
    if datatype == "ms":
        #Example Code for MS state counties map
        driver = ogr.GetDriverByName('ESRI Shapefile')
        shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')
        polygon_layer = shp.GetLayer()
        feature = polygon_layer.GetFeature(pid)
        polygon = feature.GetGeometryRef()
        
        ##PLOTTING
    if debug:
        fig = plt.figure()
        fig.set_figwidth(12)
        fig.set_figheight(12)
        ax1 = fig.add_subplot(231)
        colors=['r','b','k','g']
        poly = polygon.GetBoundary().GetPoints()
        x = [a for a,b in poly]
        y = [b for a,b in poly]    
        ax1.plot(x,y, '-',color=colors[0])

        ax1 = plt.gca()
        maxv = max(max(x),max(y))
        minv = min(min(x), min(y))
        ax1.set_xlim([minv,maxv])
        ax1.set_ylim([minv,maxv]) 

        plt.show()

    return polygon.GetBoundary().GetPoints()


def getPolygonPoints(datatype='sample', pid=0, debug = False):    
    saveFiles = True
    all_polygons = [] #list of all the polygons
    polygon = []
    #======================================================
    if datatype =="sample":
        #Manual Polygons
        wkt0 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        wkt1 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 4 1,2 0, 1 1,1 3,2 5 ))"
        wkt2 = "POLYGON ((5 3, 7 1, 10 1, 11 2,11 4,7 3,5 3 ))"
        poly0 = ogr.CreateGeometryFromWkt(wkt0)
        poly1 = ogr.CreateGeometryFromWkt(wkt1)
        poly2 = ogr.CreateGeometryFromWkt(wkt2)
        all_polygons = [poly0,poly1,poly2] #list of all the polygons
        polygon = all_polygons[pid]
    #======================================================
    if datatype == "ms":
        #Example Code for MS state counties map
        driver = ogr.GetDriverByName('ESRI Shapefile')
        shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')
        polygon_layer = shp.GetLayer()
        feature = polygon_layer.GetFeature(pid)
        polygon = feature.GetGeometryRef()
        
        ##PLOTTING
    if debug:
        fig = plt.figure()
        fig.set_figwidth(12)
        fig.set_figheight(12)
        ax1 = fig.add_subplot(231)
        colors=['r','b','k','g']
        poly = polygon.GetBoundary().GetPoints()
        x = [a for a,b in poly]
        y = [b for a,b in poly]    
        ax1.plot(x,y, '-',color=colors[0])

        ax1 = plt.gca()
        maxv = max(max(x),max(y))
        minv = min(min(x), min(y))
        ax1.set_xlim([minv,maxv])
        ax1.set_ylim([minv,maxv]) 

        plt.show()

    return polygon.GetBoundary().GetPoints()


def intersectionComputed(intersection_computed_tbl, polyid1,polyid2):
    '''Returns true if intersection is computed between polyid1 and polyid2'''
    debug=False
    if debug: print("p1,p2"), polyid1,polyid2
    if polyid1 in intersection_computed_tbl[polyid2]:
        return True
    return False

def getLinesForPolygonPointsSequences2(ptseq):
    #if close = True then constructs a closed polygon.
    plinesequenceList = []
    lineseq = []
    howmany = len(ptseq)
    for ptid in range(howmany-1):
        point1,point2 = ptseq[ptid], ptseq[ptid+1]
        line1 = ogr.Geometry(ogr.wkbLineString)
        line1.AddPoint(point1[0],point1[1])
        line1.AddPoint(point2[0],point2[1])
        lineseq +=[line1]
    return lineseq
def calcCosine(point1,midpoint,point3):
    '''Tuples of point'''
    x1,y1,z = point1
    x2,y2,z = midpoint
    x3,y3,z = point3
    mag_a = abs(math.sqrt((x2-x1)**2 + (y2-y1)**2))
    mag_b = abs(math.sqrt((x3-x2)**2 + (y3-y2)**2))
    dotAB = (x2-x3) *(x2-x1) + (y2-y3) * (y2-y1)
    try:
        theta = math.degrees(math.acos(dotAB/float(mag_a * mag_b))) #math.acos(1)
        return theta
    except:
        return 'inff'

def calcVertexAngle(geosLineString1, geosLineString2):
    
    point1, point2 = geosLineString1.GetPoints()
    point3, point4 = geosLineString2.GetPoints()
    return point2, calcCosine(point1,point2,point4)
    
def calcAngleBetweenLines(geosLineString1, geosLineString2):
    import math
    angle = 0.0
    a =  Mcs.calcSlope(geosLineString1) #slope
    b =  Mcs.calcSlope(geosLineString2) #slope
    
    try:
        angle1 = math.degrees(math.atan(a))
        if angle1 < 0: angle1 = 180 - abs(angle1)
    except:angle1 = 90  #if a == 'inf'
        
    try:
        angle2 = math.degrees(math.atan(b))
        if angle2 < 0: angle2 = 180 - abs(angle2)
    except: angle2 = 90 #if b == 'inf'
    
    print geosLineString1, geosLineString2,'\t',a,b, angle1,angle2
    
    if a >= 0 and b >= 0:
        return angle1 - angle2
    elif a >=0 and b <= 0:
        return angle2 - angle1 
    elif a <= 0 and b >= 0:
        return angle1- angle2
    elif a <= 0 and b <= 0:
        return angle2 - angle1
    return 'abc'

def removeEmptyList(listlist,inplace=True):
    '''Removes one level empty list from list  of lists.'''
    if inplace:
        for listItem in listlist:
            if not listItem:
                listlist.remove(listItem)
        return None
    else:
        return [ listItem for listItem in listlist if listItem]
def getGeosPolyObjFromPolyPoints(polygonpoints):
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x,y in polygonpoints:
        ring.AddPoint(x,y)
    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly

def countAcute(angles):
    acutecount = 0
    for angle in angles:
        if angle <= 90:
            acutecount +=1
    return acutecount
def getNeighborsT(polygons):
    '''It is more efficient version of getNeighborsD function. Returns a list of lists containing ids of polygon that are connected to a polygon with id at that index.
For example, list at index 0 contains polygons that are neighbors to polygion-0.
@polygion_layer is a osgeo layer. For example.
polygon_layer = shp.GetLayer()'''
    neighbormap = [[] for i in range(len(polygons))]
    howmany = len(polygons)
    for ffid in range(howmany): #this can be optimized
        poly =  polygons[ffid]
        for fid in range(len(polygons)):
            #escape if already computed
            if ffid != fid and (ffid not in neighbormap[fid]):
                candpoly = polygons[fid] #feature
                dist = poly.Touches(candpoly)
                if dist == True:
                    neighbormap[ffid] += [fid]
                    neighbormap[fid] += [ffid]
    return neighbormap

def getGeosPolyObjFromPolyPoints(polygonpoints):
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x,y,z in polygonpoints:
        ring.AddPoint(x,y)
    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly

def getNeighborsT2(poly_points_lists):
    polygons = []
    for polyptslist in poly_points_lists:
        polygons += [getGeosPolyObjFromPolyPoints(polyptslist)]
    return getNeighborsT(polygons)

def makeItAcute(angles):
    for i in range(0,len(angles)):
        if angles[i] >= 90:
            angles[i] = 360 - angles[i]
            
def GetAngles(pid, datatype='sample'):
    import itertools, sys,time,pickle
    ppointssequenceList = [getPolygonPoints(datatype=datatype, pid=pid, debug=False)]
    #Convert each sides into Geos LineString objects
    listLineStringObjs= []#list of Geos LineString
    
    for ppsequence in ppointssequenceList:
        lineStringObjl = getLinesForPolygonPointsSequences2(ppsequence)
        listLineStringObjs.append(lineStringObjl)
    
    removeEmptyList(listLineStringObjs) #remove [] items
    linePairs =[]
    for sideStringGeos in listLineStringObjs:
        linePairs += Tf.getConsequitivePairs(sideStringGeos,winding=True)
    angles = []
    vertices = []
    for linepair in linePairs:
        geoslinea,geoslineb = linepair
        vertex, angbetwn = calcVertexAngle(geoslinea, geoslineb)
        vertices +=[vertex]
        angles +=[angbetwn]
    return vertices, angles


def GetAnglesOfPolygon(listLineStringObjs):
    import itertools, sys,time,pickle    
    removeEmptyList(listLineStringObjs) #remove [] items
    linePairs =[]
    linePairs += Tf.getConsequitivePairs(listLineStringObjs,winding=True)
    angles = []
    vertices = []
    for linepair in linePairs:
        geoslinea,geoslineb = linepair
        vertex, angbetwn = calcVertexAngle(geoslinea, geoslineb)
        vertices +=[vertex]
        angles +=[angbetwn]
    return vertices, angles


def main(debug = False):
    #For all desired polygons, find angles formed at all vertices. 
    datatype = 'sample'
    if datatype == 'sample':howmany = 3
    else: howmany = 85
    vertexLists = [] #[[vertices for polygon-1],[vertices for polygon-2]...]
    anglesLists = [] #[[angles for polygon-1],[angles for polygon-2]...]

    #Foe each polygon, get vertices and angles at each vertices
    for pid in range(0,howmany):
        vertices, angles = GetAngles(pid, datatype=datatype)
        vertexLists.append(vertices+[vertices[0]])
        anglesLists.append(angles)

    #Find neighbors of each polygons
    neighbormap = getNeighborsT2(vertexLists)

    #find removable angles' indices for each polygon
    removeAnglesAt = [[] for i in range(howmany)]
    intersection_tbl = [[] for i in range(howmany)]
    for pid in range(howmany):
        polygonPts = vertexLists[pid]
        neighbors = neighbormap[pid]
        for each_neighbor in neighbors:
            if not intersectionComputed(intersection_tbl,pid,each_neighbor):
                neighborpolygonPts = vertexLists[each_neighbor]
                for ptIndex in range(len(polygonPts)-1):
                    if polygonPts[ptIndex] in neighborpolygonPts:
                        #Check if that point forks 
                        if polygonPts[ptIndex+1] in neighborpolygonPts and polygonPts[ptIndex-1] in neighborpolygonPts:
                            removeAnglesAt[pid].append(ptIndex)
                intersection_tbl[pid] +=[each_neighbor]
                intersection_tbl[each_neighbor] +=[pid]
    print removeAnglesAt

    #Remove all the duplicate angles and get filtered angles lists.
    from operator import itemgetter 
    filteredAngleLists = [[] for i in range(len(anglesLists))]
    for polyIdx in range(len(anglesLists)):
        remainingAnglesAt =[]
        for i in range(len(anglesLists[polyIdx])):
            if i not in removeAnglesAt[polyIdx]:
                remainingAnglesAt +=[i]   
        if remainingAnglesAt[polyIdx]:
            filteredAngleLists[polyIdx] = list(itemgetter(*remainingAnglesAt)(anglesLists[polyIdx]))

    filteredAngleLists

    #Create Table of (cnt all ang,count acute ang)
    counttable =[]
    for pid in range(howmany):
        angles = filteredAngleLists[pid]
        counttable += [(len(angles), countAcute(angles))]
    print counttable

    #Count total angle and total acute angles
    totalangles = 0
    totalact = 0
    for cntangles, cntact in counttable:
        totalangles += cntangles
        totalact += cntact
    if debug:print totalangles, totalact

def main1():
    debug = False
    dumppickle = True
    read = True
    calcang = True #False
    howmany = 84
    anglesLists = []
    if calcang:
        for pid in range(0, howmany):
            poly_pts_sequence = getPolygonPoints(datatype='ms', pid=pid, debug=False)
            poly_pts_sequenceLSobjs = getLinesForPolygonPointsSequences2(poly_pts_sequence)
            vertices,angles= GetAnglesOfPolygon(poly_pts_sequenceLSobjs)
            angles = [ round(ang, 2) if ang is not 'inff' else ang for ang in angles]
            angles = [pid]+angles
            anglesLists.append(angles)
                 
    anglefile = 'results/angles_mscounty.txt'
    if dumppickle:
        with open(anglefile, 'ab') as fp:
            pickle.dump(anglesLists, fp)
    
    if read:
            with open(anglefile, 'rb') as fp:
                ranglesLists= pickle.load(fp)
            print len(ranglesLists)
main1()


# In[ ]:



