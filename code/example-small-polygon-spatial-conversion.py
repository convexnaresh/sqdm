
# coding: utf-8

# In[14]:

import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
from descartes import PolygonPatch
import treeformation101 as Tf
import MapCommonSides101 as Mcs
import math

def removeEmptyList(listlist,inplace=True):
    '''Removes one level empty list from list  of lists.'''
    if inplace:
        for listItem in listlist:
            if not listItem:
                listlist.remove(listItem)
        return None
    else:
        return [ listItem for listItem in listlist if listItem]

def LineToSegment(line,lineid,linetype='geos'):
    #Can be customized by removing nested if-statement if linetype argument is given.
    '''@line can be 1) geos LineString geometry object 
                    2) list of two tuples representing points in a line.
                    Example: [(x1,y1),(x2,y2)]'''
    if linetype =='geos':
        if isinstance(line,ogr.Geometry):
            if line.GetGeometryName() == 'LINESTRING':
                wkbLineString = line
                point1 = wkbLineString.GetPoint(0)
                point2 = wkbLineString.GetPoint(1)
            else:
                message = "parameter should be either LineString type or list of tuples."
                raise Exception(message)

    #line may be list with tuple of two points [(a,b),(c,d)]
    if linetype=='list':
        if type(line) == list:
                point1 = line(0)
                point2 = line(1)
                
    if (point1[0]> point2[0]):
        ktuple = (point2[0],point1[0]) #forming tuple (x1,x2)
    else:
        ktuple = (point1[0],point2[0])
    segment = Tf.Segment(ktuple,lineid)
    return segment

def getLinesForPolygonPointsSequences(ptseq,close=False):
    #if close = True then constructs a closed polygon.
    plinesequenceList = []
    if close:
        ptseq.append(ptseq[0])#to construct a closed polygon
    lineseq = []
    howmany = len(ptseq)
    for ptid in range(howmany-1):
        point1,point2 = ptseq[ptid], ptseq[ptid+1]
        line1 = ogr.Geometry(ogr.wkbLineString)
        line1.AddPoint(point1[0],point1[1])
        line1.AddPoint(point2[0],point2[1])
        lineseq +=[line1]
    return lineseq

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

def getPolygonPointsForMultiPolygonIds(polygonidlist,pointstable,polygon_to_ptsindexmapping):
    series=[]
    for polygonid in polygonidlist:
        lstpts=getPolygonPointsById(polygonid,pointstable,polygon_to_ptsindexmapping,ret_allpoly=False)
        series.append(lstpts)
    return series
    
def getPolygonPointsById(polygonid,pointstable,polygon_to_ptsindexmapping,ret_allpoly=False):
    '''Returns list containing continious sequence of points in polygon with id @polygonid.
    Points are selected from @pointstable as indicated by indices in the dictionary @polygon_to_ptsindexmapping'''

    if ret_allpoly==False:
        pidxes = polygon_to_ptsindexmapping[polygonid] #return list of <i,j> tuples that indicates indices of points for polygon.
    else:
        indextuples = []
        for key,tuples in polygon_to_ptsindexmapping.items():
            if key != '0':
                indextuples +=tuples
        pidxes = list(set(indextuples)) #select unique indextuples.

    ppointssequenceList = []
    for i,j in pidxes:
        ppointssequenceList.append(pointstable[i:j+1])
    return ppointssequenceList #it is list of list like [[],[]]because to keep track of sequential lines.

def getSharedPointsOld(polygonRef,unionpoints):
    '''Returns list of shared points or common points in between points in @polygonRef and @unionpoints'''
    polygonpoints = polygonRef.GetBoundary().GetPoints()
    outerpolygon_points = []
    cont =[]
    for ppt in polygonpoints[:]:
        if ppt in unionpoints:
            cont += [True]
            outerpolygon_points +=[ppt]
        else:
            cont += [False]
            
    if cont[0] is True and cont[1] is False:
        return outerpolygon_points[1:]+[outerpolygon_points[0]]
    return outerpolygon_points

def getSharedPointsUpdated(polygonRef,unionpoints):
    '''Returns list of shared points or common points in between points in @polygonRef and @unionpoints
    @ImprovePerformance'''
    polygonpoints = polygonRef.GetBoundary().GetPoints()
    outerpolygon_points = []
    cont =[]
    for ppt in polygonpoints[:]:
        if ppt in unionpoints:
            cont += [True]
        else:
            cont += [False]
    tt =[]
    t=[]
    cont +=[False] #so that all the trues before last false are appended.
    for i in range(len(cont)):
        if cont[i] is True:
            t+=[i]
        elif len(t) > 0 :
            tt.append(t)
            t=[]
            
    return tt

def getSharedPoints(polygonRef,unionpoints):

    '''Returns list of shared points or common points in between points in @polygonRef and @unionpoints'''
    polygonpoints = polygonRef.GetBoundary().GetPoints()
    #unionpoints = unionPolygon.GetBoundary().GetPoints()
    outerpolygon_points = []

    for ppt in polygonpoints[1:]:#exclude last point which is starting point in a closed polygon.
        if ppt in unionpoints:
            outerpolygon_points +=[ppt]
    return outerpolygon_points

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

def getSidesFromMLS(intersection_geometry, point_mode=False):
    '''Returns sequence of points for a list of MultiLineString. For example, for a MultiLineString <a,b> and <b,c>, 
it returns sequence [a,b,c]'''
    debug = True
    ipoints = []
    if intersection_geometry.GetGeometryName() == 'MULTILINESTRING':
        multilinestring = intersection_geometry
        for linestr in multilinestring:
            #return single point from each linestr untill last linestr.
            ipoints += [linestr.GetPoints()[0]]
        ipoints.append(multilinestring.GetGeometryRef(multilinestring.GetGeometryCount()-1).GetPoints()[1])

    if intersection_geometry.GetGeometryName() == 'POINT':
        #if only one point is in the intersection, then return the point if point_mode=True.
        if point_mode:
            if debug: print("Single Point marked as interior point.")
            ipoints += intersection_geometry.GetPoints()
        else:
            if debug: print("Single Point marked as exterior point.")
    return ipoints

def getPolygonPointIndicesForSharedPolyIndices(polygonid,sharedpolyidx):
    '''@sharedpolyidxtuple is a list of tuples (pid1,pid2,npoints), where pid1, and pid2 are polygon ids that share
and @npoints. 
This function returns list of tuples (i,j) where points from iTh to jTh (inclusive) indices are points in polygon @polygonid.
To get items from i to jth indices in 0-indexed list, use as below:
myList[i:j+1] #total items returned is (j-1+1) '''
    ptindices = []
    for tupleidx in range(len(sharedpolyidx)):
        pid1,pid2,ptcnt = sharedpolyidx[tupleidx]
        if pid1 == polygonid or polygonid== pid2:
            previouswalk = 0
            for a,b,c in  sharedpolyidx[0:tupleidx]:
                previouswalk +=c
            ptindices +=[(previouswalk,previouswalk+ptcnt-1)]
    return ptindices

def getPolygonPointsIndicesMapping(sharedpolyidx):
    polygonids = []
    for p1,p2,npoints in sharedpolyidx:
        polygonids +=[p1,p2]
    uniquepolyids = list(sorted(set(polygonids)))
    polyidpointindexMap = {}
    for pid in uniquepolyids:
        pidxes = getPolygonPointIndicesForSharedPolyIndices(pid,sharedpolyidx)
        polyidpointindexMap[pid] = pidxes
    return polyidpointindexMap

def getPolygonPointsIndicesMapping(sharedpolyidx):
    polygonids = []
    for p1,p2,npoints in sharedpolyidx:
        polygonids +=[p1,p2]
    uniquepolyids = list(sorted(set(polygonids)))
    polyidpointindexMap = {}
    for pid in uniquepolyids:
        pidxes = getPolygonPointIndicesForSharedPolyIndices(pid,sharedpolyidx)
        polyidpointindexMap[pid] = pidxes
    return polyidpointindexMap

def getOuterBoundaryForMap(polygons,neighbormap):
    '''Returns outer boundary of a map for a @polygon_layer. 
    It is union of all the inner polygons in  @polygon_layer.'''
    howmany = len(polygons)
    unionPolygon = ogr.Geometry(ogr.wkbPolygon)
    for fid in range(howmany):
        polygon = polygons[fid]
        unionPolygon = unionPolygon.Union(polygon)
        neighbors = neighbormap[fid]  
        
    if unionPolygon.GetGeometryName() == 'MULTIPOLYGON':
        unionPolygonPoints = []
        for polygon in unionPolygon:
            unionPolygonPoints += polygon.GetBoundary().GetPoints()
        return list(set(unionPolygonPoints))
    if unionPolygon.GetGeometryName() == 'POLYGON':
        return unionPolygon.GetBoundary().GetPoints()

def getSidesFromMLS(intersection_geometry, point_mode=False):
    '''Returns sequence of points for a list of MultiLineString. For example, for a MultiLineString <a,b> and <b,c>, 
it returns sequence [a,b,c]'''
    debug = True
    ipoints = []
    if intersection_geometry.GetGeometryName() == 'MULTILINESTRING':
        multilinestring = intersection_geometry
        for linestr in multilinestring:
            #return single point from each linestr untill last linestr.
            ipoints += [linestr.GetPoints()[0]]
        ipoints.append(multilinestring.GetGeometryRef(multilinestring.GetGeometryCount()-1).GetPoints()[1])

    if intersection_geometry.GetGeometryName() == 'POINT':
        #if only one point is in the intersection, then return the point if point_mode=True.
        if point_mode:
            if debug: print("Single Point marked as interior point.")
            ipoints += intersection_geometry.GetPoints()
        else:
            if debug: print("Single Point marked as exterior point.")
    return ipoints


def intersectionComputed(intersection_computed_tbl, polyid1,polyid2):
    '''Returns true if intersection is computed between polyid1 and polyid2'''
    debug=False
    if debug: print("p1,p2"), polyid1,polyid2
    if polyid1 in intersection_computed_tbl[polyid2]:
        return True
    return False

def constructSharedPointsDatabase(polygons,neighbormap):
    '''Returns 1) list myPolyPoints of unique points (x1,y1) in a map denoted by a gis layer @layer parameter 
    and 2) mappings from polygon_id to list of <i,j> which indicates that points from index i to j 
inclusive belong to polygoin_id. It is created from mapptables that contains
indices entries <polyid1,polyid2,npoints> which means that polygons with id polyid1 and polyid2 share
npoints. Occurance of entries <polyid1,polyid2, npoints> should in proper order, it should not be messed up.'''
    debug = False
    allpoints = []
    outer_polygonpoints =[]
    howmany = len(polygons)
    intersection_tbl = [[] for i in range(howmany)]
    sharedpolyidx = [] #tuple of (pid1,pid2,howmanyptsshared) tuples containing ids of polygon, and how many such points shared from
    insert_single_shared_pt = True
    fmls = []
    for fid in range(howmany):
        polygon = polygons[fid]
        neighbors = neighbormap[fid]
        inner_polygonpoints = []
        for each_neighbor in neighbors:
            if not intersectionComputed(intersection_tbl,fid,each_neighbor):
                neighborpolygon = polygons[each_neighbor]
                intersection_mls = polygon.Intersection(neighborpolygon)
                intersection_points = getSidesFromMLS(intersection_mls,point_mode=True)
                if insert_single_shared_pt:
                    inner_polygonpoints += intersection_points
                    sharedpolyidx += [(fid, each_neighbor,len(intersection_points))] #insert fid,each_neighbor in points to mark shared points          
                elif len(intersection_points) > 1 :
                    inner_polygonpoints += intersection_points
                    sharedpolyidx += [(fid, each_neighbor,len(intersection_points))] #insert fid,each_neighbor in points to mark shared points 

                intersection_tbl[fid] +=[each_neighbor]
                intersection_tbl[each_neighbor] +=[fid]
        allpoints += inner_polygonpoints
    #find out the the points not shared by polygons
    #its slower process. Needs some optimization.
    unionpoints = getOuterBoundaryForMap(polygons,neighbormap)
    for fid in range(howmany):
        polygon = polygons[fid]
        outerpolygonpts = getSharedPointsUpdated(polygon, unionpoints)
        if outerpolygonpts:
            polygonpts = polygon.GetBoundary().GetPoints()
            for indexl in outerpolygonpts:
                outer_polygonpoints += [polygonpts[idx] for idx in indexl]
                sharedpolyidx += [(fid,'0',len(indexl))]

    allpoints +=outer_polygonpoints
    print("completed.Returns list of unique points and mapping from polygon to point indices.")
    return allpoints, getPolygonPointsIndicesMapping(sharedpolyidx)

def getRectanglePoints(xleft,xright,ybottom,ytop):
    '''Returns a list with fourt 2-tuples each represents vertices of a rectangle.
    For example:Returns [(2.0, 4.0), (3.0, 4.0), (3.0, 5.0), (2.0, 5.0)] '''
    return [(xleft,ybottom),(xright,ybottom),(xright,ytop),(xleft,ytop)]
    
def constructUnOrderedRecsPolygon(xordranges,imediateyords): 
    '''
    aks construct-unordered-rectabgular-polygon.unordered means discrete; ordered means continous.
    INPUT @lutentries (xordranges) are all unique ranges along x-cords. Example [(x1,x2),(x2,x3)...]
    lowerlevels (imediateyords) aks immediate-yordinates, are all y-ords that fall inside (x1,x2) lutentries. 
    INPUT @immediateyords = [[y1,y2,y3,y4], [y2,y3,y6,y7]..] corresponding to each range (item) in xordranges.
    OUTPUT @vcolumns = [vcolumn1,vcolumn2...] for each ith item in xordranges
    @vcolumn =[rec1,rec2,rec3...].
    @rec1 = rectangle bonded by x1,x2 along x-axis and a pair of yord in 
    '''
    debug=False
    vcolumns= []
    round_rect= False
    for i in range(len(xordranges)-1): #excluding rounding-up entry
        x1,x2 = xordranges[i][:-1] #[1.0, 2.0, 0] exclude last tuple item which is 0.
        vcolumn = []
        for j in xrange(0,len(imediateyords[i]),2): #each pairs forms a rectangle
            y1,y2 = imediateyords[i][j:j+2]
            rectangle = [(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
            vcolumn.append(rectangle)
            if round_rect:
                rectangle = [(int(round(a)),int(round(b))) for a,b in rectangle]
        vcolumns.append(vcolumn)
    vcolumns.append(['None'])
    if debug:print("len xordranges,imediateyords,vcolumns"), len(xordranges),len(imediateyords),len(vcolumns) 
    if debug:print(vcolumns)
    return vcolumns

def constructOrderedRoundedRecPolygons(xordranges,imediateyords_rpaired):
    rectangles =[]
    debug = False
    for i in range(len(xordranges)-1):
        lutent = xordranges[i]
        recs = []            
        if debug:
            print("xords-range:"),xordranges[i]
            print("\t yords-range:"), imediateyords_rpaired[i]
        
        for each_yordpair in imediateyords_rpaired[i]:
            if debug:print("\t"),each_yordpair
            recs.append(getRectanglePoints(lutent[0],lutent[1],each_yordpair[0],each_yordpair[1]))
        
        if debug:
            print('\trectangles_in_columns')
            print('\t\t'),recs
        rectangles.append(recs)
    #rectangles.append([(0,0),(0,0),(0,0),(0,0)])
    rectangles.append(['Null'])
    if debug:
        print("xords-range:"),xordranges[i+1]
        print("\tyords-range:"), imediateyords_rpaired[i+1]
        print('\t\t'),rectangles[i+1]
    assert len(imediateyords_rpaired) == (len(xordranges)), "len(vcolumns) != len(lutentries){0},{1}".format(len(lutentries)-1, len(vcolumns))
    return rectangles

def constructOrderedRoundedRecsYords(lutentries,imediateyords): 
    import copy
    '''Returns rounded yords pairs for each lutentries.
        rectangle = [(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
        vcolumn.append(rectangle)
    '''
    vcolumns= [] #set of vcolumn, vcolumn-> list of intermediate 'im' yords.
    for i in range(len(lutentries)-1):#excluding rounding-up entry
        vcolumn = []
        for j in xrange(0,len(imediateyords[i]),2): #each pairs forms a rectangle
            y1,y2 = imediateyords[i][j:j+2]
            rectangle = [tuple(sorted((y1,y2)))] #we take only lower y and upper y value
            vcolumn +=rectangle
        
        #sort items in v-column by first elements vcolumn = [(10.0, 9.0), (3.75, 4.0), (1.0, 2.0)]
        svcolumn = [i[1] for i in sorted(enumerate(vcolumn), key=lambda pairyords:pairyords[1][0],reverse=False)]
        svcolumn +=[svcolumn[0]] #to round up
        extendedsvcolumn = []
        #find other intermediate rectangles's yords
        for j in range(len(svcolumn)-1):
            newrec = (svcolumn[j][1], svcolumn[j+1][0])
            extendedsvcolumn +=[svcolumn[j],newrec]                      
        vcolumns.append(extendedsvcolumn)
    vcolumns.append(['Null'])
    assert len(vcolumns) == (len(lutentries)), "len(vcolumns) not equal to len(lutentries){0},{1}".format(len(lutentries)-1, len(vcolumns))
    return vcolumns

def getIntMediateYords(listLineStringObjs,uniquexords,len_xrange=0):
    debug = False
    polygonid = 0 ###
    imedyords = [[] for i in range(len_xrange-1)] +[['None']]
    for linesegmentObj in listLineStringObjs:
        segment = LineToSegment(linesegmentObj,polygonid,linetype='geos')
        intermediateXords = Tf.getIntermediateXords(uniquexords,segment)#include enpoint's xords
        intermediateYords = Tf.getIntermediateYords(linesegmentObj,intermediateXords)#includes endpoints' yords.
        assert len(intermediateXords) == len(intermediateYords),'im-Xords & Yords no match %r,%r'%(len(intermediateXords),len(intermediateYords))
        count = 0
        #escape vertical line.
        if Mcs.calcSlope(linesegmentObj) == 'inf':
            if debug:print(linesegmentObj),' slope inf ','escaped'
            continue
        for stxord in intermediateXords[0:-1]:
            occursindex = uniquexords.index(stxord)
            imedyords[occursindex] += intermediateYords[count:count+2]
            count +=1 
    return imedyords


# In[7]:

def main2():
    import itertools, sys,time,pickle
    
    close=True
    round_rect = False
    ret_allpoly = True
    polyid= 0
    debug = False
    pointstable,indxmap, all_polygons = getPolygonPointsMap(datatype='ms')
    if ret_allpoly:polygons= all_polygons
    else:polygons = [all_polygons[polyid]]
     
    ppointssequenceList = getPolygonPointsById(polyid,pointstable,indxmap,ret_allpoly=ret_allpoly)#True means return uniq sides.
    #Convert each sides into Geos LineString objects
    listLineStringObjs= []#list of Geos LineString
    for ppsequence in ppointssequenceList:
        lineStringObjl = getLinesForPolygonPointsSequences(ppsequence,close=False)
        listLineStringObjs +=lineStringObjl
        
    #Get all Segment objects for all Geos LineString
    segments = [] #Segment objs
    for geoslineobj in listLineStringObjs:
        segments +=[LineToSegment(geoslineobj,polyid,linetype='geos')]
        
    uniquexords = Tf.getUniqueXords(segments)#sorted unique xords
    xordranges = Tf.getLutEntries(uniquexords)#top-level lut entries

    uniqueLatvals = list(map(getXordLatVal, uniquexords))
    latvalranges = Tf.getLutEntries(uniqueLatvals)#top-level lut entries
    print uniqueLatvals[0:3]
    print
    print
    listIntMedYords = getIntMediateYords(listLineStringObjs,uniquexords,len_xrange=len(xordranges))
    #vcolumnsUnOrd= constructUnOrderedRecsPolygon(xordranges,listIntMedYords)  
    print
    print
    #Test 
    vcolumns_orderedRry = constructOrderedRoundedRecsYords(xordranges,listIntMedYords)
    '''vcolumn_orderedRry: 
        [[(0.0, 1.0), (1.0, 3.0), (3.0, 5.0), (5.0, 0.0)], 
        [(0.0, 0.5), (0.5, 4.0), (4.0, 5.0), (5.0, 5.0), (5.0, 6.25), (6.25, 0.0)], 
        [(0.5, 1.0), (1.0, 4.0), (4.0, 4.0), (4.0, 6.25), (6.25, 7.5), (7.5, 0.5)], 
        ['Null']]
    '''

    vcolumns_orderedRrecs = constructOrderedRoundedRecPolygons(xordranges,vcolumns_orderedRry)
    '''rectangles between x1,x2:
        [[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)], 
        [(1.0, 1.0), (2.0, 1.0), (2.0, 3.0), (1.0, 3.0)], 
        [(1.0, 3.0), (2.0, 3.0), (2.0, 5.0), (1.0, 5.0)], 
        [(1.0, 5.0), (2.0, 5.0), (2.0, 0.0), (1.0, 0.0)]]
    '''
    print vcolumns_orderedRrecs[0][0]
    print("mississippi point-lat/long")
    latlonlstpts = []
    for pt in vcolumns_orderedRrecs[0][0] :
        print pt
        latlonlstpts.append(ConvertCoordsPointForm(pt))
    print latlonlstpts
                            
                            
    return xordranges,uniquexords


# In[8]:

#2-D Gis points --> Lat,Long --> 32 Bit Integers
import numpy as numpy
def getXordLatVal(xord):
    point_tuple = (xord,0.0,0.0)
    return ConvertCoordsPointForm(point_tuple)[0]

def ConvertCoordsPointForm(point_tuple):
    '''Returns list of lat,long,alt values for 2-D planer cordinates passed os list or tuple of x,y,z.'''
    from osgeo import osr
    # Set up spatial reference systems
    # using EPSG 3814
    if len(point_tuple) == 3:
        x,y,z = point_tuple
    if len(point_tuple) == 2:
        x,y = point_tuple
        z = 0
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(3814)
    proj.SetTOWGS84(565.237, 50.0087, 465.658, -0.406857, 0.350733, -1.87035, 4.0812 )
    # lat/lon WGS84
    latlong = osr.SpatialReference()
    latlong.ImportFromProj4('+proj=latlong +datum=WGS84')
    #Define transform, from EPSG28992 to LatLon/WGS84
    transform = osr.CoordinateTransformation( proj, latlong  )
    (lon, lat,z) = transform.TransformPoint( x, y, z)
    return (lat,lon, z)


# In[9]:

xordranges,uniquexords = main2()
print uniquexords[0:10]
print

uniqueLatvals = list(map(getXordLatVal, uniquexords))
latvalranges = Tf.getLutEntries(uniqueLatvals)#top-level lut entries
print uniqueLatvals[0:10]
print
print latvalranges[0:10]


# In[10]:

pts = [(629496.8275247859, 0,0), (629801.2480145998, 0,0)]
xpts = [(629496.8275247859, 0,0), (629801.2480145998, 0,0)]
print list(map(ConvertCoordsPointForm, pts))
print
print list(map(ConvertCoordsPointForm, xpts))

