
# coding: utf-8

# In[12]:

# In[84]:

'''
Algorithm (linestring ls1, polygons lst):
out: all the names of polygon that share this linestring.

For a linestring LS1, find out all the polygons that shares this linestring.
    for each in polygons:
        if LS1 is in envelope of polygon-1:
            if points in LS1 are on boundary of polygon-1:
                return TRUE
        else:
            return False
        
'''
'''
For a polygon to be inserted in existing DLS
Algorithm (dls, polygon, processedtable):
    processedtable: list of linestring inserted in dls
    dls : dictionary {xord: [lineid: [yords,slope],...}
    polygon: polygon that is to be inserted in DLS
find out all linestrings from poly-1 which are not in processed tables
inserted them

call Algorithm for all polygons.

'''
'''
Algorithm (linestring ls1, polygons lst):
out: all the names of polygon that share this linestring.

For a linestring LS1, find out all the polygons that shares this linestring.
    for each in polygons:
        if LS1 is in envelope of polygon-1:
            if points in LS1 are on boundary of polygon-1:
                return TRUE
        else:
            return False
        
'''
'''
For a polygon to be inserted in existing DLS
Algorithm (dls, polygon, processedtable):
    processedtable: list of linestring inserted in dls
    dls : dictionary {xord: [lineid: [yords,slope],...}
    polygon: polygon that is to be inserted in DLS
find out all linestrings from poly-1 which are not in processed tables
inserted them

call Algorithm for all polygons.
'''


# In[2]:




# In[85]:

import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
import treeformation101 as Tf
import MapCommonSides101 as Mcs
import math
import time
import pickle
import polygons_table_builder as PTB
reload(PTB)
import LineStringCustom
from  DlsStat import DlsStat

def getPolygonPoints(datatype='sample', pid=0, debug = False):    
    saveFiles = True
    all_polygons = [] #list of all the polygons
    
    #======================================================
    if datatype =="sample":
        #Manual Polygons
        wkt0 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        wkt1 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 4 1,2 0, 1 1,1 3,2 5 ))"
        wkt2 = "POLYGON ((5 3, 7 1, 10 1, 11 2,11 4,7 3,5 3 ))"
        wkt3 = "POLYGON ((2 5, 4 6, 3 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        
        poly0 = ogr.CreateGeometryFromWkt(wkt0)
        poly1 = ogr.CreateGeometryFromWkt(wkt1)
        poly2 = ogr.CreateGeometryFromWkt(wkt2)
        poly3 = ogr.CreateGeometryFromWkt(wkt3)
        all_polygons = [poly0,poly1,poly2, poly3] #list of all the polygonsA
        polygon = all_polygons[pid]
        return all_polygons[pid] ##polygon.GetBoundary().GetPoints()
        
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

def compareGeosLineStrings(geosLS1, geoLS2):
    '''Return true if they represent same line segments.'''
    l1p1,l1p2 = geosLS1.GetPoints()
    l2p1, l2p2 = geosLS2.GetPoints()
    if comareTwoPoints(l1p1,l2p1) and comareTwoPoints(l1p1,l2p1):
        return True
    elif comareTwoPoints(l1p1,l2p2) and comareTwoPoints(l1p2,l2p1):
        return True
    return False
    
def comareTwoPoints(point1,point2):
    x1,y1 = point1
    x2, y2 = point2
    if x1 == x2 and y1 == y2:
        return True
    return False
    
def getIntermediateXordsModified(unique_xordslist,xord1,xord2):
    '''returns list of intermediate xords that falls between given xord1 and xord2.'''
    #print("inside getIntermediateXords"),segment.getSegKP()[0],segment.getSegKP()[1]
    xordmin, xordmax = sorted([xord1,xord2])
    intermediateXords = []
    for xord in unique_xordslist[0:]:
        if xord >= xordmin and xord <= xordmax:
            intermediateXords +=[xord]
    if len(intermediateXords)==1:
        return intermediateXords*2 #for vertical line with same xords in segment
    else:
        return intermediateXords


def getIntermediateYordsModified(linesegmentObj,intermediateXords):
    '''IntermediateXords contains end xords as well.Escape them before finding intermediate yords. Return list are endpoint'''
    
    if type(linesegmentObj) is LineStringCustom:
        p1 = linesegmentObj.geosls.GetPoint(0)
        p2 = linesegmentObj.geosls.GetPoint(1)
    else:
        p1 = linesegmentObj.GetPoint(0)
        p2 = linesegmentObj.GetPoint(1)
        
        
    intermediateYords = []
    for intxord in intermediateXords[1:-1]: #excluding firs and last xords.
        intermediateYords +=[((p2[1] - p1[1])*(intxord-p1[0])/float(p2[0]-p1[0])) + p1[1]]
        
    if p1[0] > p2[0]:
        return [p2[1]] + intermediateYords + [p1[1]]
    else:
        return [p1[1]] + intermediateYords + [p2[1]]
    
def updateXords(newxords):
    return newxords

def updateProcessedLS(processed_ls, newlineid):
    return processed_ls + [newlineid]

def initDLS():
    debug=True
    rnewxords = []
    processed_lsids = []
    imedyords = {}
    if debug:print("Initialized DLS for 2-D map.")
    return processed_lsids, rnewxords, imedyords


def constructDLSForPolygon(polygon_points,processed_lsids=[], oldxords=[],imedyords={}):
    poly_geos_linestrings = getLinesForPolygonPointsSequences2(polygon_points)
    return addPolygonToDLS(poly_geos_linestrings,processed_poly_goesls=[],processed_lsids=[], oldxords=[],imedyords={})

def addPolygonToDLS(new_poly_geosls,processed_poly_goesls=[],processed_lsids=[], oldxords=[],imedyords={}):
    debug = True
    #updates the existing DLS structures which are processed_lsids, oldxords, imedyords
    poly_geos_linestrings = processed_poly_goesls + new_poly_geosls
    if processed_lsids:
        min_lineid = len(processed_lsids)
    else:
        processed_lsids, oldxords,imedyords = initDLS()        
        min_lineid = 0
    
    for newlineid in range(min_lineid, min_lineid+len(new_poly_geosls),1):
        new_linestring = poly_geos_linestrings[newlineid]
        processed_lsids,oldxords,imedyords =insertNewLineStringInDLS(newlineid,new_linestring,processed_lsids,oldxords,imedyords)
               
        for proc_lsid in processed_lsids[0:-1]: #exclude recently added linestring
            if needUpdate(poly_geos_linestrings[proc_lsid], new_linestring):
                updateYordsForOldLineString(proc_lsid,poly_geos_linestrings[proc_lsid],oldxords,imedyords)
            #updateYordsForOldLineString(proc_lsid,poly_geos_linestrings[proc_lsid],oldxords,imedyords)
            
    if debug:print("Total lines added to dls:"), len(new_poly_geosls)
    return processed_lsids, oldxords,imedyords    


def insertNewLineStringInDLS(lineid,newlinestring,processed_lsids,oldxords,imedyords):
    debug = False

    point1, point2 = newlinestring.GetPoints()
    rnewxords_raw = oldxords + [point1[0],point2[0]]
    rnewxords = sorted(list(set(rnewxords_raw)))
    
    intermediateXords = getIntermediateXordsModified(rnewxords,point1[0],point2[0])
    intermediateYords =  getIntermediateYordsModified(newlinestring,intermediateXords)#includes endpoints' yords.
    count = 0
    
    slope = Mcs.calcSlope(newlinestring)
    if slope == 'inf':
        if debug:print(newlinestring),' slope inf ','escaped'
    else:
        for start_xord in intermediateXords[0:-1]:   
            if start_xord in imedyords.keys():
                imedyords[start_xord][lineid] = []
                try:
                    imedyords[start_xord][lineid] =intermediateYords[count:count+2] +[round(slope,2)] 
                except:imedyords[start_xord][lineid] = intermediateYords[count:count+2]+ [round(slope,2)] 
            else:
                imedyords[start_xord]= {}
                imedyords[start_xord][lineid] = intermediateYords[count:count+2] + [round(slope,2)]          
            count +=1
        imedyords[rnewxords[-1]] = {}
        
    processed_lsids = updateProcessedLS(processed_lsids, lineid)
    return processed_lsids,rnewxords, imedyords

def updateYordsForOldLineString(lineid,oldlinestring,rnewxords,imedyords):
    debug = False
    point1, point2 = oldlinestring.GetPoints()
    intermediateXords = getIntermediateXordsModified(rnewxords,point1[0],point2[0])
    intermediateYords =  getIntermediateYordsModified(oldlinestring,intermediateXords)#includes endpoints' yords.

    count = 0
    
    slope = Mcs.calcSlope(oldlinestring)
    if slope == 'inf':
        if debug:print(oldlinestring),' slope inf ','escaped'
    else:
        for start_xord in intermediateXords[0:-1]:   
            if start_xord in imedyords.keys():
                imedyords[start_xord][lineid] = [] #reset again.
                try:
                    imedyords[start_xord][lineid] = intermediateYords[count:count+2]+[round(slope,2)] 
                except:imedyords[start_xord][lineid] = intermediateYords[count:count+2]+[round(slope,2)] 
            else:
                imedyords[start_xord]= {}
                imedyords[start_xord][lineid] = intermediateYords[count:count+2] +[round(slope,2)]      
            count +=1    
    return imedyords

def needUpdateCustomLS(oldlinegeos, newlinegeos): #fist is perviously inserted obj, newlinegeos is recently inserted obj
    nlp1 , nlp2 = newlinegeos.geosls.GetPoints() #newline
    plp1, plp2 = oldlinegeos.geosls.GetPoints() #processed line

    minXNew,maxXNew = sorted([nlp1[0], nlp2[0]])
    minXOld, maxXOld = sorted([plp1[0], plp2[0]])
    #if (plp1[0] < nlp1[0] < plp2[0]) or (plp1[0] < nlp2[0] < plp2[0]):
    if (minXOld < minXNew < maxXOld) or (minXOld < maxXNew < maxXOld):
        return True
    else:
        return False


def needUpdate(oldlinegeos, newlinegeos): #fist is perviously inserted obj, newlinegeos is recently inserted obj
    nlp1 , nlp2 = newlinegeos.geosls.GetPoints() #newline
    plp1, plp2 = oldlinegeos.geosls.GetPoints() #processed line
    if (plp1[0] < nlp1[0] < plp2[0]) or (plp1[0] < nlp2[0] < plp2[0]):
        return True
    else:
        return False    
    
    
def constructDLSForMultiPolygon(multi_polygons):
    
    '''multi_polygon is a list of list of polygon points'''
    
    debug = False
    processed_poly_goesls = []
    new_poly_geosls =getLinesForPolygonPointsSequences2(multi_polygons[0])
    processed_lsids, oldxords,imedyords = addPolygonToDLS(new_poly_geosls,processed_poly_goesls=[],processed_lsids=[], oldxords=[],imedyords={})
    processed_poly_goesls += new_poly_geosls
    
    for new_polygon_points in multi_polygons[1:]:
        new_poly_geosls = getLinesForPolygonPointsSequences2(new_polygon_points)
        processed_lsids, oldxords,imedyords = addPolygonToDLS(new_poly_geosls,processed_poly_goesls,processed_lsids, oldxords,imedyords)
        processed_poly_goesls += new_poly_geosls
    
    return processed_lsids, oldxords,imedyords

def readDlsSavedPickle():

    dlsfile = './results/dls_ms' + '.txt'
    timefile = './results/dls_ms_timings' + '.txt'
    list_dls = []
    timings = []
    with open(dlsfile, 'rb') as fp:
        list_dls = pickle.load(fp)
    with open(timefile, 'rb') as fp:
        timings = pickle.load(fp)
    return list_dls, timings
    
def convert_mscountypolygon_to_dls():
    debug = False
    dumppickle = True
    list_dls = []
    timings = [] #list containing dls which is dictionary
    howmany = 84
    for pid in range(0,howmany):
        t0 = time.time()
        poly_pts_sequence = getPolygonPoints(datatype='ms', pid=pid, debug=False)
        poly_pts_sequenceLSobjs = getLinesForPolygonPointsSequences2(poly_pts_sequence)
        processed_lsids, oldxords,imedyords = constructDLSForPolygon(poly_pts_sequence)
        t1 = time.time()
        list_dls += [(pid,imedyords)]
        timings += [(pid,round(t1-t0,3))]
        
    if dumppickle:
        dlsfile = './results/dls_ms' + '.txt'
        with open(dlsfile, 'wb') as fp:
            pickle.dump(list_dls, fp)
        print("dls-dictionary saved in"), dlsfile
        
    timefile = './results/dls_ms_timings' + '.txt'
    if dumppickle:
        with open(timefile, 'wb') as fp:
            pickle.dump(timings, fp)
        print("timings for dls generation saved in "), timefile    
        
def test_single_polyon():
    howmany = 1
    for pid in range(0,howmany):
        poly_pts_sequence = getPolygonPoints(datatype='sample', pid=pid, debug=False)
        poly_pts_sequenceLSobjs = getLinesForPolygonPointsSequences2(poly_pts_sequence)
        processed_lsids, oldxords,imedyords = constructDLSForPolygon(poly_pts_sequence)
        print  processed_lsids, oldxords,imedyords

def test_addPolygonToDLS():
    multi_polygons = []
    for polygonId in range(0,3):
        poly_pts_sequence = getPolygonPoints(datatype='sample', pid=polygonId, debug=False)
        multi_polygons.append(poly_pts_sequence)
    firstpoly_geos =getLinesForPolygonPointsSequences2(multi_polygons[0])
    processed_lsids, oldxords,imedyords = addPolygonToDLS(firstpoly_geos,processed_poly_goesls=[],processed_lsids=[], oldxords=[],imedyords={})   
    #print results
    
def test_milti_polygon():
    multi_polygons = []
    for polygonId in range(0,3):
        poly_pts_sequence = getPolygonPoints(datatype='sample', pid=polygonId, debug=False)
        multi_polygons.append(poly_pts_sequence)
    print("len multi_polygons"),len(multi_polygons)
    processed_lsids, oldxords,imedyords = constructDLSForMultiPolygon(multi_polygons)
    #print results

#=======================================================================================================
def constructDLSForPolygonCustomLS(lscustomObjects):
    
    return addPolygonToDLSCustomLS(lscustomObjects,processed_lsids=[], oldxords=[],imedyords={})
    
        
def addPolygonToDLSCustomLS(lscustomObjects,processed_lsids=[], oldxords=[],imedyords={}):
        
    debug = True
    if processed_lsids:
        min_lineid = len(processed_lsids)
    else:
        processed_lsids, oldxords,imedyords = initDLS()        
    
    for newlineid in range(0,len(lscustomObjects),1):
	print("inserting linesegment"), newlineid, lscustomObjects[newlineid]
        custom_lsobj = lscustomObjects[newlineid]
        processed_lsids,oldxords,imedyords =insertNewLineStringInDLSCustomLS(newlineid,custom_lsobj,
                                                                             processed_lsids,oldxords,imedyords)
        for proc_lsid in processed_lsids[0:-1]: #exclude recently added linestring
            if needUpdateCustomLS(lscustomObjects[proc_lsid], custom_lsobj):
                updateYordsForOldLineStringCustomLS(proc_lsid,lscustomObjects[proc_lsid],oldxords,imedyords)
    
    if debug:print("Total lines added to dls:"), len(lscustomObjects)
    return processed_lsids, oldxords,imedyords #this is dls       
    
def insertNewLineStringInDLSCustomLS(lineid,newlinestring,processed_lsids,oldxords,imedyords):
    debug = False
    point1, point2 = newlinestring.geosls.GetPoints()
    rnewxords_raw = oldxords + [point1[0],point2[0]]
    rnewxords = sorted(list(set(rnewxords_raw)))
    
    intermediateXords = getIntermediateXordsModified(rnewxords,point1[0],point2[0])
    intermediateYords =  getIntermediateYordsModified(newlinestring.geosls,intermediateXords)#includes endpoints' yords.
    count = 0
    
    slope = str(Mcs.calcSlope(newlinestring.geosls))
    if slope == 'inf':
        if debug:print(newlinestring),' slope inf ','escaped'
    else:
        for start_xord in intermediateXords[0:-1]:   
            if start_xord in imedyords.keys():
                imedyords[start_xord][lineid] = []
                try:
                    imedyords[start_xord][lineid] =intermediateYords[count:count+2] 
                    imedyords[start_xord][lineid].append(slope+","+newlinestring.attr)
                    
                except:
                    imedyords[start_xord][lineid] = intermediateYords[count:count+2]
                    imedyords[start_xord][lineid].append(slope+","+newlinestring.attr)
                    
            else:
                imedyords[start_xord]= {}
                imedyords[start_xord][lineid] = intermediateYords[count:count+2]
                imedyords[start_xord][lineid].append(slope+","+newlinestring.attr)
            count +=1
        imedyords[rnewxords[-1]] = {}
        
    processed_lsids = updateProcessedLS(processed_lsids, lineid)
    return processed_lsids,rnewxords, imedyords
    
def updateYordsForOldLineStringCustomLS(lineid,oldlinestring,rnewxords,imedyords):
    debug = False
    point1, point2 = oldlinestring.geosls.GetPoints()
    intermediateXords = getIntermediateXordsModified(rnewxords,point1[0],point2[0])
    intermediateYords =  getIntermediateYordsModified(oldlinestring.geosls,intermediateXords)#includes endpoints' yords.

    count = 0
    
    slope = str(Mcs.calcSlope(oldlinestring.geosls))
    if slope == 'inf':
        if debug:print(oldlinestring),' slope inf ','escaped'
    else:
        for start_xord in intermediateXords[0:-1]:   
            if start_xord in imedyords.keys():
                imedyords[start_xord][lineid] = [] #reset again.
                try:
                    imedyords[start_xord][lineid] =intermediateYords[count:count+2] 
                    imedyords[start_xord][lineid].append(slope+","+oldlinestring.attr)
                    
                except:
                    imedyords[start_xord][lineid] = intermediateYords[count:count+2]
                    imedyords[start_xord][lineid].append(slope+","+oldlinestring.attr)
                    
            else:
                imedyords[start_xord]= {}
                imedyords[start_xord][lineid] = intermediateYords[count:count+2]
                imedyords[start_xord][lineid].append(slope+","+oldlinestring.attr)
            count +=1  
    return imedyords

def savePickle(object_tosaved,outfile):
    with open(outfile, 'wb') as fp:
        pickle.dump(object_tosaved, fp)
    print("Object saved in"), outfile
    return True
#===========================================================================================

def getDlsStats(dls):
    dlsStat = DlsStat(dls)
    return dlsStat.stats

def CreateSinglePolygonToDLS():

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(r'../data/gis/ms_county/stco.shp')
    layer = shp.GetLayer()
    
    datatype, howmany = 'sample', 3
    dlsfile = './results/dls_ms' + '.txt'
    dlslist = []
    polygons = []
    for pid in range(0,howmany):
        if datatype == 'ms':
            feature= layer.GetFeature(pid)
            polyref = feature.GetGeometryRef()
        if datatype == 'sample': polyref = getPolygonPoints(datatype='sample', pid=pid, debug=False)
        polygons = [polyref]

        #Construct points table and maps 
        pointstable,neighbormap,dic_ptseries = PTB.getSharedPolygonTable(polygons)
        #Get all the uniquely shared pts-sequences
        lscustomObjects = PTB.GetInsertableUniqueLineCustomSegments(pointstable,neighbormap,dic_ptseries)

        #Insert Each of the CustomLineStringObjects into DLS
        dls={}
        processed_lsids, oldxords,dls = constructDLSForPolygonCustomLS(lscustomObjects)
        dlslist +=[(pid,dls)]
    
    savePickle(dlslist,dlsfile)
    return dlslist

def CreateWholeMapToDLS():

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(r'../data/gis/ms_county/stco.shp')
    layer = shp.GetLayer()
    
    datatype, howmany = 'ms', 84
    dlsfile = './results/dls_ms_map' + '.txt'
    dlslist = []
  
    polygons = []
    for pid in range(0,howmany):
        if datatype == 'ms':
                polyref = None
                inFeature = layer.GetFeature(pid)
                geom = inFeature.GetGeometryRef()
                wkt0 = geom.ExportToWkt()
                polyref = ogr.CreateGeometryFromWkt(wkt0)
                polygons += [polyref]
                
        if datatype == 'sample': 
            polyref = getPolygonPoints(datatype='sample', pid=pid, debug=False)
            polygons += [polyref]
    #Construct points table and maps
    t0 = time.time()
    pointstable,neighbormap,dic_ptseries = PTB.getSharedPolygonTable(polygons)
    #Get all the uniquely shared pts-sequences
    lscustomObjects = PTB.GetInsertableUniqueLineCustomSegments(pointstable,neighbormap,dic_ptseries)
    print("Total line-segments to be inserted"), len(lscustomObjects)
    t1 = time.time()
    #Insert Each of the CustomLineStringObjects into DLS
    dls={}
    processed_lsids, oldxords,dls = constructDLSForPolygonCustomLS(lscustomObjects)
    dlslist +=[('map', dls)]
    t2 = time.time()
    savePickle(dlslist,dlsfile)
    print("times, sharedPolygonTable, dls-construction:"),round(t1-t0,2),round(t1-t2,2)
    return dlslist

if __name__ == '__main__':
        
    #dlslist = CreateSinglePolygonToDLS()
    dlslist = CreateWholeMapToDLS()


# In[ ]:



