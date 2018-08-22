
# coding: utf-8

# In[ ]:

import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
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
                sharedpolyidx += [(fid,'K',len(indexl))]

    allpoints +=outer_polygonpoints
    return allpoints, getPolygonPointsIndicesMapping(sharedpolyidx)

def GetInsertableUniqueLineCustomSegments(pointstable,neighbormap,dic_ptseries):
    #Construct points table and maps 
    print len(pointstable), len(neighbormap)
    
    #Get all the uniquely shared pts-sequences
    lscustomObjects = []
    for ptsidx_range, lst_pids in dic_ptseries.items():
        
        poly_pts_sequence = getPtsSequences(pointstable, ptsidx_range)
        lst_LSobjs = getLinesForPolygonPointsSequences2(poly_pts_sequence)
        
        if lst_LSobjs:
            lscustomObjects += LineStringCustom().customizeGeosLStrings(lst_LSobjs,lst_pids)
    return lscustomObjects

def getPtsSequences(polygon_table, pt_index_range):
    ''''''
    i, j = pt_index_range
    ptssequence = []
    return polygon_table[i:j+1]

def main1():
    import itertools, sys,time,pickle
    
    ppointssequenceList = getPolygonPointsById(polyid,pointstable,indxmap,ret_allpoly=ret_allpoly)#True means return uniq sides.
    
    listLineStringObjs= []#list of Geos LineString
    for ppsequence in ppointssequenceList:
        if debug: print("side:"),ppsequence
        lineStringObjl = getLinesForPolygonPointsSequences(ppsequence,close=False)
        if debug: print("\t"),[geosLS.GetGeometryName()  for geosLS in lineStringObjl]
        listLineStringObjs +=lineStringObjl

def getSharedPolygonTable(all_polygons):
    neighbors = getNeighborsT(all_polygons)
    neighbormap = {}
    howmany = len(all_polygons)
    for pid in range(0,howmany):
        neighbormap[pid] = neighbors[pid]
    pointstable,indxmap = constructSharedPointsDatabase(all_polygons,neighbormap)

    dic_ptseries = {}
    for pid,list_pt_index_range in indxmap.items():
        for pt_index_range in list_pt_index_range:
            try:
                dic_ptseries[pt_index_range] +=[pid]
            except:
                dic_ptseries[pt_index_range] = []
                dic_ptseries[pt_index_range] +=[pid]
    print("Completed building shared points table for polygons")
    return pointstable,neighbormap,dic_ptseries

class LineStringCustom(object):
    
    def __init__(self,geos_line_string=None,attr=[]):
        self.geosls =geos_line_string
        self.attr =",".join(map(str,attr))
        
    def customizeGeosLStrings(self,list_geos_ls,commonattr):
        LSCustomObjects = []
        
        for geos_ls in list_geos_ls:
            LSCustomObjects +=[LineStringCustom(geos_ls,commonattr)]
        return LSCustomObjects
    
    
    def __str__(self):
        p1,p2 = self.geosls.GetPoints()
        return str((p1,p2,self.attr))



# In[ ]:


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

if __name__ == "__main__":
    neighbors = getNeighborsT(all_polygons)
    neighbormap = {}
    howmany = len(all_polygons)
    for pid in range(0,howmany):
        neighbormap[pid] = neighbors[pid]

    print neighbormap
    print
    pointstable,indxmap = constructSharedPointsDatabase(all_polygons,neighbormap)
    print pointstable
    print("Test-3")
    dic_ptseries = {}
    for pid,list_pt_index_range in indxmap.items():
        for pt_index_range in list_pt_index_range:
            try:
                dic_ptseries[pt_index_range] +=[pid]
            except:
                dic_ptseries[pt_index_range] = []
                dic_ptseries[pt_index_range] +=[pid]
    print
    print("Test-4")
    print dic_ptseries


# In[ ]:



