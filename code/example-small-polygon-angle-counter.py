
# coding: utf-8

# In[11]:

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


# In[5]:

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


# In[10]:

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
    for sideString in listLineStringObjs:
        linePairs += Tf.getConsequitivePairs(sideString,winding=True)
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
main(True)


# In[ ]:



