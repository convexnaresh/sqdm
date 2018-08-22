import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
import numpy as np
np.set_printoptions(precision=12)
import numpy as np
from  osgeo import ogr, osr
# Import the necessary modules
'''Returns numpy array with values that represents a side of a polygon,and polygon id.'''
def GetSidesForPoly(polygonRef,polyId=0):
    import numpy as np
    geom = polygonRef.Boundary()  
    polypts = geom.GetPoints()
    sides = []
    for i in range(geom.GetPointCount()-1, -1,-1):
        a,b = polypts[i]
        c,d = polypts[i-1]
        sides += [(a,b,c,d,polyId)]
    table = np.asarray(sides)
    return table
'''Returns True if two sides represented by 4-tuple @sideTuple1 and 4-tuple @sideTuple2.'''
def isSharedSide(sideTuple1,sideTuple2):
    #sideTuple contains tuple with 4 values (x1,y1,x2,y2)
    x11,y11,x12,y12 = sideTuple1
    x21,y21,x22,y22 = sideTuple2
    if (x11 == x22 and y11 == y22) and (x12 == x21 and y12 == y21):
        return True
    elif ((x11 == x21 and y11 == y21) and (x12 == x22 and y12 == y22)):
        return True
    else:
        return False
    
'''It finds out common sides of polygons and creates an entry in table that are the polygon ids of 
the polygon that corresponds to the side.'''
def CreateSharedSidesTable(sidesTable, dimension=2):
    sides = sidesTable[:,:-1]
    polyIds = sidesTable[:,-1:]
    sortedColumn = 0
    #In plane, a side is shared by atmost two polygons
    if dimension == 2:
        #sort sidesTable by first-column-x1.
        sortedidx = [sides[:,sortedColumn].argsort()]
        pass
    #In n=3, a side is shared by atmost 3-polygons.
    if dimension == 3:
        pass
    ssides = sides[sortedidx]
    spolyIds = polyIds[sortedidx]
    #? if (a1,a2,a1,a3:p1) and (a1,a3,a1,a2:p2) are side shared by p1 & p2 polygons
    #? if (a1,a2,a3,a4:p1) and (a1,a2,a6,a7:p2) are point shared by p1 and p2
    a = np.hstack((ssides[0:,0:1],spolyIds[0:,0:1]))
    
def getNeighborsD(polygon_layer):
    '''Returns a list of lists containing ids of polygon that are connected to a polygon with id at that index.
For example, list at index 0 contains polygons that are neighbors to polygion-0.
@polygion_layer is a osgeo layer. For example.
polygon_layer = shp.GetLayer()'''
    print("...")
    neighbormap = [[] for i in range(polygon_layer.GetFeatureCount())]
    for ffid in range(polygon_layer.GetFeatureCount()): #this can be optimized
        feature =  polygon_layer.GetFeature(ffid)
        for fid in range(polygon_layer.GetFeatureCount()):
            #escape if already computed
            if ffid != fid and (ffid not in neighbormap[fid]):
                candfeature = polygon_layer.GetFeature(fid) #feature
                candpoly = candfeature.GetGeometryRef() #polygon
                dist = feature.GetGeometryRef().Distance(candpoly)
                if dist == 0:
                    neighbormap[ffid] += [fid]
                    neighbormap[fid] += [ffid]
                candfeature.Destroy()
        feature.Destroy()
    return neighbormap

def getNeighborsT(polygon_layer):
    '''T means it uses touches() function. It is more efficient version of getNeighborsD function. Returns a list of lists containing ids of polygon that are connected to a polygon with id at that index.
For example, list at index 0 contains polygons that are neighbors to polygion-0.
##############@polygion_layer is a osgeo layer. For example.'''
    neighbormap = [[] for i in range(polygon_layer.GetFeatureCount())]
    howmany = polygon_layer.GetFeatureCount()
    for ffid in range(howmany): #this can be optimized
        feature =  polygon_layer.GetFeature(ffid)
        for fid in range(polygon_layer.GetFeatureCount()):
            #escape if already computed
            if ffid != fid and (ffid not in neighbormap[fid]):
                candfeature = polygon_layer.GetFeature(fid) #feature
                candpoly = candfeature.GetGeometryRef() #polygon
                dist = feature.GetGeometryRef().Touches(candpoly)
                if dist == True:
                    neighbormap[ffid] += [fid]
                    neighbormap[fid] += [ffid]
                candfeature.Destroy()
        feature.Destroy()
    return neighbormap

def testIfTwoMapsEqual(neighborhoodmap1,neighborhoodmap2):
    '''Returns True,None,None if all the items in two maps are equal, other returns False,i,j;
i denotes match fails at list index-i, and it's jth item. '''
    for i in range(len(neighborhoodmap1)):
        for j in range(len(neighborhoodmap1[i])):
            if neighborhoodmap1[i][j] != neighborhoodmap2[i][j]:
                return False, i, j
    return True, None,None

def prunePolygonPoints(polygonRef,prunable_polygonpoints):
    '''Returns tuples (ai,bj)s with ai to bj are index of 
the points in polygon that are not in prunable_polygonpoints'''
    polypts = []
    poly= polygonRef.GetBoundary()
    polypts = poly.GetPoints()
    
    if type(prunable_polygonpoints) == ogr.Geometry :
        intersection_geometry = prunable_polygonpoints
        prunable_polygonpoints = getSidesFromMLS(intersection_geometry)

    outerpolygon_points = []
    contd = True
    contdidx = []
    inf = 0
    for ptcnt in range(len(polypts)):
        if polypts[ptcnt] not in prunable_polygonpoints:
            outerpolygon_points +=[polypts[ptcnt]]
            if contd:
                contdidx +=[(ptcnt,inf)]
                contd = False
        else:
            if not contd:
                contdidx[-1]= (contdidx[-1][0],ptcnt-1)
                contd = True         
    return contdidx,outerpolygon_points

def getSharedPoints(polygonRef,unionpoints):
    #duplicate points in polygonRef are repeated.
    
    '''Returns list of shared points or common points in between points in @polygonRef and @unionpoints'''
    polygonpoints = polygonRef.GetBoundary().GetPoints()
    #unionpoints = unionPolygon.GetBoundary().GetPoints()
    outerpolygon_points = []
    cont =[]
    for ppt in polygonpoints[:-1]:
        if ppt in unionpoints:
	    cont += [True]
            outerpolygon_points +=[ppt]
	    continue
   	cont += [False]
    if cont[0] is True and cont[1] is False:
        return outerpolygon_points[1:]+[outerpolygon_points[0]]
    return outerpolygon_points

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

def getInteriorExteriorPolygons(sharedpolyidx):
    '''@sharedpolyidxtuple is a list of tuple'''
    bpids = []
    interiorpids = []
    pids = []
    for tupleidx in range(len(sharedpolyidx)):
        pid1,pid2,ptcnt = sharedpolyidx[tupleidx]
        #Points shared by pid1 and '0' is zero.
        pids += [pid1,pid2]
        if pid2 == '0' :
            if ptcnt > 0 :
                bpids += [pid1]
            else:
                interiorpids += [pid1]
    return interiorpids, bpids

def intersectionComputed(intersection_computed_tbl, polyid1,polyid2):
    '''Returns true if intersection is computed between polyid1 and polyid2'''
    debug=False
    if debug: print("p1,p2"), polyid1,polyid2
    if polyid1 in intersection_computed_tbl[polyid2]:
        return True
    return False

def getPolygonsIds(polygon_layer):
    '''Returns a dictionary containing <feature-name : ids> where id is index of feature-name in .dbf file.'''
    howmany = polygon_layer.GetFeatureCount()
    pdict = {}
    i = 1
    for fid in range(howmany):
        feature = polygon_layer.GetFeature(fid)
        if feature.GetField('CONAME') == None:
            pdict['none_'+str(i)] = fid
            i +=1
        else:
            pdict[feature.GetField('CONAME').lower()] = fid
        feature.Destroy()
    
    return pdict

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


def getOuterBoundaryForMap(polygon_layer,neighbormap):
    '''Returns outer boundary of a map for a @polygon_layer. 
    It is union of all the inner polygons in  @polygon_layer.'''
    howmany = polygon_layer.GetFeatureCount()
    unionPolygon = ogr.Geometry(ogr.wkbPolygon)
    for fid in range(howmany):
        feature = polygon_layer.GetFeature(fid)
        polygon = feature.GetGeometryRef()
        unionPolygon = unionPolygon.Union(polygon)
        neighbors = neighbormap[fid]
        feature.Destroy()
    return unionPolygon

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

def getPolygonPointsById(polygonid,pointstable,polygon_to_ptsindexmapping,returnonlyuniquesegment=False):
    '''Returns list containing continious sequence of points in polygon with id @polygonid.
    Points are selected from @pointstable as indicated by indices in the dictionary @polygon_to_ptsindexmapping'''

    if returnonlyuniquesegment==False:
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

def constructSharedPointsDatabase(layer,neighbormap):
    '''Returns 1) list myPolyPoints of unique points (x1,y1) in a map denoted by a gis layer @layer parameter 
    and 2) mappings from polygon_id to list of <i,j> which indicates that points from index i to j 
inclusive belong to polygoin_id. It is created from mapptables that contains
indices entries <polyid1,polyid2,npoints> which means that polygons with id polyid1 and polyid2 share
npoints. Occurance of entries <polyid1,polyid2, npoints> should in proper order, it should not be messed up.'''
    debug = False
    allpoints = []
    outer_polygonpoints =[]
    polygon_layer = layer
    howmany = polygon_layer.GetFeatureCount()
    intersection_tbl = [[] for i in range(howmany)]
    sharedpolyidx = [] #tuple of (pid1,pid2,howmanyptsshared) tuples containing ids of polygon, and how many such points shared from
    insert_single_shared_pt = False
    fmls = []
    for fid in range(howmany):
        feature = polygon_layer.GetFeature(fid)
        polygon = feature.GetGeometryRef()
        neighbors = neighbormap[fid]
        inner_polygonpoints = []
        for each_neighbor in neighbors:
            if not intersectionComputed(intersection_tbl,fid,each_neighbor):
                neighborfeature = polygon_layer.GetFeature(each_neighbor)
                neighborpolygon = neighborfeature.GetGeometryRef()
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
                neighborfeature.Destroy()
        allpoints += inner_polygonpoints    
        feature.Destroy()

    #find out the the points not shared by polygons
    #its slower process. Needs some optimization.
    unionpoints = getOuterBoundaryForMap(polygon_layer,neighbormap).GetBoundary().GetPoints()    
    for fid in range(howmany):
        feature = polygon_layer.GetFeature(fid)
        polygon = feature.GetGeometryRef()
        outerpolygonpts = getSharedPoints(polygon, unionpoints)
        if insert_single_shared_pt:
            outer_polygonpoints += outerpolygonpts
            if len(outerpolygonpts) >= 1:
                sharedpolyidx += [(fid,'0',len(outerpolygonpts))] #'0' indicates that it is outer polygon
        elif len(outerpolygonpts) > 1:
            outer_polygonpoints += outerpolygonpts
            sharedpolyidx += [(fid,'0',len(outerpolygonpts))]
        feature.Destroy()
    allpoints +=outer_polygonpoints
    print("completed.Returns list of unique points and mapping from polygon to point indices.") 
    return allpoints, getPolygonPointsIndicesMapping(sharedpolyidx)

def calcSlope(line):
    '''Returns slope of a two points of wkbLineString object representing a twio point line or list of tuples
    representing two points passed as arguement'''
    
    if isinstance(line,ogr.Geometry):
        if line.GetGeometryName() == 'LINESTRING':
            wkbLineString = line
            p1 = wkbLineString.GetPoint(0)
            p2 = wkbLineString.GetPoint(1)
        else:
            message = "parameter should be either LineString type or list of tuples."
            raise Exception(message)
        
    #line may be list with tuple of two points [(a,b),(c,d)]
    if type(line) == list:
            p1 = line(0)
            p2 = line(1)
    try:
        slope = (p2[1]-p1[1])/(p2[0]-p1[0])
        pass
    except ZeroDivisionError:
        slope = 'inf'
    return slope

def getLinesForPolygonPointsSequences(ppointssequenceList,close=False):
    #if close = True then constructs a closed polygon.
    plinesequenceList = []
    
    for ptseq in ppointssequenceList:
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
	
        plinesequenceList.append(lineseq)
    return plinesequenceList
'''
#Example Code
import numpy as np
from  osgeo import ogr, osr
driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')
layer = shp.GetLayer()
neighbormap = getNeighborsT(layer)
#construct a table consisting of polygon points, and map consisting of polygion and indices of the points in the table.
pointstable, polygon_to_ptsindexmapping = constructSharedPointsDatabase(layer)
'''
