
# coding: utf-8

# In[6]:

import numpy        
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
from descartes import PolygonPatch
#%matplotlib inline  


# In[7]:

def ConvertCoords ( x, y, z):
        from osgeo import osr
        # Set up spatial reference systems
        # using EPSG 28992 \28992
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(26759)
        proj.SetTOWGS84(565.237, 50.0087, 465.658, -0.406857, 0.350733, -1.87035, 4.0812 )
        # lat/lon WGS84
        latlong = osr.SpatialReference()
        latlong.ImportFromProj4('+proj=latlong +datum=WGS84')
        #Define transform, from EPSG28992 to LatLon/WGS84
        transform = osr.CoordinateTransformation( proj, latlong  )
        (lon, lat,z) = transform.TransformPoint( x, y, z)
        return (lat,lon, z)
 
if __name__=="__main__":
        from osgeo import osr 
        from osgeo import ogr
        polygon =numpy.array ([[132817.006604435708141, 550302.852720651309937, 0.],              [131182.28895997320069, 558340.214472591876984, 0.],              [132578.61028128489852, 558748.893883707583882, 0.],              [136631.347774848196423, 553436.061539204441942, 0.],              [136631.347774848196423, 553436.061539204441942, 0.],              [132817.006604435708141, 550302.852720651309937, 0.],             [ 266385.9, 1266212.7 ,0.0]])
        # We can use two methods: convert individual points, or create an OGR feature and convert that 
        # Method 1: Use convert individual points
        result = [ConvertCoords ( polygon[i, 0], polygon[i, 1], polygon[i,2] ) for i in xrange(polygon.shape[0])]
        print result
        print
        print
        
        # Method 2: Convert whole polygon...
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(26759)
        proj.SetTOWGS84(565.237, 50.0087, 465.658, -0.406857, 0.350733, -1.87035, 4.0812 )
        latlong = osr.SpatialReference()
        latlong.ImportFromProj4('+proj=latlong +datum=WGS84')
        wkt = "POLYGON(("
        edges = ["%f %f %f,"%( polygon[i,0], polygon[i,1], polygon[i,2]) for i in xrange(polygon.shape[0])]
        wkt = wkt+"".join (edges[:-1])+edges[-1].replace(",", "")+"))"
        geom = ogr.CreateGeometryFromWkt ( wkt )
        geom.AssignSpatialReference ( proj )
        geom.TransformTo ( latlong)
        print geom.ExportToWkt()


# In[8]:

cd /home/naresh-ad/workspace/machinelrn/data/gis


# In[9]:

ls


# In[10]:

# Import the necessary modules
from  osgeo import ogr, osr
driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')

# Get Projection from layer
layer = shp.GetLayer()
spatialRef = layer.GetSpatialRef()
print spatialRef

# Get Shapefile Fields and Types
layerDefinition = layer.GetLayerDefn()
print layerDefinition.GetFieldCount()
print "Name  -  Type  Width  Precision"
print "-----------------------------------"
for i in range(layerDefinition.GetFieldCount()-100):
    fieldName =  layerDefinition.GetFieldDefn(i).GetName()
    fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
    fieldType = layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
    fieldWidth = layerDefinition.GetFieldDefn(i).GetWidth()
    GetPrecision = layerDefinition.GetFieldDefn(i).GetPrecision()
    print fieldName + " - " + fieldType+ " " + str(fieldWidth) + " " + str(GetPrecision)

# Check if rows in attribute table meet some condition
inFeature = layer.GetNextFeature()
import sys
maxcount = sys.maxint
counts = []
while inFeature:
    # get the cover attribute for the input feature
    cover = inFeature.GetField('CONAME')
    inputGeom = inFeature.GetGeometryRef()
    if cover == 'Tishomingo':
        print "Do some action..."
    #print inputGeom
    geom = inputGeom.Boundary()    
    counts += [geom.GetPointCount()]
    #print geom.GetPoint(0)
    inFeature.Destroy()
    inFeature = layer.GetNextFeature()

print("maxcount,mincount", max(counts),min(counts))


# In[11]:

# Import the necessary modules
from  osgeo import ogr, osr
driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')

# Get Projection from layer
layer = shp.GetLayer()
spatialRef = layer.GetSpatialRef()

counts = []
features = []
# get the cover attribute for the input feature
inFeature = layer.GetNextFeature()
ptcount = 0
polypts = []
dictcb = {}
totalpts = 0
while inFeature:
    cover = inFeature.GetField('CONAME')
    inputGeom = inFeature.GetGeometryRef()
    geom = inputGeom.Boundary()   
    features += [cover]
    totalpts += geom.GetPointCount()
    if cover == 'Newton' or cover == 'Scott':
        polypts = geom.GetPoints()
        dictcb [cover] = polypts
    inFeature.Destroy()
    inFeature = layer.GetNextFeature()

print totalpts


# In[12]:

import numpy as np
np.set_printoptions(precision=12)
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

'''Returns a list of lists containing ids of polygon that are connected to a polygon with id at that index.
For example, list at index 0 contains polygons that are neighbors to polygion-0.
@polygion_layer is a osgeo layer. For example.
polygon_layer = shp.GetLayer()'''
def getNeighborsD(polygon_layer):
    neighbormap = [[] for i in range(layer.GetFeatureCount())]
    for ffid in range(layer.GetFeatureCount()): #this can be optimized
        feature =  layer.GetFeature(ffid)
        for fid in range(layer.GetFeatureCount()):
            #escape if already computed
            if ffid != fid and (ffid not in neighbormap[fid]):
                candfeature = layer.GetFeature(fid) #feature
                candpoly = candfeature.GetGeometryRef() #polygon
                dist = feature.GetGeometryRef().Distance(candpoly)
                if dist == 0:
                    neighbormap[ffid] += [fid]
                    neighbormap[fid] += [ffid]
                candfeature.Destroy()
        feature.Destroy()
    return neighbormap

def getNeighborsT(polygon_layer):
    '''It is more efficient version of getNeighborsD function. Returns a list of lists containing ids of polygon that are connected to a polygon with id at that index.
For example, list at index 0 contains polygons that are neighbors to polygion-0.
@polygion_layer is a osgeo layer. For example.
polygon_layer = shp.GetLayer()'''
    neighbormap = [[] for i in range(layer.GetFeatureCount())]
    howmany = layer.GetFeatureCount()
    for ffid in range(howmany): #this can be optimized
        feature =  layer.GetFeature(ffid)
        for fid in range(layer.GetFeatureCount()):
            #escape if already computed
            if ffid != fid and (ffid not in neighbormap[fid]):
                candfeature = layer.GetFeature(fid) #feature
                candpoly = candfeature.GetGeometryRef() #polygon
                dist = feature.GetGeometryRef().Touches(candpoly)
                if dist == True:
                    neighbormap[ffid] += [fid]
                    neighbormap[fid] += [ffid]
                candfeature.Destroy()
        feature.Destroy()
    return neighbormap

def testIfTwoMapsEqual(neighborhoodmap1,neighborhoodmap2):    '''Returns True,None,None if all the items in two maps are equal, other returns False,i,j;
i denotes match fails at list index-i, and it's jth item. '''
    for i in range(len(neighborhoodmap1)):
        for j in range(len(neighborhoodmap1[i])):
            if neighborhoodmap1[i][j] != neighborhoodmap2[i][j]:
                return False, i, j
    return True, None,None


# In[29]:


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
    
    '''Returns list of shared points or common points in between points in @polygonRef and @unionpoints'''
    polygonpoints = polygonRef.GetBoundary().GetPoints()
    #unionpoints = unionPolygon.GetBoundary().GetPoints()
    outerpolygon_points = []
    for ppt in polygonpoints:
        if ppt in unionpoints:
            outerpolygon_points +=[ppt]
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


def getOuterBoundaryForMap(polygon_layer):
    
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


def constructSharedPointsDatabase():
    '''Returns 1) list myPolyPoints of unique points (x1,y1) in a map denoted by a gis layer and 2) 
mappings from polygon_id to list of <i,j> which indicates that points from index i to j 
inclusive belong to polygoin_id. It is created from mapptables that contains
indices entries <polyid1,polyid2,npoints> which means that polygons with id polyid1 and polyid2 share
npoints. Occurance of entries <polyid1,polyid2, npoints> should not be messed up.'''
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
    unionpoints = getOuterBoundaryForMap(polygon_layer).GetBoundary().GetPoints()    
    for fid in range(howmany):
        feature = polygon_layer.GetFeature(fid)
        polygon = feature.GetGeometryRef()
        outerpolygonpts = getSharedPoints(polygon, unionpoints)
        if insert_single_shared_pt:
            outer_polygonpoints += outerpolygonpts
            if len(outerpolygonpts) >= 1:
                sharedpolyidx += [(fid,'0',len(outerpolygonpts))]
        elif len(outerpolygonpts) > 1:
            outer_polygonpoints += outerpolygonpts
            sharedpolyidx += [(fid,'0',len(outerpolygonpts))]
        feature.Destroy()
    allpoints +=outer_polygonpoints
    print("completed.Returns list of unique points and mapping from polygon to point indices.") 
    return allpoints, getPolygonPointsIndicesMapping(sharedpolyidx)


# In[24]:

import numpy as np
from  osgeo import ogr, osr
driver = ogr.GetDriverByName('ESRI Shapefile')
shp = driver.Open(r'/home/naresh-ad/workspace/machinelrn/data/gis/county/county-ms/stco.shp')
layer = shp.GetLayer()
neighbormap = getNeighborsT(layer)
allpoints, sharedpolyidx = constructSharedPointsDatabase()


# In[35]:

# Get Projection from layer
import operator
pdict = getPolygonsIds(layer)
sorted_pdict = sorted(pdict.items(), key=operator.itemgetter(1))

mypoints = []
externalp = []
extppttdiffs= []
cnt = 0.0
newton = []
for pname,pid in sorted_pdict:
    feature = layer.GetFeature(pid)
    polygon = feature.GetGeometryRef()    
    pidxes = sharedpolyidx[pid]
    cnt +=1
    mypoints = []
    for i,j in pidxes:
        if i > j:
            print('i>j'), i, j, pname,pid
        mypoints += allpoints[i:j+1]
    if pname == 'lawrence':
        newton = mypoints
    externalp +=[pname]
    extppttdiffs +=[len(polygon.GetBoundary().GetPoints())-1- len(mypoints)]

print len(extppttdiffs),extppttdiffs


get_ipython().magic(u'matplotlib inline')
layer = shp.GetLayer()
f1 = layer.GetFeature(80)

f2 = layer.GetFeature(51)
p1 = f1.GetGeometryRef()
p2 = f2.GetGeometryRef()
p1points = p1.GetBoundary().GetPoints()
p2points = p2.GetBoundary().GetPoints()
print f1.GetField('CONAME'),f2.GetField('CONAME')

sequence = []
ipoints = []
ulso = ogr.Geometry(ogr.wkbLineString)
inls = p1.Intersection(p2)
idx, outerpointsp1 = prunePolygonPoints(p1, inls)
print len(p1points),len(outerpointsp1), isinstance(inls,ogr.Geometry)

print inls.GetGeometryName(),inls.GetGeometryCount(), inls.GetGeometryType(), inls.Boundary(), inls.GetPointCount(),inls.GetPointCount()
for linestr in inls:
    ulso = ulso.Union(linestr)
    ipoints += linestr.GetPoints()
    sequence += [linestr.GetPoints()[0]]

#polygon-1    
p1x = [a/1000 for a,b in p1points]
p1y = [b/1000 for a,b in p1points]   

#polygon-2
p2x = [a/1000 for a,b in p2points]
p2y = [b/1000 for a,b in p2points]   

#p1.intersection(p2), points in MLS
x = [a/1000 for a,b in ipoints]
y = [b/1000 for a,b in ipoints]     

#p1.intersection(p2), sequence
sx = [a/1000 for a,b in sequence]
sy = [b/1000 for a,b in sequence]   #use sequence.  

#p1 - intersections
opx = [a/1000 for a,b in outerpointsp1[0:1428]]
opy = [b/1000 for a,b in outerpointsp1[0:1428]]   #use sequence.  


opx = [a/1000 for a,b in outerpointsp1[0:3000]]
opy = [b/1000 for a,b in outerpointsp1[0:3000]] 


topx = [a/1000 for a,b in newton]
topy = [b/1000 for a,b in newton] 


fig = plt.figure()
fig.set_figwidth(12)
fig.set_figheight(12)

ax1 = fig.add_subplot(231)
ax1.plot(p1x,p1y, 'r-')

ax2 = fig.add_subplot(232)
ax2.plot(p2x,p2y, 'k-')

ax3 = fig.add_subplot(233)
ax3.plot(x,y, 'b-')

ax4 = fig.add_subplot(234)
ax4.plot(opx,opy, 'b-')


ax4 = fig.add_subplot(235)
ax4.plot(topx,topy, 'b-')
plt.show() 


# In[ ]:

fig = plt.figure()
ax = fig.gca() 
ax.plot(topx,topy)
plt.show()   


# In[ ]:

get_ipython().magic(u'matplotlib auto')
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
fig.set_figwidth(12)
fig.set_figheight(12)
ax = fig.gca() 
plt.ion()
xs =[]
ys = []
for i in range(2000):
    xs += [opx[i]]
    ys += [opy[i]]
    ax.plot(xs,ys)
    plt.pause(0.05)

while True:
    plt.pause(0.005)


# In[ ]:

import sys
minc = sys.maxint
minf = ''
for feat,count in dictfeat_ptcnt.items():
    if minc > count:
        minc = count
        minf = feat
print minc, minf


# In[ ]:


from osgeo import ogr
shp= "yourfile"

shapef = ogr.Open(shpor)
lyr = shapef.GetLayer()
unionc = ogr.Geometry(ogr.wkbMultiPolygon)
for feat in xrange(lyr.GetFeatureCount()):
    fit= lyr.GetFeature(feat)
    geom= fit.GetGeometryRef()
    unionc.AddGeometry(geom)
union= unionc.UnionCascaded()



from osgeo import ogr
test = ogr.Open('polygons.shp')
layer = test.GetLayer()
print layer.GetGeomType()
3 # -> polygons
# empty geometry
union_poly = ogr.Geometry(ogr.wkbPolygon)
# make the union of polygons
for feature in layer:
      geom =feature.GetGeometryRef()
      union_poly = union_poly.Union(geom)


# In[ ]:

from shapely.geometry import box, Polygon
import rtree

gridcell_shape = box(129.5, -27.0, 129.75, 27.25)
# Populate R-tree index with bounds of grid cells
for pos, cell in enumerate(gridcell_shape):
    # assuming cell is a shapely object
    idx.insert(pos, cell.bounds)
    
# Example polygon 
xy = [[130.21001, 27.200001], [129.52, 27.34], [129.45, 27.1], [130.13, 26.950001]]
polygon_shape = Polygon(xy)
# Example grid cell
gridcell_shape = box(129.5, -27.0, 129.75, 27.25)
# The intersection
polygon_shape.intersection(gridcell_shape).area


# In[ ]:

from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
from descartes import PolygonPatch
import shapely 
wkt1 = "POLYGON ((1208064.271243039 624154.6783778917, 1208064.271243039 601260.9785661874, 1231345.9998651114 601260.9785661874, 1231345.9998651114 624154.6783778917, 1208064.271243039 624154.6783778917))"
wkt2 = "POLYGON ((1199915.6662253144 633079.3410163528, 1199915.6662253144 614453.958118695, 1219317.1067437078 614453.958118695, 1219317.1067437078 633079.3410163528, 1199915.6662253144 633079.3410163528)))"

poly1 = ogr.CreateGeometryFromWkt(wkt1)
poly2 = ogr.CreateGeometryFromWkt(wkt2)


polygon1 = shapely.geometry.base.geom_from_wkt(wkt1)
polygon2 = shapely.geometry.base.geom_from_wkt(wkt2)
intersection = polygon1.intersection(polygon2)
print("polygon1"), polygon1


ls1= "LINESTRING (1208064.271243039 614453.958118695, 1208064.271243039 624154.6783778917)" #, 1219317.106743708 624154.6783778917, 1219317.106743708 614453.958118695, 1208064.271243039 614453.958118695)"

ls2= "LINESTRING (1208064.271243039 624154.6783778917, 1219317.106743708 624154.6783778917)" #, 1219317.106743708 614453.958118695, 1208064.271243039 614453.958118695)"

ls1 = shapely.geometry.base.geom_from_wkt(ls1)
ls2 = shapely.geometry.base.geom_from_wkt(ls2)
uls = ls1.union(ls2)
upoints = []
for ls in uls:
    upoints += list(ls.coords)

print ls1; print ls2; print uls
p1 = list(polygon1.exterior.coords)   
p2 = list(polygon2.exterior.coords)   

realdiff = list(set(p1) - set(upoints))
print realdiff
x = [a/1000 for a,b in p1]
y = [b/1000 for a,b in p1]

x1 = [a/1000 for a,b in p2]
y1 = [b/1000 for a,b in p2]

x2 = [a/1000 for a,b in intersections[0:3]]
y2 = [b/1000 for a,b in intersections[0:3]]

dx1 = [a/1000 for a,b in realdiff]
dy1 = [b/1000 for a,b in realdiff]

fig = plt.figure()
ax = fig.gca() 
ax.plot(x,y)
ax.plot(x1,y1)
ax.plot(dx1,dy1)

ax.axis('scaled')
plt.show()   


# In[ ]:

from shapely.geometry import LineString
a = LineString([(0, 0), (1, 1), (1,2), (2,2)])
b = LineString([(0, 0), (1, 1), (2,1), (2,2)])

pa = list(a.coords)
pb = list(b.coords)
diff1 = set(pa) - set(pb)
diff2 = set(pb) - set(pa)
print diff1, diff2
print(".........")
x = a.intersection(b)
print x
d1 = b.difference(a)
d2 = a.difference(b)
print "a.union(b)"
ms = a.union(b)
print ms
upoints = []
for ls in ms:
    upoints +=list(ls.coords)
print upoints
#print "diff"
#print d1,d2

S1 = set(upoints)
print S1


# In[ ]:

from shapely.geometry import Point
a = Point(1, 1).buffer(1.5)
b = Point(2, 1).buffer(1.5)
d = b.difference(a)
x,y = d.boundary.coords.xy
fig = plt.figure()
ax = fig.gca() 
ax.plot(x,y)

ax.axis('scaled')
plt.show()  

