"""Affine transforms, both in general and specific, named transforms."""

'''Note:https://github.com/Toblerity/Shapely/blob/master/shapely/geometry/base.py'''
from math import sin, cos, tan, pi
import os
from  osgeo import ogr, osr


def pyGetLinearRing(cords):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in cords:
        ring.AddPoint(coord[0], coord[1])
    return ring
def pyGetPolygon(cords):
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = pyGetLinearRing(cords)
    return poly.AddGeometry(ring)

def pyGetPolygonFromRings(rings):
    poly = ogr.Geometry(ogr.wkbPolygon)
    for ring in rings: poly.AddGeometry(ring)
    return poly

def affine_transform(geom, matrix,enlar=False):
    """Returns a transformed geometry using an affine transformation matrix.
    The coefficient matrix is provided as a list or tuple with 6 or 12 items
    for 2D or 3D transformations, respectively.
    For 2D affine transformations, the 6 parameter matrix is::
        [a, b, d, e, xoff, yoff]
    which represents the augmented matrix::
                            / a  b xoff \ 
        [x' y' 1] = [x y 1] | d  e yoff |
                            \ 0  0   1  /
    or the equations for the transformed coordinates::
        x' = a * x + b * y + xoff
        y' = d * x + e * y + yoff
    For 3D affine transformations, the 12 parameter matrix is::
        [a, b, c, d, e, f, g, h, i, xoff, yoff, zoff]
    which represents the augmented matrix::
                                 / a  b  c xoff \ 
        [x' y' z' 1] = [x y z 1] | d  e  f yoff |
                                 | g  h  i zoff |
                                 \ 0  0  0   1  /
    or the equations for the transformed coordinates::
        x' = a * x + b * y + c * z + xoff
        y' = d * x + e * y + f * z + yoff
        z' = g * x + h * y + i * z + zoff
    """
    enlarge=enlar
    if geom.IsEmpty():
        return geom
    if len(matrix) == 6:
        ndim = 2
        a, b, d, e, xoff, yoff = matrix
        if geom.GetCoordinateDimension() == 3:
            ndim = 3
            i = 1.0
            c = f = g = h = zoff = 0.0
            matrix = a, b, c, d, e, f, g, h, i, xoff, yoff, zoff
    elif len(matrix) == 12:
        ndim = 3
        a, b, c, d, e, f, g, h, i, xoff, yoff, zoff = matrix
        if not (geom.GetCoordinateDimension() == 3):
            ndim = 2
            matrix = a, b, d, e, xoff, yoff
    else:
        raise ValueError("'matrix' expects either 6 or 12 coefficients")

    def affine_pts(pts):
        """Internal function to yield affine transform of coordinate tuples"""
        if not enlarge:
            if ndim == 2:
                for x, y in pts:
                    xp = a * x + b * y + xoff
                    yp = d * x + e * y + yoff
                    yield (xp, yp)
            elif ndim == 3:
                for x, y, z in pts:
                    xp = a * x + b * y + c * z + xoff
                    yp = d * x + e * y + f * z + yoff
                    zp = g * x + h * y + i * z + zoff
                    yield (xp, yp, zp)
        else:
            if ndim == 2:
                for x, y in pts:
                    xp = a * x + b * y + xoff
                    yp = d * x + e * y + yoff
                    yield (int(round(enlarge*xp,0)), int(round(enlarge*yp,0)))
            elif ndim == 3:
                for x, y, z in pts:
                    xp = a * x + b * y + c * z + xoff
                    yp = d * x + e * y + f * z + yoff
                    zp = g * x + h * y + i * z + zoff
                    yield (int(round(enlarge*xp,0)), int(round(enlarge*yp,0)), int(round(enlarge*zp,0)))
                   
    # Process coordinates from each supported geometry type
    if geom.GetGeometryName() in ('POINT', 'LINESTRING', 'LINEARRING'):
        print("object type %r"%geom.GetGeometryName())
        cords=list(affine_pts(geom.GetPoints()))
        print("cords: "),cords
        if geom.GetGeometryName() == 'POINT':
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(cords[0][0],cords[0][1])
            return point
        elif geom.GetGeometryName() == 'LINEARRING':
            return pyGetLinearRing(cords)
        elif geom.GetGeometryName() =='LINEARRING':
            return pyGetLinearRing(cords)
        
    elif geom.GetGeometryName() in ('POLYGON'):
        print("object type %r"%geom.GetGeometryName())
        nr= geom.GetGeometryCount()
        poly = ogr.Geometry(ogr.wkbPolygon)
        rings =[geom.GetGeometryRef(i) for i in range(nr)]            
        for pos, ring in enumerate(rings):
            rings[pos] = pyGetLinearRing(list(affine_pts(ring.GetPoints())))
            poly.AddGeometry(rings[pos])
        return poly
        
    elif geom.GetGeometryName().startswith('MULTI') or geom.GetGeometryName() == "Collection":
        # Recursive call
        # TODO: fix GeometryCollection constructor
        print("m - object type %r"% geom.GetGeometryName())
        mpoly = ogr.Geometry(ogr.wkbMultiPolygon)
        np = geom.GetGeometryCount()
        polys=[geom.GetGeometryRef(k) for k in range(np)]
        for j in range(np):
            nr= polys[j].GetGeometryCount()
            poly = ogr.Geometry(ogr.wkbPolygon)
            rings =[polys[j].GetGeometryRef(i) for i in range(nr)]            
            for pos, ring in enumerate(rings):
                rings[pos] = pyGetLinearRing(list(affine_pts(ring.GetPoints())))
                poly.AddGeometry(rings[pos])
            mpoly.AddGeometry(poly) 
            poly=rings=nr=None #endfor
        #endfor
        polys=None
        return mpoly
    else:
        raise ValueError('Type %r not recognized' % geom.GetGeometryName())

def interpret_origin(geom, origin, ndim):
    """Returns interpreted coordinate tuple for origin parameter.
    This is a helper function for other transform functions.
    The point of origin can be a keyword 'center' for the 2D bounding box
    center, 'centroid' for the geometry's 2D centroid, a Point object or a
    coordinate tuple (x0, y0, z0).
    """
    # get coordinate tuple from 'origin' from keyword or Point type
    if origin == 'center':
        # bounding box center
        minx, miny, maxx, maxy = geom.bounds
        origin = ((maxx + minx)/2.0, (maxy + miny)/2.0)
    elif origin == 'centroid':
        origin = geom.centroid.coords[0]
    elif isinstance(origin, str):
        raise ValueError("'origin' keyword %r is not recognized" % origin)
    elif hasattr(origin, 'type') and origin.type == 'Point':
        origin = origin.coords[0]

    # origin should now be tuple-like
    if len(origin) not in (2, 3):
        raise ValueError("Expected number of items in 'origin' to be "
                         "either 2 or 3")
    if ndim == 2:
        return origin[0:2]
    else:  # 3D coordinate
        if len(origin) == 2:
            return origin + (0.0,)
        else:
            return origin


def rotate(geom, angle, origin='center', use_radians=False):
    """Returns a rotated geometry on a 2D plane.
    The angle of rotation can be specified in either degrees (default) or
    radians by setting ``use_radians=True``. Positive angles are
    counter-clockwise and negative are clockwise rotations.
    The point of origin can be a keyword 'center' for the bounding box
    center (default), 'centroid' for the geometry's centroid, a Point object
    or a coordinate tuple (x0, y0).
    The affine transformation matrix for 2D rotation is:
      / cos(r) -sin(r) xoff \ 
      | sin(r)  cos(r) yoff |
      \   0       0      1  /
    where the offsets are calculated from the origin Point(x0, y0):
        xoff = x0 - x0 * cos(r) + y0 * sin(r)
        yoff = y0 - x0 * sin(r) - y0 * cos(r)
    """
    if not use_radians:  # convert from degrees
        angle *= pi/180.0
    cosp = cos(angle)
    sinp = sin(angle)
    if abs(cosp) < 2.5e-16:
        cosp = 0.0
    if abs(sinp) < 2.5e-16:
        sinp = 0.0
    x0, y0 = interpret_origin(geom, origin, 2)

    matrix = (cosp, -sinp, 0.0,
              sinp,  cosp, 0.0,
              0.0,    0.0, 1.0,
              x0 - x0 * cosp + y0 * sinp, y0 - x0 * sinp - y0 * cosp, 0.0)
    return affine_transform(geom, matrix)


def scale(geom, xfact=1.0, yfact=1.0, zfact=1.0, origin='center'):
    """Returns a scaled geometry, scaled by factors along each dimension.
    The point of origin can be a keyword 'center' for the 2D bounding box
    center (default), 'centroid' for the geometry's 2D centroid, a Point
    object or a coordinate tuple (x0, y0, z0).
    Negative scale factors will mirror or reflect coordinates.
    The general 3D affine transformation matrix for scaling is:
        / xfact  0    0   xoff \ 
        |   0  yfact  0   yoff |
        |   0    0  zfact zoff |
        \   0    0    0     1  /
    where the offsets are calculated from the origin Point(x0, y0, z0):
        xoff = x0 - x0 * xfact
        yoff = y0 - y0 * yfact
        zoff = z0 - z0 * zfact
    """
    x0, y0, z0 = interpret_origin(geom, origin, 3)

    matrix = (xfact, 0.0, 0.0,
              0.0, yfact, 0.0,
              0.0, 0.0, zfact,
              x0 - x0 * xfact, y0 - y0 * yfact, z0 - z0 * zfact)
    return affine_transform(geom, matrix)



def skew(geom, xs=0.0, ys=0.0, origin='center', use_radians=False):
    """Returns a skewed geometry, sheared by angles along x and y dimensions.
    The shear angle can be specified in either degrees (default) or radians
    by setting ``use_radians=True``.
    The point of origin can be a keyword 'center' for the bounding box
    center (default), 'centroid' for the geometry's centroid, a Point object
    or a coordinate tuple (x0, y0).
    The general 2D affine transformation matrix for skewing is:
        /   1    tan(xs) xoff \ 
        | tan(ys)  1     yoff |
        \   0      0       1  /
    where the offsets are calculated from the origin Point(x0, y0):
        xoff = -y0 * tan(xs)
        yoff = -x0 * tan(ys)
    """
    if not use_radians:  # convert from degrees
        xs *= pi/180.0
        ys *= pi/180.0
    tanx = tan(xs)
    tany = tan(ys)
    if abs(tanx) < 2.5e-16:
        tanx = 0.0
    if abs(tany) < 2.5e-16:
        tany = 0.0
    x0, y0 = interpret_origin(geom, origin, 2)

    matrix = (1.0, tanx, 0.0,
              tany, 1.0, 0.0,
              0.0,  0.0, 1.0,
              -y0 * tanx, -x0 * tany, 0.0)
    return affine_transform(geom, matrix)


def translate(geom, xoff=0.0, yoff=0.0, zoff=0.0,enlarge=False):
    """Returns a translated geometry shifted by offsets along each dimension.
    The general 3D affine transformation matrix for translation is:
        / 1  0  0 xoff \ 
        | 0  1  0 yoff |
        | 0  0  1 zoff |
        \ 0  0  0   1  /
    """
    matrix = (1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0,
              xoff, yoff, zoff)
    return affine_transform(geom, matrix,enlarge)

#Testing.
'''
import pickle
home =r'/home/naresh-1/workspace/machinelrn/data/gis/'
tv = pickle.load(open(home+r'out/translation_vector.p','rb'))
xoff,yoff=tv['xoff'],tv['yoff']
point = ogr.Geometry(ogr.wkbPoint)
poly = ogr.Geometry(ogr.wkbPolygon)
ring = ogr.Geometry(ogr.wkbLinearRing)
mpoly = ogr.Geometry(ogr.wkbMultiPolygon)

ring.AddPoint(2,1)
ring.AddPoint(2,2)
ring.AddPoint(0,2)
ring.AddPoint(2,1)
poly.AddGeometry(ring)
mpoly.AddGeometry(poly)
mpoly.AddGeometry(poly)
point.AddPoint(2,1)

print point.GetPoints(),point.IsEmpty(),point.GetCoordinateDimension()
print ring.GetPoints(),ring.IsEmpty(),ring.GetCoordinateDimension()
print poly.GetBoundary(),poly.IsEmpty(),poly.GetCoordinateDimension()
print
pt1=translate(point,xoff,yoff)
pt2=translate(ring,xoff,yoff)
pt3=translate(poly,xoff,yoff)
pt4=translate(mpoly,xoff,yoff)
print
print("translated"),point, pt1.GetGeometryName()
print("translated"),ring, pt2.GetGeometryName()
print("translated"),poly, pt3.GetGeometryName()
print("translated"),mpoly,pt4.GetGeometryName()
'''