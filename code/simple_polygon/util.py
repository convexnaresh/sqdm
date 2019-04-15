import hashlib
import sys
import json
import os.path
import random
import pickle
from operator import itemgetter
from collections import OrderedDict
from ast import literal_eval
if sys.version_info < (3, 6):
    import sha3


#CONSTANTS
hindexi = 0
hindexj = 12
SCALE_FACTOR = 344
debug = False
SEGMENT_ID_TYPES = {'linehash': 'linehash', 'lineid': 'lineid', 'edgeid': 'edgeid'}
SEGMENT_TUPLE_KEYS=["x1","y1","x2","y2","edgeid","abv","bel"]
POLYGON_INNERID = "PHI"
POLYGON_OUTSIDEID = "GLOBE"


def hashargs(self, *args):
    inp = args
    s = hashlib.sha3_224()
    s = hashlib.sha3_256()
    for item in inp:
        s.update(str(item))
    return s.hexdigest()[hindexi:hindexj]

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    from itertools import tee, izip
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def poly_ptstoshp(polypts,outfilename,save_as_multipt=False):
    epsgdic = {'nad83': 4269, 'wgs84': 4326, 'pseudoutm': 3857, 'worldmercater': 3395}

    from shapely.geometry import LinearRing, Point, mapping
    from osgeo import ogr, osr
    import os
    save_as_multipt = save_as_multipt
    # load initial vertex order files.
    # save to shp.
    if save_as_multipt:
        outfilename +='as_points'

    outfilename += '.shp'
    driver = ogr.GetDriverByName('ESRI Shapefile')
    outprjref = osr.SpatialReference()
    outprjref.ImportFromEPSG(epsgdic["worldmercater"])

    if os.path.exists(outfilename):
        driver.DeleteDataSource(outfilename)
    outDataSet = driver.CreateDataSource(outfilename)
    if save_as_multipt:
        outLayer = outDataSet.CreateLayer("mystates", outprjref, geom_type=ogr.wkbMultiPoint)
    else:
        outLayer = outDataSet.CreateLayer("mystates", outprjref, geom_type=ogr.wkbMultiPolygon)

    outLayer.CreateField(ogr.FieldDefn('STATEFP'), ogr.OFTInteger) #create new field/column/attribute
    outLayerDefn = outLayer.GetLayerDefn()
    poly = ogr.Geometry(ogr.wkbPolygon)

    #save as point feature
    if not save_as_multipt:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x1,y1 in polypts:
            x,y = float(x1),float(y1)
            ring.AddPoint(x, y)
        poly.AddGeometry(ring)

        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(poly)

    #save as polygon feature
    if save_as_multipt:
        mp = ogr.Geometry(ogr.wkbMultiPoint)
        for x1,y1 in polypts:
            point1 = ogr.Geometry(ogr.wkbPoint)
            point1.AddPoint(float(x1),float(y1))
            mp.AddGeometry(point1)
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(mp)

    outFeature.SetField(outLayerDefn.GetFieldDefn(0).GetNameRef(), 99) #polygon id 99
    outLayer.CreateFeature(outFeature)
    # end-for
    outLayer = outFeature = poly  = ring = outLayerDefn = outLayer = outDataSet = None
    print("polygion points is saved in .shp file "),outfilename

def  non_overlaping_yspans(ytuples):
    '''
    :param ytuples: a list of Y-Span (yi,yj).
    :return: list of (yi',yj') such that every tuple forms a non-overlapping blocks
    on top of one another; construct a cover-up y-span such that the y-span's form a complete lut-entries.
    For example:
    (1,3,a) and (2,4,b) are overlapping y-spans; it this function resulsts tuples (1,4,{a,b}) and (4,1,{})
    tuples.
    '''

    #sort each tuple in ytuples.ytuples[jdx-1][1]
    for i in range(len(ytuples)):
        if ytuples[i][1] < ytuples[i][0]:
            ytuples[i] = (ytuples[i][1],ytuples[i][0]) + ytuples[i][2:]

    #sort list of tuples by a tuple's first element.
    #slower version:
    #ytuples = sorted(ytuples, key=lambda tuple: tuple[0])
    #faster version:
    ytuples = sorted(ytuples,key=itemgetter(0))
    recs = { }
    #now merge the overlapping y-spans
    idx = 0
    jdx = idx
    while idx < len(ytuples):
        y0,y1 = ytuples[idx]
        try:
            recs[y0][0].append(ytuples[idx])
        except:
            recs[y0] = []
            #recs[y0] += [ytuples[idx]]
            recs[y0].append([ytuples[idx]])
            bot,top = y0,y1 #bottom and top of y-value.
            recs[y0].append((bot, top))
        jdx = idx+1
        while jdx < len(ytuples) and ytuples[jdx][0] < top: #ytuples[jdx-1][1]:
            recs[y0][0].append(ytuples[jdx])
            #recs[y0][1]=(min(recs[y0][1][0],ytuples[jdx][0]),max(recs[y0][1][1],ytuples[jdx][1]))
            bot, top = min(bot, ytuples[jdx][0]),max(top, ytuples[jdx][1])
            recs[y0][1] = bot, top
            jdx += 1
        #endwhile
        idx = jdx
    #endwhile
    #from collections import OrderedDict
    #ordred = OrderedDict(sorted(recs.items(),key=lambda t:t[0]))
    yvalues=[]
    #extract bot/top values of each rectangles.
    for k,v in recs.items():
        ybot,ytop = v[1]
        yvalues +=[ybot,ytop]
    yvalues = sorted(list(set(yvalues)))
    return yvalues

def save_sqdm(data,outfile,dosort=True):

    newdata ={}
    for k1, dv1 in data.items():
        newdata[float(k1)] = OrderedDict()

        for k2, dv2 in dv1.items():
            newdata[k1][str(k2)] = dv2

    #delete if exists
    if os.path.exists(outfile):
        if not os.path.isdir(outfile) : os.remove(outfile)
        else: outfile += '/data-'
    outfile = outfile.replace(".json",'')
    outfile = outfile.replace(".p",'')
    outfile = outfile.replace(".pickle",'')
    with open(outfile+'.json', 'w') as fp:
        json.dump(newdata, fp, sort_keys=dosort, indent=4)

    print("json data saved in "), outfile + '.json'
def save(data,outfile,dosort=True):

    #delete if exists
    if os.path.exists(outfile):
        if not os.path.isdir(outfile) : os.remove(outfile)
        else: outfile += '/data-'
    outfile = outfile.replace(".json",'')
    outfile = outfile.replace(".p",'')
    outfile = outfile.replace(".pickle",'')
    with open(outfile+'.json', 'w') as fp:
        json.dump(data, fp, sort_keys=dosort, indent=4)

    print("json data saved in "), outfile + '.json'

def load_sqdm_from_file(infile,typecast=[]):
    import json
    data = None
    infile = infile.replace(".json",'')
    infile = infile.replace(".p",'')
    infile = infile.replace(".pickle",'')
    for ext in ['.json','.pickle']:
        try:
            with open(infile + ext, 'r') as fp:
                data = json.load(fp, object_pairs_hook=OrderedDict)
        except:
            continue

    #type conversion
    newdata = OrderedDict()
    for k1, dv1 in data.items():
        newdata[float(k1)] = OrderedDict()

        for k2, dv2 in dv1.items():
            newdata[float(k1)][literal_eval(k2)] = dv2

    return newdata

def load_from_file(infile,typecast=[]):
    import json
    data = None
    infile = infile.replace(".json",'')
    infile = infile.replace(".p",'')
    infile = infile.replace(".pickle",'')
    for ext in ['.json', '.pickle']:
        try:
            with open(infile + ext, 'r') as fp:
                data = json.load(fp, object_pairs_hook=OrderedDict)
        except Exception, e:
            print(""), str(e)
            continue

    tmp = OrderedDict()
    return data

if debug:
    yp =[(3.4000000000000004, 3.7), (3.7, 4.6), (1.0, 1.0), (4.6, 4.8999999999999995), (5.8, 5.95), (5.2, 5.8)]

    nonovr = non_overlaping_yspans(yp)
    print nonovr
    print
    pairs = pairwise(nonovr)
    for p in pairs:
        print p

    print hashargs("anders")
    print hashlib.sha256("anders").hexdigest()
