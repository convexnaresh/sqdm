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
F2_SEGMENT_TUPLE_KEYS=["x1","y1","x2","y2","edgeid","abv","bel","ischild"]
POLYGON_INNERID = "PHI"
POLYGON_OUTSIDEID = "GLOBE"
CHILD_ANNOTATION = 2
SQDM_DELEGATED_BY_ID = "USA"
SQDM_DELEGATED_TO_ID= "CALIF"
COUNTRY_DELEGATED_BY_ID = 0.0
extents={'dummy_usa_extent': [(1200, 600), (5200, 3000)]}

EQT = {}

def hashargs(*args):
    inp = args
    #s = hashlib.sha3_224()
    s = hashlib.sha3_256()
    for item in inp:
        s.update(str(item))
    return s.hexdigest()[hindexi:hindexj]

def hash_cblocks_dist(block_components):
    inp = block_components
    #s = hashlib.sha3_224()
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

def save_xun_asshp(xun, outfilename, mapextent=[]):
    epsgdic = {'nad83': 4269, 'wgs84': 4326, 'pseudoutm': 3857, 'worldmercater': 3395}

    from shapely.geometry import LinearRing, Point, mapping
    from osgeo import ogr, osr
    import os
    # load initial vertex order files.
    # save to shp.
    outfilename += '.shp'
    driver = ogr.GetDriverByName('ESRI Shapefile')
    outprjref = osr.SpatialReference()
    outprjref.ImportFromEPSG(epsgdic["worldmercater"])

    if os.path.exists(outfilename):
        #driver.DeleteDataSource(outfilename)
        pass
    outDataSet = driver.CreateDataSource(outfilename)

    outLayer = outDataSet.CreateLayer("mystates", outprjref, geom_type=ogr.wkbMultiLineString)

    outLayer.CreateField(ogr.FieldDefn('STATEFP'), ogr.OFTInteger) #create new field/column/attribute
    outLayerDefn = outLayer.GetLayerDefn()
    poly = ogr.Geometry(ogr.wkbPolygon)

    #save as polygon feature
    mls = ogr.Geometry(ogr.wkbMultiLineString)
    ymin,ymax = mapextent[0][1],mapextent[1][1]
    for x1 in xun:
        ls = ogr.Geometry(ogr.wkbLineString)
        ls.AddPoint(float(x1),float(ymin))
        ls.AddPoint(float(x1),float(ymax))
        mls.AddGeometry(ls)

    outFeature = ogr.Feature(outLayerDefn)
    outFeature.SetGeometry(mls)

    outFeature.SetField(outLayerDefn.GetFieldDefn(0).GetNameRef(), 99) #polygon id 99
    outLayer.CreateFeature(outFeature)
    # end-for
    outLayer = outFeature = poly  = ring = outLayerDefn = outLayer = outDataSet = None
    print("Vertical lines are saved as shape file in "),outfilename

def save_segments_asshp(segments, outfilename, mapextent=[]):
    epsgdic = {'nad83': 4269, 'wgs84': 4326, 'pseudoutm': 3857, 'worldmercater': 3395}

    from shapely.geometry import LinearRing, Point, mapping
    from osgeo import ogr, osr
    import os
    # load initial vertex order files.
    # save to shp.
    outfilename += '.shp'
    driver = ogr.GetDriverByName('ESRI Shapefile')
    outprjref = osr.SpatialReference()
    outprjref.ImportFromEPSG(epsgdic["worldmercater"])

    if os.path.exists(outfilename):
        driver.DeleteDataSource(outfilename)
    outDataSet = driver.CreateDataSource(outfilename)

    outLayer = outDataSet.CreateLayer("mystates", outprjref, geom_type=ogr.wkbMultiLineString)

    outLayer.CreateField(ogr.FieldDefn('STATEFP'), ogr.OFTInteger) #create new field/column/attribute
    outLayerDefn = outLayer.GetLayerDefn()
    poly = ogr.Geometry(ogr.wkbPolygon)

    #save as polygon feature
    mls = ogr.Geometry(ogr.wkbMultiLineString)
    try:
        ymin,ymax = mapextent[0][1],mapextent[1][1]
    except:
        print("Warning! Map extent not supplied.")
    for seg in segments:

        ls = ogr.Geometry(ogr.wkbLineString)
        x1, y1, x2, y2 = seg
        ls.AddPoint(float(x1),float(y1))#
        ls.AddPoint(float(x2),float(y2))#
        mls.AddGeometry(ls)

    outFeature = ogr.Feature(outLayerDefn)
    outFeature.SetGeometry(mls)

    outFeature.SetField(outLayerDefn.GetFieldDefn(0).GetNameRef(), 99) #polygon id 99
    outLayer.CreateFeature(outFeature)
    # end-for
    outLayer = outFeature = poly  = ring = outLayerDefn = outLayer = outDataSet = None
    print("Saved lines are saved as shape file in "),outfilename


def save_yun_asshp(xun, outfilename, mapextent=[]):
    epsgdic = {'nad83': 4269, 'wgs84': 4326, 'pseudoutm': 3857, 'worldmercater': 3395}

    from shapely.geometry import LinearRing, Point, mapping
    from osgeo import ogr, osr
    import os
    # load initial vertex order files.
    # save to shp.
    outfilename += '.shp'
    driver = ogr.GetDriverByName('ESRI Shapefile')
    outprjref = osr.SpatialReference()
    outprjref.ImportFromEPSG(epsgdic["worldmercater"])

    if os.path.exists(outfilename):
        driver.DeleteDataSource(outfilename)
    outDataSet = driver.CreateDataSource(outfilename)

    outLayer = outDataSet.CreateLayer("mystates", outprjref, geom_type=ogr.wkbMultiLineString)

    outLayer.CreateField(ogr.FieldDefn('STATEFP'), ogr.OFTInteger) #create new field/column/attribute
    outLayerDefn = outLayer.GetLayerDefn()
    poly = ogr.Geometry(ogr.wkbPolygon)

    #save as polygon feature
    mls = ogr.Geometry(ogr.wkbMultiLineString)
    ymin,ymax = mapextent[0][1],mapextent[1][1]
    xmin,xmax = mapextent[0][0],mapextent[1][0]
    for x1 in xun:
        ls = ogr.Geometry(ogr.wkbLineString)
        ls.AddPoint(float(xmin),float(x1))
        ls.AddPoint(float(xmax),float(x1))
        mls.AddGeometry(ls)

    outFeature = ogr.Feature(outLayerDefn)
    outFeature.SetGeometry(mls)

    outFeature.SetField(outLayerDefn.GetFieldDefn(0).GetNameRef(), 99) #polygon id 99
    outLayer.CreateFeature(outFeature)
    # end-for
    outLayer = outFeature = poly  = ring = outLayerDefn = outLayer = outDataSet = None
    print("Vertical lines are saved as shape file in "),outfilename

def geom_toshp(newgeom,outfilename,save_as_multipt=False):
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
        #driver.DeleteDataSource(outfilename)
        pass
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
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(newgeom)

    #save as polygon feature
    if save_as_multipt:
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(newgeom)

    outFeature.SetField(outLayerDefn.GetFieldDefn(0).GetNameRef(), 99) #polygon id 99
    outLayer.CreateFeature(outFeature)
    # end-for
    outLayer = outFeature = poly  = ring = outLayerDefn = outLayer = outDataSet = None
    print("polygion points is saved in .shp file "),outfilename

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
        #driver.DeleteDataSource(outfilename)
        pass
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

def Set_Equivalance(SQDM_DELEGATED_TO_ID, SQDM_DELEGATED_BY_ID):
    EQT[SQDM_DELEGATED_BY_ID] = SQDM_DELEGATED_TO_ID

def Break_Equivalance(SQDM_DELEGATED_TO_ID, SQDM_DELEGATED_BY_ID):
    if util.EQT[SQDM_DELEGATED_BY_ID]:
        util.EQT[SQDM_DELEGATED_BY_ID] = None

def merkle_root(key_val_map):
    return hashargs('1')

if __name__ == "__main__":

    nonmatchyvalues =[2005256200.0, 2003910676.0, 2006866606.0, 2003939354.0, 1996018351.0, 2057644065.0, 2057609262.0, 1999355954.0, 1999323191.0, 2057590841.0, 2006873439.0, 2006872125.0, 2005282485.0, 2006874176.0, 2003941448.0, 2003931212.0, 1993478222.0, 1993461849.0, 1999366234.0, 2004656230.0, 2006872167.0, 2003914858.0, 2005241972.0, 2006855797.0, 1993437303.0, 2005256313.0, 2024298625.0, 1998977387.0, 2006231181.0, 2005254295.0, 2006868123.0, 2003910818.0, 2006870181.0, 2004652205.0, 1985591472.0, 1993482423.0, 2003648885.0, 2003650752.0, 2003654851.0, 1999352005.0, 2006913057.0, 1936206024.0, 2005207241.0, 2003941591.0, 1993490649.0, 2006860000.0, 1996013264.0, 2003888360.0, 1985593587.0, 1993447668.0, 2005209334.0, 2005211384.0, 2003661050.0, 2004859131.0, 2006872316.0, 1999327485.0, 2003888390.0, 1999358219.0, 1993453837.0, 2003653755.0, 2006921502.0, 2003890466.0, 2003887531.0, 2006855998.0, 2005236033.0, 1999368514.0, 2003644747.0, 2003894607.0, 2006853974.0, 2057607512.0, 2029355364.0, 2029361509.0, 2005248360.0, 2006913387.0, 1985593709.0, 2003648880.0, 1999362421.0, 1999358353.0, 2006217121.0, 2029341096.0, 1999319467.0, 2003655085.0, 2003933615.0, 2015576501.0, 2003935670.0, 2005266874.0, 2003913149.0, 2057632190.0, 2003927489.0, 2003909058.0, 1962254788.0, 2006225354.0, 1993470412.0, 2006917581.0, 2003923409.0, 2006854100.0, 2006862294.0, 2006854104.0, 2005225950.0, 1999366623.0, 2034835940.0, 2029353445.0, 2023543271.0, 1999354345.0, 1993449964.0, 2004650408.0, 2003917299.0, 2006856186.0, 2005238272.0, 1936157099.0, 2003886598.0, 2006874632.0, 2003917321.0, 2003642890.0, 2001500685.0, 2003902995.0, 2023598612.0, 2006917656.0, 2003898914.0, 2003638022.0, 2003903015.0, 2003888687.0, 2006913587.0, 1993488957.0, 2013755989.0, 2003901020.0, 2057636965.0, 1962228323.0, 2003946090.0, 1993495147.0, 1993497197.0, 1993484733.0, 2003890802.0, 2003923571.0, 2006219380.0, 2005215867.0, 2005242496.0, 2006874754.0, 2003946117.0, 1999368842.0, 2005236368.0, 2003661458.0, 2006913690.0, 2006870684.0, 1999366816.0, 2004654756.0, 1999364775.0, 2006862510.0, 2057640623.0, 2057646776.0, 2006874814.0, 1993464511.0, 2003917505.0, 2003655371.0, 1998973647.0, 2003886802.0, 2003894995.0, 2006915801.0, 2003641053.0, 2056744543.0, 1993446115.0, 2057640677.0, 2003905255.0, 2057642728.0, 2003929834.0, 1993489134.0, 2003940087.0, 2006913788.0, 2029357826.0, 2003936006.0, 1956274951.0, 2003643149.0, 2005224211.0, 2005212292.0, 2003943672.0, 2014704414.0, 2003892389.0, 2005285668.0, 2006854440.0, 2003913526.0, 2005218106.0, 1985588027.0, 2003647298.0, 1985585988.0, 2006909765.0, 2003925831.0, 1993462606.0, 2006227795.0, 1999321940.0, 2005289816.0, 1993479005.0, 2003921759.0, 1996018529.0, 2029353826.0, 2003660262.0, 2006856557.0, 1999358831.0, 1993481072.0, 2005250934.0, 2006870913.0, 2005242755.0, 1980526402.0, 2003913617.0, 2003919765.0, 1975344039.0, 2003915693.0, 2005223923.0, 1975352249.0, 2013774789.0, 1993472966.0, 2003921871.0, 2005245092.0, 2056762330.0, 2006218917.0, 2003647458.0, 2005218276.0, 2006860774.0, 2003916625.0, 2003932138.0, 2057579499.0, 2003891185.0, 2003917810.0, 2003921912.0, 2006920191.0, 2003899392.0, 2003905538.0, 1985592326.0, 2006873099.0, 2003913740.0, 2005232672.0, 1999340578.0, 1985592355.0, 1993468964.0, 2003895344.0, 2006231903.0, 2003911745.0, 2003641410.0, 1934314563.0, 2001493060.0, 2003889222.0, 2005244999.0, 2003893329.0, 2006221910.0, 2056725605.0, 2005257322.0, 2005224555.0, 2006858868.0, 1966740603.0, 2006873220.0, 2006232198.0, 2003649676.0, 2006873230.0, 1993450639.0, 2003894038.0, 2003651734.0, 2006918295.0, 2003940504.0, 2003657883.0, 2003913884.0, 1999324323.0, 2003930276.0, 1996012713.0, 2003913902.0, 2029360307.0, 2004864180.0, 1998968000.0, 2024297668.0, 2003901641.0, 1993454798.0, 2006875345.0, 2006914259.0, 2004667605.0, 2003915997.0, 2003657951.0, 2034830563.0, 1998978289.0, 2003920119.0, 1993479955.0, 2005255417.0, 2005284090.0, 2003928324.0, 2003647749.0, 2029350150.0, 2056774925.0, 2003898585.0, 2003928344.0, 1999324442.0, 2001487131.0, 2005224734.0, 2003903778.0, 2003934499.0, 2003889449.0, 2005259562.0, 2006856925.0, 2005243184.0, 2057635125.0, 1999361334.0, 1999365431.0, 2006875448.0, 2005234911.0, 2003890293.0, 2003889500.0, 2003928415.0, 1999318379.0, 2006861166.0, 2003914102.0, 2006871416.0, 2003897723.0, 2057635208.0, 2056742296.0, 2003936667.0, 2003942814.0, 2005237155.0, 2001491365.0, 2003928486.0, 2003654056.0, 2003641771.0, 2004645296.0, 2024297908.0, 2005257657.0, 2056781388.0, 1999013322.0, 1993481680.0, 1999013329.0, 2006914515.0, 1999013338.0, 1999013340.0, 1999013341.0, 1999013343.0, 2001497568.0, 2003902715.0, 2057635302.0, 2001499623.0, 1934316455.0, 2057586157.0, 1993457917.0, 2003641848.0, 2006857209.0, 2003645946.0, 2015565311.0, 2003887620.0, 2003904012.0, 2003914258.0, 1985590811.0, 2003906087.0, 2003938856.0, 1993479723.0, 2003922482.0, 2003654239.0, 2003656259.0, 1999326796.0, 2006874381.0, 2001491542.0, 2003648087.0, 2003889754.0, 2003893855.0, 2029346401.0, 2005239394.0, 2006857316.0, 2003945066.0, 1993458282.0, 2003887742.0, 2003887745.0, 1996018283.0, 2003646086.0, 2057586317.0, 1993447055.0, 1993498257.0, 2003943058.0, 2003650195.0, 2029346453.0, 2003642007.0, 2006855587.0, 2003924661.0, 1999361721.0, 1993461434.0, 1993463484.0, 2003654333.0, 2005284551.0, 1993465546.0, 2014718242.0, 2056756942.0, 2006861520.0, 2005264082.0, 1993465558.0, 1999361752.0, 2006855386.0, 1985588955.0, 1993467614.0, 2003900135.0, 2003943149.0, 2005259901.0, 2003941106.0, 2004651764.0, 2006861557.0, 2006859511.0, 2006232838.0, 2003947272.0, 2006914828.0, 1999353614.0, 1999345427.0, 2003654422.0, 2005249817.0, 2005210906.0, 2005288732.0, 2003924769.0, 2003658533.0, 2056750887.0, 2003932973.0, 2029348670.0, 1993445186.0, 1966741321.0, 2057580365.0, 1993479480.0, 1993447255.0, 2006921055.0, 2057638544.0, 2005288807.0, 1980526442.0, 1998960492.0, 2003890039.0, 2006910842.0, 2003941247.0, 2057603393.0, 2003941257.0, 2006859660.0, 2003912590.0, 2005227407.0, 2004668316.0, 2057611171.0, 2057641936.0, 2005252056.0, 2003644378.0, 2057641950.0, 2003648482.0, 2003892197.0, 2003894246.0, 1956274167.0]
    save_yun_asshp(nonmatchyvalues, "Yun",[(618829602, 1930722595), (733248005, 2063037184)])

    yp =[(3.4000000000000004, 3.7), (3.7, 4.6), (1.0, 1.0), (4.6, 4.8999999999999995), (5.8, 5.95), (5.2, 5.8)]

    nonovr = non_overlaping_yspans(yp)
    print nonovr
    print
    pairs = pairwise(nonovr)
    for p in pairs:
        print p

    print hashargs("anders")
    print hashlib.sha256("anders").hexdigest()

    d = {'06ab6966910e': (618829602, 2039865788, 618831383.0, 2039877854.7153616, '3702', 777, 777),
     '81670468e2ae': (618829602, 2039865788, 618831383.0, 2039862860.5236104, '3703', 777, 777)}
    e ={}
    for key,value in d.items():
        e.update({key: [float(v) for v in value[0:4]] + [str(v) for v in value[4:]]  })

    save(e,"aaaaaa.json")