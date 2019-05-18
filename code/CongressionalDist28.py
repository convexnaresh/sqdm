import ogr, osr
from Metrics import Metri
from simple_polygon import util
#blocks for a district.
apblocks={"FEAT_IDS":{}, #block-ids
          "AREA":{},
          "PERI":{},
          "ENV":{},
          "CENTROID":{},
          "POPULATION":{},
          "BLKCTR_TO_DIST_CTR":{},
          "BLKCTR_TO_DIST_PERI":{} #length from centroid of this block to Perimeter
          }
#district properties
apstate={"FEAT_IDS":{},
         "AREA":{},
         "PERI":{},
         "ISO_PERIMETRIC_IDX":{},
         "ENV":{},
         "CENTROID":{},
         "POPULATION":{},
         "MOMENT_AREA":{},
         "MOMENT_POPU":{},
         "HULL_AREA":{},
         "OVERLAP_AREA_EQCIRCLE":{},
         "AREA_CCIRCLE":{},
         "AREA_INCIRCLE":{},
         "EXCHG_IDX":{},
         "ROHRBACH_IDX":{},
         "AREA_HULL_RATIO": {},
         "POP_HULL_RATIO":{},
         "MEAN_RADIUS_TO_DISTCTR":{},
         "DYN_RADIUS_TO_DISTCTR":{},
         "HRM_RADIUS_TO_DISTCTR":{},
         "INTERPERSONAL_DIST":{}
         }
epsgdic = {'nad83':4269,'wgs84':4326,'pseudoutm':3857,'worldmercater':3395}

home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
block_shapefile = home+"cong_2011/Cong_2011.shp"

outprjref = osr.SpatialReference()
outprjref.ImportFromEPSG(epsgdic["worldmercater"])

shapef = ogr.Open(block_shapefile)
layer = shapef.GetLayer()


#for each block, calc properties
udist_geom = ogr.Geometry(ogr.wkbMultiPolygon)
fcount = 0
for feature in layer:
    bid = feature.GetField('DISTRICT')
    blkgeo = feature.GetGeometryRef()
    udist_geom.AddGeometry(blkgeo)
    apblocks["AREA"][bid] = blkgeo.GetArea() #/1000000 #*0.386102
    apblocks["PERI"][bid] = blkgeo.GetBoundary().Length()
    apblocks["ENV"][bid] = blkgeo.GetEnvelope()
    apblocks["CENTROID"][bid] = blkgeo.Centroid().GetPoints()[0]
    apblocks["POPULATION"][bid] = 741776
    apblocks["FEAT_IDS"][bid] = None
    fcount += 1
layer = None
apblocks["POLSBY_POPPER_RATIO"]= Metri.dpolsby_popper(apblocks)
print("Blocks"), apblocks["FEAT_IDS"].keys()


#Union of all blocks is a District
did = 1234 #district
state = udist_geom.UnionCascaded() #state distg
stpts = state.GetBoundary().GetPoints()
util.poly_ptstoshp(stpts, home + "cong_2011/st", save_as_multipt=False)


#Area/Peri/Env/Centroid/Popl
apstate["FEAT_IDS"][did] = None
apstate["AREA"][did] = state.GetArea()
apstate["PERI"][did] = state.GetBoundary().Length()
apstate["ENV"][did] = state.GetEnvelope()
apstate["CENTROID"][did] = state.Centroid().GetPoints()[0]
apstate["POPULATION"][did] = 741776*fcount
apblocks["BLKCTR_TO_DIST_PERI"]= Metri.blkctr_to_dist_peri(state,apblocks)


#centroids and distance btw Blk-centroid and Dist-centroid
apblocks["BLKCTR_TO_DIST_CTR"] = Metri.dist_to_point(apblocks,apstate["CENTROID"][did])

#ISO-perimetric ratio A/P
apstate["ISO_PERIMETRIC_IDX"] = Metri.dpolsby_popper(apstate)

#Equal-perimeter-circle
apstate["DAREA_EQPC"] = Metri.darea_by_ac_eqpd(apstate)


#moments for district
Ia = Metri.moment_area(apblocks)
apstate["MOMENT_AREA"][did] = Ia
Ip = Metri.moment_popu(apblocks)
apstate["MOMENT_POPU"][did] = Ip

#Mean block-ctrroid-to-district-centroid
apstate["MEAN_RADIUS_TO_DISTCTR"][did] = Metri.mean_radius_to_distctr(apblocks)
apstate["DYN_RADIUS_TO_DISTCTR"][did] = Metri.dynamic_radius_to_distctr(apblocks)
apstate["HRM_RADIUS_TO_DISTCTR"][did] = Metri.dynamic_radius_to_distctr(apblocks)

#Interpersonal-distance
apstate["INTERPERSONAL_DIST"][did] = Metri.interpersonal_distance(apblocks)

#HULL
state_hull = state.ConvexHull()
assert state_hull.GetArea() > apstate["AREA"][did]
apstate["HULL_AREA"][did] = state_hull.GetArea()
apstate["AREA_HULL_RATIO"] = Metri.convex_hull_area_ratio(apstate)
apstate["POP_HULL_RATIO"][did] = Metri.population_polygon(state_hull,apstate["POPULATION"][did], apblocks)

#Exchange-idx
eqarea_circle = Metri.comp_equal_area_circle(apstate["AREA"][did], apstate["CENTROID"][did], None)
overlap_area_eqcircle = Metri.overlaping_area(eqarea_circle, state)
apstate["OVERLAP_AREA_EQCIRCLE"][did] = overlap_area_eqcircle
apstate["EXCHG_IDX"] = Metri.exchange_index(apstate)

#Rohrbach-index
apstate["ROHRBACH_IDX"][did] = Metri.rohrbach_index(apblocks)

print("apblocks:")
print apblocks

print
print("apstate:")
print apstate

util.save(apstate,home + "cong_2011/stprop")