from simple_polygon import util
import numpy as np
from sortedcontainers import SortedDict
import time
from osgeo import ogr, osr
import os
from Metrics import Metri
import csv
from IO import *

block_prop ={}
state_prop={}
timing_attr={'tAPEC_': 0, 'tISO_PERIMETRIC_IDX_': 0, 'SOL_A': 0, 'tAPEC_BY_B': 0, 'tPOPULATION_': 0, 'tMOMENT_AREA_': 0, 'tCOLL_BLKS_BY_DI': 0, 'tMOMENT_POPU_': 0, 'tEXCHG_IDX': 0, 'tPOP_HULL_RATIO_': 0, 'tROHRBACH_IDX': 0, 'tBLKS_UCASCADE_': 0, 'tAREA_HULL_RATIO_': 0, '_': 0, 'PROB_A': 0, 'tHULL_AREA_': 0, 't3R_TO_DISTCTR_': 0, 'tBLKCTR_TO_DIST_CTR_': 0, 'tDAREA_EQPC_': 0, 'tDBLKCTR_TO_DIST_PERI_': 0}

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
         "INTERPERSONAL_DIST":{},
         "BOUNDARY_PTS":{}
         }
epsgdic = {'nad83':4269,'wgs84':4326,'pseudoutm':3857,'worldmercater':3395}

dplan_attr={"block_to_dist_map": "block_to_dist_map",
            "nblocks": "nblocks",
            "timestamp": "timestamp",
            "checksum": "checksum",
            "merkleroot": "merkleroot",
            "sol_from_acc": "solacc",
            "prob_from_acc": "problemacc",
            "planid": "planid",
            "keyschainhash": "keyschainhash",
            "dists":"dists"}

DPROBLEM_META_DICT ={
    "attr_dists": "Names/Ids of districts to be constructed.",
    "attr_ndist": "Number of districts to be constructed",
    "attr_problemacc":"Account associated with this Districting Problem",
    "attr_keychainhash":"Accumulated Hash for sorted hashes of"
                        " census blockids=(int)<STATEFP10'+'COUNTYFP10'+'TRACTCE10'+'BLOCKCE'>.",
    "attr_trans_vec":"Translation Matrix to be used.",
    "attr_scalef": "Scaling factor to be used",
    "attr_projcs":"Projection Cordinate Sytem (PCS) to be used.",
    "attr_releasetime":"Release time for this Districting Problem Map",
    "attr_closetime":"Close time for this Districting Problem Map"
    }

dpmap_attr = {"dists": "dist",
              "ndist": "ndist",
              "problemacc": "problemacc",
              "block_to_attr_map":"block_to_attr_map",
              "nblocks": "nblocks",
              "keychainhash":"keychainhash",
              "trans_vec":"translation_vector",
              "trans_vec_size":"translation_vector_size",
              "closetime":"closetime",
              "releasetime":"releasetime"
}


DEFAULT_DISTRICT_ANNOTATIONS={"01":1, "02":2, "03":3, "04":4}

DEFAULT_PROBLEM_ACC = "dplan.com"#, util.hashargs("ms.dplan.com")

DEFAULT_SOLUTION_ACC = "cse.dplan.com"#,util.hashargs("cse.dplan.com")

PLAYTIME = 3.154e+7 #one year

class DPmap(object):
    '''
    Problem definition for a redistricting problem.
    '''

    def __init__(self,ndist=None):

        self.__dists= DEFAULT_DISTRICT_ANNOTATIONS
        self.__ndist = len(DEFAULT_DISTRICT_ANNOTATIONS)
        self.__block_to_attr_map = {}
        self.__metadict = DPROBLEM_META_DICT
        self.__problemacc= None
        self.__keyschainhash =None
        self.__trans_vec=(0, 0, 0, 0, 0, 0, 0, 0, 0)
        self.__trans_vec_size = "3x3"
        self.__scalef = 0.5
        self.__projcs = "EPSG:3398"
        self.__releasetime = None
        self.__closetime = None
        self.__nblocks = None

    def __str__(self):
        return str(self.head())

    def head(self):
        dict ={}
        dict[dpmap_attr["ndist"]]= self.__ndist
        dict[dpmap_attr["dists"]] = self.__dists
        dict[dpmap_attr["problemacc"]] = self.__problemacc
        dict[dpmap_attr["keychainhash"]] = self.get_keyschainhash()
        dict[dpmap_attr["nblocks"]] = self.__nblocks
        dict[dpmap_attr["trans_vec"]] = self.__trans_vec
        dict[dpmap_attr["trans_vec_size"]] = self.__trans_vec_size
        dict[dpmap_attr["releasetime"]] = self.__releasetime
        dict[dpmap_attr["closetime"]] = self.__closetime

        return dict

    def set_dists(self,dist_ids):
        self.__dists = dist_ids

    def set_problemacc(self,prob_accid):
        if self.__problemacc is None:
            self.__problemacc = str(prob_accid)+"."+ DEFAULT_PROBLEM_ACC

    def get_dist_keys(self):
        return self.__dists.keys()


    def get_problemacc(self):
        return self.__problemacc

    def get_keyschainhash(self):
        if self.__keyschainhash == None:
            block_keys = self.get_block_to_attr_map().keys() #sorted keys
            if block_keys:
                return util.hashargs(block_keys)
        return self.__keyschainhash

    def get_block_to_attr_map(self):
        return self.__block_to_attr_map

    def to_dict(self):
        dict=self.head()
        dict[dpmap_attr["block_to_attr_map"]] = self.__block_to_attr_map
        return dict

    def to_json(self, outfile):
        dplan_data = self.to_dict()
        util.save(dplan_data, outfile)

    def load_json_dprob(self,dplanfile):

        '''
        :param jsonfile: json file with <key:h(district_id+block_id),value=(district_id,block_id)>
        :return: a dictionary that is collection of districts and blocks created.
        '''
        dprobmap = util.load_from_file(dplanfile)
        if dprobmap:
            self.__ndist = dprobmap[dpmap_attr["ndist"]]
            self.__dists = dprobmap[dpmap_attr["dists"]]
            self.__problemacc = dprobmap[dpmap_attr["problemacc"]]
            self.__keyschainhash = dprobmap[dpmap_attr["keychainhash"]]
            self.__nblocks = dprobmap[dpmap_attr["nblocks"]]
            self.__trans_vec = dprobmap[dpmap_attr["trans_vec"]]
            self.__trans_vec_size = dprobmap[dpmap_attr["trans_vec_size"]]
            self.__releasetime = dprobmap[dpmap_attr["releasetime"]]
            self.__closetime = dprobmap[dpmap_attr["closetime"]]
            self.__block_to_attr_map = dprobmap[dpmap_attr["block_to_attr_map"]]
        else:
            raise Exception("Cannot load plan data.")

    def dproblem_from_files(self,block_shapefile):
        # open shape file.

        blk_to_attr_map = SortedDict()

        shapef = ogr.Open(block_shapefile)
        layer = shapef.GetLayer()

        #print schema
        schema = []
        ldefn = layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
        print schema

        #create maping.
        sourceSpatialRef = layer.GetSpatialRef()
        Nfeats = layer.GetFeatureCount()
        ctract ={}
        for fid in xrange(Nfeats):
            fit = layer.GetFeature(fid)
            geom = fit.GetGeometryRef()
            a = fit.GetField('STATEFP10')
            b = fit.GetField('COUNTYFP10')
            c = fit.GetField('TRACTCE10')
            d = fit.GetField('BLOCKID10')
            ctract[c] = 0
            pkey = util.hash_cblocks_dist([d])
            blk_to_attr_map[pkey] = (fit.GetField('POP10'))
        print("n(census-tracts)"), len(ctract)
        #set property
        self.__block_to_attr_map = blk_to_attr_map
        self.__keyschainhash = self.get_keyschainhash()
        self.__nblocks = len(self.__block_to_attr_map.keys())
        self.__problemacc = DEFAULT_PROBLEM_ACC
        self.__releasetime = time.time()
        self.__closetime = self.__releasetime + PLAYTIME

    def blk_feat_to_dist_mapping(self,block_shapefile,blk_to_dist_map,census_tract_file=False):
        # open block shape file.
        home = os.path.dirname(block_shapefile)
        if census_tract_file:
            idfieldname  = "GEOID10"
            popfieldname = "P0010001"
        else:
            idfieldname="BLOCKID10"
            popfieldname="POP10"

        blk_to_attr_map = SortedDict()

        shapef = ogr.Open(block_shapefile)
        layer = shapef.GetLayer()

        #print schema
        ordered_dist = []
        ldefn = layer.GetLayerDefn()

        #create maping.
        sourceSpatialRef = layer.GetSpatialRef()
        Nfeats = layer.GetFeatureCount()
        ctract ={}
        for fid in xrange(Nfeats):
            fit = layer.GetFeature(fid)
            blkgeo = fit.GetGeometryRef()
            d = fit.GetField(idfieldname) #id
            pop = fit.GetField(popfieldname) #population
            blockkey = util.hash_cblocks_dist([d])
            #get dist value for this 'd'
            ordered_dist += [blk_to_dist_map[blockkey][1:]]
        #set property

        #save as csv.
        print ordered_dist[0]
        with open(home + '/ordered_blk_eqiv.csv', 'wb') as writeFile:
            writer = csv.writer(writeFile)
            writer.writerows(ordered_dist)
        writeFile.close()

    def dproblem_from_shape_files(self,block_shapefile,census_tract_file=False):
        # open block shape file.

        if census_tract_file:
            idfieldname  = "GEOID10"
            popfieldname = "P0010001"
        else:
            idfieldname="BLOCKID10"
            popfieldname="POP10"

        blk_to_attr_map = SortedDict()

        shapef = ogr.Open(block_shapefile)
        layer = shapef.GetLayer()

        #print schema
        schema = []
        ldefn = layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
        print schema

        #create maping.
        sourceSpatialRef = layer.GetSpatialRef()
        Nfeats = layer.GetFeatureCount()
        ctract ={}
        for fid in xrange(Nfeats):
            fit = layer.GetFeature(fid)
            blkgeo = fit.GetGeometryRef()
            d = fit.GetField(idfieldname) #id
            pop = fit.GetField(popfieldname) #population
            blockkey = util.hash_cblocks_dist([d])
            # block attrs
            blk_to_attr_map[blockkey] = {"FEAT_KEY":d,"POPULATION":pop}
            blk_to_attr_map[blockkey]["AREA"] = blkgeo.GetArea()
            blk_to_attr_map[blockkey]["PERI"] = blkgeo.GetBoundary().Length()
            blk_to_attr_map[blockkey]["ENV"] = blkgeo.GetEnvelope()
            blk_to_attr_map[blockkey]["CENTROID"] = blkgeo.Centroid().GetPoints()[0]

        #set property
        self.__block_to_attr_map = blk_to_attr_map
        self.__keyschainhash = self.get_keyschainhash()
        self.__nblocks = len(self.__block_to_attr_map.keys())
        if self.__problemacc is None:
            self.__problemacc= DEFAULT_PROBLEM_ACC
        self.__releasetime = time.time()
        self.__closetime = self.__releasetime + PLAYTIME

    @classmethod
    def geom_by_fid(self,blockids,block_shapefile):
        # open shape file.
        shapef = ogr.Open(block_shapefile)
        layer = shapef.GetLayer()
        for feature in layer:
            if [int(feature.GetField('STATEFP10')),
                int(feature.GetField('COUNTYFP10')),
                int(feature.GetField('TRACTCE10')),
                int(feature.GetField('BLOCKCE'))] == blockids:

                geom = feature.GetGeometryRef()
                print geom.GetGeometryName()
                return geom

    @classmethod
    def test_dpmap(self):
        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
        block_shapefile = home+"census_blocks_ms_2015/tabblock2010_28_pophu.shp"

        dpmap = DPmap()
        dpmap.dproblem_from_files(block_shapefile)
        print("dpmap")
        print dpmap.head()

        dpmap.to_json(home+"census_blocks_ms_2015/dpmap.json")

        return dpmap

class Dplan(object):

    def __init__(self):

        self.__problemacc = None
        self.__solacc = None #solution acc
        self.__block_to_dist_map = {} #block eqiv dict.<blockid,distid>
        self.__merkleroot = None
        self.__checksum = None
        self.__timestamp = None
        self.__nblocks = None
        self.__planid = None
        self.__keyschainhash=None
        self.__dists = {} #list of districts


    def get_dists(self):
        return self.__dists

    def get_planid(self):

        if self.__planid == None:
            return util.hashargs(self.get_problemacc(),self.get_solacc(),self.get_merkleroot())
        else:
            return self.__planid

    def get_keyschainhash(self):
        if self.__keyschainhash == None:
            block_keys = self.get_block_to_attr_map().keys() #sorted keys
            if block_keys:
                return util.hashargs(block_keys)
        return self.__keyschainhash

    def set_problemacc(self,problemacc):
        self.__problemacc = problemacc

    def get_problemacc(self):
        return self.__problemacc

    def get_solacc(self):
        return self.__solacc

    def get_block_to_attr_map(self):
        return self.__block_to_dist_map

    def get_timestamp(self):
        return self.__timestamp

    def get_nblocks(self):
        return self.__nblocks

    def get_checksum(self):
        return self.__checksum

    def get_merkleroot(self):
        return self.__merkleroot

    def get_dists(self):
        return self.__dists

    def reset(self):
        self.__problemacc = None
        self.__solacc = None
        self.__block_to_dist_map = {} #block eqiv dict.
        self.__merkleroot = None
        self.__checksum = None
        self.__timestamp = None
        self.__nblocks = None

    def header(self):
        dplan_data ={}
        dplan_data[dplan_attr['nblocks']] = self.get_nblocks()
        dplan_data[dplan_attr['timestamp']] = self.get_timestamp()
        dplan_data[dplan_attr['merkleroot']] = self.get_merkleroot()
        dplan_data[dplan_attr['prob_from_acc']] = self.get_problemacc()
        dplan_data[dplan_attr['checksum']] = self.get_checksum()
        dplan_data[dplan_attr['sol_from_acc']] = self.get_solacc()
        dplan_data[dplan_attr['planid']] = self.get_planid()
        dplan_data[dplan_attr['keyschainhash']] = self.get_keyschainhash()

        return dplan_data

    def load_json_dplan(self,dplanfile):
        '''

        :param jsonfile: json file with <key:h(district_id+block_id),value=(district_id,block_id)>
        :return: a dictionary that is collection of districts and blocks created.
        '''
        plandata = util.load_from_file(dplanfile)
        if plandata:
            self.__block_to_dist_map = plandata[dplan_attr["block_to_dist_map"]]
            self.__nblocks = plandata[dplan_attr["nblocks"]]
            self.__timestamp = plandata[dplan_attr["timestamp"]]
            self.__solacc = plandata[dplan_attr["sol_from_acc"]]
            self.__problemacc = plandata[dplan_attr["prob_from_acc"]]
            self.__checksum = plandata[dplan_attr["checksum"]]
            self.__merkleroot = plandata[dplan_attr["merkleroot"]]
            self.__keyschainhash = plandata[dplan_attr["keyschainhash"]]
        else:
            raise Exception("Cannot load plan data.")

    def build_from_blkequiv_csv(self, blk_equivfile):

        dists =SortedDict()
        dict = SortedDict()

        with open(blk_equivfile) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')

            for blkvals in csv_reader:
                pkey = util.hash_cblocks_dist([blkvals[0]]) #key:blockid
                dict[pkey] = tuple(blkvals)
                dists[blkvals[-1]] = 0

        #
        self.__block_to_dist_map = dict
        self.__nblocks = len(list(dict.keys()))
        self.__merkleroot = util.merkle_root(self.__block_to_dist_map)
        self.__problemacc = DEFAULT_PROBLEM_ACC
        self.__solacc = DEFAULT_SOLUTION_ACC
        self.__planid = self.get_planid()
        self.__keyschainhash = self.get_keyschainhash()
        self.__timestamp = time.time()
        self.__dists = dists

    def build_from_blkequiv_csvbak(self, blk_equivfile):

        arr = np.loadtxt(blk_equivfile,delimiter=',')
        print arr.shape
        # if <blockid is a not single string.
        if arr.shape[1] >2:
            arr = arr.astype(str)

        dists =SortedDict()
        dict = SortedDict()
        if arr.shape[1] > 2:
            for blkvals in arr:
                pkey = util.hash_cblocks_dist(tuple(blkvals[0:-1])) #key:blockid
                dict[pkey] = tuple(blkvals)
                dists[tuple(blkvals)[-1]] = 0

        #
        self.__block_to_dist_map = dict
        self.__nblocks = len(list(dict.keys()))
        self.__merkleroot = util.merkle_root(self.__block_to_dist_map)
        self.__problemacc = DEFAULT_PROBLEM_ACC
        self.__solacc = DEFAULT_SOLUTION_ACC
        self.__planid = self.get_planid()
        self.__keyschainhash = self.get_keyschainhash()
        self.__timestamp = time.time()
        self.__dists = dists

    def to_json(self, outfile):

        dplan_data ={}
        dplan_data[dplan_attr['block_to_dist_map']] = self.get_block_to_attr_map()
        dplan_data[dplan_attr['nblocks']] = self.get_nblocks()
        dplan_data[dplan_attr['timestamp']] = time.time()
        dplan_data[dplan_attr['merkleroot']] = self.get_merkleroot()
        dplan_data[dplan_attr['prob_from_acc']] = self.get_problemacc()
        dplan_data[dplan_attr['checksum']] = self.get_checksum()
        dplan_data[dplan_attr['sol_from_acc']] = self.get_solacc()
        dplan_data[dplan_attr['planid']] = self.get_planid()
        dplan_data[dplan_attr['keyschainhash']] = self.get_keyschainhash()
        #dplan_data[dplan_attr['dists']] = self.get_dists()
        util.save(dplan_data, outfile)

    def validate_dplan(self):

        cond = False
        #check if solacc exists
        #check if merkleroot correct
        #check if probacc exists
        #check if all blocks in problem defn are assigned at exactly once to a dist
        #check if required no of district in problem defn is created
        #check if signature for this plan is valid.
        #check if planid(problemid+solby+merkleroot) is not already processed i.e it is not in MetricChain

        return cond

    @classmethod
    def test_dpaln(self):

        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
        #blk_equivfile = "congress_blk_equiv_coma.txt" #with 5 cols
        #blk_equivfile = "congress_blk_equiv_coma.txt"  # with 5 cols
        blk_equivfile = "tr_equiv.csv" #with 2 cols #tract equiv file

        dplan = Dplan()
        dplan.build_from_blkequiv_csv(home +blk_equivfile)
        print dplan.header()

        dplan.to_json(home+"dplan")
        return dplan

    @classmethod
    def test_load_dpaln(self):
        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
        dplan = Dplan()
        print("Load dplan from json.")
        dplan.load_json_dplan(home+"dplan.json")

        print dplan.header()
        return dplan

    @classmethod
    def test_metric2(self):
        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/cd116/"
        blk_equivfile = "blk_eqiv28.csv" #with 5 cols

        dplan = Dplan()
        dplan.build_from_blkequiv_csv(home +blk_equivfile)
        print dplan.header()

        dplan.to_json(home+"dplan")
        return dplan

    def block_id_components(self,strblockid):
        '''

        :param strblockid: example:28 121 020402 3035 returns
        :return:

        >>STATEFP10,C,2	; COUNTYFP10,C,3; TRACTCE10,C,6 ; BLOCKCE,C,4
        >> 28	001	000900	2042

        '''
        stid = strblockid[0:2]
        cid  = strblockid[2:5]
        tid  = strblockid[5:11]
        bid =  strblockid[11:]

        return (stid,cid,tid,bid)

    @classmethod
    def test_metric_ctract_integration(self):
        timing={}
        state_prop = {}

        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
        ctract_shapefile = home+"Tracts10PopnHou28/Tracts10PopnHou.shp"
        blk_equivfile = "tr_equiv28.csv" #with 2 cols #tract equiv file

        #save and get a plan
        dplan = Dplan()
        dplan.build_from_blkequiv_csv(home +blk_equivfile)
        dplan.to_json(home+"ctdplan28")
        del dplan

        dplan = Dplan()
        dplan.load_json_dplan(home+"ctdplan28.json")
        IO().insert_plan(dplan)
        timing["SOL_ACC"] =dplan.get_solacc()

        #save and get a dist problem
        t0 = time.clock()
        dpmap = DPmap()
        dpmap.dproblem_from_shape_files(ctract_shapefile,census_tract_file=True)
        dpmap.to_json(home+"ctdpmap28.json")
        del dpmap

        dpmap = DPmap()
        dpmap.load_json_dprob(home + "ctdpmap28.json")
        timing["tAPEC_BY_BLK"] = round(time.clock() - t0,6)
        timing["PROB_ACC"] = dpmap.get_problemacc()

        #for which prob, which solution
        state_prop["PROB_ACC"] = dpmap.get_problemacc()
        state_prop["SOL_ACC"] = dplan.get_solacc()
        state_prop["DIST_KEYS"] = dpmap.get_dist_keys()

        # fill up block's prop.
        block_prop = dpmap.get_block_to_attr_map() #transient properties
        util.save(block_prop, home + "blkprop28-")
        IO().insert_block_prop(dpmap)
        del dpmap

        t0 = time.clock()
        # collection of blocks by dist.
        districts_blocks = {}
        for ctractkey, idvalues in dplan.get_block_to_attr_map().items():
            #TODO:
            #this blockkey must be in dpmap.
            distid = idvalues[-1]
            try:
                districts_blocks[distid] += [ctractkey]
            except:
                districts_blocks[distid] = [ctractkey]

        del dplan
        timing["tCOLL_BLKS_BY_DIST"] = round(time.clock() - t0,6)
        ###

        #Validate:
        ##dist_keys in plan must be in problem-maps's dist keys.
        ##all of the input blocks must be mapped to one of the district.

        for distkey  in sorted(state_prop["DIST_KEYS"])[:]:
            t0 = time.clock()
            state_prop[distkey] ={}
            ctractkeys = districts_blocks[distkey][:]
            print("# of blocks in this dist:"),len(ctractkeys)
            udist_geom = ogr.Geometry(ogr.wkbMultiPolygon)

            #get each block's shape to create a district's map.
            fcount = 0
            shapef = ogr.Open(ctract_shapefile)
            layer = shapef.GetLayer()
            #
            for feature in layer:
                bkey = util.hash_cblocks_dist([feature.GetField('GEOID10')])
                blkgeo = feature.GetGeometryRef()
                if bkey in ctractkeys:
                    if blkgeo.GetGeometryName() == "LINEARRING":
                        print("\t-.-"), "linearring"
                        poly = ogr.Geometry(ogr.wkbPolygon)
                        poly.AddGeometry(blkgeo)

                    if blkgeo.GetGeometryName() == 'MULTIPOLYGON':
                        print("\t-.-"), "multipolygon"
                        for geom in blkgeo:
                            udist_geom.AddGeometry(geom)

                    elif blkgeo.GetGeometryName() == 'POLYGON':
                        udist_geom.AddGeometry(blkgeo)
                    else:
                        print("\t\t-.-"),blkgeo.GetGeometryName()
                        udist_geom.AddGeometry(blkgeo)

            #end-for
            layer = None

            #take union of blocks.
            print("--"),udist_geom.GetGeometryName(), udist_geom.GetGeometryCount(), fcount
            udist_geom = udist_geom.UnionCascaded()
            count = 0
            #for each ring geom
            for geom in udist_geom:
                print("\t"),count,geom.GetGeometryName()
                if geom.GetGeometryName() == 'LINEARRING':
                    state_prop[distkey]["BOUNDARY_PTS"] =geom.GetPoints()
                    #poly = ogr.Geometry(ogr.wkbPolygon)
                    #poly.AddGeometry(geom)
                    #util.geom_toshp(poly, home + "dists/st-" + str(distkey) + "-" + str(count), save_as_multipt=False)
                else:
                    #util.geom_toshp(geom,home + "dists/st-"+str(distkey)+"-"+str(count), save_as_multipt=False)
                    pass
                count +=1
            print("Completed union of blks for a district")
            timing["tUNION_CASCADE_BLOCKS_BY_DIST_"+str(distkey)] = round(time.clock() - t0,6)

        util.save(state_prop, home + "temp_distprop28-")
        print("Completed all district construction")

        #compute state's APEC properties
        t0 = time.clock()
        for distkey in state_prop["DIST_KEYS"]:
            polypts = state_prop[distkey]["BOUNDARY_PTS"]
            poly = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in polypts:
                x, y = float(x), float(y)
                ring.AddPoint(x, y)
            poly.AddGeometry(ring)

            state_prop[distkey]["AREA"]= poly.GetArea()
            state_prop[distkey]["PERI"] = poly.GetBoundary().Length()
            state_prop[distkey]["ENV"] = poly.GetEnvelope()
            state_prop[distkey]["CENTROID"] = poly.Centroid().GetPoints()[0]
            timing["tAPEC_BY_DIST"+str(distkey)] = round(time.clock() - t0,6)

        #for each block in district, update block's properties relative to comprising district.
        for distkey in state_prop["DIST_KEYS"]:
            st_centroid=state_prop[distkey]["CENTROID"]

            t0 = time.clock()
            Metri.bctr_to_dctr(st_centroid,districts_blocks[distkey],block_prop)
            timing["tBLKCTR_TO_DIST_CTR_"+distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            Metri.blkctr_to_dist_peri(state_prop[distkey]["BOUNDARY_PTS"],
                                      districts_blocks[distkey],
                                      block_prop)
            timing["tDBLKCTR_TO_DIST_PERI_"+distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            #compute distr population
            state_prop[distkey]["POPULATION"]=Metri.dist_population(districts_blocks[distkey], block_prop)
            timing["tPOPULATION_" + distkey] = round(time.clock() - t0,6)


            t0 = time.clock()
            # moments for district
            Ia = Metri.moment_area(districts_blocks[distkey],block_prop)
            state_prop[distkey]["MOMENT_AREA"] = Ia
            timing["tMOMENT_AREA_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            Ip = Metri.moment_popu(districts_blocks[distkey],block_prop)
            state_prop[distkey]["MOMENT_POPU"] = Ip
            timing["tMOMENT_POPU_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # ISO-perimetric ratio A/P
            Metri.dpolsby_popper2(distkey,state_prop)
            timing["tISO_PERIMETRIC_IDX_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # Equal-perimeter-circle
            Metri.darea_by_ac_eqpd2(distkey,state_prop)
            timing["tDAREA_EQPC_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # Mean block-ctrroid-to-district-centroid
            state_prop[distkey]["MEAN_RADIUS_TO_DISTCTR"] = Metri.mean_radius_to_distctr(districts_blocks[distkey],
                                                                                  block_prop)
            state_prop[distkey]["DYN_RADIUS_TO_DISTCTR"]= Metri.dynamic_radius_to_distctr(districts_blocks[distkey],
                                                                                    block_prop)
            state_prop[distkey]["HRM_RADIUS_TO_DISTCTR"] = Metri.harmonic_radius_to_distctr(districts_blocks[distkey],
                                                                                    block_prop)
            timing["t3R_TO_DISTCTR_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # Interpersonal-distance
            state_prop[distkey]["INTERPERSONAL_DIST"] = Metri.interpersonal_distance(districts_blocks[distkey],
                                                                                     block_prop)
            timing["tINTERPERSONAL_DIST_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # HULL
            state_hull = Metri.get_convexhull(state_prop[distkey]["BOUNDARY_PTS"])
            state_prop[distkey]["HULL_AREA"]= state_hull.GetArea()
            timing["tHULL_AREA_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()

            state_prop[distkey]["AREA_HULL_RATIO"] = Metri.convex_hull_area_ratio(distkey,
                                                                                  state_prop)
            timing["tAREA_HULL_RATIO_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            state_prop[distkey]["POP_HULL_RATIO"]= Metri.population_polygon(state_hull,
                                                                            state_prop[distkey]["POPULATION"],
                                                                            districts_blocks[distkey],
                                                                            block_prop)

            timing["tPOP_HULL_RATIO_" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # Exchange-idx
            eqarea_circle = Metri.comp_equal_area_circle(state_prop[distkey]["AREA"], state_prop[distkey]["CENTROID"], None)
            overlap_area_eqcircle = Metri.overlaping_area(eqarea_circle, state_prop[distkey]["BOUNDARY_PTS"])
            state_prop[distkey]["OVERLAP_AREA_EQCIRCLE"] = overlap_area_eqcircle
            state_prop[distkey]["EXCHG_IDX"] = Metri.exchange_index(distkey,state_prop)
            timing["tEXCHG_IDX" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            # Rohrbach-index
            state_prop[distkey]["ROHRBACH_IDX"]= Metri.rohrbach_index(districts_blocks[distkey], block_prop)
            timing["tROHRBACH_IDX" + distkey] = round(time.clock() - t0,6)

            t0 = time.clock()
            print("Completed prop for "), distkey
        #update

        IO().insert_timings(timing)
        del timing
        IO().insert_dist_metric(state_prop)
        del state_prop

    def create_problem_solution_object(self,state_id="28"):

        data_info ={}
        home_cong= "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/cd116/"
        home_blkpop = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/census_blocks_by_states/"+"tabblock2010_"+state_id+"_pophu/"
        ctract_shapefile = home_blkpop+"prjtabblock2010_"+state_id+"_pophu.shp"
        blk_equivfile = home_cong+"blk_eqiv"+state_id+".csv" #with 2 cols #tract equiv file
        dist_idsfile = home_cong+"dist_ids"+state_id+".csv" #with 2 cols #tract equiv file
        nbeqivf = home_cong + "National_CD116.csv"
        data_info["state_id"] = state_id

        census_tract_file = False
        if census_tract_file:
            idfieldname  = "GEOID10"
            dist_id_field = "STID3" #dist-id
        else:
            idfieldname="BLOCKID10"
            dist_id_field="STID3"

        prefix="blk" #integration-type


        dplanproblemfile= prefix+"dpmap"+state_id
        dplandatafile = prefix + "dplan" + state_id
        data_info["size_shp"] = util.get_size(os.path.dirname(ctract_shapefile))
        data_info["dists"] ={}

        #create dist-ids for this state.
        dist_ids = util.extract_dist_ids_for_state(nbeqivf,state_id=state_id)

        #save and get a dist problem
        t0 = time.clock()
        dpmap = DPmap()
        dpmap.set_problemacc(state_id)
        dpmap.set_dists(dist_ids) #or (util.load_dist_ids_by_state(dist_idsfile))
        dpmap.dproblem_from_shape_files(ctract_shapefile,census_tract_file=False)
        dpmap.to_json(home_blkpop+dplanproblemfile)
        problemacc = dpmap.get_problemacc()
        data_info["storsize_dprob"] = util.get_size(home_blkpop+dplanproblemfile+'.json')
        del dpmap

        # district plan
        util.extract_blk_eq_by_state(nbeqivf,state_id=state_id) #national block equiv file.

        dplan = Dplan()
        dplan.build_from_blkequiv_csv(blk_equivfile)
        dplan.set_problemacc(problemacc)
        dplan.to_json(home_cong+dplandatafile)
        #IO().insert_plan(dplan)
        data_info["storsize_dplan"] = util.get_size(home_cong+dplandatafile+'.json')
        IO().insert_data_info(data_info)
        del dplan

    def test_metric_block_integration(self,state_id="28",filter_dists=[]):
        timing={}
        state_prop = {}
        data_info ={}
        home_cong= "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/cd116/"
        home_blkpop = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/census_blocks_by_states/"+"tabblock2010_"+state_id+"_pophu/"
        ctract_shapefile = home_blkpop+"prjtabblock2010_"+state_id+"_pophu.shp"
        blk_equivfile = home_cong+"blk_eqiv"+state_id+".csv" #with 2 cols #tract equiv file
        dist_idsfile = home_cong+"dist_ids"+state_id+".csv" #with 2 cols #tract equiv file

        census_tract_file = False
        if census_tract_file:
            idfieldname  = "GEOID10"
            dist_id_field = "STID3" #dist-id
        else:
            idfieldname="BLOCKID10"
            dist_id_field="STID3"

        prefix="blk" #integration-type

        dplandatafile =prefix+"dplan"+state_id
        dplanproblemfile= prefix+"dpmap"+state_id
        dplan_block_propfile= prefix+"dpblock_prop"+state_id
        data_info["size_shp"] = util.get_size(os.path.dirname(ctract_shapefile))
        data_info["dists"] ={}

        #get a dist problem
        t0 = time.clock()
        dpmap = DPmap()
        dpmap.load_json_dprob(home_blkpop + dplanproblemfile)
        timing["tAPEC_BY_BLK"] = round(time.clock() - t0,6)
        timing["PROB_ACC"] = dpmap.get_problemacc()

        #get a plan
        dplan = Dplan()
        dplan.load_json_dplan(home_cong+dplandatafile)
        timing["SOL_ACC"] =dplan.get_solacc()

        #for which prob, which solution
        state_prop["PROB_ACC"] = dpmap.get_problemacc()
        state_prop["SOL_ACC"] = dplan.get_solacc()
        state_prop["DIST_KEYS"] = sorted(dpmap.get_dist_keys())

        #filter by district ids.
        for k in filter_dists:
            try:state_prop["DIST_KEYS"].remove(k)
            except:
                continue

        # fill up block's prop.
        block_prop = dpmap.get_block_to_attr_map() #transient properties
        #IO().insert_block_prop(dpmap) #TODO: isssue with unicode keys in dict
        del dpmap

        t0 = time.clock()
        # collection of blocks by dist.
        districts_blocks = {}
        for ctractkey, idvalues in dplan.get_block_to_attr_map().items():
            #TODO:
            #this blockkey must be in dpmap.
            distid = idvalues[-1]
            try:
                districts_blocks[distid] += [ctractkey]
            except:
                districts_blocks[distid] = [ctractkey]

        timing["tCOLL_BLKS_BY_DIST"] = round(time.clock() - t0,6)
        print "tCOLL_BLKS_BY_DIST",timing["tCOLL_BLKS_BY_DIST"]
        ###

        #Validate:
        ##dist_keys in plan must be in problem-maps's dist keys.
        ##all of the input blocks must be mapped to one of the district.
        print("# of Districts :"),len(state_prop["DIST_KEYS"])
        for distkey  in state_prop["DIST_KEYS"]:

            t0 = time.clock()
            state_prop[distkey] ={}
            print("# of blocks in this dist:"),distkey, len(districts_blocks[distkey])
            udist_geom = ogr.Geometry(ogr.wkbMultiPolygon)
            cd_container = {distkey:udist_geom}
            #get each block's shape to create a district's map.

            fcount = 0
            shapef = ogr.Open(ctract_shapefile)
            layer = shapef.GetLayer()
            #
            t0 = time.clock()
            t2 = time.time()
            for feature in layer:
                bkey = util.hash_cblocks_dist([feature.GetField(idfieldname)])
                blkgeo = feature.GetGeometryRef()
                #distid = feature.GetField(dist_id_field) #or
                distid = dplan.get_block_to_attr_map()[bkey][1]
                try:
                    cd_container[distid].AddGeometry(blkgeo)
                except:
                    continue

            #end-for
            layer = None
            print("time:"), time.clock() -t0, time.time()-t2
            print("-udist_geom count-"), udist_geom.GetGeometryCount()

            #take union of blocks.
            t0 = time.clock()
            udist_geomc = udist_geom.UnionCascaded()
            timing["tBLKS_UCASCADE_" + str(distkey)] = round(time.clock() - t0, 6)

            count = 0
            #for each ring geom, search for a largest ring.
            for geom in udist_geomc:
                if geom.GetGeometryName() == 'LINEARRING':
                    state_prop[distkey]["BOUNDARY_PTS"] = geom.GetPoints()
                    poly = ogr.Geometry(ogr.wkbPolygon)
                    poly.AddGeometry(geom)
                    util.geom_toshp(poly, home_cong + "dists/"+str(state_id)+"-"+"dist-" + str(distkey) + "-" + str(count), save_as_multipt=False)
                    break
                elif geom.GetGeometryName() == "POLYGON":
                    for geo in geom:
                        state_prop[distkey]["BOUNDARY_PTS"] = geo.GetPoints()
                        poly = ogr.Geometry(ogr.wkbPolygon)
                        poly.AddGeometry(geo)
                        util.geom_toshp(poly,home_cong + "dists/" + str(state_id) + "-" + "dist-" + str(distkey) + "-" + str(count), save_as_multipt=False)
                        break
                count += 1
            print("\t\t dist-boundary completed"),len(state_prop[distkey]["BOUNDARY_PTS"])
            #TODO: cascadedunion of these linearring must be performed.
        ##
        #util.save(state_prop, home_blkpop + "temp_distprop28-")
        print("Completed block integration to construct a district")

        #compute state's APEC properties
        t0 = time.clock()
        for distkey in state_prop["DIST_KEYS"][:]:

            polypts = state_prop[distkey]["BOUNDARY_PTS"]
            poly = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in polypts:
                x, y = float(x), float(y)
                ring.AddPoint(x, y)
            poly.AddGeometry(ring)

            state_prop[distkey]["AREA"]= poly.GetArea()
            state_prop[distkey]["PERI"] = poly.GetBoundary().Length()
            state_prop[distkey]["ENV"] = poly.GetEnvelope()
            state_prop[distkey]["CENTROID"] = poly.Centroid().GetPoints()[0]
            timing["tAPEC_"+str(distkey)] = round(time.clock() - t0,6)
            state_prop[distkey]["FAILED"] =[]
            print("Completed APEC for district "),distkey

        #for each block in district, update block's properties relative to comprising district.
        for distkey in state_prop["DIST_KEYS"][:]:
            st_centroid=state_prop[distkey]["CENTROID"]
            data_info["dists"]["#"+prefix+"-"+distkey]=len(districts_blocks[distkey])
            t0 = time.clock()
            try:
                Metri.bctr_to_dctr(st_centroid,districts_blocks[distkey],block_prop)
                timing["tBLKCTR_TO_DIST_CTR_"+distkey] = round(time.clock() - t0,6)
            except Exception, e:
                state_prop[distkey]["FAILED"] = ["BLKCTR_TO_DIST_CTR_"]
                print("Exception:"), str(e)
            print("1")

            t0 = time.clock()
            try:
                Metri.blkctr_to_dist_peri_np(state_prop[distkey]["BOUNDARY_PTS"],districts_blocks[distkey],block_prop)
                timing["tDBLKCTR_TO_DIST_PERI_"+distkey] = round(time.clock() - t0,6)
                pass
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["DBLKCTR_TO_DIST_PERI_"]
                print("Exception--"), str(e)
            print timing["tDBLKCTR_TO_DIST_PERI_" + distkey]
            print("2")

            t0 = time.clock()
            #compute distr population
            state_prop[distkey]["POPULATION"]=Metri.dist_population(districts_blocks[distkey], block_prop)
            timing["tPOPULATION_" + distkey] = round(time.clock() - t0,6)
            print("3")

            t0 = time.clock()
            # moments for district
            try:
                Ia = Metri.moment_area(districts_blocks[distkey],block_prop)
                state_prop[distkey]["MOMENT_AREA"] = Ia
                timing["tMOMENT_AREA_" + distkey] = round(time.clock() - t0,6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["MOMENT_AREA"]
                print("Exception"), str(e)

            print("4")


            t0 = time.clock()
            try:
                Ip = Metri.moment_popu(districts_blocks[distkey],block_prop)
                state_prop[distkey]["MOMENT_POPU"] = Ip
                timing["tMOMENT_POPU_" + distkey] = round(time.clock() - t0,6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["MOMENT_POPU"]
                print("Exception"), str(e)

            print("5")

            t0 = time.clock()
            # ISO-perimetric ratio A/P
            try:
                Metri.dpolsby_popper2(distkey,state_prop)
                timing["tISO_PERIMETRIC_IDX_" + distkey] = round(time.clock() - t0,6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["ISO_PERIMETRIC_IDX_"]
                print("Exception"), str(e)

            print("6")

            t0 = time.clock()
            # Equal-perimeter-circle
            Metri.darea_by_ac_eqpd2(distkey,state_prop)
            timing["tDAREA_EQPC_" + distkey] = round(time.clock() - t0,6)

            # Mean block-ctrroid-to-district-centroid
            t0 = time.clock()
            try:
                state_prop[distkey]["MEAN_RADIUS_TO_DISTCTR"] = Metri.mean_radius_to_distctr(districts_blocks[distkey],
                                                                        block_prop)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["MEAN_RADIUS_TO_DISTCTR"]
                print("Exception"), str(e)
            print("7")
            try:
                state_prop[distkey]["DYN_RADIUS_TO_DISTCTR"]= Metri.dynamic_radius_to_distctr(districts_blocks[distkey],
                                                                                    block_prop)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["DYN_RADIUS_TO_DISTCTR"]
                print("Exception"), str(e)
            print("8")
            try:
                state_prop[distkey]["HRM_RADIUS_TO_DISTCTR"] = Metri.harmonic_radius_to_distctr(districts_blocks[distkey],
                                                                                    block_prop)
                timing["t3R_TO_DISTCTR_" + distkey] = round(time.clock() - t0, 6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["HRM_RADIUS_TO_DISTCTR"]
                print("Exception"), str(e)
            print("9")

            # Interpersonal-distance
            t0 = time.clock()
            try:
                #state_prop[distkey]["INTERPERSONAL_DIST"] = Metri.interpersonal_distance(districts_blocks[distkey],
                                                                            # block_prop)
                #timing["tINTERPERSONAL_DIST_" + distkey] = round(time.clock() - t0,6)
                pass
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["INTERPERSONAL_DIST"]
                print("Exception"), str(e)
            print("10")

            t0 = time.clock()
            # HULL
            state_hull = Metri.get_convexhull(state_prop[distkey]["BOUNDARY_PTS"])
            state_prop[distkey]["HULL_AREA"]= state_hull.GetArea()
            timing["tHULL_AREA_" + distkey] = round(time.clock() - t0,6)
            print("11")

            t0 = time.clock()
            try:
                state_prop[distkey]["AREA_HULL_RATIO"] = Metri.convex_hull_area_ratio(distkey,
                                                                                  state_prop)
                timing["tAREA_HULL_RATIO_" + distkey] = round(time.clock() - t0, 6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["AREA_HULL_RATIO"]
                print(str(e))

            t0 = time.clock()
            try:
                state_prop[distkey]["POP_HULL_RATIO"]= Metri.population_polygon(state_hull,
                                                                            state_prop[distkey]["POPULATION"],
                                                                            districts_blocks[distkey],
                                                                            block_prop)
                timing["tPOP_HULL_RATIO_" + distkey] = round(time.clock() - t0, 6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["POP_HULL_RATIO"]
                print("Exception"), str(e)
            print("13")


            t0 = time.clock()
            try:
                #Exchange-idx
                eqarea_circle = Metri.comp_equal_area_circle(state_prop[distkey]["AREA"], state_prop[distkey]["CENTROID"], None)
                overlap_area_eqcircle = Metri.overlaping_area(eqarea_circle, state_prop[distkey]["BOUNDARY_PTS"])
                state_prop[distkey]["OVERLAP_AREA_EQCIRCLE"] = overlap_area_eqcircle
                state_prop[distkey]["EXCHG_IDX"] = Metri.exchange_index(distkey,state_prop)
                timing["tEXCHG_IDX_" + distkey] = round(time.clock() - t0, 6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["EXCHG_IDX"]
                print("Exception"), str(e)
            print("14")

            t0 = time.clock()
            try:
                #Rohrbach-index
                state_prop[distkey]["ROHRBACH_IDX"]= Metri.rohrbach_index(districts_blocks[distkey], block_prop)
                timing["tROHRBACH_IDX_" + distkey] = round(time.clock() - t0, 6)
            except Exception, e:
                state_prop[distkey]["FAILED"] += ["ROHRBACH_IDX"]
                print("Exception"), str(e)

            print("15")
            print("Completed prop for "), distkey
        #update

        IO().insert_timings(timing)
        del timing
        IO().insert_dist_metric(state_prop)
        #del state_prop
        data_info["integration_type"] = prefix
        data_info["state_id"] = state_id
        IO().insert_data_info(data_info)
        print("Completed Metric for State:"),state_id
        print
        print

import csv
states_code = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/states-codes.csv"
completed_states="D:/workspace/sqdm-repo/sqdm/out/tmp/redist/completed_states.csv"

state_dic ={}
with open(states_code) as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        state_dic[str(row[1])] = [row[0],row[2]]



filter_states ={"00":'memory error'}
for stid in filter_states.keys():
    try:
        state_dic.pop(stid) #del state_dic[stid]
    except:
        continue



#load record for completed states.
completed_dic ={}
with open(completed_states) as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        completed_dic[str(row[0])] = row[1:]

print("Cd"), completed_dic
'''
#do for ms and al
for state_id in sorted(state_dic.keys()):
    #create problem
    #create solution
    if state_id not in completed_dic:
        print("state-id"), state_id
        try:
            Dplan().create_problem_solution_object(state_id=state_id)
            blk_dist =[[state_id,"P","S"]]
            with open(completed_states, 'ab') as writeFile:
                writer = csv.writer(writeFile)
                writer.writerows(blk_dist)
            writeFile.close()
        except Exception, e:
            print("Exception "),str(e), state_id
'''

f = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/census_blocks_by_states/tabblock2010_06_pophu/blkdpmap06.json"
import json
from collections import OrderedDict
import json

print("Calculating Metric")
#metric computation for states
for state_id in sorted(state_dic.keys()):
    #metric computation
    if state_id in completed_dic and len(completed_dic[state_id])<3:
        print("state-id"), state_id,state_dic[state_id]
        try:
            if state_id == "06":
                filter_dists = ["01"]
            else:
                filter_dists = ["ZZ"]

            Dplan().test_metric_block_integration(state_id=state_id,filter_dists=filter_dists)
            blk_dist =[[state_id,"P","S","M"]]
            with open(completed_states, 'ab') as writeFile:
                writer = csv.writer(writeFile)
                writer.writerows(blk_dist)
            writeFile.close()
        except Exception, e:
            print("Exception "),str(e), state_id



















