import collections
import sys

class IO:

    def __init__(self):
        from pymongo import MongoClient
        client = MongoClient()
        self.db = client.local
        self.cdplan = self.db["dplan"]
        self.cblock_prop = self.db["block_prop"]
        self.ctiming=self.db['timing']
        self.cdistmetric = self.db["dist_metric"]

    def insert_plan(self,plan):
        d = dict(plan.get_block_to_attr_map())
        key={"_id":plan.get_solacc()}
        try:
            #self.cdplan.insert_one(d).inserted_id
            self.cdplan.update(key,d, upsert=True)
        except:
            print("Exists plan in DB.")

    def insert_block_prop(self,block_prop={}):
        '''
            :param block_prop:
            :return:
        '''
        d = dict(block_prop.get_block_to_attr_map())
        try:
            for k,dict_val in d.items()[0:]:

                key = {"_id": block_prop.head()["problemacc"]+k}
                self.cblock_prop.update(key,self.convert(dict(dict_val)), upsert=True)


        except Exception as e:
            print("Exists plan in DB."), sys.getsizeof(block_prop)/1000,str(e)

    def insert_dist_metric(self, dist_metric):
        #remove boundary pts
        for distkey in dist_metric["DIST_KEYS"]:
            dist_metric[distkey]["LEN_BOUNDARY_PTS"] = len(dist_metric[distkey]["BOUNDARY_PTS"])
            dist_metric[distkey]["BOUNDARY_PTS"] = 0

        key={"_id": dist_metric["SOL_ACC"]+"/"+dist_metric["PROB_ACC"]}
        try:
            self.cdistmetric.update(key,dist_metric, upsert=True)
        except:
            print("Exists metric in DB.")

    def insert_timings(self,timing):
        d = timing
        key = {"_id":d["SOL_ACC"]+"/"+d["PROB_ACC"]}
        self.ctiming.update(key, d, upsert=True)

    def convert(self, data):
        if isinstance(data, basestring):
            return str(data)
        elif isinstance(data, collections.Mapping):
            return dict(map(self.convert, data.iteritems()))
        elif isinstance(data, collections.Iterable):
            return type(data)(map(self.convert, data))
        else:
            return data

DATA = { u'spam': u'eggs', u'foo':[u'Gah!'], u'bar': { u'baz': 97 },
         u'list': [u'list', (True, u'Maybe'),[u'and', u'a', u'set', 1]]}
d2 = {u'CENTROID': [-89.27711427042178, 33.82729807770301], u'ENV': [-89.27742599999999, -89.27676, 33.821625, 33.831120999999996], u'AREA': 3.4796230000515294e-06, u'PERI': 0.019328495060411812, u'FEAT_KEY': u'280139503003004', u'POPULATION': 1}
