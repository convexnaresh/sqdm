import collections
import sys
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

class IO:

    def __init__(self):
        from pymongo import MongoClient
        client =MongoClient('localhost', 27017)
        self.db = client.local
        self.cdplan = self.db["dplan"]
        self.cblock_prop = self.db["block_prop"]
        self.ctiming=self.db['timing']
        self.cdistmetric = self.db["dist_metric"]
        self.cdatainfo = self.db["data_info"]

    def insert_plan(self,plan):
        d = dict(plan.get_block_to_attr_map())
        print len(d)
        key={"_id":plan.get_solacc()+plan.get_problemacc()}
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
            try:
                dist_metric[distkey]["LEN_BOUNDARY_PTS"] = len(dist_metric[distkey]["BOUNDARY_PTS"])
                dist_metric[distkey]["BOUNDARY_PTS"] = 0
            except:
                continue

        key={"_id": dist_metric["SOL_ACC"]+"/"+dist_metric["PROB_ACC"]}
        try:
            self.cdistmetric.update(key,dist_metric, upsert=True)
        except:
            print("Exists metric in DB.")

    def insert_timings(self,timing):
        d = timing
        key = {"_id":d["SOL_ACC"]+"/"+d["PROB_ACC"]}
        self.ctiming.update(key, d, upsert=True)

    def insert_data_info(self,data_info):

        key = {"_id":data_info["state_id"]}
        self.cdatainfo.update(key, data_info, upsert=True)

    def convert(self, data):
        if isinstance(data, basestring):
            return str(data)
        elif isinstance(data, collections.Mapping):
            return dict(map(self.convert, data.iteritems()))
        elif isinstance(data, collections.Iterable):
            return type(data)(map(self.convert, data))
        else:
            return data

    def get_all_metrics(self):
        return self.cdistmetric.find()

    def get_data_info_by_key(self,search_key, order_by="keys"):
        curser = self.cdatainfo.find()
        keys = []
        vals = []
        for doc in curser:
            for k,v in doc.items():

                if k.find(search_key) > -1:
                    if k == "dists":
                        for k1,v1 in v.items():
                            keys += [str(doc["_id"])+str(k1)]
                            vals += [v1]

                    if k == "size_shp":
                        keys += [str(doc["_id"])]
                        vals += [v]
        if order_by == "keys":
            return collections.OrderedDict(sorted(dict(zip(keys,vals)).items()))
        else:
            d = dict(zip(keys, vals))
            od = collections.OrderedDict(sorted(d.items(), key =lambda kv:(kv[1], kv[0])))
            return od

    def get_timings_by_key(self,metric_name,avg_over_doc=False):
        curser = self.ctiming.find()
        keys=[]
        vals=[]
        dtime={}
        if avg_over_doc:
            for doc in curser:
                n =0
                avg = 0
                for k,v in doc.items():
                    if k.find(metric_name) > -1:
                        avg += v
                        n +=1
                #if a key exists.
                if n > 0:
                    keys += [doc["_id"].split("/")[1].split(".")[0]] #retrieve st id.
                    vals +=[round(avg/float(n),6)]
        else:
            for doc in curser:
                stid= doc["_id"].split("/")[1].split(".")[0]
                for k,v in doc.items():
                    if k.find(metric_name) > -1:
                        keys += [stid+"-"+k[-2:]]
                        vals +=[v]

        return collections.OrderedDict(sorted(dict(zip(keys,vals)).items()))
    def plot_data_info(self):

        search_key = "dists"
        od = IO().get_data_info_by_key(search_key, order_by="values")

        y_pos = np.arange(len(od.keys()))
        plt.bar(y_pos, od.values(), align='center', alpha=0.5)
        plt.xticks(y_pos, od.keys())

        if search_key == "dists":
            plt.ylabel('#blocks')
            plt.title('Distribution of Blocks/District (State wide)')

        if search_key == "size_shp":
            plt.ylabel('Size(MB)')
            plt.title('Spatial Data Size by State')
        plt.show()

    def plot_metric_timings(self):
        io = IO()
        metric_names = ["ROHRBACH_IDX","BLKCTR_TO_DIST_CTR_","ISO_PERIMETRIC_IDX_","DBLKCTR_TO_DIST_PERI_","BLKS_UCASCADE_"]
        #plot
        from math import log

        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.ticker import NullFormatter  # useful for `logit` scale
        # plot with various axes scales
        plt.figure(1)

        # linear
        plt.subplot(211)
        for metric_name in metric_names:
            dk = io.get_timings_by_key(metric_name, avg_over_doc=False)
            plt.plot(dk.keys(), [log(v,2) for v in dk.values()],label=metric_name)

        plt.yscale('linear')
        plt.ylabel("Time(s)")
        plt.xlabel("States/Districts")
        plt.title('Runtime for Differnt Metrics by Districts')
        plt.legend(loc='upper left')
        plt.grid(True)

        #average over state
        plt.subplot(212)
        for metric_name in metric_names:
            dkavg = io.get_timings_by_key(metric_name, avg_over_doc=True)

            plt.plot(dkavg.keys(),[log(v,2) for v in dkavg.values()] ,label="avg:"+metric_name)

        plt.yscale('linear')
        plt.ylabel("Time(s)")
        plt.xlabel("States/Districts")
        plt.title('Avg. Runtime (Over States) for Metric:'+metric_name)
        plt.grid(True)
        plt.legend(loc='upper left')
        plt.tight_layout()
        plt.show()
        print("Complete plot.")

#IO().plot_metric_timings()
#IO().plot_data_info()
import numpy as np
import matplotlib.pyplot as plt
