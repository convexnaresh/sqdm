from simple_polygon import util
import numpy as np
from sortedcontainers import SortedDict
import time
dplan_attr={"map":"map",
            "nblocks":"nblocks",
            "timestamp":"timestamp",
            "checksum":"checksum",
            "merkleroot":"merkleroot",
            "sol_from_acc":"solacc",
            "prob_from_acc":"problemacc"}

PROBLEM_ACC = util.hashargs("ms.dplan.com")
SOLUTION_ACC = util.hashargs("cse.dplan.com")

class DPmap(object):

    def __init__(self,ndist=None):
        self._ndist = ndist
        self.blocks ={}
        self.bvertices ={}
        self.problemacc= None

    def dproblem_from_files(self):
        pass

class Dplan(object):

    def __init__(self):

        self.__problemacc = None
        self.__solacc = None
        self.__map = {} #block eqiv dict.
        self.__merkleroot = None
        self.__checksum = None
        self.__timestamp = None
        self.__nblocks = None

    def get_problemacc(self):
        return self.__problemacc

    def get_solacc(self):
        return self.__solacc

    def get_map(self):
        return self.__map

    def get_timestamp(self):
        return self.__timestamp

    def get_nblocks(self):
        return self.__nblocks

    def get_checksum(self):
        return self.__checksum

    def get_merkleroot(self):
        return self.__merkleroot

    def reset(self):
        self.__problemacc = None
        self.__solacc = None
        self.__map = {} #block eqiv dict.
        self.__merkleroot = None
        self.__checksum = None
        self.__timestamp = None
        self.__nblocks = None

    def load_json_dplan(self,dplanfile):
        '''

        :param jsonfile: json file with <key:h(district_id+block_id),value=(district_id,block_id)>
        :return: a dictionary that is collection of districts and blocks created.
        '''
        plandata = util.load_from_file(dplanfile)
        if plandata:
            self.__map = plandata[dplan_attr["map"]]
            self.__nblocks = plandata[dplan_attr["nblocks"]]
            self.__timestamp = plandata[dplan_attr["timestamp"]]
            self.__solacc = plandata[dplan_attr["sol_from_acc"]]
            self.__problemacc = plandata[dplan_attr["prob_from_acc"]]
            self.__checksum = plandata[dplan_attr["checksum"]]
            self.__merkleroot = plandata[dplan_attr["merkleroot"]]

        else:
            raise Exception("Cannot load plan data.")

    def build_from_blkequiv_csv(self, blk_equivfile):

        arr = np.loadtxt(blk_equivfile,delimiter=',')
        print arr.shape
        # if <blockid is a not single string.
        if arr.shape[1] >2:
            arr = arr.astype(int)

        dict = SortedDict()
        if arr.shape[1] > 2:
            for blkvals in arr[0:3]:

                pkey = util.hash_cblocks_dist(tuple(blkvals[0:]))
                dict[pkey] = tuple(blkvals)

        self.__map = dict
        self.__nblocks = len(list(dict.keys()))
        self.__problemacc = PROBLEM_ACC

    def to_json(self,outfile):
        dplan_data ={}
        dplan_data[dplan_attr['map']] = self.get_map()
        dplan_data[dplan_attr['nblocks']] = self.get_nblocks()
        dplan_data[dplan_attr['timestamp']] = time.time()
        dplan_data[dplan_attr['merkleroot']] = self.get_merkleroot()
        dplan_data[dplan_attr['prob_from_acc']] = self.get_problemacc()
        dplan_data[dplan_attr['checksum']] = self.get_checksum()
        util.save(dplan_data,outfile)

        dplan_data[dplan_attr['sol_from_acc']] = SOLUTION_ACC

home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"
blk_equivfile = "congress_blk_equiv_coma.txt" #with 5 cols
#blk_equivfile = "congress_blk_equiv.txt" #with 2 cols

dplan = Dplan()
dplan.build_from_blkequiv_csv(home +blk_equivfile)
for k,v in dplan.get_map().items():
    print k,v

dplan.to_json(home+"dplan")

dplan.reset()
print("Load dplan from json.")
dplan.load_json_dplan(home+"dplan.json")

for k,v in dplan.get_map().items():
    print k,v













