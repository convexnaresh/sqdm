'''This is method-1 for delegation.'''


from SQDM import SQDM, Segment
from collections import OrderedDict
from simple_polygon import util
from ast import literal_eval
#CONSTANTS
SQDM_DELEGATED_BY_ID = 777
SQDM_DELEGATED_TO_ID= 444
debug = False

class Delegation:
    def __init__(self,delegated_dic={}):

        self.__sqdm_dic = delegated_dic
        self.__delegated_dic =delegated_dic #this dic to be delegated

    def tosegsdict(self, segments, segment_id_type=util.SEGMENT_ID_TYPES['linehash']):
        '''returns a set of items of this form for each segment:
        (hash(x1,y1,x2,y2),(x1,y1,x2,y2,AreaAbove,AreaBelow,etc.))'''
        seg_dictionary = OrderedDict()

        for seg in segments:
            if segment_id_type == util.SEGMENT_ID_TYPES['linehash']:
                #entry = {seg.dictentry_linehashkey(): seg.dictentry_value()}
                seg_dictionary[seg.dictentry_linehashkey()] = seg.dictentry_value()

            if segment_id_type in [util.SEGMENT_ID_TYPES['lineid'], util.SEGMENT_ID_TYPES['lineid']]:
                seg_dictionary[seg.attr['edgeid']] =seg.dictentry_value()

        return seg_dictionary

    def mesh_vertical_sweeplines(self,mesh_segents):
        '''
        :return: a dictionary whose key is xvalue, and value is
        'index' of xvalue in sorted array.
        '''
        import collections

        if mesh_segents:
            xunik =[]
            for seg in mesh_segents:
                xunik +=[seg.getLeftPoint().getX(),seg.getRightPoint().getX()]
            xunik = sorted(set(xunik))
        else:
            return {}

        xunik_idx = enumerate(xunik)
        orderedxunik =OrderedDict()
        for idx,xval in xunik_idx:
            orderedxunik[xval] = idx

        if debug:
            print("vertical_sweeplines \n"),orderedxunik
            print("\n")
        return orderedxunik


    def segment_mapping_to_yblocks(self,segments, template_sqdm):
        '''
        :param segments: a list of Segment objects to be mapped to template_sqdm
        :param empty_sqdm: is an SQDM, a nested dictionary of <xval,<(y1-y2),<seghash:value>>>
        :return: returns the input parameter by mapping segments in self.sides.
        '''
        if not template_sqdm:
            return {}
        for seg in segments:
            x1,y1,x2,y2 = seg.co_ordinates()
            #map to x-column
            xcolumn_yblocks = template_sqdm[x1]
            for yblock in xcolumn_yblocks.keys():
                ylow, yhigh = literal_eval(yblock)
                if y1 < y2:
                    if ylow <= y1 and y2 <= yhigh:
                        #TODO: make an ordered dictionary inside.
                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1,y1,x2,y2]
                else:
                    if ylow <= y2 and y1 <= yhigh:
                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1,y1,x2,y2]

        return template_sqdm

    def split_mesh_sides_at_x(self,mesh_segents, xunikdic):

        '''splits a given random set of segments at xunikdic values.
        xunikdic must be an ordred dictionary of x-values with xkey:index, index is the
        index of xkey in sorted list.
        sides:Type Segment list
        xunikdic: Type Dictionary of sorted <x-value,index>
        retun: Type Segment list'''

        mesh_splits = []
        for seg in mesh_segents:
            # split seg
            # 1)get x-values at which split must take place
            li = xunikdic[seg.getLeftPoint().getX()]
            hi = xunikdic[seg.getRightPoint().getX()]
            splits=seg.split_at_multiple_x(xunikdic.keys()[li + 1:hi])
            mesh_splits +=splits
        return mesh_splits

    def slice_sqdm(self,lowkey,highkey,sqdm_dic):
        '''returns all dictionary keys between lower key (inclusive) and higher key (exclusive).
        This is a interval overlapping search problem. One can implement range search using RangeTrees.
        lowkey: value to search for the lowest key which is greater than this value.
        highkey: value to search for the highest key which is less than this value.
        sqdm_dic is a key-sorted dictionary.'''
        if highkey < lowkey:
            lowkey,highkey = highkey, lowkey

        print("low,high"),lowkey,highkey
        slice_dic ={}
        i = 0
        just_higherkey=None
        for key in sqdm_dic.keys():
            if float(key) >= lowkey:
                if float(key) < highkey:
                    slice_dic[float(key)] = sqdm_dic[key]
                else:
                    just_higherkey = key
                    break
            else:
                i +=1 #this key at index i just crosses this lowkey.
                continue

        if float(sqdm_dic.keys()[i]) > lowkey: #this key just past the lowkey, which means capture previous key
            lowboundary_key = sqdm_dic.keys()[i-1]
            slice_dic[float(lowboundary_key)] =sqdm_dic[lowboundary_key]
        if just_higherkey:
            slice_dic[float(just_higherkey)] ={}
        else:
            slice_dic[float(key)] = {}
        return slice_dic

    def xcolumns_yblocks(self,union_segments, unionkeys_dictionary):
        '''
        Returns an sqdm nested dictionary
        :param union_splits: list of Segment objects
        :param unionkeys_dictionary: dictionary of <x-value,index>
        :return: sqdm, a nexted dictionary <x-value,<(yi,y2),<list of segment-hash>>
        '''

        #collect all x:[(y1,y2),...] for each segments
        xcolumn_yblock_dic = OrderedDict()
        for xkey in unionkeys_dictionary.keys():
            xcolumn_yblock_dic[xkey] =[]

        for seg in union_segments:
            x1, y1, x2, y2 = seg.co_ordinates()
            if y1 < y2:
                try:
                    xcolumn_yblock_dic[x1].append((y1,y2))
                except:
                    xcolumn_yblock_dic[x1]= []
                    xcolumn_yblock_dic[x1].append((y1, y2))
            else:
                try:
                    xcolumn_yblock_dic[x1].append((y2, y1))
                except:
                    xcolumn_yblock_dic[x1] = []
                    xcolumn_yblock_dic[x1].append((y2, y1))


        #construct nested dictionary
        #  {xkey:{(y1,y2):{segid:segvalue,...}...(y1,y2):{segid:segvalue}..}, ..}
        for xkey, list_yblocks in xcolumn_yblock_dic.items():
            print list_yblocks
            ypairs_dic = OrderedDict()
            list_yspans_pairs =util.pairwise(util.non_overlaping_yspans(list_yblocks))

            count=0
            for ypair in list_yspans_pairs:
                ypairs_dic[str(ypair)] ={}
                count+=1

            xcolumn_yblock_dic[xkey] =ypairs_dic

        return xcolumn_yblock_dic

    @classmethod
    def test_delegation(self):

        #driver
        infile_Psqdm = "../out/tmp/Psqdm.json" #sqdm for P #parent polygon sqdm.
        infile_Pssegs = "../out/tmp/Ps.json" #Split of P segs #parent polygon
        infile_Cssegs = "../out/tmp/Cs.json" #split of C #child polygon

        Psqdm = util.load_sqdm_from_file(infile_Psqdm) #Parent sqdm
        Pssegs = util.load_from_file(infile_Pssegs) #parent segment
        Cssegs = util.load_from_file(infile_Cssegs) #child polygon's original segments
        Csqdm = SQDM.test3() #child sqdm


        #get lower xkey and higher xkey from delegation_sqdm
        delegation = Delegation()
        lowkey, highkey = Csqdm.keys()[0], Csqdm.keys()[-1]

        #find the vcolumns intersected by delegation_sqdm.
        Pslice_sqdm = delegation.slice_sqdm(lowkey, highkey, Psqdm)

        #find union of two set of x-keys
        PuCkeys = Pslice_sqdm.keys() + Csqdm.keys()
        PuCkeys = SQDM().get_ordered_keys(PuCkeys)

        #A)
        #split slice_dict at x={dsqdm.keys()}
        # 1)get segments for Pslice_sqdm
        Pslice_seg_keys = SQDM().sqdm_segments(Pslice_sqdm)
        Pslice_segments = [Segment().ntuples_to_seg(Pssegs[segkey],util.SEGMENT_TUPLE_KEYS) for segkey in Pslice_seg_keys]
        # 2)split segments UnionKeys
        P_slice_splitsPuC = delegation.split_mesh_sides_at_x(Pslice_segments, SQDM().get_ordered_keys(PuCkeys))
        print len(P_slice_splitsPuC)

        #B)
        #split Csqdm at x={isqdm.keys()}
        C_seg_keys = SQDM().sqdm_segments(Csqdm)
        C_segments = [Segment().ntuples_to_seg(Cssegs[segkey],util.SEGMENT_TUPLE_KEYS) for segkey in C_seg_keys]
        # 2)split segments at UnionKeys
        C_splitsPuC = delegation.split_mesh_sides_at_x(C_segments, SQDM().get_ordered_keys(PuCkeys))
        print len(C_segments), len(C_splitsPuC)

        UnionSplits = C_splitsPuC + P_slice_splitsPuC

        templateMerged_Sqdm = delegation.xcolumns_yblocks(UnionSplits,PuCkeys) #
        print templateMerged_Sqdm
        print("Construction of merged-sqdm completed:")
        for xkey,yblocks in templateMerged_Sqdm.items():
            print xkey
            for yblock,segments in yblocks.items():
                print("\t"), yblock, segments

        print

        #mapping the segments to merged sqdm.
        # Extract assignment Vector.
        #1) Tag region above and region below by Owner's Tags for All segments in C_splitsPuC
        print("Changing Child Region codes to Original Owner.")
        AssignmentVector =[]
        for C_seg in C_splitsPuC:
            abv,bel = C_seg.getAttrByName('abv'),C_seg.getAttrByName('bel')
            C_seg.setAttrByName('abv',SQDM_DELEGATED_BY_ID)
            C_seg.setAttrByName("bel",SQDM_DELEGATED_BY_ID)
            AssignmentVector.append((abv,bel))

        print AssignmentVector
        #1) Delegate segments by Tagging segments in C_sqdm  individually according to an assignment vector.
        print
        seg_index = 0
        for C_seg in C_splitsPuC:
            abv,bel = AssignmentVector[seg_index]
            if abv == util.POLYGON_INNERID:
                C_seg.setAttrByName("abv", SQDM_DELEGATED_TO_ID)
            else:
                C_seg.setAttrByName("bel", SQDM_DELEGATED_TO_ID)
            seg_index +=1

        print("Segment Delegation Completed. Save these segments.")
        delegated_seg_dic = delegation.tosegsdict(C_splitsPuC)
        util.save(delegated_seg_dic,"../out/tmp/Cd.json")

        ## 1. Map Delegated Segments to templateMerged_Sqdm
        ## 2. Save delegated SQDM
        ## 3. Remove all y-blocks outside boundary for Child
        templateMerged_Sqdm = delegation.segment_mapping_to_yblocks(C_splitsPuC,templateMerged_Sqdm)
        print("Child Segments Mapped templateMerged_Sqdm completed:")
        util.save(templateMerged_Sqdm,"../out/tmp/Cdsqdm.json")

        # Map original splits from Parent's Slice to templateMerged_SQDM
        ##Merge this parital sdm to the main Psqdm
        templateMerged_Sqdm = delegation.segment_mapping_to_yblocks(P_slice_splitsPuC,templateMerged_Sqdm)
        print("Sliced Parent Segments Mapped templateMerged_Sqdm completed:")

        ##Merge this templateMerged_Sqdm to the main Psqdm
        #1. Paritally Lock the Psqdm from (lowkey,highkey(exclusive))
            #which means that, query to Psqdm cannot respond to queries
            #for all x>=lowkey, and x<highkey.
            #locking must be done for certain time to avoid permanent locking.
        #2. Insert new or Update old xkey from templateMerged_Sqdm into main sqdm Psqdm
        #3. Remove

        for xkey in templateMerged_Sqdm.keys()[:-1]: #skip the highest key.
            Psqdm[xkey] =templateMerged_Sqdm[xkey]
        print("Insert/Update original sqdm after delegation completed.")

        #save the updated Psqdm
        util.save(Psqdm,"../out/tmp/Psqdm-1.1.json")

Delegation.test_delegation()











