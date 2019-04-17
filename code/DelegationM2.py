'''This is method-2 for delegation.'''


from SQDM import SQDM, Segment, Point
from collections import OrderedDict
from simple_polygon import util
from ast import literal_eval
from sortedcontainers import SortedDict
import copy
#CONSTANTS
SQDM_DELEGATED_BY_ID = 777
SQDM_DELEGATED_TO_ID= 444
debug = True

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

    def containing_box(self,segment_obj,sqdm):
        x1, y1, x2, y2 = segment_obj.co_ordinates()
        xkyes = list(sqdm.keys())
        # search along vertical slabs (x_i -- x_j)
        xsearchidx = self.xrange_search(xkyes, x1, x2)
        if xsearchidx < 0:
            return None, None
        else:
            # search along yblocks (y_i -- y_j)
            ykeys_tuples = sqdm[xkyes[xsearchidx]].keys()
            if y1 > y2:
                y1,y2 = y2,y1
            ysearchidx = self.yrange_search(ykeys_tuples, y1, y2)
            if ysearchidx < 0:
                return None, None
            else:
                xkey = xkyes[xsearchidx]
                xnkey = xkyes[xsearchidx+1]
                ykeytuple = ykeys_tuples[ysearchidx]
                # return singular_sqdm
                # singular_sqdm =SortedDict()
                #singular_sqdm[xkey] = {}
                #singular_sqdm[xkey][ykeytuple] = sqdm[xkey][ykeytuple]
                #singular_sqdm[xnkey] ={}

                return ((xkey,xnkey),ykeytuple) #,singular_sqdm #(xkyes[xsearchidx],xkyes[xsearchidx+1]),ykeys_tuples[ysearchidx]

    def xrange_search(self,arr, search_key1,search_key2):
        '''
        :arr : sorted array of items.
        :param search_key: is a search-key
        :return: lower index i of the item in arr such that
        arr[i]<= search_key1 and arr[i+1] >=search_key1 and
        arr[i]<= search_key2 and arr[i+1] >=search_key2
        not found returns -1
        '''
        # range search
        arr += arr[0:1]
        l = 0
        r = len(arr) - 2

        while l <= r:
            mid = l + (r - l) / 2;
            # Check if x is present at mid
            # print mid

            if (arr[mid] <= search_key1 and search_key1 <= arr[mid + 1]) and \
                    (arr[mid] <= search_key2 and search_key2 <= arr[mid + 1]):
                return mid

            # If x is greater, ignore left half
            elif arr[mid] < search_key1:
                l = mid + 1

            # If x is smaller, ignore right half
            else:
                r = mid - 1

        # If we reach here, then the element was not present
        return -1

    def yrange_search(self,arr, search_key1,search_key2):
        '''
        :param arr: sorted arrays of tuples like (ya,yb)
        :param search_key1: is a search-key
        :param search_key2: is a search-key
        :return: lower index i of the item in arr such that
        arr[i][0]<= search_key1 and arr[i][1] >=search_key1 and
        arr[i][0]<= search_key2 and arr[i][1] >=search_key2
        not found returns -1
        '''
        if search_key1 > search_key2:
            search_key1, search_key2 = search_key2, search_key1

        # range search
        l = 0
        r = len(arr) - 1

        while l <= r:
            mid = l + (r - l) / 2;
            # Check if x is present at mid
            # print mid
            if (arr[mid][0] <= search_key1 and search_key1 <= arr[mid][1]) and \
                    (arr[mid][0] <= search_key2 and search_key2 <= arr[mid][1]):
                return mid

            # If x is greater, ignore left half
            elif arr[mid][0] < search_key1:
                l = mid + 1

            # If x is smaller, ignore right half
            else:
                r = mid - 1

        # If we reach here, then the element was not present
        return -1


    def split_mesh_segs_at_x(self,mesh_segments, xunikdic,doround=False):

        '''splits a given random set of segments at xunikdic values.
        xunikdic must be an ordred dictionary of x-values with xkey:index, index is the
        index of xkey in sorted list.
        sides:Type Segment list
        xunikdic: Type Dictionary of sorted <x-value,index>
        retun: Type Segment list'''

        mesh_splits = []
        for seg in mesh_segments:
            # split seg
            # 1)get x-values at which split must take place
            li = xunikdic[seg.getLeftPoint().getX()]
            hi = xunikdic[seg.getRightPoint().getX()]
            splits = seg.split_at_multiple_x(xunikdic.keys()[li + 1:hi],doround)
            mesh_splits +=splits
        return mesh_splits

    def split_mesh_side_at_y(self,seg, yunikdict,doround=False):

        '''splits a given random set of segments at xunikdic values.
        yunikdict must be an ordred dictionary of y-values with xkey:index, index is the
        index of xkey in sorted list.
        sides:Type Segment list
        yunikdict: Type Dictionary of sorted <x-value,index>
        retun: Type Segment object seg'''
        mesh_splits=[]
        # split seg
        # 1)get y-values at which split must take place
        y1 = seg.getLeftPoint().getY()
        y2 = seg.getRightPoint().getY()
        if y1 > y2:
            y1,y2 = y2, y1

        li = yunikdict[y1]
        hi = yunikdict[y2]

        splits = seg.split_at_multiple_y(yunikdict.keys()[li + 1:hi],doround)
        mesh_splits +=splits
        return mesh_splits
        #return yunikdict.keys()[li + 1:hi]

    def slice_sqdm(self,lowkey,highkey,sqdm_dic):
        '''returns all dictionary keys between lower key (inclusive) and higher key (exclusive).
        This is a interval overlapping search problem. One can implement range search using RangeTrees.
        lowkey: value to search for the lowest key which is greater than this value.
        highkey: value to search for the highest key which is less than this value.
        sqdm_dic is a key-sorted dictionary.'''
        if highkey < lowkey:
            lowkey,highkey = highkey, lowkey

        slice_dic = SortedDict()
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
            print xkey, list_yblocks
            ypairs_dic = OrderedDict()
            list_yspans_pairs =util.pairwise(util.non_overlaping_yspans(list_yblocks))

            count=0
            for ypair in list_yspans_pairs:
                ypairs_dic[ypair] ={}
                count+=1

            xcolumn_yblock_dic[xkey] =ypairs_dic

        return xcolumn_yblock_dic

    def ravel_sqdm(self,oldsqdm):
        '''

        :param oldsqdm: is a nested dictionary {x1:{yblock1:{line1,line2..},...}, x2:{yblock1:{}}}
        :return: (x1,x2,y1,y2,[line1,line2,..]),...
        '''
        rect_tuples = []
        xkeys = oldsqdm.keys() #xkeys
        print("\t \t xkeys:"), xkeys
        for i in xrange(len(xkeys)-1):
            xkey = xkeys[i]
            nxkey = xkeys[i+1]
            for yb,lines in oldsqdm[xkey].items():
                rect_tuples += [(xkey,nxkey) + yb + (lines,)]

        #wrap up space after xn in sqdm is (xn,x0)
        xkey = xkeys[-1] #last key
        nxkey = xkeys[0] #first item
        for yb,lines in oldsqdm[xkey].items():
            rect_tuples += [(xkey,nxkey)+ yb + (lines,)] #yb is tuple
        else:
            #TODO:
            #rect_tuples += [(xkey, nxkey, 0,0, {})]
            pass

        if debug:
            for tup in rect_tuples:
                print("\t rect:"),tup
        return  rect_tuples

    def search_containing_rect(self, segment,list_rectangles):

        if len(list_rectangles) ==0:
            return None
        if len(list_rectangles[0]) < 4:
            print("Invalid rectangle entry")
            return None
        for rect in list_rectangles:
            xmin,xmax,ymin,ymax = rect[0:4]
            if xmin > xmax:
                xmin, xmax = xmax,xmin
            if ymin > ymax:
                ymin,ymax = ymax,ymin

            x1,y1,x2,y2 = segment.co_ordinates()
            if xmin <= x1 <= xmax and xmin <= x2 <= xmax:
                if y1 < y2:
                    if ymin <= y1 <= ymax and ymin <= y2 <= ymax:
                        return rect
                else:
                    if ymin <= y2 <= ymax and ymin <= y1 <= ymax:
                        return rect

        return None


    def islpsegment(self,new_segment, singular_sqdm):
        '''
        legalparition segmment

        :param new_segment: new segment object
        :param singular_sqdm: it is a singular_sqdm dictionary with (xkey,nxkey):{(ykey,nykey):{segment_ids}}
        :return: if this new_segment is a legal partition segment for a given singular_sqdm. It must test
        legality against all the existing segments in this singular_sqdm.
        '''
        #A,B = new_segment.p1, new_segment.p2
        for lines in singular_sqdm.values():
            #TODO: check the validity of the new_segment against all the existing segents in this yblock
            #get line L in this box
            #test A and B on L, return True
            #test A on L and Above(L)==BY and B.Above(L)
            #test B on L and Above(L)==By and A.Above(L)
            #test A.Above(L) and B.Above(L) and Above(L)==BY
            #singular_sqdm[xkey][ybkey[0]][new_segment.dictentry_linehashkey()] = list(new_segment.co_ordinates())
            pass

        return True

    def sqdm_storageid(self,xkey,yblock):
        return str(xkey) + "-" + str(yblock[0])+"-"+str(yblock[1]) + ".json"

    @classmethod
    def pre_process_dummy_polygon(self):
        from SQDM import Polygon
        delegation = Delegation()
        subsqdm_file = "../out/tmp/"
        infile_Psqdm = "../out/tmp/aPsqdm.json" #sqdm for P #parent polygon sqdm.
        infile_Pssegs = "../out/tmp/aPs.json" #Split of P segs #parent polygon
        Psqdm = util.load_sqdm_from_file(infile_Psqdm) #Parent sqdm
        Pssegs_dictionary = util.load_from_file(infile_Pssegs)  # parent segment

        #Ca = [(3.6, 4.8), (3.8, 4.5),(3.6,4.2), (3.2,4),
        #(3.6,3.8),(3.8,3.6), (3.2,3.2),(3,3.6),(2.6,3.4),(2.8,3.8),(2.6,4), (3,4.2),(3.2,4.5)]
        Ca = [(3.6, 4.8), (3.8, 4.6), (3.6, 4.2), (3.4, 4), (3.6, 3.8),
                (3.8, 3.6), (3.4, 3.2),
                (2.8, 3.6), (2.6, 3.4),
                (2.8, 3.8), (2.6, 4)]
        Cpoly = Polygon(Ca)
        AV = Cpoly.label_bottom_up(util.POLYGON_INNERID)
        Cseg_dictionary = Cpoly.tosegsdict()

        #1. VerifySimplePolygon(CPoly) : it must be done  prior to adding segments into an existing sqdm.
        #it is because, an intersecting segment forces to check it's intersection with all other segments in
        #the same polygon#.

        #since Cpoly is completely inside Polygon P, region abv and bel is
        #both assigned to be same.

        for C_seg in Cpoly.sides():
            C_seg.setAttrByName('abv',SQDM_DELEGATED_BY_ID)
            C_seg.setAttrByName("bel",SQDM_DELEGATED_BY_ID)

        #find the part of Psqdm intersected by C .
        Cxkeys = Cpoly.vertical_sweeplines().keys()
        lowkey= Cxkeys[0]
        highkey = Cxkeys[-1]
        Pslice_sqdm = delegation.slice_sqdm(lowkey, highkey, Psqdm)
        psqdmlowkey = Pslice_sqdm.keys()[0]
        print("Completed Slicing sqdm for Cpoly..\n")

        ##
        print("C.xmin, C.xmax"), lowkey, highkey
        print("Pslice.xmin,Pslice.xmax"), Pslice_sqdm.keys()[0], Pslice_sqdm.keys()[-1]



        ##new code.
        ##new segments are split at x-values unique to CPoly and that
        # is overlapped by CPoly on Psqdm.
        #(Optimal: only at x-values unique to Psqdm.x values.)
        PuCxkeys = Cpoly.vertical_sweeplines().keys() + list(Pslice_sqdm.keys())
        Csegs_splits_atx = delegation.split_mesh_segs_at_x(Cpoly.sides(),
                                                           SQDM().get_ordered_keys(PuCxkeys),
                                                           doround=True)
        for seg in Cpoly.sides():
            print("..."), seg.co_ordinates()
        for seg in Csegs_splits_atx:
            print("..........."), seg.co_ordinates()
        Csegs_splits_dict = OrderedDict(enumerate(Csegs_splits_atx))

        #these splits from new segments are split at y-values on each of the vertical slabs.
        #1. collect new-splits for each vertical columns in Slice_sqdm
        seg_buckes =OrderedDict()
        for xkey in list(Pslice_sqdm.keys()):
            seg_buckes[xkey] = []
        print seg_buckes

        for seg_id,new_split in Csegs_splits_dict.items():
            #search vcolumn for each seg.
            #add this new-split to the bucket of vcolumn
            xlow,xhigh = new_split.getLeftPoint().getX(),new_split.getRightPoint().getX()
            xindex = delegation.xrange_search(list(Pslice_sqdm.keys()),xlow,xhigh)
            xkey = list(Pslice_sqdm.keys())[xindex]
            seg_buckes[xkey].append(seg_id)
        #2. for each segments in vertical buckets/columns, try splitting at y-values
        print("------------")
        CPolySides = []
        for xkey,seg_ids in seg_buckes.items()[:]:
            y_blocks = Pslice_sqdm[xkey].keys()
            y_values = list(set( yval for tup in y_blocks for yval in tup))

            for seg_id in seg_ids[:]:
                seg =Csegs_splits_dict[seg_id]
                print("seg"),seg.co_ordinates()
                y1,y2 = seg.getLeftPoint().getY(),seg.getRightPoint().getY()

                y_values += [y1,y2]
                ordered_ydict = SQDM().get_ordered_keys(y_values)
                splits = delegation.split_mesh_side_at_y(seg,ordered_ydict,doround=True)
                CPolySides += splits
                for sp in splits:
                    print("\t"),xkey, sp.co_ordinates()
                #delete last two elements
                del y_values[-1]
                del y_values[-1]
            del y_values
        print
        print
        #Reconstruct the dictionary for segents.
        ftCPoly_splitted = Polygon([])
        ftCPoly_splitted.setSides(CPolySides)
        ftCseg_dictionary = ftCPoly_splitted.tosegsdict()
        print("length after splitting at x and at y."),len(CPolySides)
        for k, v in ftCseg_dictionary.items():
            print k, v
        print("length(fitable-segments in Cpoly)"),len(CPolySides)

        Pslice_segs_dictionary = {}

        for xkey, yblocks in Pslice_sqdm.items():
            for yb, lines in yblocks.items():
                for linehash in lines.keys():
                    Pslice_segs_dictionary[linehash] = Pssegs_dictionary[linehash]
        del Psqdm
        del Pssegs_dictionary

        util.save(Pslice_segs_dictionary, "../out/tmp/aPslicesegs.json")  # segment dictionary
        util.save_sqdm(Pslice_sqdm, "../out/tmp/aPslicesqdm.json")
        util.save(ftCseg_dictionary, "../out/tmp/aftCsegs.json")  # segment dictionary

        #return Pslice_segs_dictionary, Pslice_sqdm, ftCseg_dictionary

    @classmethod
    def test_dummy_poly_delegation(self):
        import time
        delegation = Delegation()
        t0 = time.time()
        Pslice_segs_dictionary = util.load_from_file("../out/tmp/aPslicesegs.json")
        Pslice_sqdm = util.load_sqdm_from_file("../out/tmp/aPslicesqdm.json")
        ftCseg_dictionary = util.load_from_file("../out/tmp/aftCsegs.json")

        print("time:"), time.time() - t0
        ##
        '''
        Find a yblock in current sqdm containing a new-segment; and test if a segment is legal to draw.
        Collect all segents in the box (x1,x2)x(y1,y2)
        1. for each s_i in Cpoly:
                #get containing box:(x1,x2)x(y1,y2) from Pslice_sqdm
                #test if s_i is a legal segment that can be inserted in the box; 
                    #s_i does not intersect any segment in the box
                    #s_i completely fall inside the box or 
                    #s_i is on boundary of the box
                #if legal: a) add to a dictionary[x1][(y1,y2)]
                           b) add all segments to a dictionary dictionary[x1][(y1,y2)]     
        '''
        #Map every segments in ftCseg_dictionary to PSlice_sqdm box.
        fixables_xkey_yblocks={} #it contains valid partition segments, old segments for each box (x1,x2)x(y1,y2)
        for segkey in ftCseg_dictionary.keys():
            seg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
            #get box (x1,x2)x(y1,y2) containing the segment seg.
            xpair,yblock = delegation.containing_box(seg,Pslice_sqdm)
            if xpair == None:
                xpair, yblock = delegation.containing_box(seg, Pslice_sqdm)
            xkey,xnkey = xpair

            legalpartition = delegation.islpsegment(seg,Pslice_sqdm[xkey][yblock])
            #any illegal partition line found, return False
            if not legalpartition:
                return False

            Pslice_sqdm[xkey][yblock][seg.dictentry_linehashkey()] = [seg.co_ordinates()]

            #collect all the segments in the containing box for further divisions.
            if xkey in fixables_xkey_yblocks:
                if yblock not in fixables_xkey_yblocks[xkey]:
                    fixables_xkey_yblocks[xkey][yblock]={}
            else:
                fixables_xkey_yblocks[xkey]={}
                fixables_xkey_yblocks[xkey][yblock]={}

        print("Completed checking legal paritions and collecting all legal partiion in fixable-xkey-yblocks. \n")

        #TODO:Change new segment's abv/bel labels to DELEGATED_TO and Establish equivalent ..
        #TODO: ..relation between DELEGATED_TO and DELEGATED_BY relation.

        #
        for xkey, yblocks in Pslice_sqdm.items():
            print xkey
            for yblock, seg_dict in yblocks.items():
                print("\t yblock:"),yblock
                for lh,lc in seg_dict.items(): #linehash, line cord.
                    print("\t -- \t"),lh,lc
        print
        print ("fixables:")

        for xkey,yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:"), yblock, "n(segments)",len(Pslice_sqdm[xkey][yblock])
        '''
        For segments (new segment + existing segments) in its containing box, develop a template_subsqdm.
        1. for each segment_mesh = fixables_xkey_yblocks[x_i][yblock_j]:
            a) compute a new mesh_sqdm
            b) update Psqdm[x_i][yblock_j] by a reference to new mesh_sqdm
        '''
        #iterate through fixable y-blocks and split all segments inside it.
        print("splitting new segments.")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:"), yblock
                yblock_segs_keys = Pslice_sqdm[xkey][yblock].keys()
                yblock_segs = []

                #iterate through segment's keys.
                for segkey in yblock_segs_keys:
                    try:
                        yblock_segs += [Segment().ntuples_to_seg(ftCseg_dictionary[segkey],
                                                                 util.SEGMENT_TUPLE_KEYS)]
                    except:
                        yblock_segs += [Segment().ntuples_to_seg(Pslice_segs_dictionary[segkey],
                                                                 util.SEGMENT_TUPLE_KEYS)]

                yblock_xkeys = delegation.mesh_vertical_sweeplines(yblock_segs)
                yblock_split_keys = {} #5f7cb0df8e21
                for seg in yblock_segs:
                    print("\t\t-"), seg.co_ordinates()
                    # split seg
                    # 1)get x-values at which split must take place
                    li = yblock_xkeys[seg.getLeftPoint().getX()]
                    hi = yblock_xkeys[seg.getRightPoint().getX()]
                    splits = seg.split_at_multiple_x(yblock_xkeys.keys()[li + 1:hi], doround=True)
                    for splitseg in splits:
                        print("\t\t\t splits:"), splitseg.dictentry_linehashkey(),splitseg.co_ordinates()
                        yblock_split_keys[splitseg.dictentry_linehashkey()] = splitseg.dictentry_value()

                Pslice_sqdm[xkey][yblock] = yblock_split_keys
                print("len splits in yblock"), len(yblock_split_keys),len(Pslice_sqdm[xkey][yblock])

            print

                #TODO: 1) issue: if a parital vertical slab has only one horizontal line, then x: yblock: is empty.
                #singular_template_sqdm = delegation.xcolumns_yblocks(mesh_splits, mesh_xkeys)
                #rectangles_tuples = delegation.ravel_sqdm(singular_template_sqdm)
                #mappable_splits_rectangles += [(mesh_splits,rectangles_tuples)]

        print("Completed splitting segs in fixable y-blocks")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:"), yblock, "n(segments)"
                for segk,val in Pslice_sqdm[xkey][yblock].items():
                    print("\t\t\t"), segk,val
                print
        print
        ##
        ##Construct template sqdm for each fixable xkey,yblock
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print("xkey:--"),xkey
            for yblock in yblocks.keys():
                print("\t\tyblock:-"),yblock
                yblock_segs_keys = Pslice_sqdm[xkey][yblock].keys()
                yblock_segs = []
                for segkey in yblock_segs_keys:
                    segobj = Segment().ntuples_to_seg(Pslice_sqdm[xkey][yblock][segkey],
                                                      util.SEGMENT_TUPLE_KEYS)
                    yblock_segs += [segobj]
                    print("\t\t"),segkey, segobj.co_ordinates()

                yblock_xkeys = delegation.mesh_vertical_sweeplines(yblock_segs)
                yblock_template_sqdm = delegation.xcolumns_yblocks(yblock_segs, yblock_xkeys)
                fixables_xkey_yblocks[xkey][yblock] = yblock_template_sqdm #
            print

        print("Completed constructing template-sqdm for fixable-xkey-yblock")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:--"), yblock
                for tkey in fixables_xkey_yblocks[xkey][yblock].keys():
                    print("--\t\t"), tkey
                    for tyblock in fixables_xkey_yblocks[xkey][yblock][tkey].keys():
                        print("--\t\t\t"), tyblock

        print

        print("--------------------------------------------------------------------")
        print
        #map each of the split segents in PSlice_sqdm fixable xkey,yblock
        for fxkey, fyblocks in fixables_xkey_yblocks.items():
            print fxkey
            for fyblock in fyblocks.keys():
                #get template sqdm for this yblock
                yblock_template_sqdm = fixables_xkey_yblocks[fxkey][fyblock]

                #get fixable segments from PSlice_sqdm
                yblock_segs_keys = Pslice_sqdm[fxkey][fyblock].keys()
                print("keys:"),yblock_segs_keys
                for segkey in yblock_segs_keys:
                    print("..."), segkey, fxkey, yblock
                    seg = Segment().ntuples_to_seg(Pslice_sqdm[fxkey][fyblock][segkey],
                                                   util.SEGMENT_TUPLE_KEYS)
                    #map this segobj to yblock_template_sqdm
                    #find a containing box
                    x1, y1, x2, y2 = seg.co_ordinates()
                    # map to x-column
                    xcolumn_yblocks = yblock_template_sqdm[x1]
                    for yblock in xcolumn_yblocks.keys():
                        ylow, yhigh = yblock
                        if y1 < y2:
                            if ylow <= y1 and y2 <= yhigh:
                                # TODO: make an ordered dictionary inside.
                                try:
                                    xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1, y1, x2, y2]
                                except:
                                    print("Exception in matching."), seg.co_ordinates()
                        else:
                            if ylow <= y2 and y1 <= yhigh:
                                try:
                                    xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1, y1, x2, y2]
                                except:
                                    print("Exception in matching"), seg.co_ordinates()
        print
        print("Completed mapping segments to template-sqdm for fixable-xkey-yblock")
        for fxkey, fyblocks in fixables_xkey_yblocks.items():
            print("fxkey"),fxkey
            for fyblock in fyblocks.keys():
                print("\t\t fyblock:--"), fyblock
                fixed_sqdm = fixables_xkey_yblocks[fxkey][fyblock]
                Delegation.print_sqdm(fixed_sqdm)

        return fixables_xkey_yblocks

        #merge fixed_sqdm to Psqdm.

    @classmethod
    def print_sqdm(self,sqdm):
        for xkey, yblocks in sqdm.items():
            print("\t \t \t xkey:"), xkey
            for yb, lines in yblocks.items():
                print("\t\t\t\t"), yb
                for lhash, lval in lines.items():
                    print("\t\t\t\t\t"), lhash, lval

    @classmethod
    def preprocess_us_sqdm(self):
        from SQDM import Polygon
        delegation = Delegation()

        subsqdm_file = "../out/tmp/"
        infile_Psqdm = "../out/tmp/2USA_sqdm.json" #sqdm for P #parent polygon sqdm.
        infile_Pssegs = "../out/tmp/2USAs.json" #Split of P segs #parent polygon

        Psqdm = util.load_sqdm_from_file(infile_Psqdm) #Parent sqdm
        Cpoly = SQDM.get_usa_state_boundary_by_name("CALIFORNIA")
        Cpoly = Cpoly.slice_poly(500)
        print("len(Pssegs),len(Csegs)"), Cpoly.nv
        Cxkeys = Cpoly.vertical_sweeplines().keys()

        #find the part of Psqdm intersected by C .
        lowkey, highkey = Cxkeys[0], Cxkeys[-1]
        Pslice_sqdm = delegation.slice_sqdm(lowkey, highkey, Psqdm)
        del Psqdm
        print("Completed Slicing sqdm for Cpoly. len(Pslice_sqdm)"), len(Pslice_sqdm)

        #AV = Cpoly.label_bottom_up(util.POLYGON_INNERID)
        #util.poly_ptstoshp(Cpoly.get_vertices(), "../out/tmp/2Calif-500")

        def verifySimplePolygon():
            #1. VerifySimplePolygon(CPoly) : it must be done  prior to adding segments into
            # an existing sqdm.
            #it is because, an intersecting segment forces to check it's intersection
            # with all other segments in
            #the same polygon# .

            #since Cpoly is completely inside Polygon P, region abv and bel is
            #both assigned to be same
            pass
        #endf

        def annotate_del_by(Cpoly):
            for C_seg in Cpoly.sides():
                C_seg.setAttrByName('abv',SQDM_DELEGATED_BY_ID)
                C_seg.setAttrByName("bel",SQDM_DELEGATED_BY_ID)
        #endf
        annotate_del_by(Cpoly)

        xun = sorted(set(Pslice_sqdm.keys()))
        #plot only random 50%
        util.save_xun_asshp(xun,"../out/tmp/Xun",Cpoly.polygon_extent(Cpoly.sides()))
        #print Cpoly.polygon_extent(Cpoly.sides())

        ##new code.
        ##new segments are split at x-values unique to CPoly and that is overlapped by CPoly on Psqdm.
        #(Optimal: only at x-values unique to Pslice_sqdm.x values.)
        def split_C_at_CxPx(Cxkeys,Pslice_sqdm):
            PuCxkeys = Cxkeys + list(Pslice_sqdm.keys())
            Csegs_splits_atx = delegation.split_mesh_segs_at_x(Cpoly.sides(),
                                                               SQDM().get_ordered_keys(PuCxkeys),
                                                               doround=True)
            Csegs_splits_dict = OrderedDict(enumerate(Csegs_splits_atx))
            del Csegs_splits_atx
            del Cxkeys
            return Csegs_splits_dict

        Csegs_splits_dict = split_C_at_CxPx(Cxkeys,Pslice_sqdm)

        def fittable_segs(Pslice_sqdm,Csegs_splits_dict):

            #these splits from new segments are split at y-values on each of the vertical slabs.
            #1. collect new-splits for each vertical columns in Slice_sqdm
            seg_buckes =OrderedDict()
            for xkey in list(Pslice_sqdm.keys()):
                seg_buckes[xkey] = []
            print len(seg_buckes)

            #partition segments split at x into vertical columns/buckets.
            for seg_id,new_split in Csegs_splits_dict.items():
                #search vcolumn for each seg.
                #add this new-split to the bucket of vcolumn
                xlow,xhigh = new_split.getLeftPoint().getX(),new_split.getRightPoint().getX()
                xindex = delegation.xrange_search(list(Pslice_sqdm.keys()),xlow,xhigh)
                xkey = list(Pslice_sqdm.keys())[xindex]
                seg_buckes[xkey].append(seg_id)

            #2. for each segments in vertical buckets/columns, try splitting at y-values
            ftCPolySides = []
            for xkey,seg_ids in seg_buckes.items()[:]:
                y_blocks = Pslice_sqdm[xkey].keys()
                y_values = list(set( yval for tup in y_blocks for yval in tup))

                for seg_id in seg_ids[:]:
                    seg =Csegs_splits_dict[seg_id]

                    y1,y2 = seg.getLeftPoint().getY(),seg.getRightPoint().getY()

                    y_values += [y1,y2]
                    ordered_ydict = SQDM().get_ordered_keys(y_values)
                    splits = delegation.split_mesh_side_at_y(seg, ordered_ydict, doround=True)
                    ftCPolySides += splits
                    #delete last two elements
                    del y_values[-1]
                    del y_values[-1]
                del y_values
            print("length(fitable-segments in Cpoly)"), len(ftCPolySides)

            # Reconstruct the dictionary for segments.
            ftCPoly_splitted = Polygon([])
            ftCPoly_splitted.setSides(ftCPolySides)
            ftCseg_dictionary = ftCPoly_splitted.tosegsdict()
            ftCseg_dictionary_ = {}
            for key, value in ftCseg_dictionary.items():
                ftCseg_dictionary_.update({key: [float(v) for v in value[0:4]] + [str(v) for v in value[4:]]})
            util.save(ftCseg_dictionary_, "../out/tmp/ftCsegs.json")  # segment dictionary
        #end
        fittable_segs(Pslice_sqdm,Csegs_splits_dict)

        #collect segments in Pslice_sqdm
        def get_Pslice_sqdm_segs(Pslice_sqdm):
            Pssegs_dictionary = util.load_from_file(infile_Pssegs)  # parent segment
            Pslice_segs_dictionary={}
            for xkey, yblocks in Pslice_sqdm.items():
                for yb, lines in yblocks.items():
                    for linehash in lines.keys():
                        Pslice_segs_dictionary[linehash] = Pssegs_dictionary[linehash]
            return Pslice_segs_dictionary

        util.save_sqdm(Pslice_sqdm, "../out/tmp/Pslicesqdm.json")
        Pslice_segs_dictionary = get_Pslice_sqdm_segs(Pslice_sqdm)
        util.save(Pslice_segs_dictionary, "../out/tmp/Pslicesegs.json")  # segment dictionary
        #return Pslice_segs_dictionary,Pslice_sqdm,ftCseg_dictionary

    @classmethod
    def test_usa_state_delegation(self):

        delegation = Delegation()
        Pslice_segs_dictionary = util.load_from_file("../out/tmp/Pslicesegs.json")
        Pslice_sqdm = util.load_sqdm_from_file("../out/tmp/Pslicesqdm.json")
        ftCseg_dictionary = util.load_from_file("../out/tmp/ftCsegs.json")

        print len(Pslice_sqdm), len(Pslice_segs_dictionary), len(ftCseg_dictionary)
        ##
        '''
        Find a yblock in current sqdm containing a new-segment; and test if a segment is legal to draw.
        Collect all segents in the box (x1,x2)x(y1,y2)
        1. for each s_i in Cpoly:
                #get containing box:(x1,x2)x(y1,y2) from Pslice_sqdm
                #test if s_i is a legal segment that can be inserted in the box; 
                    #s_i does not intersect any segment in the box
                    #s_i completely fall inside the box or 
                    #s_i is on boundary of the box
                #if legal: a) add to a dictionary[x1][(y1,y2)]
                           b) add all segments to a dictionary dictionary[x1][(y1,y2)]     
        '''
        fixables_xkey_yblocks={} #it contains valid partition segments, old segments for each box (x1,x2)x(y1,y2)
        cntnomatch =0
        nomatchsegs = []

        for segkey in ftCseg_dictionary.keys():
            seg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
            #get box (x1,x2)x(y1,y2) containing the segment seg.
            xpair,yblock = delegation.containing_box(seg,Pslice_sqdm)
            if xpair == None:
                nomatchsegs +=[segkey]
                xpair, yblock = delegation.containing_box(seg, Pslice_sqdm)
                cntnomatch +=1
                continue
            xkey,xnkey = xpair

            legalpartition = delegation.islpsegment(seg,Pslice_sqdm[xkey][yblock])
            #any illegal partition line found, return False
            if not legalpartition:
                return False


            Pslice_sqdm[xkey][yblock][segkey] = [seg.co_ordinates()]
            #collect all the segments in the containing box for further divisions.
            if xkey in fixables_xkey_yblocks:
                if yblock not in fixables_xkey_yblocks[xkey]:
                    fixables_xkey_yblocks[xkey][yblock]={}
            else:
                fixables_xkey_yblocks[xkey]={}
                fixables_xkey_yblocks[xkey][yblock]={}

        def non_matching():
            print("cnt-nomatch"), cntnomatch
            nonmatchyvalues =[]
            for segkey in nomatchsegs:
                seg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                x1,y1,x2,y2 = seg.co_ordinates()
                nonmatchyvalues += [y1,y2]

            print("ylen-non matching"),len(nonmatchyvalues)
            nonmatchyvalues = list(set(nonmatchyvalues))
            print nonmatchyvalues


        print("Completed mapping fitable segments to Pslice_sqdm. \n")
        #TODO:Change new segment's abv/bel labels to DELEGATED_TO and Establish equivalent ..
        #TODO: ..relation between DELEGATED_TO and DELEGATED_BY relation.

        cntyblocks =0
        for xkey,yblocks in fixables_xkey_yblocks.items():
            cntyblocks += len(yblocks.keys())
        print("Len(fixable-xkeys),len(fixable-yblocks)"), len(fixables_xkey_yblocks),cntyblocks


        '''
        For segments (new segment + existing segments) in its containing box, develop a template_subsqdm.
        1. for each segment_mesh = fixables_xkey_yblocks[x_i][yblock_j]:
            a) compute a new mesh_sqdm
            b) update Psqdm[x_i][yblock_j] by a reference to new mesh_sqdm
        '''
        #iterate through fixable y-blocks and split all segments inside it.
        print("splitting new segments.")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print yblock
                yblock_segs_keys = Pslice_sqdm[xkey][yblock].keys()
                yblock_segs = []
                print("\t"),yblock_segs_keys

                #iterate through segment's keys.
                for segkey in yblock_segs_keys:
                    try:
                        yblock_segs += [Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)]
                    except:
                        yblock_segs += [Segment().ntuples_to_seg(Pslice_segs_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)]

                yblock_xkeys = delegation.mesh_vertical_sweeplines(yblock_segs)
                yblock_split_keys = {}

                for seg in yblock_segs:
                    # split seg
                    # 1)get x-values at which split must take place
                    li = yblock_xkeys[seg.getLeftPoint().getX()]
                    hi = yblock_xkeys[seg.getRightPoint().getX()]
                    splits = seg.split_at_multiple_x(yblock_xkeys.keys()[li + 1:hi], doround=True)

                    for splitseg in splits:
                        yblock_split_keys[splitseg.dictentry_linehashkey()] = splitseg.dictentry_value()

                Pslice_sqdm[xkey][yblock] = yblock_split_keys
                print("len splits in yblock"), len(yblock_split_keys),len(Pslice_sqdm[xkey][yblock])

            print

        ftCseg_dictionary = None

        #TODO: 1) issue: if a parital vertical slab has only one horizontal line, then x: yblock: is empty.
        #singular_template_sqdm = delegation.xcolumns_yblocks(mesh_splits, mesh_xkeys)
        #rectangles_tuples = delegation.ravel_sqdm(singular_template_sqdm)
        #mappable_splits_rectangles += [(mesh_splits,rectangles_tuples)]

        print("Completed splitting segs in fixable y-blocks")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:"), yblock, "n(segments)"
                for segk,val in Pslice_sqdm[xkey][yblock].items():
                    print("\t\t\t"), segk,val
                print
        print
        ##
        ##Construct template sqdm for each fixable xkey,yblock
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print("xkey:--"),xkey
            for yblock in yblocks.keys():
                print("\t\tyblock:-"),yblock
                yblock_segs_keys = Pslice_sqdm[xkey][yblock].keys()
                yblock_segs = []
                for segkey in yblock_segs_keys:
                    segobj = Segment().ntuples_to_seg(Pslice_sqdm[xkey][yblock][segkey], util.SEGMENT_TUPLE_KEYS)
                    yblock_segs += [segobj]
                    print("\t\t"),segkey, segobj.co_ordinates()

                yblock_xkeys = delegation.mesh_vertical_sweeplines(yblock_segs)
                yblock_template_sqdm = delegation.xcolumns_yblocks(yblock_segs, yblock_xkeys)
                fixables_xkey_yblocks[xkey][yblock] = yblock_template_sqdm #
            print

        print("Completed constructing template-sqdm for fixable-xkey-yblock")
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print xkey
            for yblock in yblocks.keys():
                print("\t\t yblock:--"), yblock
                for tkey in fixables_xkey_yblocks[xkey][yblock].keys():
                    print("--\t\t"), tkey
                    for tyblock in fixables_xkey_yblocks[xkey][yblock][tkey].keys():
                        print("--\t\t\t"), tyblock

        print

        print("--------------------------------------------------------------------")
        print
        #map each of the split segents in PSlice_sqdm fixable xkey,yblock
        for fxkey, fyblocks in fixables_xkey_yblocks.items():
            print fxkey
            for fyblock in fyblocks.keys():
                #get template sqdm for this yblock
                yblock_template_sqdm = fixables_xkey_yblocks[fxkey][fyblock]

                #get fixable segments from PSlice_sqdm
                yblock_segs_keys = Pslice_sqdm[fxkey][fyblock].keys()
                print("keys:"),yblock_segs_keys
                for segkey in yblock_segs_keys:
                    print("..."), segkey, fxkey, yblock
                    seg = Segment().ntuples_to_seg(Pslice_sqdm[fxkey][fyblock][segkey], util.SEGMENT_TUPLE_KEYS)
                    #map this segobj to yblock_template_sqdm
                    #find a containing box
                    x1, y1, x2, y2 = seg.co_ordinates()
                    # map to x-column
                    xcolumn_yblocks = yblock_template_sqdm[x1]
                    for yblock in xcolumn_yblocks.keys():
                        ylow, yhigh = yblock
                        if y1 < y2:
                            if ylow <= y1 and y2 <= yhigh:
                                # TODO: make an ordered dictionary inside.
                                try:
                                    xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1, y1, x2, y2]
                                except:
                                    print("Exception in matching."), seg.co_ordinates()
                        else:
                            if ylow <= y2 and y1 <= yhigh:
                                try:
                                    xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = [x1, y1, x2, y2]
                                except:
                                    print("Exception in matching"), seg.co_ordinates()
        print
        print("Completed mapping segments to template-sqdm for fixable-xkey-yblock")
        for fxkey, fyblocks in fixables_xkey_yblocks.items():
            print("fxkey"),fxkey
            for fyblock in fyblocks.keys():
                print("\t\t fyblock:--"), fyblock
                fixed_sqdm = fixables_xkey_yblocks[fxkey][fyblock]
                Delegation.print_sqdm(fixed_sqdm)

        return fixables_xkey_yblocks

        #merge fixed_sqdm to Psqdm.

#Delegation.pre_process_dummy_polygon()
#Delegation.test_dummy_poly_delegation()
print
print
Delegation.preprocess_us_sqdm()
#Delegation.test_usa_state_delegation()

