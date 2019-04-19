'''This is method-2 for delegation.'''


from SQDM import SQDM, Segment, Point
from collections import OrderedDict
from simple_polygon import util
from ast import literal_eval
from sortedcontainers import SortedDict
import copy
#CONSTANTS

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

    @classmethod
    def test_dummy_poly_delegation(self):
        import time
        delegation = Delegation()
        t0 = time.time()
        Pslice_segs_dictionary = util.load_from_file("../out/tmp/aPslicesegs.json")
        Pslice_sqdm = util.load_sqdm_from_file("../out/tmp/aPslicesqdm.json")
        ftCseg_dictionary = util.load_from_file("../out/tmp/aftCsegs.json")

    @classmethod
    def print_sqdm(self,sqdm, segment_key_dictionary):
        for xkey, yblocks in sqdm.items():
            print("\t \t \t xkey:"), xkey
            for yb, lines in yblocks.items():
                print("\t\t\t\t"), yb
                for lhash, lval in lines.items():
                    print("\t\t\t\t\t"), segment_key_dictionary[lhash]

    @classmethod
    def clip_polygon(self):
        from SQDM import Polygon
        delegation = Delegation()

        Cpoly = SQDM.get_usa_state_boundary_by_name("CALIFORNIA")

        Nsegs = 3110
        Cpoly1 = Cpoly.slice_poly(Nsegs)
        util.poly_ptstoshp(Cpoly1.get_vertices(), "../out/tmp/calif-1000/3Calif-" + str(Nsegs))

        util.save_segments_asshp([Cpoly1.sides()[0].co_ordinates(),Cpoly1.sides()[-1].co_ordinates()], "../out/tmp/" + str(Nsegs)+'--')

        return Cpoly1


    @classmethod
    def preprocess_us_sqdm(self):
        from SQDM import Polygon
        delegation = Delegation()

        subsqdm_file = "../out/tmp/"
        infile_Psqdm = "../out/tmp/2USA_sqdm.json" #sqdm for P #parent polygon sqdm.
        infile_Pssegs = "../out/tmp/2USAs.json" #Split of P segs #parent polygon
        Nsegs = 3110
        Psqdm = util.load_sqdm_from_file(infile_Psqdm) #Parent sqdm
        Cpoly = SQDM.get_usa_state_boundary_by_name("CALIFORNIA")
        Cpoly = Cpoly.slice_poly(Nsegs)
        extent = Cpoly.polygon_extent(Cpoly.sides())
        util.poly_ptstoshp(Cpoly.get_vertices(), "../out/tmp/calif-1000/2Calif-" + str(Nsegs))
        abv_bel_labels = Cpoly.label_bottom_up(util.SQDM_DELEGATED_TO_ID,
                                                 util.SQDM_DELEGATED_BY_ID)  # inner region is delegated.

        Cpoly_dictionary = Cpoly.tosegsdict()
        Cxkeys = Cpoly.vertical_sweeplines().keys()

        #find the part of Psqdm intersected by C
        Pslice_sqdm = delegation.slice_sqdm(Cxkeys[0], Cxkeys[-1], Psqdm)
        #plot only random 50%
        util.save_xun_asshp(sorted(set(Pslice_sqdm.keys())),"../out/tmp/Xun"+str(Nsegs),extent)

        del Psqdm
        del Cpoly

        print("Completed Slicing sqdm for Cpoly. len(Pslice_sqdm)"), len(Pslice_sqdm)
        #save the lengths of the child-segments

        def get_csegs_length(Cpoly_dictionary):
            csegs_lengths = {}
            for csegkey in Cpoly_dictionary.keys():
                segobj = Segment().ntuples_to_seg(Cpoly_dictionary[csegkey], util.SEGMENT_TUPLE_KEYS)
                csegs_lengths[segobj.getAttrByName('edgeid')]= segobj.length()

            util.save(csegs_lengths, "../out/tmp/cseg_lengths" + str(Nsegs) + ".json")  # segment dictionary
            return csegs_lengths

        csegs_lengths =get_csegs_length(Cpoly_dictionary)

        #create assignment vector.
        def create_assignmentvector(Cpoly_dictionary,abv_bel_labels):
            AssignmentVector={}
            for csegkey in Cpoly_dictionary.keys():
                segobj = Segment().ntuples_to_seg(Cpoly_dictionary[csegkey], util.SEGMENT_TUPLE_KEYS)
                AssignmentVector[csegkey] = abv_bel_labels[int(segobj.getAttrByName('edgeid'))]
            return AssignmentVector

        AssignmentVector = create_assignmentvector(Cpoly_dictionary,abv_bel_labels)

        #assign the segment's above and below region.
        def assign_segment_regions(Cpoly_dictionary,AssignmentVector):
            util.Set_Equivalance(util.SQDM_DELEGATED_TO_ID, util.SQDM_DELEGATED_BY_ID)
            for csegkey in Cpoly_dictionary.keys():
                segobj = Segment().ntuples_to_seg(Cpoly_dictionary[csegkey], util.SEGMENT_TUPLE_KEYS)
                abv_label,bel_label = AssignmentVector[csegkey]
                segobj.setAttrByName("abv", abv_label)
                segobj.setAttrByName("bel", bel_label)
                Cpoly_dictionary[csegkey] =segobj.dictentry_value()

        assign_segment_regions(Cpoly_dictionary,AssignmentVector)

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

        ##new code.
        ##new segments are split at x-values unique to CPoly and that is overlapped by CPoly on Psqdm.
        #(Optimal: only at x-values unique to Pslice_sqdm.x values.)

        def split_C_at_CxPx(Cpoly_dictionary,Cxkeys,Pslice_sqdm):
            cseg_mesh = []
            for csegkey in Cpoly_dictionary.keys():
                segobj = Segment().ntuples_to_seg(Cpoly_dictionary[csegkey], util.SEGMENT_TUPLE_KEYS)
                cseg_mesh += [segobj]


            PuCxkeys = Cxkeys + list(Pslice_sqdm.keys())
            Csegs_splits_atx = delegation.split_mesh_segs_at_x(cseg_mesh,
                                                               SQDM().get_ordered_keys(PuCxkeys),
                                                               doround=True)

            Csegs_splits_atx = [seg.dictentry_value() for seg in Csegs_splits_atx]
            Csegs_splits_dict = OrderedDict(enumerate(Csegs_splits_atx))
            #TODO: think about other way to put these splits @x into a dictionary so
            #TODO: that it is easy to collect into buckets while splitting @y.
            #collect segment's hashkey.
            del Csegs_splits_atx
            del Cxkeys
            return Csegs_splits_dict

        Csegs_splits_dict = split_C_at_CxPx(Cpoly_dictionary,Cxkeys,Pslice_sqdm) #split @x

        def fittable_segs(Pslice_sqdm,Csegs_splits_dict):

            #these splits from new segments are split at y-values on each of the vertical slabs.
            #1. collect new-splits for each vertical columns in Slice_sqdm
            seg_buckes =OrderedDict()
            for xkey in list(Pslice_sqdm.keys()):
                seg_buckes[xkey] = []
            print len(seg_buckes)

            #2. partition segments split at x into vertical columns/buckets.
            for seg_id in Csegs_splits_dict.keys():
                #search vcolumn for each seg.
                #add this new-split to the bucket of vcolumn
                segobj = Segment().ntuples_to_seg(Csegs_splits_dict[seg_id], util.SEGMENT_TUPLE_KEYS)
                xlow,xhigh = segobj.getLeftPoint().getX(),segobj.getRightPoint().getX()
                xindex = delegation.xrange_search(list(Pslice_sqdm.keys()),xlow,xhigh)
                xkey = list(Pslice_sqdm.keys())[xindex]
                seg_buckes[xkey].append(seg_id)

            #3. for each segments in vertical buckets/columns, try splitting at y-values
            ftCPolySides = []
            for xkey,seg_ids in seg_buckes.items()[:]:
                y_blocks = Pslice_sqdm[xkey].keys()
                y_values = list(set( yval for tup in y_blocks for yval in tup))

                for seg_id in seg_ids[:]:
                    seg = Segment().ntuples_to_seg(Csegs_splits_dict[seg_id], util.SEGMENT_TUPLE_KEYS)
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
            util.save(ftCseg_dictionary_, "../out/tmp/ftCsegs"+str(Nsegs)+".json")  # segment dictionary
            return ftCseg_dictionary_

        #end
        ftCseg_dictionary= fittable_segs(Pslice_sqdm,Csegs_splits_dict) #split @y

        def get_ftcsegs_length(ftCseg_dictionary):
            ftCseg_lengths ={}
            for segkey in ftCseg_dictionary.keys():
                ftseg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                if ftseg.getAttrByName('edgeid') in ftCseg_lengths:
                    ftCseg_lengths[ftseg.getAttrByName('edgeid')] += ftseg.length()
                else:
                    ftCseg_lengths[ftseg.getAttrByName('edgeid')] = ftseg.length()

            util.save(ftCseg_lengths, "../out/tmp/ftCseg_lengths" + str(Nsegs) + ".json")  # segment dictionary

            return ftCseg_lengths

        get_ftcsegs_length(ftCseg_dictionary)

        #collect segments in Pslice_sqdm
        def get_Pslice_sqdm_segs(Pslice_sqdm,infile_Pssegs):
            Pssegs_dictionary = util.load_from_file(infile_Pssegs)  # parent segment
            Pslice_segs_dictionary={}
            for xkey, yblocks in Pslice_sqdm.items():
                for yb, lines in yblocks.items():
                    for linehash in lines.keys():
                        Pslice_segs_dictionary[linehash] = Pssegs_dictionary[linehash]
            return Pslice_segs_dictionary

        util.save_sqdm(Pslice_sqdm, "../out/tmp/Pslicesqdm"+str(Nsegs)+".json")


        Pslice_segs_dictionary = get_Pslice_sqdm_segs(Pslice_sqdm,infile_Pssegs)
        util.save(Pslice_segs_dictionary, "../out/tmp/Pslicesegs"+str(Nsegs)+".json")  # segment dictionary
        #return Pslice_segs_dictionary,Pslice_sqdm,ftCseg_dictionary

    @classmethod
    def test_usa_state_delegation(self):
        Nsegs=3110
        delegation = Delegation()
        Pslice_segs_dictionary = util.load_from_file("../out/tmp/Pslicesegs"+str(Nsegs)+".json")
        Pslice_sqdm = util.load_sqdm_from_file("../out/tmp/Pslicesqdm"+str(Nsegs)+".json")
        ftCseg_dictionary = util.load_from_file("../out/tmp/ftCsegs"+str(Nsegs)+".json")

        print("len(pslice_sqdm,len(pslice_segs_dic),len(ftcsegs_dic"), len(Pslice_sqdm), \
            len(Pslice_segs_dictionary), \
            len(ftCseg_dictionary)
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
         #it contains valid partition segments, old segments for each box (x1,x2)x(y1,y2)

        #common = set( ftCseg_dictionary.keys()).intersection(set(Pslice_segs_dictionary.keys()))
        #print("common between P and C"), len(common), len(set(ftCseg_dictionary.keys()))

        #function
        def ft_cseg_mapping(ftCseg_dictionary,Pslice_sqdm):
            i =0
            j = 0
            nomatchsegs = []
            fixables_xkey_yblocks = {}
            for segkey in ftCseg_dictionary.keys():
                seg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                #get box (x1,x2)x(y1,y2) containing the segment seg.
                xpair,yblock = delegation.containing_box(seg,Pslice_sqdm)

                if xpair == None:
                    nomatchsegs +=[segkey]
                    if segkey in ["06f8461feae2","36e410aabf5b","56eda1ac83f2"]:
                        xpair, yblock = delegation.containing_box(seg, Pslice_sqdm)
                    continue
                xkey,xnkey = xpair

                legalpartition = delegation.islpsegment(seg, Pslice_sqdm[xkey][yblock])
                #any illegal partition line found, return False
                if not legalpartition:
                    return False

                #? happens if there are duplicate/overlapping segments between P and C.
                if segkey in Pslice_sqdm[xkey][yblock]:
                    existing_seg = Segment().ntuples_to_seg(Pslice_segs_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                    new_child_seg= Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                    #swap the region code.
                    area_abv,area_bel = existing_seg.getAttrByName('abv'), existing_seg.getAttrByName('bel')
                    if area_abv in [util.SQDM_DELEGATED_BY_ID, util.COUNTRY_DELEGATED_BY_ID]:
                        existing_seg.setAttrByName('abv', new_child_seg.getAttrByName('abv'))
                    if area_bel in [util.SQDM_DELEGATED_BY_ID, util.COUNTRY_DELEGATED_BY_ID]:
                        existing_seg.setAttrByName('bel', new_child_seg.getAttrByName('bel'))
                    ftCseg_dictionary[segkey] = existing_seg.dictentry_value()
                    #update
                    Pslice_sqdm[xkey][yblock][segkey] = util.CHILD_ANNOTATION #mapped as child
                    i +=1
                #collect all the segments in the containing box for further divisions.
                else:
                    #insert
                    Pslice_sqdm[xkey][yblock][segkey] = util.CHILD_ANNOTATION
                    j +=1
                if xkey in fixables_xkey_yblocks:
                    if yblock not in fixables_xkey_yblocks[xkey]:
                        fixables_xkey_yblocks[xkey][yblock]={}
                else:
                    fixables_xkey_yblocks[xkey]={}
                    fixables_xkey_yblocks[xkey][yblock]={}

            non_matching(nomatchsegs)
            print("Unique segs in C"),j
            print("Common segs between C and P"),i
            return fixables_xkey_yblocks

        def non_matching(nomatchsegs):
            print("non-matching c-seg key"),nomatchsegs
            nonmatchyvalues =[]
            for segkey in nomatchsegs:
                seg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                x1,y1,x2,y2 = seg.co_ordinates()
                nonmatchyvalues += [y1,y2]

            print("non-matching csegs"),len(nomatchsegs)
            nonmatchyvalues = list(set(nonmatchyvalues))
        #end

        fixables_xkey_yblocks = ft_cseg_mapping(ftCseg_dictionary, Pslice_sqdm)

        print("Completed mapping fitable segments to Pslice_sqdm. \n")

        #TODO:Change new segment's abv/bel labels to DELEGATED_TO and Establish equivalent ..
        #TODO: ..relation between DELEGATED_TO and DELEGATED_BY relation.

        #Function
        def count_fixable_blocsk(fixables_xkey_yblocks):
            cntyblocks = 0
            for xkey,yblocks in fixables_xkey_yblocks.items():
                cntyblocks += len(yblocks.keys())
            print("Len(fixable-xkeys),len(fixable-yblocks)"), len(fixables_xkey_yblocks),cntyblocks

        count_fixable_blocsk(fixables_xkey_yblocks)

        '''
            For segments (new segment + existing segments) in its containing box, develop a template_subsqdm.
            1. for each segment_mesh = fixables_xkey_yblocks[x_i][yblock_j]:
                a) compute a new mesh_sqdm
                b) update Psqdm[x_i][yblock_j] by a reference to new mesh_sqdm
        '''

        #function
        #iterate through fixable y-blocks and split all segments inside it.
        def in_place_split_fixable_sqdm(fixables_xkey_yblocks, Pslice_sqdm, ftCseg_dictionary, Pslice_segs_dictionary):
            cnt = 0 #new segs
            cnt2 = 0 #old segs
            global_splits_dictionary = {}
            for xkey, yblocks in fixables_xkey_yblocks.items():

                for yblock in yblocks.keys():
                    yblock_segs_keys_values = Pslice_sqdm[xkey][yblock].items()
                    yblock_segs = []
                    #iterate through segment's keys.
                    for segkey,segval in yblock_segs_keys_values:
                        if segval == util.CHILD_ANNOTATION:
                            cnt +=1
                            childseg = Segment().ntuples_to_seg(ftCseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)
                            childseg.setAttrByName("ischild",True)
                            yblock_segs += [childseg]
                        else:
                            cnt2 +=1
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
                            yblock_split_keys[splitseg.dictentry_linehashkey()] = 1
                            global_splits_dictionary[splitseg.dictentry_linehashkey()] = splitseg.dictentry_value()

                    Pslice_sqdm[xkey][yblock] = yblock_split_keys
            #end for
            print("Completed splitting segs in fixable y-blocks")
            print("\t n(ftCsegs), n(splitted/mapped-ftCsegs)"), len(ftCseg_dictionary), cnt
            print("\t splited n(old-segments"), cnt2
            print("\t total splits"),len(global_splits_dictionary)
            return global_splits_dictionary

        global_splits_dictionary = in_place_split_fixable_sqdm(fixables_xkey_yblocks,
                                                               Pslice_sqdm,
                                                               ftCseg_dictionary,
                                                               Pslice_segs_dictionary)


        ftCseg_dictionary = None
        Pslice_segs_dictionary = None

        #TODO: 1) issue: if a parital vertical slab has only one horizontal line, then x: yblock: is empty.
        #singular_template_sqdm = delegation.xcolumns_yblocks(mesh_splits, mesh_xkeys)
        #rectangles_tuples = delegation.ravel_sqdm(singular_template_sqdm)
        #mappable_splits_rectangles += [(mesh_splits,rectangles_tuples)]

        def print_stats_fixables(fixables_xkey_yblocks,Pslice_sqdm):
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
        def template_sqdm_each_fixables(fixables_xkey_yblocks,Pslice_sqdm,global_splits_dictionary):

            for xkey, yblocks in fixables_xkey_yblocks.items():
                for yblock in yblocks.keys():
                    yblock_segs_keys = Pslice_sqdm[xkey][yblock].keys()
                    yblock_segs = []

                    for segkey in yblock_segs_keys:
                        segobj = Segment().ntuples_to_seg(global_splits_dictionary[segkey], util.F2_SEGMENT_TUPLE_KEYS)
                        yblock_segs += [segobj]

                    yblock_xkeys = delegation.mesh_vertical_sweeplines(yblock_segs)
                    yblock_template_sqdm = delegation.xcolumns_yblocks(yblock_segs, yblock_xkeys)
                    fixables_xkey_yblocks[xkey][yblock] = yblock_template_sqdm #
            return  fixables_xkey_yblocks

        fixables_xkey_yblocks = template_sqdm_each_fixables(fixables_xkey_yblocks, Pslice_sqdm, global_splits_dictionary)
        print("Completed constructing template-sqdm for fixable-xkey-yblock")

        #function
        def printf_fixables_xkey_yblocks():
            for xkey, yblocks in fixables_xkey_yblocks.items():
                print xkey
                for yblock in yblocks.keys():
                    print("\t\t yblock:--"), yblock
                    for tkey in fixables_xkey_yblocks[xkey][yblock].keys():
                        print("--\t\t"), tkey
                        for tyblock in fixables_xkey_yblocks[xkey][yblock][tkey].keys():
                            print("--\t\t\t"), tyblock

        #map each of the split segents in global_splits_dictionary corresponding templates in fixable xkey, yblock
        def map_global_splits_to_templatesqdm(fixables_xkey_yblocks,Pslice_sqdm, global_splits_dictionary):
            for fxkey, fyblocks in fixables_xkey_yblocks.items():

                for fyblock in fyblocks.keys():
                    #get template sqdm for this yblock
                    yblock_template_sqdm = fixables_xkey_yblocks[fxkey][fyblock]

                    #get fixable segments from PSlice_sqdm
                    yblock_segs_keys = Pslice_sqdm[fxkey][fyblock].keys()

                    for segkey in yblock_segs_keys:
                        seg = Segment().ntuples_to_seg(global_splits_dictionary[segkey], util.F2_SEGMENT_TUPLE_KEYS)
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
                                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = 1 #mapped
                                    except:
                                        print("Exception in matching."), seg.co_ordinates()
                            else:
                                if ylow <= y2 and y1 <= yhigh:
                                    try:
                                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = 1 #mapped.
                                    except:
                                        print("Exception in matching"), seg.co_ordinates()
            print

        map_global_splits_to_templatesqdm(fixables_xkey_yblocks, Pslice_sqdm, global_splits_dictionary)
        print("Completed mapping segments to template-sqdm for fixable-xkey-yblock")

        def print_after_mapping():
            for fxkey, fyblocks in fixables_xkey_yblocks.items():
                print("fxkey"),fxkey
                for fyblock in fyblocks.keys():
                    print("\t\t fyblock:--"), fyblock
                    fixed_sqdm = fixables_xkey_yblocks[fxkey][fyblock]
                    Delegation.print_sqdm(fixed_sqdm,global_splits_dictionary)
                print

        #print_after_mapping()
        #
        def delegation_test_by_reconstruction(fixables_xkey_yblocks):
            '''
            Given a parent polygon, P, and a child polygon C, this method delegated C out of P.

            extract segments mapped to this set of sqdm's after delegation, and plot
            the graph to match the original child-polygon.
            :param fixables_xkey_yblocks:
            :return:
            '''
            tdelegated_segs =[]
            for fxkey, fyblocks in fixables_xkey_yblocks.items():
                for fyblock in fyblocks.keys():
                    fixed_sqdm = fixables_xkey_yblocks[fxkey][fyblock]
                    for xkey, yblocks in fixed_sqdm.items():
                        for yb, lines in yblocks.items():
                            for lhash, lval in lines.items():
                                tdelegated_segs += [global_splits_dictionary[lhash][0:4]]

            print("Colllected segments # from delegated Pslice sqdm."), len(tdelegated_segs)
            util.save_segments_asshp(tdelegated_segs,"../out/tmp/"+str(Nsegs))

        def delegation_test_distance_equivalance(fixables_xkey_yblocks):
            '''
            Given a parent polygon, P, and a child polygon C, this method delegated C out of P.

            extract segments mapped to this set of sqdm's after delegation, and plot
            the graph to match the original child-polygon.
            :param fixables_xkey_yblocks:
            :return:
            '''
            delegated_csegs_lengths ={} #<segment-id,segment-distance>
            for fxkey, fyblocks in fixables_xkey_yblocks.items():
                for fyblock in fyblocks.keys():
                    fixed_sqdm = fixables_xkey_yblocks[fxkey][fyblock]
                    for xkey, yblocks in fixed_sqdm.items():
                        for yb, lines in yblocks.items():
                            for segkey, lval in lines.items():

                                delegated_seg = Segment().ntuples_to_seg(global_splits_dictionary[segkey],
                                                                         util.F2_SEGMENT_TUPLE_KEYS)

                                #check if this segment is child.
                                if delegated_seg.getAttrByName('ischild') != '':
                                    if delegated_seg.getAttrByName('edgeid') in delegated_csegs_lengths:
                                        delegated_csegs_lengths[delegated_seg.getAttrByName('edgeid')] += delegated_seg.length()
                                    else:
                                        delegated_csegs_lengths[delegated_seg.getAttrByName('edgeid')] = delegated_seg.length()

            util.save(delegated_csegs_lengths, "../out/tmp/delegated_csegs_lengths" + str(Nsegs) + ".json")  # segment dictionary
            return delegated_csegs_lengths

            print("Colllected segments # from delegated Pslice sqdm."), len(tdelegated_segs)
            util.save_segments_asshp(tdelegated_segs,"../out/tmp/"+str(Nsegs))
        delegation_test_by_reconstruction(fixables_xkey_yblocks)
        delegation_test_distance_equivalance(fixables_xkey_yblocks)

        return fixables_xkey_yblocks

        #merge fixed_sqdm to Psqdm.

#Delegation.pre_process_dummy_polygon()
#Delegation.test_dummy_poly_delegation()
print
print
#Delegation.clip_polygon()

Delegation.preprocess_us_sqdm()
Delegation.test_usa_state_delegation()

