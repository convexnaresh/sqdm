from SQDM import SQDM, Segment
from collections import OrderedDict
from simple_polygon import util
from ast import literal_eval
from sortedcontainers import SortedDict
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

    def containing_box(self,segment_obj,sqdm):
        x1, y1, x2, y2 = segment_obj.co_ordinates()
        xkyes = list(sqdm.keys())
        # search along vertical slabs (x_i -- x_j)
        xsearchidx = self.xrange_search(xkyes, x1, x2)
        if xsearchidx < 0:
            return {}
        else:
            # search along yblocks (y_i -- y_j)
            ykeys_tuples = sqdm[xkyes[xsearchidx]].keys()
            ysearchidx = self.yrange_search(ykeys_tuples, y1, y2)
            if ysearchidx < 0:
                return {}
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
        # range search
        arr += arr[0:1]
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
            elif arr[mid] < search_key1:
                l = mid + 1

            # If x is smaller, ignore right half
            else:
                r = mid - 1

        # If we reach here, then the element was not present
        return -1


    def split_mesh_segs_at_x(self,mesh_segments, xunikdic):

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
            print xkey, list_yblocks
            ypairs_dic = OrderedDict()
            list_yspans_pairs =util.pairwise(util.non_overlaping_yspans(list_yblocks))

            count=0
            for ypair in list_yspans_pairs:
                ypairs_dic[str(ypair)] ={}
                count+=1

            xcolumn_yblock_dic[xkey] =ypairs_dic

        return xcolumn_yblock_dic

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
    def test_delegation(self):
        from SQDM import Polygon
        delegation = Delegation()
        subsqdm_file = "../out/tmp/"
        infile_Psqdm = "../out/tmp/Psqdm.json" #sqdm for P #parent polygon sqdm.
        infile_Pssegs = "../out/tmp/Ps.json" #Split of P segs #parent polygon
        Psqdm = util.load_sqdm_from_file(infile_Psqdm) #Parent sqdm
        Pssegs_dictionary = util.load_from_file(infile_Pssegs)  # parent segment

        Cl = [(3.6, 4.8), (3.8, 4.5),(3.6,4.2), (3.2,4),
                (3.6,3.8),(3.8,3.6), (3.2,3.2),(3,3.6),(2.6,3.4),(2.8,3.8),(2.6,4), (3,4.2),(3.2,4.5)]

        Cpoly = Polygon(Cl)
        AV=Cpoly.label_bottom_up(util.POLYGON_INNERID)
        Cseg_dictionary = Cpoly.tosegsdict()

        #1. VerifySimplePolygon(CPoly) : it must be done  prior to adding segments into an existing sqdm.
        #it is because, an intersecting segment forces to check it's intersection with all other segments.

        #since Cpoly is completely inside Polygon P, region abv and bel is
        #both assigned to be same

        for C_seg in Cpoly.sides():
            C_seg.setAttrByName('abv',SQDM_DELEGATED_BY_ID)
            C_seg.setAttrByName("bel",SQDM_DELEGATED_BY_ID)

        #find the part of Psqdm intersected by C .
        Cxkeys = Cpoly.vertical_sweeplines().keys()
        lowkey= Cxkeys[0]
        highkey = Cxkeys[-1]
        Pslice_sqdm = delegation.slice_sqdm(lowkey, highkey, Psqdm)
        print("Completed Slicing sqdm for Cpoly..\n")

        '''
        Find a containing yblock in current sqdm and test if a segment is legal to draw.
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
        for seg in Cpoly.sides()[:]:
            #get box (x1,x2)x(y1,y2) containing the segment seg.
            xpair,yblock = delegation.containing_box(seg,Pslice_sqdm)
            xkey,xnkey = xpair

            legalpartition = delegation.islpsegment(seg,Pslice_sqdm[xkey][yblock])
            #any illegal partition line found, return False
            if not legalpartition:
                return False

            #collect all the segments in the containing box for further divisions.
            if xkey in fixables_xkey_yblocks:
                if yblock in fixables_xkey_yblocks[xkey]:
                    fixables_xkey_yblocks[xkey][yblock][seg.dictentry_linehashkey()]=[seg.co_ordinates()]
                else:
                    fixables_xkey_yblocks[xkey][yblock]={}
                    fixables_xkey_yblocks[xkey][yblock][seg.dictentry_linehashkey()] = [seg.co_ordinates()]
            else:
                fixables_xkey_yblocks[xkey]={}
                fixables_xkey_yblocks[xkey][yblock]={}
                fixables_xkey_yblocks[xkey][yblock][seg.dictentry_linehashkey()] = [seg.co_ordinates()]

            #collect segments from Pslice_sqdm for further splitting.
            fixables_xkey_yblocks[xkey][yblock].update(Pslice_sqdm[xkey][yblock])

        print("Completed constructing legal-fixable-xkey-yblocks. \n")
        '''
        For segments (new segment + existing segments) in its containing box, develop a template_subsqdm.
        1. for each segment_mesh = fixables_xkey_yblocks[x_i][yblock_j]:
            a) compute a new mesh_sqdm
            b) update Psqdm[x_i][yblock_j] by a reference to new mesh_sqdm
        '''
        for xkey, yblocks in fixables_xkey_yblocks.items():
            print("-xkey"), xkey,
            for yb, segs in yblocks.items():
                print("\t -yb"), yb
                mesh_segs_keys = segs.keys()
                singular_mesh_segs = []
                for segkey in mesh_segs_keys:
                    try:
                        singular_mesh_segs += [Segment().ntuples_to_seg(Cseg_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)]
                    except:
                        singular_mesh_segs += [Segment().ntuples_to_seg(Pssegs_dictionary[segkey], util.SEGMENT_TUPLE_KEYS)]

                # a) construct a sqdm out of these segments.

                mesh_xkeys = delegation.mesh_vertical_sweeplines(singular_mesh_segs)

                mesh_splits = delegation.split_mesh_segs_at_x(singular_mesh_segs, mesh_xkeys)

                #TODO: 1) issue: if a parital vertical slab has only one horizontal line, then x: yblock: is empty.
                singular_template_sqdm = delegation.xcolumns_yblocks(mesh_splits, mesh_xkeys)
                print singular_template_sqdm
                #b)
                sqdm_id = delegation.sqdm_storageid(xkey,yb)

                Pslice_sqdm[xkey][yb] = sqdm_id #singular_template_sqdm #put id
                util.save(singular_template_sqdm, subsqdm_file+sqdm_id)
                print
        print

        #TODO: Segment's mapping to singlular_template_sqdm for a box.

        #save slice sqdm
        sqdm_id = delegation.sqdm_storageid(0,(0,0))
        util.save_sqdm(Pslice_sqdm, subsqdm_file + sqdm_id)

        for xkey, yblocks in Pslice_sqdm.items():
            print(":=xkey"), xkey
            if yblocks:
                for yb, sqdm_segs in yblocks.items():
                    try:
                        print("\t"), yb, sqdm_segs.keys()
                    except:
                        sqdm_id = delegation.sqdm_storageid(xkey,yb)
                        print("\t"), util.load_from_file(subsqdm_file+sqdm_id)


Delegation.test_delegation()











