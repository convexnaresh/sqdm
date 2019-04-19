import sys
import math
from ast import literal_eval
from simple_polygon import util
from datatypes import S, STYPE, EVENT_TYPE
from collections import  OrderedDict

#CONSTANTS

Left = 0
Right = 1
EVENTNAME = {0: "Start", 1: "End"}
SEGMENT_ID_TYPES = {'linehash': 'linehash', 'lineid': 'lineid', 'edgeid': 'edgeid'}
debug = False

class Point:
    def __init__(self, tup):
        self.x, self.y = tup

    def __eq__(self, other):
        return (self.x,self.y) == (other.x, other.y)

    def __add__(self, other):
        return (self.x,self.y) + (other.x, other.y)

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def isLeft(self, ostartpoint, oendpoint):
        # if self point is left of line joing two points
        # ostartpoint and oendpoint
        return (ostartpoint.x - self.x) * (oendpoint.y - self.y) - \
               (oendpoint.x - self.x) * (ostartpoint.y - self.y)

    def areatriangle(self, op1, op2):
        return self.isLeft(op1, op2)

    def xyorder(self, other):
        # returns if this point is left to other point
        if self.x > other.x:
            return 1
        if self.x < other.x:
            return -1
        if self.y > other.y:
            return 1
        if self.y < other.y:
            return -1
        return 0  # for all cases.

    def __str__(self):
        return str((self.x, self.y))

    def totuple(self):
        return self.x, self.y

class Segment(object):

    def __init__(self, endpoint1=None, endpoint2=None, attr={}, dolrorder=True):
        self.lp = endpoint1
        self.rp = endpoint2
        self.attr = attr #attributes for this segment is a dictionary
        self.attr['isswapped'] = 0
        if dolrorder:
            try:
                if self.orderleftright():
                    self.attr['isswapped'] = 1 #original dir changed.
            except:
                pass
        self.attr['islrordered'] = dolrorder

    def __str__(self):
        s = ''
        s += str(self.lp) + "--" + str(self.rp) + str(self.attr)
        return s

    def ntuples_to_seg(self,ntuplevalues,tuplekeys):
        if len(ntuplevalues)< 4 or (len(tuplekeys) < len(ntuplevalues)):
            return None

        attr = {}
        lpoint = Point(ntuplevalues[0:2])
        rpoint = Point(ntuplevalues[2:4])
        id = 0

        for i in range(4,len(ntuplevalues)):
            key = tuplekeys[i]
            attr[key] = ntuplevalues[i]

        segObj = Segment(rpoint, lpoint, attr)
        return segObj

    ##getters
    def getEdgeid(self):
        return self.getAttrByName('edgeid')

    def getLeftPoint(self):
        return self.lp

    def getRightPoint(self):
        return self.rp

    def getAttrByName(self,name):
        try:
            return self.attr[name]
        except:
            return ''

    def setAttrByName(self,name,value):
        self.attr[name] = value

    def delx(self):
        return self.getRightPoint().getX() - self.getLeftPoint().getX()

    def dely(self):
        return self.getRightPoint().getY() - self.getLeftPoint().getY()

    def length(self):
        return math.hypot(self.delx(),self.dely())

    def setAttr(self,name,value):
        self.attr[name] = value

    def isLRordered(self):
        return self.islrordered

    def isleftright(self):
        '''returns -1 if this segment is represented as left to right.
         Else 1; else 0 if two points coincides.'''
        x1,y1,x2,y2 = self.co_ordinates()
        if x1 < x2:
            return -1
        elif x1 == x2:
            if y1 < y2:
                return -1
            elif y1 == y2:
                return 0
        return 1

    def orderleftright(self):
        if self.isleftright() == 1:
            self.lp, self.rp = self.rp, self.lp
            return True
        return False

    def origco_ordinates(self):
        if self.getAttrByName('isswapped'):
            return self.rp, self.lp
        return self.lp, self.rp

    def co_ordinates(self):
        return self.lp.x, self.lp.y, self.rp.x, self.rp.y

    # get m for a given seg in (p1, p2) form
    def get_slope(self):
        x0, y0, x1, y1 = self.co_ordinates()
        if (x1 - x0) == 0:
            return None
        else:
            return float(y1 - y0) / (x1 - x0)

    # given a point p, return the point on s that shares p's y-val
    def get_x_at(self, p, doround=False):

        m = self.get_slope()

        # for now that it would have been deleted already if not
        if m == 0:  # horizontal segment
            return None
        # ditto; should check if y-val on seg
        if m is None:  # vertical segment
            return Point((self.lp.x, p[1]))

        if doround:
            x1 = int(round(self.lp.x - (self.lp.y - p[1]) / m))
        else:
            x1 = self.lp.x - (self.lp.y - p[1]) / m

        # this should check if p's x-val is actually on seg; we're assuming
        if self.lp.x <= x1 <= self.rp.x:
            return Point((x1, p[1]))
        return None

    # given a point p, return the point on s that shares p's x-val
    def get_y_at_x(self, x1,doround=False):
        '''this is wrong. correct it.'''
        m = self.get_slope()

        # ditto; should check if y-val on seg
        if m is None:  # vertical segment
            return None

        y1 = m * (x1 - self.rp.x) + self.rp.y  # yn = float(y2 - y1) / (x2 - x1) * (xn - x1) + y1
        if doround:
            y1 = int(round(y1))

        return Point((x1,y1))

    def split_at_x(self,xvalue,doround=False):
        import copy
        splitpt = self.get_y_at_x(xvalue,doround)

        #split point is none or split point is one of the end point of a segment.
        if (xvalue <= self.getLeftPoint().getX()) or (xvalue >= self.getRightPoint().getX()):
            return None,None

        if splitpt is None or (splitpt == self.lp or splitpt == self.rp):
            return None,None


        #left split
        lseg = copy.deepcopy(self)
        lseg.rp = splitpt
        lseg.setAttr('edgeid',lseg.getAttrByName('edgeid'))#+"."+str(1))
        #right split
        rseg = copy.deepcopy(self)
        rseg.lp = splitpt
        rseg.setAttr('edgeid', rseg.getAttrByName('edgeid')) #+ "." + str(2))
        return lseg,rseg

    def split_at_y(self,ypoint,doround=False):
        import copy
        splitpt = self.get_x_at(ypoint,doround)
        if splitpt is None or (splitpt == self.lp or splitpt == self.rp):
            return None,None

        #left split
        lseg = copy.deepcopy(self)
        lseg.rp = splitpt
        lseg.setAttr('edgeid',str(lseg.getAttrByName('edgeid')))#+"."+str(1))
        #right split
        rseg = copy.deepcopy(self)
        rseg.lp = splitpt
        rseg.setAttr('edgeid', str(rseg.getAttrByName('edgeid'))) #+ "." + str(2))
        return lseg,rseg

    def split_at_multiple_x(self,xvalue_list,doround=False):

        splits =[]
        curseg = self
        for xval in xvalue_list:
            segl,segr = curseg.split_at_x(xval,doround)
            if segl is not None:
                splits.append(segl)
                curseg = segr
        #append last part.
        splits.append(curseg)

        if debug:
            print("splitting "), str(self)
            for sp in splits:
                print("\t split"),str(sp)

        #reverse lists if necessary.
        if self.getAttrByName('isswapped') ==1:
            splits.reverse()
        return splits

    def split_at_multiple_y(self,yvalue_list,doround=False):

        splits =[]
        curseg = self
        for yval in yvalue_list:
            segl,segr = curseg.split_at_y((0,yval),doround)
            if segl is not None:
                splits.append(segl)
                curseg = segr
        #append last part.
        splits.append(curseg)

        if debug:
            print("splitting "), str(self)
            for sp in splits:
                print("\t split"),str(sp)

        #reverse lists if necessary.
        if self.getAttrByName('isswapped') ==1:
            splits.reverse()
        return splits

    def dictentry_linehashkey(self):
        # hashes for (1,11)--(11,1) is same as (11,11)-(1,1) for two different segments
        #TODO : handle hash for duplicating points as shown above
        x1, y1, x2, y2 = self.co_ordinates()
        key = util.hashargs(str(int(x1)), str(int(y1)), str(int(x2)), str(int(y2)))
        return key

    def dictentry_value(self):
        t = self.getLeftPoint() + self.getRightPoint() + (self.attr['edgeid'],)
        try:
            t+= (self.attr['abv'],self.attr['bel'])
        except:
            pass

        try:
            t+= (self.attr['ischild'],)
        except:
            pass

        return t

    def tupleToSegment(self, tup):
        '''tup is 4 tuple segment'''
        x1, y1, x2, y2 = tup
        p1 = Point((x1, y1))
        p2 = Point((x2, y2))
        if p1.xyorder(p2) < 0:
            self.lp = p1
            self.rp = p2
        else:
            self.lp = p2
            self.rp = p1
        return self

class Polygon(object):
    import hashlib
    '''A polygon object can be constructed in three ways:
    1) po = Polygon([(x1,y1),(x2,y2)...(xn,yn)]), where the list of tuples represents a points (x,y) in the polygon.
    2) po = Polygon([po1,po2,po3..pon]), where list of Point objects represents points in the polygon. 
    3) po = Polygon([])'''
    def __init__(self, list_pts=[],pid=None):
        # polygon is set of Points
        if len(list_pts) >=3:
            #polygon with Point obj as vertices
            if isinstance(list_pts[0], Point):
                self.vertices = list_pts
            #polygon with Segment obj as sides.
            if isinstance(list_pts[0],Segment):
                self._sides = list_pts

            if isinstance(list_pts[0], tuple):
                #point tuple (x,y)
                if len(list_pts[0]) == 2:
                    self.vertices = [Point(tup) for tup in list_pts]

                #segment tuple (x1,y1,x2,y2,other-attributes)
                if len(list_pts[0]) >= 4:
                    self.vertices = []
                    for idx in xrange(0, len(list_pts), 2):
                        x1, y1, x2, y2 = list_pts[idx][0:4]
                        self.vertices += [Point((x1, y1)), Point((x2, y2))]
                    if len(list_pts) % 2 == 0:  # n is even
                        self.vertices += [self.vertices[0]]
            self.makeVertexLoop()
        else:
            self.vertices =[]
        self.nv = len(self.vertices) - 1
        #properties
        self.properties={}
        self._sides =[]
        self._sides_attr=[]

    def __len__(self):
        if self.nv < 0:
            return len(self._sides)

    def __str__(self):
        st = ''
        if self.vertices:
            for ptobj in self.vertices:
                st += str(ptobj)
            return "Polygon:" + st

        if self.getSides():
            for side in self.getSides():
                st += str(side)
            return "Polygon:" + st

    def get_vertices(self):
        '''
        :return: a list of vertices ordered in the same order of segments/vertices in this polygon.
        '''
        if self.vertices:
            return [(p.getX(),p.getY()) for p in self.vertices]
        if self.sides():
            vertices =[]
            segments = self.sides()
            for idx in xrange(0,len(segments),2):
                x1, y1, x2, y2 = segments[idx].co_ordinates()
                if segments[idx].getAttrByName("isswapped"):
                    vertices +=[(x2,y2),(x1,y1)]
                else:
                    vertices +=[(x1,y1),(x2,y2)]
            if len(segments) % 2 == 0: #len is even
                vertices +=[vertices[0]]
            return vertices

        return []

    def setSidesAttr(self,list_attr_dict):
        '''one to one correspoindance from side of this polygon to attribute dictionary in the
        list list_attr_dict'''
        self._sides_attr = list_attr_dict

    def getSidesAttr(self):
        return self._sides_attr

    def addSegment(self,segmentObj):
        self._sides += [segmentObj]

    def addVertex(self,vertex_point):
        self.vertices += [vertex_point]

    def setSides(self,sides_seg_objects):
        self._sides = sides_seg_objects
        self.setPropertyByName('n(sides)', len(sides_seg_objects))

    def getSides(self):
        return self._sides

    def polygon_extent(self,segs=None):
        xvals = []
        yvals = []
        if segs == None:
            segs = self.sides()
        for t in segs:
            x1,y1,x2,y2 = t.co_ordinates() #xvals += [t[0], t[2]]
            xvals +=[x1,x2]
            yvals += [y1,y2]

        xmax, xmin, ymax, ymin = max(xvals), min(xvals), max(yvals), min(yvals)
        return [(xmin,ymin),(xmax,ymax)]

    def slice_poly(self,nsegments):
        '''get only nsegments of this polygon.'''
        slicePoly = Polygon([])
        vertices = self.get_vertices()[290:290+nsegments+1]
        return Polygon(vertices)


    def compute_properties(self):
        minx,maxx = self.minx_maxx()
        self.setPropertyByName('minx',minx)
        self.setPropertyByName('maxx',maxx)
        self.setPropertyByName('mainy',0)
        self.setPropertyByName('maxy',0)

        if self.vertices:
            self.properties['n(v)'] = len(self.vertices) - 1
            #other properties
            for i in range(0, len(self.vertices) - 1):
                seg = Segment(self.vertices[i], self.vertices[i+1])
                self.setPropertyByName('seg_delx', round(seg.delx(),2))
                self.setPropertyByName('seg_dely', round(seg.dely(),2))
                self.setPropertyByName('seg_lengths', round(seg.length(),2))
        elif self.getSides():
            self.properties['n(v)'] = len(self.getSides())
            for seg in self.getSides():
                self.setPropertyByName('seg_delx', round(seg.delx(),2))
                self.setPropertyByName('seg_dely', round(seg.dely(),2))
                self.setPropertyByName('seg_lengths', round(seg.length(),2))

        return self.properties

    def minx_maxx(self):
        try:
            xvalues = self.vertical_sweeplines().keys()
            return xvalues[0],xvalues[-1]
        except:
            return None,None

    def sides(self):
        '''returns a set of Segment objects which are sides of this polygon.'''

        #if this polygon has sides.
        if self.getSides():
            #update with attributes
            return self._sides

        #if this polygon has vertices
        seglist = []
        for i in range(0, len(self.vertices) - 1):
            dattr = {'edgeid':i}
            try:
                dattr= self.getSidesAttr()[i]

            except:
                pass
            seg = Segment(self.vertices[i], self.vertices[i+1],dattr)
            seglist.append(seg)

        if debug:
            print("sides:")
            for s in seglist:
                print str(s)
            print("\n")
        return seglist

    def setPropertyByName(self,name, value):

        if name in ['seg_lengths','split_counts', 'seg_delx','seg_dely']: #number of splits for each segment
            try:
                self.properties[name] += [value]
            except:
                self.properties[name]=[]
                if type(value) == list:
                    self.properties[name] += value
                else:
                    self.properties[name] += [value]

        else:
            self.properties[name] = value

    def getPropertyByName(self,name):
        try:
            return self.properties[name]
        except:
            return []

    def setProperties(self,properties):
        self.properties = properties

    def getProperties(self):
        return self.properties

    def isVertexLoop(self):

        if self.vertices[-1] == self.vertices[0]:
            return True
        return False

    def makeVertexLoop(self):
        if not self.isVertexLoop():
            self.vertices.append(self.vertices[0])

    def isSimple(self):
        return True

    def tosegsdict(self, segment_id_type=SEGMENT_ID_TYPES['linehash']):
        '''returns a set of items of this form for each segment:
        (hash(x1,y1,x2,y2),(x1,y1,x2,y2,AreaAbove,AreaBelow,etc.))'''
        seg_dictionary = OrderedDict()
        if self.vertices:
            for i in range(0, len(self.vertices) - 1):
                p1, p2 = self.vertices[i], self.vertices[i + 1]
                try:
                    dattr = self.getSidesAttr()[i]
                except:
                    dattr ={}
                dattr['edgeid'] = i

                seg = Segment(p1, p2,dattr)#attr for this side.

                if segment_id_type == SEGMENT_ID_TYPES['linehash']:
                    #entry = {seg.dictentry_linehashkey(): seg.dictentry_value()}
                    try:
                        #checking segment collision.
                        seg_dictionary[seg.dictentry_linehashkey()]
                    except:
                        seg_dictionary[seg.dictentry_linehashkey()] = seg.dictentry_value()

                if segment_id_type in [SEGMENT_ID_TYPES['lineid'], SEGMENT_ID_TYPES['lineid']]:

                    #entry = {seg.attr['edgeid'] : seg.dictentry_value()}
                    seg_dictionary[seg.attr['edgeid']] = seg.dictentry_value()

        elif self.sides:
            for seg in self.sides():
                if segment_id_type == SEGMENT_ID_TYPES['linehash']:
                    #entry = {seg.dictentry_linehashkey(): seg.dictentry_value()}
                    try:
                        #checking segment collision.
                        seg_dictionary[seg.dictentry_linehashkey()]
                    except:
                        seg_dictionary[seg.dictentry_linehashkey()] = seg.dictentry_value()

                if segment_id_type in [SEGMENT_ID_TYPES['lineid'], SEGMENT_ID_TYPES['lineid']]:
                    seg_dictionary[seg.attr['edgeid']] =seg.dictentry_value()

        return seg_dictionary

    def sides_as_collection(self, segment_id_type=SEGMENT_ID_TYPES['linehash']):
        collection = []

        #insert first entry
        p1, p2 = self.vertices[0], self.vertices[1]
        seg = Segment(p1, p2, {'edgeid': 0})
        struct_s = S(stype=STYPE['DIC'])
        struct_s.setK(k=seg.dictentry_linehashkey())
        struct_s.setNk(None)
        struct_s.setValue(seg.dictentry_value())
        collection.append(struct_s)

        for i in range(1, len(self.vertices) - 1):
            p1, p2 = self.vertices[i], self.vertices[i + 1]
            seg = Segment(p1, p2, {'edgeid': i})

            struct_s = S(stype=STYPE['DIC'])
            struct_s.setK(k=seg.dictentry_linehashkey())
            struct_s.setNk(None)
            struct_s.setValue(seg.dictentry_value())

            #set next-key for previous structure.
            collection[i-1].setNk(struct_s.getK())
            collection.append(struct_s)

        #last item's next-key is the first item's key
        collection[-1].setNk(collection[0].getK())
        #
        if debug:
            for d in collection:
                print("S"),str(d)

        return collection

    def from_segment_objects(self,sides):
        '''given a sequence of segment objects, construct a polygon object.
        the segment objects must be in the order of their occurances in the polygon.'''
        '''
        self.vertices = []
        for idx in xrange(0, len(sides), 2):
            pointa, pointb = sides[idx].origco_ordinates()
            self.vertices += [pointa,pointb]
        if len(sides) % 2 == 0:  # n is even
            self.vertices += [self.vertices[0]]
        
        newpoly =Polygon(self.vertices)
        return newpoly
        '''
        newpoly = Polygon([])
        if sides:
            newpoly.setSides(sides)
            self.setPropertyByName('n(sides)', len(sides))
            self.setPropertyByName("n(vertices)",len(sides))
        return newpoly

    def convert_sides_to_collection(self,sides_list):

        '''given a sequence of sides, returns a list of object of type S
        whose are of form (k, nk, value)'''
        collection = []

        #insert first entry
        seg = sides_list[0]
        struct_s = S(stype=STYPE['DIC'])
        struct_s.setK(k=seg.dictentry_linehashkey())
        struct_s.setNk(None)
        struct_s.setValue(seg.dictentry_value())
        collection.append(struct_s)

        for i in range(1, len(sides_list)):
            seg = sides_list[i]
            struct_s = S(stype=STYPE['DIC'])
            struct_s.setK(k=seg.dictentry_linehashkey())
            struct_s.setNk(None)
            struct_s.setValue(seg.dictentry_value())

            #set next-key for previous structure.
            collection[i-1].setNk(struct_s.getK())
            collection.append(struct_s)

        #last item's next-key is the first item's key
        collection[-1].setNk(collection[0].getK())

        if debug:
            for d in collection:
                print("S"),str(d)
        return collection


    def all_segpairs(self):
        seglist = self.sides()
        pairs = []
        n = len(seglist)
        for i in range(0, n):
            for j in range(i + 1, n):
                pair = seglist[i], seglist[j]
                pairs.append(pair)
        return pairs

    def shoelace_area(self):

        polygon_segs = self.vertices
        '''polygon-segs is a sequence of connecting point in a polygon'''
        if len(polygon_segs)<3:
            return 0
        if polygon_segs[0] != polygon_segs[-1]:
            polygon_segs +=[polygon_segs[0]] #the first point must be last point in the series.

        lra = 0 #left-to-right product
        rla = 0
        for i in range(len(polygon_segs)-1):
            lra += polygon_segs[i].getX()* polygon_segs[i+1].getY()
            rla += polygon_segs[i].getY()* polygon_segs[i+1].getX()
        return abs(lra-rla)/float(2)

    def vertical_sweeplines(self):
        '''
        :return: a dictionary whose key is xvalue, and value is
        'index' of xvalue in sorted array.
        '''
        import collections
        if self.vertices:
            xunik = sorted(set([p.getX() for p in self.vertices]))
        elif self.getSides():
            xunik =[]
            for seg in self.getSides():
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


    def split_mesh_sides_at_x(self,xunikdic):

        xunikdic = self.vertical_sweeplines()
        sides = self.sides()
        all_splits = []
        for seg in sides:
            # split seg
            # 1)get x-values at which split must take place
            li = xunikdic[seg.getLeftPoint().getX()]
            hi = xunikdic[seg.getRightPoint().getX()]

            splits=seg.split_at_multiple_x(xunikdic.keys()[li + 1:hi])
            all_splits +=splits
            temppoly.setPropertyByName('split_counts',len(splits)) #segment's split count

        newpoly = temppoly.from_segment_objects(all_splits)
        newpoly.setPropertyByName('split_counts',temppoly.getPropertyByName('split_counts'))
        newpoly.setPropertyByName('vcolumns_count', len(xunikdic)) #number of unique x-values in this poly.
        newpoly.setPropertyByName("issplitpoly", True) #number of unique x-values in this poly.

        xvalues = sorted(xunikdic.keys())
        newpoly.setPropertyByName('minx',xvalues[0])
        newpoly.setPropertyByName('maxx',xvalues[-1])

        #segment's properties.
        for seg in all_splits:
            newpoly.setPropertyByName('seg_delx', round(seg.delx(),2))
            newpoly.setPropertyByName('seg_dely', round(seg.dely(),2))
            newpoly.setPropertyByName('seg_lengths', round(seg.length(),2))

        return newpoly


    def split_sides_at_x(self,doround=False):
        temppoly = Polygon([])
        xunikdic = self.vertical_sweeplines()
        sides = self.sides()
        all_splits = []
        for seg in sides:
            # split seg
            # 1)get x-values at which split must take place
            li = xunikdic[seg.getLeftPoint().getX()]
            hi = xunikdic[seg.getRightPoint().getX()]

            splits = seg.split_at_multiple_x(xunikdic.keys()[li + 1:hi],doround)
            all_splits += splits
            temppoly.setPropertyByName('split_counts',len(splits)) #segment's split count

        newpoly = temppoly.from_segment_objects(all_splits)
        newpoly.setPropertyByName('split_counts',temppoly.getPropertyByName('split_counts'))
        newpoly.setPropertyByName('vcolumns_count', len(xunikdic)) #number of unique x-values in this poly.
        newpoly.setPropertyByName("issplitpoly", True) #number of unique x-values in this poly.

        xvalues = sorted(xunikdic.keys())
        newpoly.setPropertyByName('minx',xvalues[0])
        newpoly.setPropertyByName('maxx',xvalues[-1])

        #segment's properties.
        for seg in all_splits:
            newpoly.setPropertyByName('seg_delx', round(seg.delx(),2))
            newpoly.setPropertyByName('seg_dely', round(seg.dely(),2))
            newpoly.setPropertyByName('seg_lengths', round(seg.length(),2))

        return newpoly

    def  non_overlaping_yspans(self,ytuples):

        '''given list of (yi,yj,li) that means line li has y-span yi-yj, construct a coverup y-span
        such that the y-span's form a complete lut-entries. For instance:
        (1,3,a) and (2,4,b) are overlapping y-spans; it must result a tuple (1,4,{a,b}) and (4,1,{})
        tuples.'''
        #sort each tuple in ytuples.ytuples[jdx-1][1]

        for i in range(len(ytuples)):
            if ytuples[i][1] < ytuples[i][0]:
                ytuples[i] = (ytuples[i][1],ytuples[i][0]) + ytuples[i][2:]

        #sort list of tuples by a tuple's first element.
        from operator import itemgetter
        #slower version:
        #ytuples = sorted(ytuples, key=lambda tuple: tuple[0])
        #faster version:
        ytuples = sorted(ytuples,key=itemgetter(0))
        recs ={}
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
                recs[y0].append((bot,top))
            jdx = idx+1
            while jdx < len(ytuples) and ytuples[jdx][0] < top: #ytuples[jdx-1][1]:
                recs[y0][0].append(ytuples[jdx])
                #recs[y0][1]=(min(recs[y0][1][0],ytuples[jdx][0]),max(recs[y0][1][1],ytuples[jdx][1]))
                bot,top = min(bot,ytuples[jdx][0]),max(top,ytuples[jdx][1])
                recs[y0][1] = bot,top
                jdx +=1
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

    def xcolumns_yblocks(self):
        #check if this poly's segments are split
        if not self.getPropertyByName("issplitpoly"):
            print("Exception: Excepted a polygon whose segements are split at unique x-values.")
            return {}
        #collect all x:[(y1,y2),...] for each segments

        xcolumn_yblock_dic = OrderedDict()
        for xkey in self.vertical_sweeplines().keys():
            xcolumn_yblock_dic[xkey] =[]

        segments = self.sides()
        for seg in segments:
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


        count_vcols_nsplits_yblocks={} #x:(count-number-of-splits in x, count-number-of-yblocks)

        #construct nested dictionary for sqdm.{xkey:{(y1,y2):{segid:segvalue},(y1,y2):{segid:segvalue},(y1,y2):{segid:segvalue}}}
        for xkey, list_yblocks in xcolumn_yblock_dic.items():
            ypairs_dic = {}
            list_yspans_pairs =util.pairwise(self.non_overlaping_yspans(list_yblocks))
            count=0
            for ypair in list_yspans_pairs:
                ypairs_dic[str(ypair)] ={}
                count+=1

            xcolumn_yblock_dic[xkey] =ypairs_dic
            #count of splits, yblocks in each of vcolumns (x1-x2)
            count_vcols_nsplits_yblocks[xkey] = (len(list_yblocks),count)

            if debug:
                print("xkey:yblocks"), xkey,xcolumn_yblock_dic[xkey]

        self.setPropertyByName("count_xcolumns_nsplits_yblocks",count_vcols_nsplits_yblocks)
        return xcolumn_yblock_dic

    def segment_mapping_to_yblocks(self,xcolumn_yblock_dic):
        '''

        :param xcolumn_yblock_dic: is an SQDM, a nested dictionary of <xval,<(y1-y2),<seghash:value>>>
        :return: returns the input parameter by mapping segments in self.sides.
        '''
        if not xcolumn_yblock_dic:
            return {}
        for seg in self.sides():
            x1,y1,x2,y2 = seg.co_ordinates()
            #map to x-column
            xcolumn_yblocks = xcolumn_yblock_dic[x1]
            for yblock in xcolumn_yblocks.keys():
                ylow, yhigh = literal_eval(yblock)
                if y1 < y2:
                    if ylow <= y1 and y2 <= yhigh:
                        #TODO: make an ordered dictionary inside.
                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = 1
                else:
                    if ylow <= y2 and y1 <= yhigh:
                        xcolumn_yblocks[yblock][seg.dictentry_linehashkey()] = 1

        xkey_max_seg_holding={}
        for x,dyb in xcolumn_yblock_dic.items():
            max_holding = -sys.maxint
            for yb, sd in dyb.items():
                max_holding = max(max_holding,len(sd))
                xkey_max_seg_holding[x] = max_holding
        self.setPropertyByName("xkey_max_seg_holding",xkey_max_seg_holding)

        return xcolumn_yblock_dic

    def label_bottom_up(self,pid,outid=None):
        '''give a polygon segments in counter clockwise directions.'''

        L =[]
        for seg in self.sides():
            p1,p2 = seg.origco_ordinates()
            L.append([pid] + list(p1.totuple()) + list(p2.totuple()) + [0,0,0])

        if outid == None:
            INF = util.POLYGON_OUTSIDEID
        else:
            INF = outid

        for ti in range(len(L)):  # each tuple-index
            t = L[ti]  # [rid,x1,y1,x2,y2,iswp?,rid-left,rid-right,isvisited?]
            polyid, s1x1, s1y1, s1x2, s1y2, s1isswapped = t[0:6]
            # check is swapped.
            if s1isswapped:
                x1, y1, x2, y2 = s1x2, s1y2, s1x1, s1y1  # get original arrangement of points for the seg.
            else:
                x1, y1, x2, y2 = s1x1, s1y1, s1x2, s1y2

            if t[-1] == False:  # this segment is no shared by any other rings.
                t[-1] = True  # now visti
                if x1 < x2:
                    t[6] = polyid  # top
                    t[7] = INF  # bot
                elif x1 > x2:
                    t[6] = INF  # top
                    t[7] = polyid  # bot
                elif x1 == x2:
                    if y1 < y2:  # vup
                        t[6] = polyid  # top
                        t[7] = INF  # bot
                    elif y1 > y2:  # vdown
                        t[6] = INF  # top
                        t[7] = polyid  # bot
        attrlist = [(l[6],l[7]) for l in L] #return (above,bel)
        #self.setSidesAttr([{'abv':t[0], 'bel':t[1]} for t in attrlist])
        return attrlist



class SQDM:

    def __init__(self,ikeys=None):
        if ikeys:
            self.keys = OrderedDict()
            for idx,key in enumerate(ikeys):
                self.keys[idx] = k

    def get_ordered_keys(self,xvals):
        xvals = sorted(set(xvals))
        dictxvals = OrderedDict()
        for idx,key in enumerate(xvals):
            dictxvals[key] = idx
        return dictxvals

    @classmethod
    def save(self,data,outfile,dosort=True):
        import json,sys,os.path,random, pickle
        #delete if exists
        if os.path.exists(outfile):
            if not os.path.isdir(outfile): os.remove(outfile)
            else:outfile += '/data-'
        outfile = outfile.replace(".json",'')
        outfile = outfile.replace(".p",'')
        outfile = outfile.replace(".pickle",'')
        with open(outfile+'.json', 'w') as fp:
            json.dump(data, fp,sort_keys=dosort,indent=4)

        print("json data saved in "),outfile + '.json'


    def sqdm_segments(self,inputsqdm):
        segments_keys =[]
        for xkey, yblocks in inputsqdm.items():
            for yb,segs in yblocks.items():
                segments_keys += segs.keys()

        return segments_keys

    @classmethod
    def test1(self):
        #Driver to Test Polygon
        doround=False
        poly24 = [(6.81,5.05), (3.98,3.95), (4.98,2.02), (2.00,0.97), (0.58,3.88), (2.00,3.86), (3.00,4.96), (3.00,6.26)] #collinear, revisited twice.

        pobj2 = [(3.0, 6.0), (4.0, 5.0), (2.0, 4.0), (1.0, 4.0),
                 (2.0, 1.0), (5.0, 1.0), (4.0, 4.0), (7.0, 5.0)]

        #pobj2 = [(7, 5), (7, 3), (7, 1), (2, 1.00), (1, 4), (1, 5), (1, 6), (3, 6)] #with vertical lines.
        #with v.lines
        #pobj2 = [(2, 8), (6, 6), (4, 4), (2, 0), (4, 2), (6, 4), (10, 8), (10, 6), (8, 4), (8, 0), (10, 2), (16, 0),
        #(14, 2), (12, 4), (14, 8), (18, 6), (16, 4), (16, 2), (18, 0), (20, 2), (20, 4),
                 #(20, 6), (18, 10), (16, 8), (12, 10), (14, 14), (12, 12), (10, 14), (6, 14), (6, 12), (4, 10)]

        '''pobj2 =[(1002451148, 2157155400), (1002985974, 2156979255), (1002949639, 2157027380),
                (1002934867, 2157037231), (1002863289, 2157053914), (1002817837, 2157080433),
                (1002809521, 2157109224), (1002785398, 2157118993), (1002737353, 2157112525),
                (1002747795, 2157163692), (1002682339, 2157195986), (1002698135, 2157230886),
                (1002705605, 2157276072), (1002646360, 2157300939), (1002625521, 2157338066),
                (1002621358, 2157383904), (1002604682, 2157424835), (1002537612, 2157424835),
                (1002530978, 2157400047), (1002502135, 2157405008), (1002451585, 2157428700),
                (1002425447, 2157429364), (1002433796, 2157386559), (1002409172, 2157344104),
                (1002380752, 2157356613), (1002281845, 2157408492), (1002228801, 2157405821),
                (1002171160, 2157503183),
                (1002161263, 2157544235), (1002129448, 2157555469), (1002049833, 2157598249)]'''

        #pobj2.reverse() #anti-clockwise
        pobj2 = Polygon(pobj2) #anti-clockwise
        pobj2.label_bottom_up(777)

        #save polygon
        seg_dict = pobj2.tosegsdict()
        self.save(seg_dict,"../out/tmp/aPo.json") #original

        ##split polygon and make 2D grids.
        split_poly = pobj2.split_sides_at_x(doround=True)

        print("verties of split poly"), split_poly.get_vertices()
        #save polygon
        seg_dict = split_poly.tosegsdict()
        self.save(seg_dict,"../out/tmp/aPs.json",dosort=True) #splits do not sort by keys.

        xcolumns_yblocks_dic = split_poly.xcolumns_yblocks()
        mapped_segments_xcolumn_yblock_dic = split_poly.segment_mapping_to_yblocks(xcolumns_yblocks_dic)
        '''
        if debug:
            for xkey,yblocks in mapped_segments_xcolumn_yblock_dic.items():
                print xkey
                for yblock in yblocks:
                    print("\t"),"key:",yblock, "values:",yblocks[yblock]
        '''

        #save dictionary sqdm.
        print type(split_poly.getProperties())
        util.poly_ptstoshp(pobj2.get_vertices(), "../out/tmp/aPo")
        if doround:
            util.poly_ptstoshp(split_poly.get_vertices(),"../out/tmp/aPs.int")
        else:
            util.poly_ptstoshp(split_poly.get_vertices(), "../out/tmp/aPs")
        self.save(split_poly.getProperties(),"../out/tmp/aPsqdm_prop.json")
        self.save(mapped_segments_xcolumn_yblock_dic,"../out/tmp/aPsqdm.json")

        return mapped_segments_xcolumn_yblock_dic

    @classmethod
    def test3(self):
        pobj2 = [ (3,4.5),(3.5,4),(3,3.7),(2,3.5),(3,4.5)]
        pobj2 = [(2.6, 3.5), (2.8, 3.2), (3.8, 3.7), (3.2, 4.6)]

        pobj2 = Polygon(pobj2) #anti-clockwise
        pobj2.label_bottom_up(util.POLYGON_INNERID)
        seg_dict = pobj2.tosegsdict()


        util.save(seg_dict,"../out/tmp/C.json")
        print

        ##split polygon and make 2D grids.
        split_poly = pobj2.split_sides_at_x()
        seg_dict = split_poly.tosegsdict()

        util.save(seg_dict,"../out/tmp/Cs.json",dosort=False)
        print

        xcolumns_yblocks_dic = split_poly.xcolumns_yblocks()
        mapped_segments_xcolumn_yblock_dic = split_poly.segment_mapping_to_yblocks(xcolumns_yblocks_dic)
        return mapped_segments_xcolumn_yblock_dic

    @classmethod
    def test_usa_sqdm(self):
        ##use usa polygon
        import numpy as np
        home = "D:/workspace/sqdm-repo/sqdm/out/tmp"
        infile1 = home + "/usa.prj.lbl.txn.int.txt"

        def get_usa_polygon(infile1):
            arr = np.loadtxt(infile1)
            arr = arr[:, :]#.astype(np.dtype('int64'))
            print arr.shape
            segments = tuple(arr)

            segment_tuples = [t[1:] for t in segments] #remove first value which is polygon id.
            usa_polygon = Polygon([])
            edgeid = 0
            for seg in segment_tuples:
                attr = {}
                lpoint = Point(seg[0:2])
                rpoint = Point(seg[2:4])
                isswp = seg[4]
                attr['abv'] = seg[5]
                attr['bel'] = seg[6]
                attr['edgeid'] = edgeid

                if isswp:
                    segObj = Segment(rpoint,lpoint,attr)
                else:
                    segObj = Segment(lpoint, rpoint, attr)

                usa_polygon.addSegment(segObj)
                edgeid +=1
            util.poly_ptstoshp(usa_polygon.get_vertices(), "../out/tmp/2USA")

            return usa_polygon

        def split_usa_polygon_at_x(usa_polygon):
            ##split polygon and make 2D grids.
            usa_split_poly = usa_polygon.split_sides_at_x(doround=True)
            seg_dict = usa_split_poly.tosegsdict()


            util.poly_ptstoshp(usa_split_poly.get_vertices(), "../out/tmp/2USA.int")
            util.save(seg_dict, "../out/tmp/2USAs.json") #segment dictionary
            return usa_split_poly

        usa_polygon = get_usa_polygon(infile1)
        usa_split_poly = split_usa_polygon_at_x(usa_polygon)
        del usa_polygon

        print("completed splitting us map.")
        print("minx"),usa_split_poly.getPropertyByName("minx")
        print("maxx"),usa_split_poly.getPropertyByName("maxx")
        print("vcolumns-count"),usa_split_poly.getPropertyByName("vcolumns_count")

        xcolumns_yblocks_dic = usa_split_poly.xcolumns_yblocks()
        segment_mapped_xcolumn_yblock_dic = usa_split_poly.segment_mapping_to_yblocks(xcolumns_yblocks_dic)

        for xkey,yblocks in segment_mapped_xcolumn_yblock_dic.items()[0:2]:
            print xkey
            for yblock in yblocks:
                print("\t"),"key:",yblock, "values:",yblocks[yblock]


        util.save(usa_split_poly.getProperties(),"../out/tmp/2USA_sqdm_prop.json")
        util.save(segment_mapped_xcolumn_yblock_dic, "../out/tmp/2USA_sqdm.json")
    @classmethod
    def get_usa_state_boundary_by_name(self,state_name):
        ##use usa polygon
        import numpy as np
        home = "D:/workspace/sqdm-repo/sqdm/"
        infile1 = home + "/states.prj.lbl.txn.int.txt"
        arr = np.genfromtxt(infile1,dtype='str')
        arr = arr[:, :]#.astype(np.dtype('int64'))
        print arr.shape
        segments = tuple(arr)
        state_polygon = Polygon([])
        edgeid = 0
        for seg in segments:
            if str(seg[0]) == state_name:
                attr = {}
                lpoint = Point(seg[1:3].astype(np.dtype('int64')))
                rpoint = Point(seg[3:5].astype(np.dtype('int64')))
                isswp = int(seg[5])
                attr['edgeid'] = edgeid

                if isswp:
                    segObj = Segment(rpoint,lpoint,attr)
                else:
                    segObj = Segment(lpoint, rpoint, attr)

                state_polygon.addSegment(segObj)
                edgeid +=1
            #else continue

        print("len state-polygon"),len(state_polygon)
        #util.poly_ptstoshp(state_polygon.get_vertices(), "../out/tmp/"+state_name+".int")
        #util.save(seg_dict, "../out/tmp/"+state_name+".json") #segment dictionary

        return state_polygon

if __name__ == "__main__":

    #SQDM.test1()
    #SQDM.test_usa_sqdm()
    #SQDM.test3()
    #SQDM.get_usa_state_boundary_by_name("CALIFORNIA")
    pass