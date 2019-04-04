Left = 0
Right = 1
EVENTNAME = {0:"Start",1:"End"}
SEGMENT_ID_TYPES={'linehash':'linehash','lineid':'lineid','edgeid':'edgeid'}

class Point:
    def __init__(self,tup):
        self.x, self.y = tup

    def isLeft(self,ostartpoint,oendpoint):
        #if self point is left of line joing two points
        #ostartpoint and oendpoint
        return (ostartpoint.x - self.x)*(oendpoint.y - self.y) -\
               (oendpoint.x-self.x)*(ostartpoint.y - self.y)

    def areatriangle(self,op1,op2):
        return self.isLeft(op1,op2)

    def xyorder(self,other):
        #returns if this point is left to other point
        if self.x > other.x:
            return True
        if self.x < other.x:
            return -1
        if self.y > other.y:
            return True
        if self.y < other.y:
            return -1
        return 0 #for all cases.
    
    def __str__(self):
        return str((self.x,self.y))

    def totuple(self):
        return self.x,self.y

class SlSegment(object):
    def __init__(self, edgeid=None, leftpoint=None,rightpoint=None,segabv=None,segbel=None):
        self.edgeid = edgeid
        self.lp = leftpoint
        self.rp = rightpoint
        self.segabv = segabv
        self.segbel = segbel
        if self.lp:
            self.key = str(self.lp.y)+str(self.edgeid)
        else:
            self.key=None

    def co_ordinates(self):
        return self.lp.x,self.lp.y, self.rp.x, self.rp.y

    # get m for a given seg in (p1, p2) form
    def get_slope(self):
        x0, y0, x1,y1= self.co_ordinates()
        if (x1 - x0) == 0:
            return None
        else:
            return float(y1 - y0) / (x1 - x0)

    # given a point p, return the point on s that shares p's y-val
    def get_x_at(self,p):

        m = self.get_slope()

        # for now that it would have been deleted already if not
        if m == 0:  # horizontal segment
            return p
        # ditto; should check if y-val on seg
        if m is None:  # vertical segment
            return (self.lp.x, p[1])
        x1 = self.lp.x - (self.lp.y - p[1]) / m
        # this should check if p's x-val is actually on seg; we're assuming
        if self.lp.x <= x1 <= self.rp.x:
            return (x1, p[1])
        return None

    # given a point p, return the point on s that shares p's x-val
    def get_y_at(self,p):
        '''this is wrong. correct it.'''
        m = self.get_slope()

        # ditto; should check if y-val on seg
        if m is None:  # vertical segment
            return (self.lp.x, p[1])
        y1 = m * (p[0] - self.rp.x) + self.rp.y #yn = float(y2 - y1) / (x2 - x1) * (xn - x1) + y1
        # this should check if p's y-val is actually on seg; we're assuming
        if self.lp.y <= y1 <= self.rp.y:
            return (x1, p[1])
        return None


    def __str__(self):
        if not self.key:
            self.key=''

        st = ''
        st += "key:"+self.key+","+str(self.edgeid) + ","+str(self.lp) + "--"+str(self.rp)
        st += " abv:"+ str(self.segabv) + str(" bel:") +str(self.segbel)
        return st

    def dictentry_linehashkey(self):
        x1,y1 = self.lp.totuple()
        x2,y2 = self.rp.totuple()
        key = util.hashargs(x1,y1,x2,y2)
        val = x1,y1,x2,y2,self.edgeid,self.segabv,self.segbel
        return (key,val)

    def dictentry_linidkey(self):
        x1,y1 = self.lp.totuple()
        x2,y2 = self.rp.totuple()
        val = x1,y1,x2,y2,self.edgeid,self.segabv,self.segbel
        return (self.edgeid,val)

    def tupleToSeg(self,index,pts):
        '''pts is 4 tuple segment'''
        x1, y1, x2, y2 = pts
        p1 = Point((x1, y1))
        p2 = Point((x2, y2))
        self.edgeid = index
        if p1.xyorder(p2) < 0:
            self.lp = p1
            self.rp = p2
        else:
            self.lp = p2
            self.rp = p1


import util
class Polygon(object):
    import hashlib

    def __init__(self,list_pts):
        #polygon is set of Points
        self.vertices = [Point(tup) for tup in list_pts]
        self.nv = len(list_pts)-1

    def xorderedsegs(self):
        seglist = []
        for i in range(0, len(self.vertices) - 1):
            p1, p2 = self.vertices[i], self.vertices[i + 1]

            if p1.xyorder(p2) < 0:
                seg = SlSegment(p1, p2)
            else:
                seg = SlSegment(p2, p1)

            seglist.append(seg)
        for seg in seglist:
            print seg
        return seglist

    def xorderednamedsegs(self):
        seglist = []
        for i in range(0, len(self.vertices) - 1):
            p1, p2 = self.vertices[i], self.vertices[i + 1]

            if p1.xyorder(p2) < 0:
                seg = SlSegment(i,p1, p2)
            else:
                seg = SlSegment(i,p2, p1)

            seglist.append(seg)
        for seg in seglist:
            print seg
        return seglist

    def tupleToSegs(self,index,pts):
        '''pts is 4 tuple segment'''
        x1, y1, x2, y2 = pts
        p1 = Point((x1, y1))
        p2 = Point((x2, y2))
        if p1.xyorder(p2) < 0:
            seg = SlSegment(index,p1,p2)
        else:
            seg = SlSegment(index,p2,p1)

        return seg


    def tosegsdict(self,segment_id_type=SEGMENT_ID_TYPES['linehash']):
        '''returns a set of items of this form for each segment:
        (hash(x1,y1,x2,y2),(x1,y1,x2,y2,AreaAbove,AreaBelow,etc.))'''
        seglist =[]
        for i in range(0, len(self.vertices) - 1):
            p1, p2 = self.vertices[i], self.vertices[i + 1]

            if p1.xyorder(p2) < 0:
                seg = SlSegment(edgeid=i, leftpoint=p1,rightpoint=p2,segabv='A',segbel='B')
            else:
                seg = SlSegment(edgeid=i, leftpoint=p2,rightpoint=p1,segabv='A',segbel='B')

            if segment_id_type == SEGMENT_ID_TYPES['linehash']:
                segentry = seg.dictentry_linehashkey()

            if segment_id_type in [SEGMENT_ID_TYPES['lineid'],SEGMENT_ID_TYPES['lineid']]:
                segentry = seg.dictentry_linidkey()

            seglist.append(segentry)
        return seglist

    def segs(self):
        seglist = []
        for i in range(0,len(self.vertices)-1):
            p1,p2 = self.vertices[i],self.vertices[i+1]

            if p1.xyorder(p2) < 0:
                seg =SlSegment(i,p1,p2)
            else:
                seg =SlSegment(i,p2,p1)

            seglist.append(seg)

        return seglist
    def all_segpairs(self):
        seglist = self.segs()
        pairs =[]
        n = len(seglist)
        for i in range(0,n):
            print("\t",i)
            for j in range(i+1,n):
                pair = seglist[i],seglist[j]
                print("pair",str(pair[0]),str(pair[1]))
                pairs.append(pair)
            print
        
        print("completed pairs")
        return pairs
    def __str__(self):
        st = ''
        for ptobj in self.vertices:
            st += str(ptobj)
        return "Polygon:"+st

class Event:
    def __init__(self,edgeid=None,etype=None,evpt=None):
        self.edgeid = edgeid
        self.eventtype = etype #int repr
        self.eventvertex = evpt
        self.key = (self.eventvertex.x ,self.eventvertex.y)
    def __str__(self):
        return "["+EVENTNAME[self.eventtype]+","+str(self.eventvertex)+","+str(self.edgeid)+","+"key:"+str(self.key)+"]"

class EventQ(object):
    def __init__(self,polyobj=None):
        self.Nv = polyobj.nv
        self.ne = 2* self.Nv
        polypts = polyobj.vertices
        self.events = []
        self.eventidx = 0
        for i in range(0,self.Nv):
            p1 = polypts[i]
            p2 = polypts[i+1]
            edgeid = i
            #find which is left and which is right event.
            if p1.xyorder(p2) < 0:
                event1 = Event(edgeid,Left,p1)
                event2 = Event(edgeid,Right,p2)
            else:
                event1 = Event(edgeid,Right,p1)
                event2 = Event(edgeid,Left,p2)
            self.events.append(event1)
            self.events.append(event2)
        #end-for
        #sort eventd depending upon their
        #self.event_sort_by_x()

    def __iter__(self):
        return self

    def __next__(self):
        self.eventidx +=1
        try:
            return self.events[self.eventidx-1]
        except IndexError:
            self.eventidx = 0
            raise StopIteration


    def __str__(self):
        st = ''
        for ev in self.events:
            st += str(ev) +'\n'
        return st

    def event_sort_by_x(self):
        self.events = sorted(self.events, key=lambda event: event.eventvertex.x, reverse=False)

        print("Sorted Events in Increasing X-value")
        return True

    def event_sort_by_key(self):
        self.events = sorted(self.events, key=lambda event: event.key, reverse=False)
        print("Sorted Events in Increasing key")
        return True
    def next(self):
        if self.eventidx >= self.ne:
            self.eventidx = 0
            raise StopIteration
        nxtevent = self.events[self.eventidx]
        self.eventidx+=1
        return  nxtevent

    def anext(self):
        self.eventidx+=1
        if self.eventidx-1 >= self.ne:
            self.eventidx =0
            return None
        return self.events[self.eventidx-1]

    def iterate(self):
        event = Event()
        print self.eventidx
        event = self.next()
        while(event):
            print(event)
            print self.eventidx
            event = self.next()
        self.eventidx = 0   

class SweepLine:
    def __init__(self,polyObj=None,bstroot=[]):
        self.nv = polyObj.nv #number of vertices in poly
        self.poly = polyObj #polygon object
        self.bstroot = bstroot  #balanced binary tree root

    def __str__(self):
        st = "\n\n"+"input Polygon"+ str(self.poly)+"\n"
        st += "nv:"+str(self.nv) +"\n"
        st += "slstatus:"+str(self.bstroot)+"\n\n"
        return st

    '''Adds a segment to a sweep line arrary'''
    def addEventKey(self,newevent):
        #self is a sweepline obj, it has poly object as a reference.
         p1,p2 = self.poly.vertices[newevent.edgeid],self.poly.vertices[newevent.edgeid+1]
         if p1.xyorder(p2) < 0:
             lp = p1
             rp = p2
         else:
             lp = p2
             rp = p1
         #(self,edgeid = None, leftpoint = None, rightpoint = None, segabv = None, segbel = None)
         slseg = SlSegment(newevent.edgeid,lp,rp) #create slseg obj with other prop. none

         #add this segment as node in BST
         #print("Pre:append slseg"),str(slseg)
         self.bstroot.append(slseg)
         self.bstroot = sorted(self.bstroot, key=lambda seg: seg.key, reverse=False) #sort by y-value
         print("\t sweep-status:")
         for item in self.bstroot:
             print("\t "),item
         #find slseg's above and below slsegments
         foundidx,curidx = 0,0
         for each_in in self.bstroot:
             #print("sl"),str(each_in)
             if slseg.edgeid == each_in.edgeid:
                 foundidx = curidx
                 break
             curidx +=1
         #print("\tfoundidx:",slseg.edgeid,foundidx)
         if foundidx+1 <= len(self.bstroot)-1:
             segabv = self.bstroot[foundidx+1]
         else:
             segabv = None
         if foundidx-1 >=0:
             segbel = self.bstroot[foundidx-1]
         else:
             segbel = None

         #slseg.segabv = segabv
         #slseg.segbel = segbel
         #print("\tappend slseg"),str(slseg)
         return slseg,segabv, segbel#added sweep line segement
    def addEvent(self,newevent):
         slseg = SlSegment(newevent.edgeid) #create slseg obj with other prop. none
         edgeid = newevent.edgeid

         p1,p2 = self.poly.vertices[newevent.edgeid],self.poly.vertices[newevent.edgeid+1]
         if p1.xyorder(p2) < 0:
             slseg.lp = p1
             slseg.rp = p2
         else:
             slseg.lp = p2
             slseg.rp = p1

         #add this segment as node in BST
         #print("Pre:append slseg"),str(slseg)
         self.bstroot.append(slseg)
         self.bstroot = sorted(self.bstroot, key=lambda seg: seg.lp.y, reverse=False) #sort by y-value   
         print("\t sweep-status:")
         for item in self.bstroot:
             print("\t "),item
         #find slseg's above and below slsegments
         foundidx,curidx = 0,0
         for each_in in self.bstroot:
             #print("sl"),str(each_in)
             if slseg.edgeid == each_in.edgeid:
                 foundidx = curidx
                 break
             curidx +=1
         #print("\tfoundidx:",slseg.edgeid,foundidx)
         if foundidx+1 <= len(self.bstroot)-1:
             segabv = self.bstroot[foundidx+1]
         else:
             segabv = None
         if foundidx-1 >=0:
             segbel = self.bstroot[foundidx-1]
         else:
             segbel = None

         #slseg.segabv = segabv
         #slseg.segbel = segbel
         #print("\tappend slseg"),str(slseg)
         return slseg,segabv, segbel#added sweep line segement

    def searchAbvBel(self,searchevent):
         #find searchevent's above and below slsegments
         foundidx,curidx = 0,0
         foundsegl = None
         for each_in in self.bstroot:
             print("sl"),str(each_in)
             if searchevent.edgeid == each_in.edgeid:
                 foundidx = curidx
                 foundsegl = each_in
                 break
             curidx +=1
         print("\t foundidx:"),foundidx
         if foundidx+1 <= len(self.bstroot)-1:
             segabv = self.bstroot[foundidx+1]
         else:
             segabv = None
         if foundidx-1 >=0:
             segbel = self.bstroot[foundidx-1]
         else:
             segbel = None

         return foundsegl,segabv,segbel
        
    def findEvent(self,searchevent):
        #find a segment in a tree.
        print("\t find and delete:"),searchevent
        idx =0
        for each_slseg in self.bstroot:
            if each_slseg.edgeid == searchevent.edgeid:
                return each_slseg
            idx +=1
        return None
        
    def intersect(self,firstseg, secondseg):
        #test if two segments intersects        
        #either of the segment doesn't exist
        if firstseg == None or secondseg == None:
            return False

        #test if firstseg and secondseg are adjoining segments
        edgeid1, edgeid2 = firstseg.edgeid, secondseg.edgeid
        if ((edgeid1+1) % self.nv) == edgeid2 or ((edgeid2+1)%self.nv)==edgeid1:
            print("-.-",edgeid1,edgeid2,False)
            return False

        #test intersection
        lsign = firstseg.lp.isLeft(secondseg.lp,secondseg.rp)
        rsign = firstseg.rp.isLeft(secondseg.lp,secondseg.rp)
        if lsign*rsign > 0: #firstseg endpoints have same sign relative to second segment
            print("--",edgeid1,edgeid2,False)
            return False #onsame side non intersection
        #
        lsign = secondseg.lp.isLeft(firstseg.lp,firstseg.rp)
        rsign = secondseg.rp.isLeft(firstseg.lp,firstseg.rp)
        if lsign*rsign > 0:
            print("--",edgeid1,edgeid2,False)
            return False
        #the segments s1 and s2 intersect
        print("-*-",edgeid1,edgeid2,True)
        return True

    def removeEvent(self,removeseg):
        #rslseg = self.findEvent(removeseg)
        #self.bstroot.remove(rslseg)
        try:
            self.bstroot.remove(removeseg)
        except:
            pass

    def issimplepolygon(self, eq):
        for event in eq:
            if event.eventtype == Left:
                print("event"), str(event)
                addedseg, abv, bel = self.addEvent(event)
                #test intersection abv
                if self.intersect(addedseg, abv):
                    print("intersected")
                    return False
                if self.intersect(addedseg, bel):
                    print("intersected")
                    return False
                #test intersection bel
            else:
                #find it in sweepline
                print
                print("event"),str(event)
                removseg, segabv, segbel= self.searchAbvBel(event)
                #test intersection
                if self.intersect(segabv, segbel):
                    print("intersected")
                    return False
                try:
                    self.removeEvent(removseg)
                except:
                    pass
        return True

    def issimplepolygonkey(self, eq):
        for event in eq:
            if event.eventtype == Left:
                print("event"), str(event)
                addedseg, abv, bel = self.addEventKey(event)
                #test intersection abv
                print("\t abv")
                if self.intersect(addedseg, abv):
                    print("intersected")
                    return False
                print("\t bel")
                if self.intersect(addedseg, bel):
                    print("intersected")
                    return False
                #test intersection bel
            else:
                #find it in sweepline
                print("event"), str(event)
                removseg, segabv, segbel= self.searchAbvBel(event)
                #test intersection
                if self.intersect(segabv, segbel):
                    print("intersected")
                    return False
                self.removeEvent(removseg)
        return True

def main1(poly):
    '''this algorithm does not always work because events are only ordered by x-values, not by (x,y) tuples.'''
    polyObj = Polygon(poly)
    print polyObj
    print polyObj.segs()
    print("all pairs")
    segpairs = polyObj.all_segpairs()

    print("\n Testing EventQ")
    eq = EventQ(polyObj)
    print eq
    print eq.event_sort_by_x()
    print eq
    #print eq.iterate()

    print
    for ev in eq:
        print("ev",str(ev))
    print("end EventQ")

    #driver
    sweepline = SweepLine(polyObj)
    print("sweepline init"), sweepline
    #use brute-force way
    decision=False
    for segpair in segpairs:
        s1,s2 = segpair
        decision = sweepline.intersect(s1,s2)
        if decision:break

    print("Testing sweepline algorithm")
    decision2 = sweepline.issimplepolygon(eq)


    print("By algorithm, this poly is simpe:"),decision2
    print("By brute force, this poly is simple:"),not decision


def main2(poly):
    '''This is correct implementation of Shamos-hoey algorithm for determining a self-intersecting polygon'''
    poly = poly
    polyObj = Polygon(poly)
    print polyObj
    #print polyObj.segs()

    print("\n Testing EventQ")
    eq = EventQ(polyObj)

    print eq.event_sort_by_key() #iomt implementation.
    print("EQ sortyed by key:")
    print eq
    sweepline = SweepLine(polyObj)
    print("Testing sweepline algorithm")
    decision2 = sweepline.issimplepolygonkey(eq)
    return decision2

if __name__ == '__main__':
    ##Driver
    print("Testing Polygon")
    ##https://www.mathsisfun.com/geometry/polygons-interactive.html
    poly0 = [(6,4.0), (4,2.00), (3,5),(6,4)]
    poly0 = [(6,4.0), (6,2.00), (3,5),(6,4)] #have vertical line so throws error.(does not work)
    poly0 = [(6,4.3), (6,2.00), (3,4),(6,4)] #have horizontal line so throws error.(does not work)
    poly0 =[(4,4), (2.00,3.0), (3.00,1.0), (6.00,1.0),(4,4)] #horizontal line with 4 sides.works
    poly0 =[(4,4), (6.00,3.0), (3.00,1.0), (4.00,1.0),(4,4)] #horizontal line with 4 sides.works
    poly1 = [(0,2),(1,4),(6,1),(3,5),(0,2)]
    poly2 = [(1,2),(4,3),(8,2),(6,6),(3,7),(1,2)]
    poly3 = [(1,2),(9,5),(8,2),(6,6),(3,7),(1,2)]
    poly4 = [(6,8),(5,4),(8,2),(6,6),(3,7),(6,8)]
    poly5 = [(4,3),(2,1),(4,0),(5,2),(6.5,1),(6,3),(8,3.5),(7,2)
             ,(8,5),(1,6),(3,7),(4,3)] #intersects
    poly6 = [(4,3),(2,1),(4,0),(5,2),(6.5,1),(6,3),(8,3.5),(5,4)
             ,(8,5),(1,6),(2,4),(4,3)] #do not intersect
    poly7 =[(0.00,2.00), (1.00,4.00), (6.00,1.00), (2.00,2.00), (1.00,1.00),(0,2)] #do not intersect

    poly8 = [(5.09,5.80), (3.54,5.97), (2.12,5.34),
             (1.20,4.09), (4.22,4.03), (1.66,1.12), (2.91,0.20),
             (4.46,0.03), (5.88,0.66), (2.65,3.13), (6.97,3.46), (6.34,4.88),(5.09,5.80)] #intersect
    poly9 = [(5.09,5.80), (3.54,5.97), (2.12,5.34), (1.20,4.09), (4.22,4.03),
             (1.66,1.12), (2.91,0.20),
             (4.46,0.03), (5.88,0.66), (3.99,2.99), (6.97,3.46), (6.34,4.88),(5.09,5.80)]

    poly10 = [(5.09,5.80), (2.95,4.98), (7.08,5.54),
              (1.03,2.54), (2.02,0.74), (3.93,0.00), (5.88,0.66), (7.00,2.00), (6.57,3.89),
              (5.09,5.80)]
    poly11 = [(5.09,5.80), (2.95,4.98), (1.99,6.03), (1.03,2.54), (2.02,0.74), (7.35,3.34), (5.88,0.66), (7.00,2.00), (-0.22,2.97), (5.09,5.80)]

    poly12 =[(5.09,5.80), (2.95,4.98), (1.99,6.03), (1.03,2.54), (2.02,0.74), (7.35,3.34), (5.88,0.66), (7.00,2.00), (8.48,3.86), (5.09,5.80)]

    poly13 = [(3.58,2.54), (2.95,4.98), (1.99,6.03), (1.03,2.54), (2.02,0.74), (5.03,3.02), (5.88,0.66),
    (6.19,2.24), (8.48,3.86), (2.69,3.01)]
    poly13 =[(7,5.00), (8.00,3.00), (2.00,1), (3.00,5.00), (6.00,6.00),(7,5)]
    poly14 = [(5.00,4.50), (6.00,6.00), (4.00,4.00), (2.00,5.00), (1.00,3.00),
              (3.00,2.00), (4.00,3.00), (3.00,1.00), (5.00,0.50), (7.00,2.00), (5.00,3.00), (7.00,4.00),(5,4.5)]

    poly15 =[(5.00,4.50), (6.00,6.00), (4.00,4.00), (2.00,5.00),
             (1.00,3.00), (3.00,2.00),(6.00,3.00), (3.00,1.00),
             (5.00,0.50), (7.00,2.00),(5.00,3.00), (7.00,4.00),(5,4.5)] #intersects
    poly16 =[(5.00,4.50), (6.00,6.00), (4.00,4.00), (2.00,5.00), (1.00,3.00), (3.00,2.00),
             (4.00,3.00), (3.00,1.00), (6.00,3.), (7.00,2.00), (5.00,3.00), (7.00,4.00),(5.0,4.5)]

    poly16 = [(5.00, 4.50), (6.00, 6.00), (4.00, 4.00), (2.00, 5.00), (1.00, 3.00), (3.00, 2.00),
              (4.00, 3.00), (3.00, 1.00), (6.00, 3), (7.00, 2.00), (5.00, 3.00), (7.00, 4.00), (5.0, 4.5)]
    poly17 = [(5.00,6), (2.00,3.00), (5.00,3.00),(5,6)] #triangle
    poly18 = [(5,6), (1.68,4.90), (1.00,3.00), (4.76,0.10), (7.02,3.00),(5,6)]
    poly19 = [(5.00,6.0), (2.00,6.00), (2.00,3.00), (6,3.00), (5.00,1.00), (7.00,1.00), (7.00,6.00),(5,6)] #is simple
    '''For this polygon19, this algorithm which orderes events only by x-values results that this polygon is Simple.
    But it is not. Because it sorts events only by x-values. must use this case to give correct results:
    '''
    poly20a =[(1.00,5.00), (5.00,2.00), (3.00,1.00), (3.00,6.00),(1,5)] #not simple
    poly20b =[(1.00,5.00), (5.00,2.00), (6.00,4.00), (3.00,6.00),(1,5)] #simple
    poly21 =[(0.00,0.00), (10.00,0.00), (20.00,0.00), (15,0), (20.00,0.00),(20,20),(0,20),(0,0)] #simple
    poly22=[(0,0), (7,8), (20,0), (24,7), (10,5), (20,20), (3,19),(0,0)]
    #two vertical line segments.
    v1=(0, 0);
    v2=(0, 1);
    v3=(0, 3);
    v4=(0, 4);
    v5=(2, 2);
    v6=(4, 2);
    poly23 =[v1,v3,v5,v3,v4,v6,v1] #simple
    poly23 = [v1, v2, v6, v3, v1] #collinear, revisited twice.
    poly=poly23
    print poly23
    #print main1(poly)
    print main2(poly)
