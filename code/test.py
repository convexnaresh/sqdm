##link https://www.geogebra.org/graphing

ctypes = {'nodeverify': 0, 'equivalentcert': 1}
class TestCertificate(object):
    a = "asdfaddsd"
    def __init__(self,ctype,*args):

        self.certtype = ctype
        self.arg = args
    def __str__(self):
        s =''
        print("...dd",a)
        self.a ="llllllllllll"
        for it in self.arg:
            s+=str(it)
        return s

def enclosedby(sk, snk, k):
    # key,nextkey,new-key
    c1 = (sk < k < snk)
    c2 = (k > sk >= snk)
    c3 = (sk >= snk > k)
    print c1,c2,c3
    if (sk < k < snk) or (k > sk >= snk) or (sk >= snk > k):
        return True
    return False

print("enclosed"),enclosedby(1,1,15)

class A:
    def __init__(self,a=[]):
        self.a =a
        self.mod()
    def mod(self):
        if not self.a:
            self.a.append("a")
    def __str__(self):
        return str(self.a)

obja = A()
print(obja)

objb = A(['xx'])
print objb

from collections import OrderedDict
from itertools import tee, izip

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

print list(pairwise((1,2,3,4,5,6)))
import matplotlib.pyplot as plt
import numpy as np
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    from itertools import tee, izip
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

pairs = pairwise(['1'])
for p in pairs:
    print p
print ("completed.")


import collections





def load_data(infile):
    import json
    data = None
    infile = infile.replace(".json", '')
    infile = infile.replace(".p", '')
    infile = infile.replace(".pickle", '')
    for ext in ['.json', '.pickle']:
        try:
            with open(infile + ext, 'r') as fp:
                data = json.load(fp)
        except:
            continue
    return data


x1,y1,x2,y2 = 2,4,4,5
xn = 3.2
yn = float(y2 - y1) / (x2 - x1) * (xn - x1) + y1
print yn

d =OrderedDict({1:1,1:5,3:22})

d2 = OrderedDict(enumerate(d.keys()))
print d2

print set([1,2] + [1, 2,3,4])

for key in d2.keys():
    print key
print (key)

a = [int, float, str]
print float == a[0]
##



pobj2 = [(3.2,4.6),(3.8,3.7),(2.8,3.2),(2.6,3.5)]
pobj2.reverse()
print pobj2

d ={1:22,4:55,6:44}
od = OrderedDict(d)
for item in od.items():
    print item
print
od[1] =33
for item in od.items():
    print item


class Util:

    def __init__(self):
        self.i = 0

    @classmethod
    def getx(self,j=None):
        return str(j)

u = Util()
print u.getx('c')

u2 = Util()
print Util().getx('d')

#------------------------------------------
import bisect
# initializing list
test_list = [1,2,3,4,5]
# printing original list
print ("The original list is : " + str(test_list))
# insert element
k = 5
# using bisect.insort()
# insertion in sorted list
# using naive method
bisect.insort(test_list, k)

# printing result
print ("The list after insertion is : " + str(test_list))

#------------------------------------------
from sortedcontainers import SortedDict
d = SortedDict({'z':1,'a':3})
print d
d['c'] = 45
print d
odd = OrderedDict(d)
print odd
'''
A. Input: P_sqdm, Cs, DELEGATED_BY, DELEGATED_TO

1.1. lk,hk = LowHighXKeys(Cs)
    #lk,hk = min(Cs.x-vals), max(Cs.x-vals) 

1.2. IsLoop(Cs)
    n=Cs.SegentsCount
    #for each si (i<n),s(i+1), si.B == S(i+1).A #A is start point, B is end point
    # sn.B == s1.A

1.3. IsSimple(Cs)
    #bool=ShamosHoey(Cs)
    #if bool==False, Exit

1.4. AV = AssignmentVector(Cs)
    #bool = Cs.CounterClockwise()
    #if bool==False, Exit
    #Walk CounterClockwise along points, for each segments note down 
    (region-left,region-right), depending on turning right or left turn.
    
    AV=[(0,0) for i in range(Cs.SegmentsCount)]
    set THIS_POLYID, OUT_SPACE_ID
    for each s_i in Cs.Segments:
        if s_i.A.x < s_i.B.x:
            AV[i] = THIS_POLYID, OUT_SPACE_ID
        elif x1 > x2:
            AV[i] = OUT_SPACE_ID,THIS_POLYID
        elif x1 == x2:
            if y1 < y2:  # vup
                AV[i] = THIS_POLYID, OUT_SPACE_ID
            elif y1 > y2:  # vdown
                AV[i] = OUT_SPACE_ID, THIS_POLYID                   
        
    
1.5. PartitionVerification(Cs, DELEGATED_BY,P_Sqdm)
    #for each s_i in Cs, bool=find r=(xi,xj,yi,yj,D) in P_Sqdm : that
    # s_i completely lies inside 'r'.
    # IF not not found, Exit.

2. sP_sqdm, sPs = Slice(lk,hk,P_sqdm) #slice SQDM, segments
    #find all (xkey,next-x-key) from P_sqdm that is intersected by (lk,0) and (hk,0)
    
3. UxK = Cs.Xkeys + sP_sqdm.XKeys
4. sorted(UxK)
5. sPsUxk= splited(sPs, UxK)
    #for each s_i in sPs, find x_k in UxK : x_k intersects s_i vertically.
    #split at each x=x_k, find (x_k, y_k) on s_i.
    
6. CsUxK = splited(Cs, UxK)
6.1 UnionSplits= sPsUxK + CsUxk
7. tempSqdm = xcolumns_yblocks(UnionSplits, UxK)
    #for each in CxK, collect (yi,yj) blocks from each in UnionSplits
    (x1,x2,y1,y2,0),...,
8. UniformAssignment(CsUxK,DELEGATED_BY) 
    #for each in CsUxk, update ABV := BEL := DELEGATED_BY
9. Delegation(CsUxK, AV)
    #for each s_i in CsUxk
    #update ABV, BEL according to av_i in AV
10. Mapping(CsUxK, tempSqdm)
    #for each s_i in CsUxk, find r=(xi,xj,yi,yj,D_i) in tempSQDM : that 
    # 'r' contains segment s_i.
    # add s_i into dictionary D_i
11. Mapping(sPsUxk, tempSqdm)
12. UpdateOriginalSqdm(P_sqdm, tempSqdm)
12.1 SET-DELEGATION_LOCK ON P_sqdm
    #for each xKey in tempSqdm, 
        #Update corresponding key in P_sqdm by value for xkey in tempSqdm.
    
12.1 RELEASE-DELEGATION_LOCK ON P_sqdm



# Papers
1. Reliable Point Location in Blockchain
2. Executing Some Computational Geometry Problems in Blockchain
3. Trustworthy Regional Delegation on Blockchain Infrastructures.
4.

Thesis work:
1. SQDM in TM
2. Reliable PL in BC
3. Reliable Line Intersection Problem / Simple Polygon Test in BC/TM
4. Trustworty Regional Delegation in Blockchain Infrastructure
    4.1. Method -1, Naresh's method
    4.2. Method -2, Dr. Ramkumar's


'''


def range_search(arr, x):
    #range search
    arr += arr[0:1]
    l = 0
    r = len(arr) -2

    while l <= r:
        mid = l + (r - l) / 2;
        # Check if x is present at mid
        #print mid
        if arr[mid] <= x and x <= arr[mid+1]:
            return arr[mid]

        # If x is greater, ignore left half
        elif arr[mid] < x:
            l = mid + 1

        # If x is smaller, ignore right half
        else:
            r = mid - 1

    # If we reach here, then the element was not present
    return -1

print("range search as binary search")
ar =[1,2,6,7,8]
print range_search(ar,1.5)


# given a point p, return the point on s that shares p's y-val
def get_x_at(self, p):
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