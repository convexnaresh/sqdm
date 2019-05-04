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


# Python code to implement iterative Binary
# Search.

# It returns location of x in given array arr
# if present, else returns -1
def binarySearch(arr, l, r, x):
    while l <= r:

        mid = l + (r - l) / 2;

        # Check if x is present at mid
        if arr[mid] == x:
            return mid

        # If x is greater, ignore left half
        elif arr[mid] < x:
            l = mid + 1

        # If x is smaller, ignore right half
        else:
            r = mid - 1

    # If we reach here, then the element
    # was not present
    return -1


# Test array
arr = [1,2,3,10]
x = 10

# Function call
result = binarySearch(arr, 0, len(arr) - 1, x)

if result != -1:
    print "Element is present at index % d" % result
else:
    print "Element is not present in array"


def yrange_search(arr, search_key1,search_key2):
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

print("------")
arr= [(1993404080.6895585, 1993406082.342383), (1993406082.342383, 1999368516.2586746), (1999368516.2586746, 1999372733.7514358), (1999372733.7514358, 1999544026.563476), (1999544026.563476, 1999545826.5626845), (1999545826.5626845, 1999769654.502416), (1999769654.502416, 1999771449.6001568), (1999771449.6001568, 2000993400.823986), (2000993400.823986, 2000994573.8145158), (2000994573.8145158, 2001208160.0685017), (2001208160.0685017, 2001209315.9731112), (2001209315.9731112, 2001288822.3559847), (2001288822.3559847, 2001289610.7142859), (2001289610.7142859, 2001295074.3167095), (2001295074.3167095, 2001295420.331931), (2001295420.331931, 2001316356.373553), (2001316356.373553, 2001317976.350815), (2001317976.350815, 2001326718.317555), (2001326718.317555, 2001330525.0418687), (2001330525.0418687, 2003619352.891878), (2003619352.891878, 2003620556.353051), (2003620556.353051, 2004656897.6541672), (2004656897.6541672, 2004658709.8808515), (2004658709.8808515, 2005022361.9470072), (2005022361.9470072, 2005041613.5880013), (2005041613.5880013, 2005118064.7964506), (2005118064.7964506, 2005136648.5187082), (2005136648.5187082, 2006875343.2067804), (2006875343.2067804, 2006876662.5516841), (2006876662.5516841, 2146202637.04653), (2146202637.04653, 2146203950.999924), (2146203950.999924, 2148055422.028144), (2148055422.028144, 2148056258.7551184), (2148056258.7551184, 2148239163.322128), (2148239163.322128, 2148242848.0), (2148242848.0, 2148261139.7779994), (2148261139.7779994, 2148472110.568072), (2148472110.568072, 2148476562.3490524), (2148476562.3490524, 2148776796.0918427), (2148776796.0918427, 2148782309.8709307), (2148782309.8709307, 2150306529.069954), (2150306529.069954, 2150307312.6930666), (2150307312.6930666, 2150941248.549057), (2150941248.549057, 2150945225.7185845), (2150945225.7185845, 2154294645.5904617), (2154294645.5904617, 2154298991.284665), (2154298991.284665, 2158083676.5801497), (2158083676.5801497, 2158084416.8564086), (2158084416.8564086, 2158378034.098157), (2158378034.098157, 2158385663.319034), (2158385663.319034, 2159784066.718593), (2159784066.718593, 2159799323.551218), (2159799323.551218, 2159949323.8395), (2159949323.8395, 2159956019.6727138), (2159956019.6727138, 2160317488.510711), (2160317488.510711, 2160320147.422396), (2160320147.422396, 2160426358.360983), (2160426358.360983, 2160434206.2532525), (2160434206.2532525, 2160449350.2683377), (2160449350.2683377, 2160458581.399523), (2160458581.399523, 2173955605.3637404), (2173955605.3637404, 2173955613.248858)]
arr =[(2005218277, 2005223924), (2005223924, 2143484807), (2143484807, 2143495596), (2143495596, 2143544922), (2143544922, 2143559482), (2143559482, 2144587976), (2144587976, 2144593427), (2144593427, 2144792124), (2144792124, 2144822537.0), (2144822537.0, 2145375622), (2145375622, 2145388920), (2145388920, 2151176089L), (2151176089L, 2151194366L), (2151194366L, 2152881614L), (2152881614L, 2152919399L), (2152919399L, 2153295523L), (2153295523L, 2153308081L), (2153308081L, 2153928300L), (2153928300L, 2153929157L), (2153929157L, 2157078748L), (2157078748L, 2157081010L), (2157081010L, 2157396092L), (2157396092L, 2157409981L), (2157409981L, 2157498846L), (2157498846L, 2157517297L), (2157517297L, 2158273041L), (2158273041L, 2158283023L), (2158283023L, 2159153298L), (2159153298L, 2159164663L)]

for t in arr:
    print t
print
print("."), len(arr)
search_key1 = 2005218276.0
search_key2 =  2005223923.0

findex = yrange_search(arr, search_key1, search_key2)
if findex >= 0:
    print arr[findex]

for t in arr:
    if (t[0] <= search_key1 and search_key1 <= t[1]) and \
            (t[0] <= search_key2 and search_key2 <= t[1]):
        print("found"), t

def yrange_search(arr, search_key1,search_key2):
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
def xrange_search(arr, search_key1, search_key2):
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

        if (arr[mid] <= search_key1 and search_key1 < arr[mid + 1]) and \
                (arr[mid] <= search_key2 and search_key2 < arr[mid + 1]):
            return mid

        # If x is greater, ignore left half
        elif arr[mid] < search_key1:
            l = mid + 1

        # If x is smaller, ignore right half
        else:
            r = mid - 1

    # If we reach here, then the element was not present
    return -1

print
print
from simple_polygon import util
a1 =[1,5,6,9,79,89,900,8949,499948]

a = list(util.pairwise(a1))
print a
k1 = 499948
k2 = 499948
print yrange_search(a,k1,k2), xrange_search(a1,k1,k2)
