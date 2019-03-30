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
'''
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

sqdm_prop = load_data("../out/tmp/usa_sqdm_properties.json")
#sqdm_prop = load_data("../out/tmp/fsqdm_properties.json")

print sqdm_prop.keys()
print

a =sqdm_prop["count_xcolumns_nsplits_yblocks"]
dic_count_nsplit_nyblocks = OrderedDict(sorted(a.items())) #sort from lowkey to high key.

seg_delx = sqdm_prop['seg_delx']
seg_dely = sqdm_prop['seg_dely']
unsignedseg_dely = [abs(int(dely)) for dely in seg_dely]
seg_lengths=sqdm_prop['seg_lengths']
split_count = sqdm_prop['split_counts']
slopes =[]
for dely in seg_dely:
    if dely <0:
        slopes +=[2]
    elif dely > 0:
        slopes +=[4]
    else:
        slopes += [16]
##
cnyb = []
cns = []
for xkey,ns_nyb in dic_count_nsplit_nyblocks.items():
    ns, nyb = ns_nyb
    cnyb +=[nyb]
    cns +=[ns]

print max(cns),min(cns), max(cnyb), min(cnyb)

x = cns

figtitle="Histogram-Segment's split counts"+"(N="+str(len(x))+")"
figtitle="Histogram-Segment's X-spans"+"(N="+str(len(x))+")"
figtitle="Histogram-Segment Lengths"+"(N="+str(len(x))+")"
figtitle="Histogram-Number Of Y-Blocks in Vertical Columns"+"(N="+str(len(x))+")"
figtitle="Histogram-Number Of Splits Segments in Vertical Columns"+"(N="+str(len(x))+")"
xvarname = "Segment's Length"
xvarname = "Number of Y-Blocks"
xvarname = "Number of Splits"
print("max,min,len"), len(x), max(x), min(x)
##
nbins=10
#x = [math.log((v+1),2) for v in x]

N, bins, patches = plt.hist(x, color='#0504aa', alpha=0.7, rwidth=1, bins=nbins)
fracs = N / N.max()
# we need to normalize the data to 0..1 for the full range of the colormap
norm = colors.Normalize(fracs.min(), fracs.max())
# Now, we'll loop through our objects and set the color of each accordingly
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.viridis(norm(thisfrac))
    thispatch.set_facecolor(color)
plt.ylabel('Frequency')
#plt.xlabel('$Log_2$'+'('+xvarname+')')
plt.xlabel(xvarname)
plt.title(figtitle)
#plt.savefig("../out/tmp/"+figtitle+".png")
#plt.show()

import numpy as np
import matplotlib.pyplot as plt

n=100
N = len(cns[0:n])

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)


rects1 = ax.bar(ind, cns[0:n], width, color='royalblue')
rects2 = ax.bar(ind+width, cnyb[0:n], width, color='seagreen')

# add some
ax.set_ylabel('Vertical Columns at X')
ax.set_title('Number of Segments and Y-Blocks in each Vertical Columns')
#ax.set_xticks(ind + width / 2)
#ax.set_xticklabels( [str(k)[0:3] for k in dic_count_nsplit_nyblocks.keys()[0:n]] )

ax.legend( (rects1[0], rects2[0]), ('#SegmentSplits', '#Y-Blocks') )

#plt.show()
def stacked_bar(data, series_labels, category_labels=None,
                show_values=False, value_format="{}", y_label=None,
                grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        axes.append(plt.bar(ind, row_data, bottom=cum_size,
                            label=series_labels[i]))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)

    if y_label:
        plt.ylabel(y_label)

    plt.legend()

    #if grid:
    #    plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2,
                         value_format.format(h), ha="center",
                         va="center")
import numpy as np
import matplotlib.pyplot as plt
plt.figure(figsize=(20, 20))

series_labels = ['# Splits', '# Y-Blocks']

data = [
    cns[0:60000],
    cnyb[0:60000]
]

category_labels = [] #['Cat A', 'Cat B', 'Cat C', 'Cat D']

stacked_bar(
    data,
    series_labels,
    category_labels=category_labels,
    show_values=False,
    value_format="{:.1f}",
    y_label="Numbers/Counts"
)
plt.title("Bar-Diagram: # of Splits & # of Y-Blocks in Vertical Columns")
plt.xlabel("Vertical Columns of X-Values")
plt.savefig('../out/tmp/Bar-Diagram- # of Splits & # of Y-Blocks.png')
plt.show()
'''

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