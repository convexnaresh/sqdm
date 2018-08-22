import random
from  osgeo import ogr, osr
class Segment(object):
    
    def __init__(self, segk,segname='inf'):
        """"""
        self.segk = segk     #tuple
        self.segname = str(segname)
        
    def getSegKP(self):
        return self.segk
    
    def getSegN(self):
        return self.segname
    
    def __str__(self):
        #return str(self.getSegKP()) + ','+ self.getSegN()
        return str((self.getSegKP()[0],self.getSegKP()[1],self.getSegN()))
        
    def overlaps(self,othersegment):
        A= ' this segment internally cuts other segment at ??'
        a,b = othersegment.getSegKP() #segment A
        c,d = self.getSegKP() #segment C
        if c >= a and c <= b:
            if d <= b and d > a:
                return (A),[c,d]
            else:
                return (A), [c]
        if  c < a:
            if d<=b and d > a:
                return (A), [d]
            else:
                return None,None
        return None,None
        
def initializeTree(segment):
    print("initialization complete")
    tree = []
    tree += [segment]
    return tree

from copy import deepcopy
def sortSegments(segments,sort_by='startvalue'):
    sorted_segments = []
    sort_byvalues=[]
    if sort_by == 'startvalue':
        for segment in segments:
            sort_byvalues +=[segment.getSegKP()[0]]   
    if sort_by == 'endvalue':
        for segment in segments:
            sort_byvalues +=[segment.getSegKP()[1]] 
    if sort_by == 'property_name':
        for segment in segments:
            sort_byvalues +=[segment.getSegN()] 
    sorted_indices = sorted(range(len(sort_byvalues)), key=lambda k: sort_byvalues[k])   
    for sorted_index in sorted_indices:
        sorted_segments +=[segments[sorted_index]]
    return sorted_segments
    
def updateLSP(tree,segment):
    LSP = []
    for seg in tree:
        LSP += list(seg.getSegKP())
    LSP += list(segment.getSegKP())
    sortedLSP = sorted(list(set(LSP))) #remove duplicate items, sort in increasing order.
    return sortedLSP

def initializeTree(segment):
    print("completion inititalization")
    tree = []
    tree += [segment]
    return tree
    
def addToTreeDynamic(tree,segment): #list of segments
    debug = True
    if debug:print("segment"),segment
    if tree:
	print("inside tree")
        fragments = []
        LSP = updateLSP(tree, segment)
        if debug: print("SLP"),LSP, len(tree)
        #Do divide all the existing lines in tree.
        for seg in tree:
            reqdLSP = []
            for point in LSP:
                if point >= seg.getSegKP()[0] and point <= seg.getSegKP()[1]:
                    reqdLSP +=[point]
            if debug:print("\t reqdLSP:"),reqdLSP
            #remove this segment and replace by fragments
            tree.remove(seg)
            for i in range(0,len(reqdLSP)-1):
                fragments +=[Segment((reqdLSP[i],reqdLSP[i+1]),seg.getSegN())]
        
        #Do divide new-segment
        reqdLSP4segment = []
        for point in LSP:
            if point >= segment.getSegKP()[0] and point <= segment.getSegKP()[1]:
                reqdLSP4segment +=[point] 
        for i in range(0,len(reqdLSP4segment)-1):
            fragments +=[Segment((reqdLSP4segment[i],reqdLSP4segment[i+1]),segment.getSegN())]        
        tree +=fragments
	if debug: print("inside tree2",len(tree))
        return tree   
    else:
        print("Initializing tree.")
        return initializeTree(segment)

def getLSP(segments):
    LSP = []
    for seg in segments:
        LSP += list(seg.getSegKP())
    sortedLSP = sorted(list(set(LSP))) #remove duplicate items, sort in increasing order.
    return sortedLSP

def getUniqueXords(listsegments):
    listuniquexords = []
    for seg in listsegments:
        listuniquexords += list(seg.getSegKP())
    return sorted(list(set(listuniquexords)))

def updateUniqueXords(olduniquexords=[],segments=[]):
    uniquexords = []
    xordsfromsegments =[]
    for segment in segments:
        xordsfromsegment += list(segment.getSegKP())
    newxords = olduniquexords + xordsfromsegments
    sortedNewXords= sorted(list(set(newxords)))
    return sortedNewXords

def getIntermediateXords(unique_xordslist,segment):
    '''returns list of intermediate xords that falls between given segment.'''
    #print("inside getIntermediateXords"),segment.getSegKP()[0],segment.getSegKP()[1]
    intermediateXords = []
    for xord in unique_xordslist[0:]:
	if xord >= segment.getSegKP()[0] and xord <= segment.getSegKP()[1]:
	    intermediateXords +=[xord]
    if len(intermediateXords)==1:
	return intermediateXords*2 #for vertical line with same xords in segment
    else:
	return intermediateXords

def calcSlope(line):
    '''Returns slope of a two points of wkbLineString object representing a twio point line or list of tuplesrepresenting two points passed as arguement'''

    if isinstance(line,ogr.Geometry):
        if line.GetGeometryName() == 'LINESTRING':
	    wkbLineString = line
	    p1 = wkbLineString.GetPoint(0)
	    p2 = wkbLineString.GetPoint(1)
	else:
	    message = "parameter should be either LineString type or list of tuples."
	    raise Exception(message)
	#line may be list with tuple of two points [(a,b),(c,d)]
	if type(line) == list:
	    p1 = line(0)
	    p2 = line(1)
	try:
	    slope = (p2[1]-p1[1])/(p2[0]-p1[0])
	    pass
	except ZeroDivisionError:
	    slope = 'inf'
    return slope

def getLutEntries(uniquexords):
    '''Returns ranges for a discrete series passed in @uniquexords arguments'''
    lutentries = []
    for i in range(len(uniquexords)-1):
        lutentries += [[uniquexords[i],uniquexords[i+1],0]]
    lutentries +=[[uniquexords[-1],uniquexords[0],0]]#winding entry
    return lutentries

def getConsequitivePairs(itemList,winding=False):
    '''Returns ranges for a discrete series passed in @itemList arguments'''
    pairs = []
    for i in range(len(itemList)-1):
        pairs += [(itemList[i],itemList[i+1])]
    if winding and itemList:
        pairs +=[[itemList[-1],itemList[0]]]#winding entry
    return pairs

def getIntermediateYords(linesegmentObj,intermediateXords):
    '''IntermediateXords contains end xords as well.Escape them before finding intermediate yords. Return list are endpoint'''
    p1 = linesegmentObj.GetPoint(0)
    p2 = linesegmentObj.GetPoint(1)
    intermediateYords = []
    for intxord in intermediateXords[1:-1]:
        intermediateYords +=[((p2[1] - p1[1])*(intxord-p1[0])/float(p2[0]-p1[0])) + p1[1]]
    if p1[0] > p2[0]:
        return [p2[1]] + intermediateYords + [p1[1]]
    else:
	return [p1[1]] + intermediateYords + [p2[1]]
    #return intermediateYords

def addToTreeDynamic2R(segments):
    LSP = getLSP(segments)
    keypairs = []
    attributes = [[] for i in range(0, len(LSP)-1)]
    for i in range(0,len(LSP) -1):
        keypairs += [(LSP[i], LSP[i+1])]
        for segment in segments:
            if segment.getSegKP()[0] <= LSP[i]  and segment.getSegKP()[1] >= LSP[i+1] :
                attributes[i].append(segment.getSegN())
    
    return keypairs, attributes

def getDicts(keypairs, attributes):
    
    t = {}
    for index in range(len(attributes)):
        for key in attributes[index]:
            try:
                t[key] += [keypairs[index]]
            except:
                t[key] = [keypairs[index]]
    return t
    
##Construct larger segments.
'''
segments = []
for i in range(0,10):
    x1 = (random.randint(0,90))      
    segments +=[Segment((x1,random.randint(x1+3,100)), i)]    #(a,b), a<b

for segment in segments:
    print(segment)          

import sys
##Construct LUTs [(x1,y1,'segment-id'),(x2,y2,'segment-id')..]; x1, y1 are start and end segments.
ttree = tree =[]
for segment in segments:
    ttree = addToTreeDynamic(ttree,segment)

sttree = sortSegments(ttree,'property_name')
print("Structure-1")
for segment in sttree:
    print segment

print("size of tree using method 1:"), len(sttree), sys.getsizeof(sttree)
keys,attlist = addToTreeDynamic2R(segments)
print
print("Structure-2")
print("Keys"), keys
print
print("Attributes"), attlist
print(len(keys),',',len(attlist), sys.getsizeof(keys)+sys.getsizeof(attlist))
print
print("Structure-3")
sum([len(litem) for litem in attlist])

print getDicts(keys, attlist)
print(sys.getsizeof(getDicts(keys, attlist)))
'''
