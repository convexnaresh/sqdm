import itertools
from osgeo import ogr
import shapely
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
import math
import time
import pickle
import numpy as np
import sys
import time
import datetime
chuncked_lines = []
unique_xordslist = []

def calcSlope(four_tupline):
    p1 = (four_tupline[0],four_tupline[1])
    p2 = (four_tupline[2],four_tupline[3])
    try:
        slope = (p2[1]-p1[1])/(p2[0]-p1[0])
    except ZeroDivisionError:
        slope = 'inf'
    return slope
                  
def getChunksOfLines(linesarr):
    global unique_xordslist
    global chuncked_lines
    x1,y1 = linesarr[0], linesarr[1]
    x2,y2 =linesarr[2], linesarr[3]
    intermediateXords = getIntermediateXordsModified(unique_xordslist,x1,x2)
    intermediateYords = getIntermediateYordsModified((x1,y1),(x2,y2),intermediateXords)
   
    try:
        chuncked_lines += getLineSegementsFromPtSeq(zip(intermediateXords, intermediateYords))
    except Exception, e:
        print("Exception when chunking."), str(e)
        chuncked_lines +=[]
    return True
    
def getIntermediateXordsModified(unique_xordslist,xord1,xord2):
    '''returns list of intermediate xords that falls between given xord1 and xord2.'''
    #print("inside getIntermediateXords"),segment.getSegKP()[0],segment.getSegKP()[1]
    xordmin, xordmax = sorted([xord1,xord2])
    intermediateXords = []
    for xord in unique_xordslist[0:]:
        if xord >= xordmin and xord <= xordmax:
            intermediateXords +=[xord]
    if len(intermediateXords)==1:
        return intermediateXords*2 #for vertical line with same xords in segment
    else:
        return intermediateXords

def getLineSegementsFromPtSeq(ptseq):
    linesegments = []
    for ptid in range(len(ptseq)-1):
        point1,point2 = ptseq[ptid], ptseq[ptid+1]
        lsTuple = (point1[0],point1[1],point2[0],point2[1])
        linesegments.append(lsTuple)
    return linesegments

def getIntermediateYordsModified(p1,p2,intermediateXords): #p1 and p2 are points in line 
    '''IntermediateXords contains end xords as well.Escape them before finding intermediate yords. Return list are endpoint'''
    global ROUNDAT
    intermediateYords = []
    for intxord in intermediateXords[1:-1]: #excluding firs and last xords.
        intermediateYords +=[round(((p2[1] - p1[1])*(intxord-p1[0])/float(p2[0]-p1[0])), ROUNDAT) + p1[1]]
        
    if p1[0] > p2[0]:
        return [p2[1]] + intermediateYords + [p1[1]]
    else:
        return [p1[1]] + intermediateYords + [p2[1]]

def reverse(lsTuple):
    x1,y1,x2,y2 = lsTuple
    return (x2,y2,x1,y1)

def getTranslationMatrix(tx,ty):
    #T = np.zeros((3,3),int)
    #np.fill_diagonal(T,1)
    #V = np.asarray([tx,ty,1]).reshape((3,1))
    #T = np.hstack((T[0:,0:2],V))   
    return np.asarray([tx,ty,tx,ty])
    
def divide_lines_into_chunks(all_polygons):    
    global chuncked_lines
    global unique_xordslist
    global FACTORING_NUMBER
    global INTERGER_CONVERSION
    chuncked_lines = []
    unique_xordslist = []
    linesegments =[]
    
    for eachPoly in all_polygons:
        try:
            ptseq = eachPoly.GetBoundary().GetPoints() #if geom obj cannot have getbondary
        except:
            ptseq = []
        for ptid in range(len(ptseq)-1):
            point1,point2 = ptseq[ptid], ptseq[ptid+1]
            lsTuple = (point1[0],point1[1],point2[0],point2[1])
            if lsTuple not in linesegments and reverse(lsTuple) not in linesegments :
                linesegments.append(lsTuple)
    
    #find max value in linesegment's tuples, divide by 1 Billion, get factor.
    linesegmentarr = np.asarray(linesegments)
    linesegmentarr = linesegmentarr.astype(np.uint32)
    #set param.
    MAP_MAX_VALUE= linesegmentarr.max() #max(list(linesegmentarr.max(axis=0)))
    MAP_MIN_VALUE = linesegmentarr.min()
    if MAP_MIN_VALUE < 0:
        MINX = min(linesegmentarr[0:,0].min(), linesegmentarr[0:,2].min())
        MINY = min(linesegmentarr[0:,1].min(), linesegmentarr[0:,3].min())
        TM = getTranslationMatrix(abs(MINX), abs(MINY))
        linesegmentarr = linesegmentarr + TM
    #print MAP_MIN_VALUE,linesegmentarr 
    #MAKE THEM INTEGERS
    if INTERGER_CONVERSION:
        linesegmentarr = linesegmentarr * FACTORING_NUMBER
        linesegmentarr = np.round(linesegmentarr, ROUNDAT)
    #make unique_xords
    unique_xordslist = sorted(np.unique(np.hstack((linesegmentarr[0:,0:1],linesegmentarr[0:,2:3]))))
    
    #For all chunked lines
    try:
        np.apply_along_axis(getChunksOfLines, 1, linesegmentarr) #apply along row.
        unique_xordslist = []
    except Exception,e:
        print("Exception in dividing_lines_into_chunks. Total lines Nill")
        return None,None, None
    
    return chuncked_lines, len(linesegmentarr)

def env(polygons):
    maxx = 0
    maxy = 0
    minx = sys.maxint
    miny = sys.maxint

    for p in polygons:
        e = p.GetEnvelope()
        maxx = max(maxx,e[0],e[2])
        maxy = max(maxy,e[1],e[3])
        minx = min(minx,e[0],e[2])
        miny = min(miny,e[1],e[3])
    print minx,miny, maxx,maxy    
    
def constructDls(chuncked_lines):
    dls = {}
    for chunked_line  in chuncked_lines:
        x1,y1,x2,y2 = chunked_line
        if x1 != x2:
            try:
                key = min(x1,x2)
                dls[key] = dls[key] + [(y1,y2)]
            except Exception, e:
                dls[key] = [(y1,y2)]
        else:
            pass
    return dls

def savePickle(object_tosaved,outfile):
    try:
	os.remove(outfile)
    except:
	pass
    with open(outfile, 'wb') as fp:
        pickle.dump(object_tosaved, fp)
    print("Object saved in"), outfile
    return True

def getPolygons(datatype='datatype'):
    #sample polygon           
    if datatype =="sample":
        #Manual Polygons
        wkt0 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 7 3,11 4,11 2, 13 2, 13 7, 10 10,6 10,2 5 ))"
        wkt1 = "POLYGON ((2 5, 3 4, 4 4, 5 3, 4 1,2 0, 1 1,1 3,2 5 ))"
        wkt2 = "POLYGON ((5 3, 7 1, 10 1, 11 2,11 4,7 3,5 3 ))"
        poly0 = ogr.CreateGeometryFromWkt(wkt0)
        poly1 = ogr.CreateGeometryFromWkt(wkt1)
        poly2 = ogr.CreateGeometryFromWkt(wkt2)
        polyref = ogr.CreateGeometryFromWkt(wkt0)
        all_polygons = [poly0,poly1,poly2] #list of all the polygonsA
    return all_polygons

if __name__ == "__main__":
    import DataSources as mydatasource
    reload(mydatasource)
    #Constants
    FACTORING_NUMBER = 100 #values divide this number to rescale to integer.
    cchuncked_lines = []
    ROUNDAT = 0 #to convert to int, make it 0
    INTERGER_CONVERSION = True
    
    #Load Data.
    dlslist = []
    datatype = 'mscounty'
    t0,t1= 0,0
    polygons, inputspatialref = mydatasource.getPolygons(datatype=datatype)
    print len(polygons)
    for pid in [0]:#,63,64,65,66,67,68,69]:#range(len(polygons)):
        print("pid:"), pid
        polygon = [polygons[pid]]
        try:
            t0 = time.time()
            cchuncked_lines, cntoriglines = divide_lines_into_chunks(polygon) 
            t1=time.time()
        except:
            continue
        dls = constructDls(cchuncked_lines)
        dlslist =[(datatype+','+str(pid)+","+str(cntoriglines), dls)]
        print("Constructed a dictionary type DLS.time:div,dlsconst:"),round(t1-t0, 3), ',',round(time.time()-t1, 3)
        print
        dlsfile = './results/' +datatype+ str(pid) + '_newex.txt'
	savePickle(dlslist,dlsfile)
    print
    '''
    #Create DLS for whole map.|
    dlslist = []
    cchuncked_lines = []
    unique_xordslist= []
    try:
        t0= time.time()
        cchuncked_lines,cntoriglines  = divide_lines_into_chunks(polygons)
        t1 = time.time()
    except Exception, e:
        print("Exception in map dls."), str(e)
    
    dls = constructDls(chuncked_lines)
    dlslist +=[(str(pid)+","+str(cntoriglines),dls)]
    dlsfile = './results/' + datatype+'_map' + '.txt'
    print("Constructed a dictionary type DLS."),round(t1-t0, 3), ',',round(time.time()-t1, 3)
    
    savePickle(dlslist,dlsfile)
    t1 = time.time()
    print("Map max-value, min-value"), MAP_MAX_VALUE,MAP_MIN_VALUE
    '''

