
# coding: utf-8

# In[6]:




# In[27]:

import pickle
import sys
class DlsStat(object):
    
    def __init__(self,dls={},id_='id'):
        self.dlsid = str(id_)
        self.stats = self.getDlsAttributes2(dls)
        
    def setDls(self,dls):
        self.stats = self.getDlsAttributes(dls)
    
    def getDlsAttributes1(self,dls={}):
        attrtable = []
        yords_cols = []
        count_boundary_recs = 0
        for xord,line_id_yords_dic in dls.items():
            vcolumn = []
            count_boundary_recs += len(line_id_yords_dic)
            for lineid_ask, yord_list in line_id_yords_dic.items():
                #lineid, slope,polyid1,polyid2
                attrstr= str(lineid_ask) + ',' + yord_list[-1]
                attrtup = tuple(yord_list[:-1]) + (attrstr,) #this comma is required.
                attrtable += [attrstr]
                vcolumn.append(attrtup)
            yords_cols.append(vcolumn)

        count_neg_diag = 0
        for attrstr in attrtable:
            if float(attrstr.split(",")[1]) < 0: #count slope < 0
                count_neg_diag +=1
        return stats
    
    def getDlsAttributes2(self,dls={}):
        total_shared_rect,total_pure_rect, total_overlapping_rect = self.getRectTypeAndCount2(dls)
        
        count_boundary_recs = 0
        for xord,yords_pairs in dls.items():
            count_boundary_recs += len(yords_pairs)
        
        count_lutentries= count_boundary_recs +total_pure_rect +(len(dls.keys())-1)
        stats = {}
        stats['shape_name'] = self.dlsid.split(',')[0]
	stats['polyid'] = self.dlsid.split(',')[1]
        try:
            stats['#origlines'] = self.dlsid.split(',')[2]
        except:
            stats['#origlines'] = 'NA'
        stats['attribute_table'] = []
        stats['#neg_diag'] = 'NA'
        stats['#pos_diag'] = 'NA'
        stats['#pure_recs'] = 'NA'
        stats['#boundary_recs'] = count_boundary_recs
        stats['#xords'] = len(dls.keys())
        stats['dls_bytes'] = sys.getsizeof(dls)
        stats['lutentries_bytes'] = sys.getsizeof(dls.values()[0][0])*count_lutentries
        stats['#lutentries'] = count_lutentries
        stats['#shared_recs'] = total_shared_rect
        stats['#pure_rect'] = total_pure_rect
        stats['#overlapping_rect'] = total_overlapping_rect
        return stats

    def getRectTypeAndCount1(self,dls):
        '''
        dls = [(0, {2.0: {7: [2.0, 2.5, '-1.0,1,K'],
                          9: [2.0, 3.0, '2.0,1,K']},
                    4.0: {3: [1.0, 2.0, '-0.0,0,K'],
                          11: [1.0, 2.5, '0.25,0,2'],
                          14: [5.0, 5.5, '0.0,2,K'],
                          15: [0.5, 0.5, '0.0,2,K']},

                    6.0: {2: [4.5, 5.0, '-1.0,0,K'],
                          17: [0.5, 0.5, '0.0,2,K']},
                    8.0: {0: [4.5, 5.5, '0.0,0,K'],
                          19: [0.5, 0.5, '0.0,2,K']
                          }})]
        '''
        for key, value in dls.iteritems():
            def getEachRectTypeAndCount(lines):
                listOfLines = []
                # inner dict for lines
                for key, value in lines.iteritems():
                    listOfLines.append([value[0], value[1]])
                rectCountList = []
                shared_rect = 0
                pure_rect = 0
                overlapping_rect = 0

                # pure rectangles
                same_y_coord = []
                for line in listOfLines:

                    if(line[0] == line[1]):
                        same_y_coord.append(line)
                        pure_rect += 1

                for y in same_y_coord:
                    listOfLines.remove(y)

                diff_y_coord = []
                for line in listOfLines:
                    isOverlapping = False
                    for l in listOfLines:
                        import copy
                        temp_line = copy.deepcopy(line)
                        joined_l = list(set(temp_line+l))
                        # overlapping rectangles and shared rectangles
                        if len(joined_l) == 3:
                            joined_l.sort()
                            y_mid = joined_l[1]
                            if (temp_line+l).count(y_mid) != 2:
                                listOfLines.remove(l)
                                isOverlapping = True
                        # pure rectangles
                        elif len(joined_l) == 4:
                            temp_line.sort()
                            l.sort()
                            def hasBetween(a, b):
                                for lin in listOfLines:
                                    #lin.sort()
                                    #print lin[0],'>',a, 'or', lin[1],'<',b
                                    if lin[0]>=a and lin[1]<=b:
                                        return True
                                return False
                            if temp_line[1] < l[0] and not hasBetween(temp_line[1], l[0]):
                                #print temp_line, l
                                pure_rect += 1
                    diff_y_coord.append(line)
                    if isOverlapping == True:
                        overlapping_rect += 1
                    else:
                        shared_rect += 1
                for y in diff_y_coord:
                    listOfLines.remove(y)
                return shared_rect, pure_rect, overlapping_rect

            shared_rect, pure_rect, overlapping_rect                 = getEachRectTypeAndCount(value)
            listOfRectCount.append([key, shared_rect, pure_rect, overlapping_rect])

        total_shared_rect = 0
        total_pure_rect = 0
        total_overlapping_rect = 0

        for i in listOfRectCount:
            total_shared_rect += i[1]
            total_pure_rect += i[2]
            total_overlapping_rect += i[3]

        return total_shared_rect,total_pure_rect, total_overlapping_rect    
    
    def getRectTypeAndCount2(self,dls):
        '''
        #Handles dictionary of following type.
            dls = [(0, {2.0: [[2.0, 2.5, '-1.0,1,K'],[2.0, 3.0, '2.0,1,K']],
            4.0: [[1.0, 2.0, '-0.0,0,K'],
                  [1.0, 2.5, '0.25,0,2'],
                  [5.0, 5.5, '0.0,2,K'],
                  [0.5, 0.5, '0.0,2,K']],
            6.0: [[4.5, 5.0, '-1.0,0,K'],
                  [0.5, 0.5, '0.0,2,K']],
            8.0: [[4.5, 5.5, '0.0,0,K'],
                  [0.5, 0.5, '0.0,2,K']]})]
        '''
        listOfRectCount = []
        for key, value in dls.iteritems():
            def getEachRectTypeAndCount(values):
                listOfLines = []
                # inner dict for lines
                for value in values:
                    listOfLines.append([value[0], value[1]])
                rectCountList = []
                shared_rect = 0
                pure_rect = 0
                overlapping_rect = 0

                # pure rectangles
                same_y_coord = []
                for line in listOfLines:

                    if(line[0] == line[1]):
                        same_y_coord.append(line)
                        pure_rect += 1

                for y in same_y_coord:
                    listOfLines.remove(y)

                diff_y_coord = []
                for line in listOfLines:
                    isOverlapping = False
                    for l in listOfLines:
                        import copy
                        temp_line = copy.deepcopy(line)
                        joined_l = list(set(temp_line+l))
                        # overlapping rectangles and shared rectangles
                        if len(joined_l) == 3:
                            joined_l.sort()
                            y_mid = joined_l[1]
                            if (temp_line+l).count(y_mid) != 2:
                                listOfLines.remove(l)
                                isOverlapping = True
                        # pure rectangles
                        elif len(joined_l) == 4:
                            temp_line.sort()
                            l.sort()
                            def hasBetween(a, b):
                                for lin in listOfLines:
                                    #lin.sort()
                                    #print lin[0],'>',a, 'or', lin[1],'<',b
                                    if lin[0]>=a and lin[1]<=b:
                                        return True
                                return False
                            if temp_line[1] < l[0] and not hasBetween(temp_line[1], l[0]):
                                #print temp_line, l
                                pure_rect += 1
                    diff_y_coord.append(line)
                    if isOverlapping == True:
                        overlapping_rect += 1
                    else:
                        shared_rect += 1
                for y in diff_y_coord:
                    listOfLines.remove(y)
                return shared_rect, pure_rect, overlapping_rect

            shared_rect, pure_rect, overlapping_rect                 = getEachRectTypeAndCount(value)
            listOfRectCount.append([key, shared_rect, pure_rect, overlapping_rect])

        total_shared_rect = 0
        total_pure_rect = 0
        total_overlapping_rect = 0

        for i in listOfRectCount:
            total_shared_rect += i[1]
            total_pure_rect += i[2]
            total_overlapping_rect += i[3]
        return total_shared_rect,total_pure_rect, total_overlapping_rect

    def getSummary(self, attlist):
        summary = {}
        valstr = ''
        for attr in attlist:
            try:
                valstr = valstr +','+ str(self.stats[attr])
            except:
                valstr = valstr +','+ str('NA')
        return valstr[1:] #remove initial ,

    def __add__(self,other):
        
        sumDlsStat = DlsStat()
        return sumDlsStat
    
    def readDlsSavedPickle(self,dlsfile):
        with open(dlsfile, 'rb') as fp:
            dls = pickle.load(fp)
        if type(dls) == list:
            self.setDls(dls[0][1]) #reads first dls.
        else:
            self.setDls(dls)
            
    
    def test(self):
        dls = [(0,
                {1.0: {7: [1.0, 0.0, '-1.0,1,K'], 9: [3.0, 5.0, '2.0,1,K']},
                 2.0: {4: [5.0, 6.25, '1.25,0,K'],
                      6: [0.0, 0.5, '0.5,1,K'],
                      16: [5.0, 4.0, '-1.0,0,1']},
                 3.0: {4: [6.25, 7.5, '1.25,0,K'],
                      6: [0.5, 1.0, '0.5,1,K'],
                      17: [4.0, 4.0, '0.0,0,1']},
                 4.0: {4: [7.5, 8.75, '1.25,0,K'],
                      5: [1.0, 3.0, '2.0,1,K'],
                      18: [4.0, 3.0, '-1.0,0,1']},
                 5.0: {4: [8.75, 10.0, '1.25,0,K'],
                      10: [3.0, 3.0, '0.0,0,2'],
                      13: [3.0, 2.0, '-1.0,2,K']},
                 6.0: {3: [10.0, 10.0, '-0.0,0,K'],
                      10: [3.0, 3.0, '0.0,0,2'],
                      13: [2.0, 1.0, '-1.0,2,K']},
                 7.0: {3: [10.0, 10.0, '-0.0,0,K'],
                      11: [3.0, 3.75, '0.25,0,2'],
                      14: [1.0, 1.0, '0.0,2,K']},
                 10.0: {2: [10.0, 9.0, '-1.0,0,K'],
                      11: [3.75, 4.0, '0.25,0,2'],
                      15: [1.0, 2.0, '1.0,2,K']},
                 11.0: {0: [2.0, 2.0, '0.0,0,K'], 
                        2: [9.0, 7.0, '-1.0,0,K']},
                 13.0: {}})]
        d = DlsStat()
        print d.stats

def getdls():
        dls = [(0,
                {1.0: [[1.0, 0.0, '-1.0,1,K'],[3.0, 5.0, '2.0,1,K']],
                 2.0: [[5.0, 6.25, '1.25,0,K'],
                      [0.0, 0.5, '0.5,1,K'],
                      [5.0, 4.0, '-1.0,0,1']],
                 3.0: [[6.25, 7.5, '1.25,0,K'],
                      [0.5, 1.0, '0.5,1,K'],
                      [4.0, 4.0, '0.0,0,1']],
                 4.0: [[7.5, 8.75, '1.25,0,K'],
                      [1.0, 3.0, '2.0,1,K'],
                      [4.0, 3.0, '-1.0,0,1']],
                 5.0: [[8.75, 10.0, '1.25,0,K'],
                      [3.0, 3.0, '0.0,0,2'],
                      [3.0, 2.0, '-1.0,2,K']],
                 6.0: [[10.0, 10.0, '-0.0,0,K'],
                      [3.0, 3.0, '0.0,0,2'],
                      [2.0, 1.0, '-1.0,2,K']],
                 7.0: [[10.0, 10.0, '-0.0,0,K'],
                      [3.0, 3.75, '0.25,0,2'],
                      [1.0, 1.0, '0.0,2,K']],
                 10.0: [[10.0, 9.0, '-1.0,0,K'],
                      [3.75, 4.0, '0.25,0,2'],
                      [1.0, 2.0, '1.0,2,K']],
                 11.0: [[2.0, 2.0, '0.0,0,K'], 
                        [9.0, 7.0, '-1.0,0,K']],
                 13.0: []})]    
        return dls
        
if __name__ == "__main__":
    
    attlist = ['shape_name','polyid','dls_bytes','#origlines','#xords','#boundary_recs','#shared_recs','#pure_rect','#overlapping_rect','#lutentries']
    print ['filename','filesize'] +attlist
    print     
    import os,sys
    files = os.listdir('./results')
    files = sorted([os.path.join('./results',f) for f in files if '_newex' in f])
    
    for dlsfile in files[0:]:
        dlstup =[]
        with open(dlsfile, 'rb') as fp:
            dlstup = pickle.load(fp)
        for tup in dlstup:
            pid,dls = tup[0], tup[1]
	    try:
                dlsStat = DlsStat(dls,pid)
                print dlsfile +","+ str(os.path.getsize(dlsfile))+'\t'+ dlsStat.getSummary(attlist)
	    except Exception, e:
		print("Exception"), str(e),dlsfile,pid
		continue
    
