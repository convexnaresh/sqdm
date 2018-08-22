
# coding: utf-8

# In[1]:

#find out all y-ords that cuts the line AB
def getAllYordsForLine(dls,lineid,linestring_geos):
    '''
    Example:
        print dls
        lineid = 0
        poly_pts_sequence = getPolygonPoints(datatype='sample', pid=0, debug=False)
        linestring_geos = getLinesForPolygonPointsSequences2(poly_pts_sequence)[lineid]
        getAllYordsForLine(dls,lineid,linestring_geos)  
    '''
    point1,point2 = linestring_geos.GetPoints()
    xord1,xord2 = [point1[0], point2[0]]
    unique_xords = dls.keys()
    intermediate_xords = getIntermediateXordsModified(unique_xords,xord1,xord2)
    yordsList = []
    for xord in intermediate_xords[:-1]:#exclude last x-ord
        line_id_yords_dict = dls[xord]
        yords = line_id_yords_dict[lineid]
        yordsList.append(yords)
    return intermediate_xords[0:-1],yordsList

def getYordsForXords(dls):
    '''
    Example:
        xords,yords = getYordsForXords(dls)
        print xords
        print yords
    '''
    vcolumns = []
    for kxord, vdic_lid_yords in dls.items():
        vcol = []
        for lineid,yords in vdic_lid_yords.items():
            vcol +=[(tuple(yords + [lineid]))]
        vcolumns.append(vcol)
    return dls.keys(),vcolumns
        


# In[ ]:



