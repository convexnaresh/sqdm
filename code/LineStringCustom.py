
# coding: utf-8

# In[ ]:

class LineStringCustom(object):
    
    def __init__(self,geos_line_string=None,attr=[]):
        self.geosls =geos_line_string
        self.attr =attr
        
    def customizeGeosLStrings(self,list_geos_ls,commonattr):
        LSCustomObjects = []
        
        for geos_ls in list_geos_ls:
            LSCustomObjects +=[LineStringCustom(geos_ls,commonattr)]
        return LSCustomObjects
    
    
    def __str__(self):
        p1,p2 = self.geosls.GetPoints()
        return str((p1,p2,self.attr))

