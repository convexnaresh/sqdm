class St:

    def __init__(self,key=None,nkey=None,val=None,type_=None,posx=None,posy=None,colid=None,dic={}):
        '''A leaf leval entity is has a key, next-key, value, type as main components.
        When a hash table is generated for a collection of leaf level entities, the hashes generated for an
        entity has a position of (x,y) in a hash table. 
        x-denotes breadth position in a binary tree of hashes.
        y-denotes depth positon in a binary tree of hashes.
        
        When ever an entry splited by an insertion, two hashes are generated at positions (2*x,y+1) and (2*x+1,y+1).
        The oritinal hash is updated by single hash of these two hashes.
        '''
        self.key = key
        self.nkey = nkey
        self.val = val
        self.type_= type_
        self.x = posx
        self.y = posy
        self.colid=colid #collection id
        if dic:
            self.key = dic["key"]
            self.nkey = dic["nkey"]
            self.val = dic["val"]
            self.type_= dic["type_"]
            self.x = dic["posx"]
            self.y = dic["posy"]
            self.colid=dic["colid"]
        
    def castToDict(self):
        return {"key":self.key,"nkey":self.nkey,"val":self.val,"type_":self.type_,"posx":self.x,"posy":self.y
                ,"colid":self.colid
               }
    
    def __str__(self):
        return str(self.castToDict())
if __name__ == '__main__':
    pass