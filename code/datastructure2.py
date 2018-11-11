import hashlib
import sha3
from collections import OrderedDict

from Hasher import Hasher
class St:

    def __init__(self, key=None, nkey=None, ykey=None,ynkey=None,val={}, type_=None, posx=None, posy=None, colid=None, dic={}):
        '''A leaf leval entity has a key, next-key, value, type as main components.
        When a hash table is generated for a collection of leaf level entities, the hashes generated for an
        entity has a position of (x,y) in a hash table.
        x-denotes breadth position in a binary tree of hashes.
        y-denotes depth positon in a binary tree of hashes.

        When ever an entry splited by an insertion, two hashes are generated at positions (2*x,y+1) and (2*x+1,y+1).
        The oritinal hash is updated by single hash of these two hashes.
        '''
        try:
            if dic:
                self.xkey = dic["xkey"]
                self.xnkey = dic["xnkey"]
                self.val = dic["val"]
                self.type_ = dic["type_"]
                self.x = dic["posx"]
                self.y = dic["posy"]
                self.colid = dic["colid"]
                self.ykey = dic["ykey"]
                self.ynkey = dic["ynkey"]

            else:
                self.xkey = key
                self.xnkey = nkey
                self.ykey = ykey
                self.ynkey = ynkey
                self.val = val
                self.type_ = type_
                self.x = posx
                self.y = posy
                self.colid = colid  # collection id
        except Exception,e:
            print("Exception:"),str(e)

    def castToDict(self):
        return {"xkey": self.xkey,
                "xnkey": self.xnkey,#key for x
                "ykey": self.ykey,#next-key for y
                "ynkey": self.ynkey,
                "val": self.val,
                "type_": self.type_,
                "posx": self.x,
                "posy": self.y,
                "colid": self.colid
                }

    def __str__(self):
        s ='{'
        s +='xkey: '+str(self.xkey)+', '
        s += 'xnkey: ' + str(self.xnkey)+ ', '
        s += 'ykey: ' + str(self.ykey)+ ', '
        s += 'ynkey: ' + str(self.ynkey) + ', '
        s += 'val: ' + str(self.val)+ ', '
        s += 'type_: ' + str(self.type_)+ ', '
        s += 'posx: ' + str(self.x) + ', '
        s += 'posy: ' + str(self.y) + ', '
        s += 'colid: ' + str(self.colid)
        s += '}'
        return s

    def h(self,str1=None):
        s = Hasher().getHasher()
        if str1:
            args =[str1]
        else:
            args = [self.xkey,self.xnkey,self.ykey, self.ynkey,self.val,self.type_]

        for p in args:
            if type(p) is dict:
                p = self.serialize(p)
            hasher = Hasher().getHasher()
            hasher.update(str(p))
            s.update(hasher.hexdigest()[0:12])
        return s.hexdigest()[0:12]


    def serialize(self,obj):
        if len(obj) > 0:
            ordd =OrderedDict(obj)
            s = '{'
            for k,v in ordd.items():
                s +=str(k)+":"+str(v)+","
            return '}'+s[-1] #exclude last ','
        return 0

if __name__ == '__main__':

    pass