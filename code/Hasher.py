import hashlib
import sha3
from collections import OrderedDict

class Helper(object):
    """Helper Class for Keeping track of Insert Order"""
    def __init__(self, arg):
        super(Helper, self).__init__()

    dictContainer = dict()
    ordering = list()

    @staticmethod
    def addItem(dictItem):
        for key,value in dictItem.iteritems():
            #print key,value
            Helper.ordering.append(key)
            Helper.dictContainer[key] = value

    @staticmethod
    def getPreviousValue(key):
        try:
            index = (Helper.ordering.index(key)-1)
            return Helper.dictContainer[Helper.ordering[index]]
        except:
            return None

    @staticmethod
    def getPreviousKey(key):
        try:
            index = (Helper.ordering.index(key)-1)
            return Helper.ordering[index]
        except:
            return None

    @staticmethod
    def getNextKey(key):
        try:
            index = (Helper.ordering.index(key)+1)
            return Helper.ordering[index]
        except:
            return None

#Your unordered dictionary
d = {'aaaa': 'a', 'bbbb':'b', 'cccc':'c', 'dddd':'d', 'eeee':'e', 'ffff':'f'}

#Create Order over keys
ordered = OrderedDict(sorted(d.items(), key=lambda t: t[0]))

#Push your ordered list to your Helper class
Helper.addItem(ordered)

print ordered
#Get Previous of
print Helper.getPreviousValue('eeee')
print Helper.getNextKey('eeee')

class Hasher(object):

    def __init__(self):
        self.hasher= sha3.sha3_224()

    def getHasher(self):
        return hashlib.sha3_224()

    def getHasherSha333(self):
        return sha3.sha3_224()

    def hash_segment(self,segment):
        s = self.getHasher()
        for e in segment:
            h = Hasher().getHasher()
            h.update(str(e))
            s.update(h.hexdigest()[0:12])
        return s.hexdigest()[0:12]
