import hashlib
import sys
import copy
STYPE = {'DIC':'DIC','LUT':'LUT'} #structure types.
EVENT_TYPE={'START':'START','STOP':'END',"X":"X"}
class S:
    #class variables.
    hindexi = 0  # these two takes a regular hash value (224 characters for eg.) from character index hindexi to hindexj
    hindexj = 12
    def __init__(self,stype=STYPE['LUT'],k=None,nk=None,val=None,x=None,y=None,collid=None,pcollid=None):
        '''By default this Structure S type os 0 which means to be a dictionary.'''
        self.k = k
        self.nk = nk
        self.val = val
        self.x = x
        self.y = y
        self.stype = stype
        self.downtree = None #it is a reference to an OMT.
        self.collid = collid
        self.pcollid= pcollid

    def setK(self,k):
        self.k = k

    def getK(self):
        return self.k

    def setNk(self,nk):
        self.nk = nk

    def getNk(self):
        return self.nk

    def setValue(self,newvalue):
        self.val = newvalue
        return self
    def getValue(self):
        return self.val

    def __str__(self):
        return str((self.stype,self.k,self.nk,self.val,self.x,self.y))

    def updatevalue(self,newvalue):
        self.val = newvalue
        return self

    def iscopy(self,other):
        c = self.k == other.k and self.nk == other.nk and self.x == other.x and self.y == other.y
        return c and self.stype == other.stype and self.collid == other.collid

    def __eq__(self, other):
        return self.k == other.k and self.val == other.val

    def todic(self):
        d = dict()
        d.update({"stype":self.stype})
        d.update({"key": self.nk})
        d.update({"nkey":self.val})
        d.update({"val":self.stype})
        return d

    def hashIt(self):
        inp = [str(self.k),str(self.nk),str(self.val),str(self.stype)]
        s = hashlib.sha3_224()
        #ToDO this type of hashing only string values of k,nk,val can generate same hash for different k,nk,val.
        #TODO for example, (1,11,2) and (11,1,2) can give same hash.

        #this method concat hashes of k,nk & v and then takes hash of the string.
        for item in inp:
            hasher =hashlib.sha3_224()
            hasher.update(str(item))
            s.update(hasher.hexdigest()[self.hindexi:self.hindexj])
        return s.hexdigest()[self.hindexi:self.hindexj]

    def toHashNode(self):
        return HashNode(self.collid,self.pcollid,self.x,self.y,self.hashIt())

class HashNode(object):
    #class variables.
    hindexi = 0
    hindexj = 12

    def __init__(self,collid=None,pcollid=None,x=None,y=None,hashv=None):
        self.treeid = 0
        self.x = x #xposition in tree
        self.y = y #yposition in tree.
        self.hashv = hashv
        self.collid = collid
        self.pcollid = pcollid #parent collection id

    def gethashvalue(self):
        return self.hashv

    def hashIt(self):
        inp = [str(self.pcollid),str(self.x),str(self.y)]
        s = hashlib.sha3_224()
        for item in inp:
            hasher =hashlib.sha3_224()
            hasher.update(str(item))
            s.update(hasher.hexdigest()[self.hindexi:self.hindexj])
            s.update(str(s))
        return s.hexdigest()[self.hindexi:self.hindexj]

    def __str__(self):
        s = str((self.collid,self.x,self.y,self.hashv))
        return s
    def todict(self):
        key = self.hashIt()
        value = self.hashv
        return {key: value}

    def updatehv(self, newhashvalue):
        '''updates hash value of this node'''
        self.hashv = newhashvalue

    def hashargs(self, *args):
        inp = args
        s = hashlib.sha3_224()
        for item in inp:
            s.update(str(item))
        return s.hexdigest()[self.hindexi:self.hindexj]

    def __eq__(self, other):
        if self.collid == other.collid and self.x == other.x and self.y == other.y and self.hashv == other.hashv:
            return True
        return False

    def issibling(self, other):
        '''return True is self is sibling hashnode of 'other'''
        if abs(self.x - other.x) == 1 and (self.y == other.y):
            return True
        return False

    def parenthash(self, rhnode):
        lhv = self.hashv
        rhv = rhnode.hashv
        return self.hashargs(lhv, rhv)

    def parentHashNode(self, rhnode):
        lhv = self.hashv
        rhv = rhnode.hashv
        #self.hashargs(lhv,rhv)
        if abs(self.x - rhnode.x) == 1 and (self.y == rhnode.y):
            return HashNode(self.collid, self.pcollid, max(self.x, rhnode.x)//2, self.y-1,self.hashargs(lhv,rhv))
        return None

class HashTree:
    hindexi =0
    hindexj =12

    def __init__(self):
        self.hashes = []
        #self.hashes =dict() #you can make it a dictionary

    def getrootnode(self):
        for hnode in self.hashes:
            if hnode.x == 0 and hnode.y == 0:
                return hnode

    def insert(self,hnode):
        self.hashes.append(hnode)

    def __str__(self):
        s = ''
        for hnode in self.hashes:
            s += str(hnode)+"\n"
        return s

    def isempty(self):
        return self.hashes == []

    def gethashnode(self,collid,x,y):
        '''returns a hashnode at position x,y for a collection with id: collid'''
        for hnode in self.hashes:
            if hnode.x == x and hnode.y == y and hnode.collid == collid:
                return hnode

    def updatehashes(self, collid, x, y, oldhash, newhash):
        # in collection cid, hash x, y is modified from old to new
        # this functions updates all parent hashes accordingly
        oldrootnode = self.gethashnode(collid, 0, 0) # this is the root of the collection that will be
        xc = x;  # current values of x and y as we move up the tree.
        yc = y;
        ch = HashNode(collid=collid,x=xc,y=yc,hashv=newhash)  # current hash. we will keep hash extending with complementary nodes
        while yc > 0:
            if (xc % 2 == 1):  # xc is odd, the sibiling hash is to the left)
                shnode = self.gethashnode(collid, xc - 1, yc);  # sh is siblihg hash
                phv= shnode.parenthash(ch)  # current hash ch, left or right information
            else:
                # sc is even so sibling is to the right
                shnode = self.gethashnode(collid, xc + 1, yc)
                phv = ch.parenthash(shnode) #parent hash value
            # end el
            yc = yc - 1  # climb up
            xc = xc >> 1  # dividing by 2, move left/right
            # @TODO: collect all (sibling hashes) complementary hashes (comp.nodes) for creating certs
            ##also need to add things to collect all complementary hashes to help create certificates
            ch = self.gethashnode(collid,xc,yc)
            ch.updatehv(phv)

        # endwhile
        # when wile loop is done, the root should be set
        return oldrootnode, ch  # ch is the new root.

    def updatehashtree(self,updatedhashnode):
        '''Updates hashes up the hashtree from this updatedhashnode.'''
        oldrootnode = self.gethashnode(updatedhashnode.collid, 0, 0) # this is the root of the collection that will be
        xc = updatedhashnode.x;  # current values of x and y as we move up the tree.
        yc = updatedhashnode.y;
        ch = updatedhashnode  # current hash. we will keep hash extending with complementary nodes
        while yc > 0:
            if (xc % 2 == 1):  # xc is odd, the sibiling hash is to the left)
                shnode = self.gethashnode(updatedhashnode.collid, xc - 1, yc);  # sh is siblihg hash
                phv= shnode.parenthash(ch)  # current hash ch, left or right information
            else:
                # sc is even so sibling is to the right
                shnode = self.gethashnode(updatedhashnode.collid, xc + 1, yc)
                phv = ch.parenthash(shnode) #parent hash value
            # end el
            yc = yc - 1  # climb up
            xc = xc >> 1  # dividing by 2, move left/right
            ch = self.gethashnode(updatedhashnode.collid,xc,yc)
            ch.updatehv(phv)
        # endwhile
        # when wile loop is done, the root should be set
        return oldrootnode, ch  # ch is the new root.

class OMT:
    def __init__(self,omtref='omt', collection=[]):
        self.collection =collection
        self.hashtree = HashTree()
        self.omtstats =dict()
        self.omtstats['depth'] = 0
        self.omtstats['maxkey'] = -sys.maxint
        self.omtstats['minkey'] = sys.maxint
        self.__omtref =omtref

    def getomtref(self):
        return self.__omtref
    def setomtref(self,refname):
        self.__omtref = refname

    def __str__(self):
        s ='--'+self.getomtref()+':collection--\n'
        for snode in self.collection:
            s += str(snode)+"\n"
        #s +='--'+self.getomtref()+':hash tree--\n'
        #s += str(self.hashtree)
        return s

    def extract_min(self):
        return self.omtstats['minkey']
    def extract_max(self):
        return self.omtstats['maxkey']

    def stats(self):
        self.omtstats['len(collection)'] = len(self.collection)
        self.omtstats['len(hashtree)'] = len(self.hashtree.hashes)
        for k,v in self.omtstats.items():
            print k,":",v

    def omtroot(self):
        return self.hashtree.getrootnode()

    def verifywithomtroot(self,secondroot):
        '''Compares secondroot with omt hash root.'''
        if type(secondroot) is str:
            return self.hashtree.getrootnode().hashv == secondroot
        if isinstance(secondroot,HashNode):
            return self.hashtree.getrootnode().hashv == secondroot.hashv

    def rootofcomphashes(self,compnodes,snode):

        '''rootOfComplementaryHashes, Computes a resultant hash from given complementary nodes. If given
        complete complementary hashes, it computes a root hash of an OMT.
        Treat compnodes as a stack to compute root.'''

        nodehash =snode.toHashNode() #snode is a leaf data node, its complementary nodes are in compnodes.
        roothnode = nodehash
        for i in range(0, len(compnodes)):
            if compnodes[i][0] ==1: #right
                rooth=roothnode.parenthash(compnodes[i][1])
                roothnode = HashNode(hashv=rooth)
            else:
                rooth=compnodes[i][1].parenthash(roothnode)
                roothnode = HashNode(hashv=rooth)
        del compnodes
        return roothnode

    def complementaries(self,snode):
        '''return list of commitment nodes for a leaf node with key searchkey'''
        if not isinstance(snode,S):
            return None,None

        if snode:
            collid = snode.collid
            comphnodes =[] #complementary hashes.
            curx = snode.x
            cury = snode.y
            while cury > 0:
                if curx % 2 == 1:
                    #left child is sibling.
                    siblinghashnode = self.gethashnode(collid,curx-1,cury)
                    comphnodes +=[(0, siblinghashnode)] #left 0
                else:
                    #right child is sibling
                    siblinghashnode = self.gethashnode(collid, curx+1, cury)
                    comphnodes += [(1, siblinghashnode)] #right 1
                curx = curx >> 1
                cury =cury - 1
            return comphnodes

    def complementariesxy(self,searchkey):
        '''return list of commitment nodes for a leaf snode: corresponding to a key searchkey.
        It returns three tuple entry like (x,y,siblinghashnode), where x and y are the x position
        along the breath of a binary treee, and y is the level-y of the binary tree. If x%2==1 then this
        siblinghashnode is the right sibling, else left sibling.'''
        if not isinstance(searchkey,S):
            snode,found = self.serchsnode(searchkey)
        else:
            snode,found = self.serchsnode(searchkey.k)

        if snode:
            collid = snode.collid
            comphnodes =[] #complementary hashes.
            curx = snode.x
            cury = snode.y
            while cury > 0:
                if curx % 2 == 1:
                    #left child is sibling.
                    siblinghashnode = self.gethashnode(collid,curx-1,cury)
                    comphnodes +=[(curx-1,cury, siblinghashnode)] #left 0
                else:
                    #right child is sibling
                    siblinghashnode = self.gethashnode(collid, curx+1, cury)
                    comphnodes += [(curx+1,cury, siblinghashnode)] #right 1
                curx = curx >> 1
                cury =cury - 1
            return comphnodes

    def extractminmax(self):
        '''return a snode (k,nk,v) such that k>nk.'''
        for snode in self.collection:
            if snode.k >= snode.nk :
                return snode
        return None

    def serchsnode(self,searchkey): #snode for searchkey
        '''Depending upon searchkey position left,eq, in or out of (k,kn), it returns a flag
        1,2,3,4. (.) position of search key.
        ....(3)...k=(1)----------(2)------------nk...(4)...'''
        for snode in self.collection:
            if snode.k < searchkey < snode.nk:
                return snode,2
            elif snode.k == searchkey:
                return snode,1
            elif snode.k >= snode.nk > searchkey:
                return snode, 3
            elif (searchkey > snode.k >= snode.nk):
                return snode,4

        return None,0

    def serchsnodebykey(self,searchkey,findprev=False): #snode for searchkey
        '''returns a snode whose key is equal to searchkey.'''
        #TODO: Write a common search function that searches in a binary search tree (BST) for efficiencly.
        if findprev:
            for snode in self.collection:
                if snode.nk == searchkey:
                    return snode

        #else find the key
        for snode in self.collection:
            if snode.k == searchkey:
                return snode
        return None

    def enclosedby(self,sk, snk, k):
        # key,nextkey,new-key
        if (sk < k < snk) or (k > sk >= snk) or (sk >= snk > k):
            return True
        return False

    def getbalancingnode(self,splitsnode):
        '''return the snode whose vertical level difference with splitsnode is  >= 1.'''
        if not splitsnode:
            return None
        for balancingnode in self.collection:
            if (splitsnode.y - balancingnode.y) >= 1:
                return balancingnode
        return None

    def getsplitsnode(self,searchkey,searchkind=STYPE['LUT']):
        #search for keys to insert searchkey in a collection.
        #ToDO : This search can be made in Balanced Binary Search Tree. It reduces search to O(LogN).
        if searchkind == STYPE['LUT']:
            for snode in self.collection:
                if self.enclosedby(snode.k, snode.nk, searchkey):
                        return snode
        else:
            for snode in self.collection:
                if snode.k == searchkey:
                    return snode
        return None

    def isempty(self):
        return self.collection == [] and self.hashtree.isempty()

    def gethashnode(self, collid, x, y):
        return self.hashtree.gethashnode(collid,x,y)

    def init(self):
        ihnode = HashNode()
        ihnode.x = 0
        ihnode.y = 0
        ihnode.pcollid = 0
        ihnode.hashv = 0
        self.hashtree.insert(ihnode)

    def insertfirstkey(self,k,val,collid,pcollid=0,stype=STYPE['LUT']):
        '''inserts the first S nodes and HashNode in an OMT'''
        '''insert an item of type S into self.collection.'''
        #if OMT empty, then insert init keys,
        #else, get a split key and split the sk to insert k,val.
        if self.isempty():
            snode = S()# data node
            snode.k = k
            snode.nk = k
            snode.collid = collid
            snode.pcollid=pcollid
            snode.stype = stype
            if stype == STYPE['DIC']:
                snode.val = 0 #do nodeupdate to update this snode's value
            elif stype == STYPE['LUT']:
                snode.val = val #you can also do nodeupdate to update this snode's value later on.

            snode.x = 0
            snode.y = 0
            hnode = snode.toHashNode() #hashnode for this snode.
            self.collection.append(snode)
            self.hashtree.insert(hnode)

            self.omtstats['minkey'] = min(self.omtstats['minkey'],snode.k)
            self.omtstats['maxkey'] = max(self.omtstats['maxkey'], snode.k)
            self.omtstats['depth'] = max(snode.y, self.omtstats['depth'])
            return snode

    def nodeupdate(self,updatesnode,updatevalue,forcetonull=False):
        '''Updates a S-node by a new update value.'''

        if updatesnode:
            if updatevalue == 0 or updatevalue == None:
                if forcetonull:
                    updatesnode.val = 0
                else:
                    return False
            else:
                updatesnode.val = updatevalue

            oldhash = self.gethashnode(updatesnode.collid, updatesnode.x, updatesnode.y)
            oldhashval = oldhash.hashval
            newhashv = updatesnode.hashIt()
            self.hashtree.updatehashes(updatesnode.collid, updatesnode.x, updatesnode.y, oldhash, newhashv)
            return True

        # ToDO: Recursive update of hashes.
        return False

    from copy import copy
    def insert(self,splitsnode, knew,vnew):
        '''inserts S nodes and HashNode in an OMT'''

        newsnode = copy.deepcopy(splitsnode)
        newsnode.k = knew  # new S
        collid= splitsnode.collid
        if newsnode.k is None:
            return False

        if splitsnode.stype == STYPE['DIC']:
            newsnode.val = 0 #do nodeupdate to update this snode's value
        elif splitsnode.stype == STYPE['LUT'] :
            # S.kn next to the split key
            if not ((splitsnode.k < newsnode.k < splitsnode.nk) or (newsnode.k > splitsnode.k >= splitsnode.nk) or (
                    splitsnode.k >= splitsnode.nk > newsnode.k)):
                return None# new key k is not enclosed by sk, sk ......(k < skp <= sk) or (sk <= skp < k):
            newsnode.val = vnew #you can also do nodeupdate to update this snode's value.

        old_nk = splitsnode.nk;# #back-up
        splitsnode.nk = newsnode.k;
        newsnode.nk = old_nk;
        ox = splitsnode.x;
        oy = splitsnode.y; # old x, y values of item sk, it's going to be pulled down.

        #xpulling down both items
        splitsnode.x = 2 * ox;
        newsnode.x = 2 * ox + 1;
        splitsnode.y = newsnode.y = oy + 1

        oldhashnode = self.gethashnode(collid, ox, oy) #old hash corresponding to split snode.
        oldhashv = oldhashnode.hashv
        shashnode = splitsnode.toHashNode()  # hash of split-key
        newhashnode = newsnode.toHashNode()  # self.h(S SN.key, SN.nkey, SN.val,N.type_)  # hash of new item

        # compute parent hash value # do not need to delete splitsnode because it's values are modified already.
        newphv = shashnode.parenthash(newhashnode)  # newparent hash value
        oldhashnode.updatehv(newphv)

        # write back newsnode and  store new hashes
        self.collection.append(newsnode)
        self.hashtree.insert(shashnode)
        self.hashtree.insert(newhashnode)

        [oldroot, newroot] = self.hashtree.updatehashes(collid, ox, oy, oldhashv, newphv) #oldhashv was replaced by newphv
        #print oldroot,newroot
        #[p_cid, p_key] = self.get(collectid)
        #if p_cid != 0:
        #    self.Recursive(p_cid, p_key, oldroot, newroot)
        #return S, SN
        self.omtstats['depth'] = max(newsnode.y, self.omtstats['depth'])
        self.omtstats['minkey'] = min(self.omtstats['minkey'], newsnode.k)
        self.omtstats['maxkey'] = max(self.omtstats['maxkey'], newsnode.k)
        return True

    def insertkeybalanced(self,splitsnode,balancingnode,knew,vnew):
        import copy
        collid = splitsnode.collid
        if balancingnode is None:
            if splitsnode:
                self.insert(splitsnode, knew, vnew)
                return True
            return False
        #splitsnodeis the key that is split; knew is the inserted key; balancingnode is the key pulled down to become sibling of k
        newsnode = copy.deepcopy(splitsnode)
        newsnode.k = knew  # new S
        # First check k enclosed by sk and sk
        if splitsnode.stype == STYPE['DIC']:
            newsnode.val = 0 #do nodeupdate to update this snode's value
        elif splitsnode.stype == STYPE['LUT'] :
            # S.kn next to the split key
            if not ((splitsnode.k < newsnode.k < splitsnode.nk) or (newsnode.k > splitsnode.k >= splitsnode.nk) or (
                    splitsnode.k >= splitsnode.nk > newsnode.k)):
                return False  # new key k is not enclosed by sk, sk ......(k < skp <= sk) or (sk <= skp < k):
            newsnode.val = vnew #you can also do nodeupdate to update this snode's value.

        old_nk = splitsnode.nk;# #back-up
        splitsnode.nk = newsnode.k;
        newsnode.nk = old_nk;

        #if (newsnode.stype == STYPE['DIC']):
        #    newsnode.val = 0
        oldhashnode = self.gethashnode(collid,splitsnode.x,splitsnode.y)
        oldhashv, newhashv = oldhashnode.gethashvalue(), splitsnode.hashIt()
        oldhashnode.updatehv(newhashv)

        [oldroot, newroot] = self.hashtree.updatehashes(collid, splitsnode.x, splitsnode.y, oldhashv, newhashv)

        # fix moving newnode to make sibling of balancingnode.
        ox = balancingnode.x;
        oy = balancingnode.y;
        oldbhashnode = self.gethashnode(collid, ox, oy)
        oldblhashv = oldbhashnode.gethashvalue()

        balancingnode.x = balancingnode.x << 1;  # double x
        newsnode.y = balancingnode.y = balancingnode.y + 1;
        newsnode.x = balancingnode.x + 1;  # make newsnode and balancingnode siblings

        newhashnode = newsnode.toHashNode()
        self.collection.append(newsnode) #append newsnode.
        self.hashtree.insert(newhashnode) #insert new hashnode
        self.hashtree.insert(balancingnode.toHashNode()) #insert balancing node's hashnode.

        newphv = oldbhashnode.parenthash(newhashnode)  # newparent hash value
        oldbhashnode.updatehv(newphv)

        [newroot, newnewroot] = self.hashtree.updatehashes(collid, ox, oy, oldblhashv, newphv)
        self.omtstats['minkey'] = min(self.omtstats['minkey'], newsnode.k)
        self.omtstats['maxkey'] = max(self.omtstats['maxkey'], newsnode.k)
        #[p_cid, p_key] = self.get(cid);
        #if p_cid != 0:
        #    self.Recursive(p_cid, p_key, old, newnew);
        return True

    #method-1
    def deletesnode(self,deletesnode):
        '''#method-1:deletes an item of type S in self.collection'''

        prevtodelnode = self.serchsnodebykey(deletesnode.k,findprev=True)
        prevtodelnode.nk = deletesnode.nk
        #get hashnode for this prev snode
        poldhashnode = self.gethashnode(prevtodelnode.collid,prevtodelnode.x, prevtodelnode.y)
        #update hashvalue for this oldhashnode
        poldhashnode.updatehv(poldhashnode.hashIt())
        #update up the tree from this hashnode
        self.hashtree.updatehashtree(poldhashnode)
        #delete the leaf snode.
        self.collection.remove(deletesnode)
    #method-2
    def nullifysnode(self,deletesnode):
        '''#method-1:set's deletesnode's value to 0/NULL in self.collection'''
        self.nodeupdate(deletesnode,0,forcetonull=True)

    '''Inserting a collection of type LUT.'''
    def bulkinsert(self,collection_dict,balancetree=0):
        '''given a list of objects of type S, it inserts into self.collection.'''
        cnt=0
        for k,v in collection_dict.items()[0:]:
            if cnt==0:
                print("---first-key--"),k
                self.insertfirstkey(k,v,collid=1,stype=STYPE['LUT']) #inserted node returned
                print self
                cnt+=1
            else:
                print("---inserting--"),k
                #get enclosing snode for k,v
                splitsnode = self.getsplitsnode(k)
                #splitsnode is not found. because duplicate key already found.
                if not splitsnode:
                    continue
                if balancetree:
                    balancingnode = self.getbalancingnode(splitsnode)
                    self.insertkeybalanced(splitsnode, balancingnode, k, v)
                else:
                    self.insert(splitsnode,k,v)
        pass

    '''Inserting a collection of type DIC'''
    def bulkinsertsegs(self,namedsegments,balancetree=0):
        '''given a list of objects of type S, it inserts into self.collection.
        namedsegments is a tuple of form(k,value) for a line segment.'''
        cnt=0
        for i in range(len(namedsegments)): # in namedsegments:
            k,v = namedsegments[i]
            if cnt==0:
                print("---Inserting first-key--"),k
                result = self.insertfirstkey(k,v,collid=1,stype=STYPE['DIC']) #inserted node returned
                if not result:
                    raise ("First key insertion failed! key:"),k
                cnt+=1
            else:
                print("---Inserting key--"), k
                prevkey,prevval=namedsegments[i-1] #do not use previous key if this item is LUT.
                splitsnode = self.getsplitsnode(prevkey,searchkind=STYPE['DIC'])
                #no splitnode found for the prevkey beacuase a node with that key might already exist. In this case execute nodeupdate func.
                if not splitsnode:
                    continue
                if balancetree:
                    balancingnode = self.getbalancingnode(splitsnode)
                    result = self.insertkeybalanced(splitsnode, balancingnode, k, v)
                    if not result:
                        raise ("Could not insert key:"),k
                else:
                    self.insert(splitsnode,k,v)

    def update_snode(self,oldsnode,newvalue):
        '''Updates oldsnode's value by newvalue and updates hashes up the tree.'''
        if newvalue:
            newsnode = oldsnode.updatevalue(newvalue)
            oldhashnode = self.gethashnode(oldsnode.collid, oldsnode.x, oldsnode.y)  # old hash corresponding to split snode.
            olhhashval,newhashvalue =oldhashnode.hashv, oldsnode.hashIt()
            oldhashnode.updatehv(newhashvalue)
            self.hashtree.updatehashtree(oldhashnode) #
            #Recursive Updates.

    def searchsnodebyvalue(self,searchvalue,findprev=False):
        '''returns a snode whose value is equal to searchvalue.'''
        #TODO: Write a common search function that searches in a binary search tree (BST) for efficiencly.
        if findprev:
            for snode in self.collection:
                if snode.val[:-1] == searchvalue[:-1]:
                    return snode

        #else find the key
        for snode in self.collection:
            if snode.val[:-1] == searchvalue[:-1]:
                return snode
        return None




class Transaction(object):
    def __init__(self,trantype,udi,cdi=None):
        self.udi = udi
        self.cdi = cdi
        self.trantype = trantype

    def getUdi(self):
        return self.udi

class TCB:
    def __init__(self):

        self.rootdic={'segomt':0,'eventomt':0, 'activeomt':0}
        self.mindic={'segomt':sys.maxint,'eventomt':sys.maxint}
        self.maxdic={'segomt':-sys.maxint,'eventomt':-sys.maxint}

    def ICurrentState(self):
        print("----------TCB-STATE-----------")
        print("roots:"),self.rootdic
        print("min dic:"),self.mindic
        print("max dic:"),self.maxdic


    def processTransac(self,transaction):
        ##
        if transaction.trantype == 'initprocessevent': #initialize to process an event.
            return True
            udiProcessEvent = transaction.getUdi()
            initeventsnode =udiProcessEvent['event']
            complementaries =udiProcessEvent['eventvos'] #verification objects for this event
            activeomtref = udiProcessEvent['activeomtref'] # = activeomt.getomtref()
            eventomtref = udiProcessEvent['eventomtref']

            #Pre-Condition
            #-------------
            #initevent exists in event tree
            comproot = OMT().rootofcomphashes(complementaries, initeventsnode)
            #initevent is the minimum key in the event tree
            comproot.gethashvalue() == self.rootdic[eventomtref]
            self.rootdic[activeomtref] == 0
            self.mindic[eventomtref] == initeventsnode.k #must be a minimum key
            eventtype = initeventsnode.val[-1]  # start or end
            #Post-Condition
            # -------------
            #set minimum key to currentevent.nextkey
            self.mindic[eventomtref] = initeventsnode.nk  # next event expected.
            #insert a node in activesegomt
            if initeventsnode.k == self.maxdic[eventomtref]:
                #return result or terminate the process
                return True
        ##
        if transaction.trantype == 'furhterprocessevent': #initialize to process an event.
            udiProcessEvent = transaction.getUdi()
            currenteventsnode =udiProcessEvent['event']
            eventomtref = udiProcessEvent['eventomtref']
            #ENDING CASE
            if currenteventsnode.k == self.maxdic[eventomtref]:
                #return result or terminate the process
                return True
            #OR
            if currenteventsnode.nk < currenteventsnode.k:
                #return result or terminate the process.
                return True
            #

        ##
        if transaction.trantype == 'nodeupdate':
            #verify the transaction and return a certificate for the validity.
            #calculate the new root and set new root for this tree
            #UDI contains following items.

            udi = transaction.getUdi()
            proofsnode = udi['updatesnode'] #= updatesnode
            searchkey = udi['updatekey'] #= searchkey
            newval = udi['updatevalue'] # = newval
            omtref = udi['omtref']
            complementaries = udi['vos'] #= segsomt.complementaries(updatesnode)

            #Pre-Conditions
            comproot = OMT().rootofcomphashes(complementaries, proofsnode)
            if comproot.gethashvalue() == self.rootdic[omtref] and proofsnode.k == searchkey:#exists
                updatedsnode = proofsnode.setValue(newval) #update proofsnode
                newcomproot = OMT().rootofcomphashes(complementaries, updatedsnode) #find new root after updating value.
                #Post-Conditions
                self.rootdic[omtref] = newcomproot.gethashvalue() #updated omtroot
                return True
        ##
        if transaction.trantype == 'eventnodeupdate':
            '''1) takes each snode from segomt, 2) Creates two event nodes for this segment'''
            udiCreateEvent = transaction.getUdi()
            updatesnode = udiCreateEvent['updatesnode']
            searchkey = udiCreateEvent['updatekey']
            newval = udiCreateEvent['updatevalue']
            eventtype = udiCreateEvent['eventtype']
            voseventomt = udiCreateEvent['voseventomt']
            eventomtref = udiCreateEvent['eventomtref']

            segsnode = udiCreateEvent['segmentsnode']
            segomtref = udiCreateEvent['segomtref']
            vossegomt = udiCreateEvent['vossegomt']

            #Pre-Conditions
            eventomtroot = OMT().rootofcomphashes(voseventomt,updatesnode)
            c1 = eventomtroot.gethashvalue() == self.rootdic[eventomtref]

            segomtroot = OMT().rootofcomphashes(vossegomt, segsnode)
            c2 = segomtroot.gethashvalue() == self.rootdic[segomtref]

            x1, y1, x2, y2 = segsnode.val[0], segsnode.val[1], segsnode.val[2], segsnode.val[3]
            c3 = (x1 < x2 and eventtype == EVENT_TYPE['START']) or (x2 > x1 and eventtype == EVENT_TYPE['STOP'])

            c4 = x1 == updatesnode.x or x2 == updatesnode.x #value x1 or x2 is a key in updatesnode in event-omt.

            if c1 and c2 and c3:
                updatedsnode = updatesnode.setValue(newval)  # update proofsnode
                newcomproot = OMT().rootofcomphashes(voseventomt, updatedsnode)
                #Post-Conditions
                self.rootdic[eventomtref] = newcomproot.gethashvalue() #updated omtroot
                self.mindic[eventomtref] = min(self.mindic[eventomtref],x1,x2)
                self.maxdic[eventomtref] = max(self.maxdic[eventomtref], x1, x2)
                return True

        if transaction.trantype == 'omtinit':
            if transaction.getUdi()['omtref'] in ['segomt','eventomt']:
                omtref =transaction.getUdi()['omtref']
                if self.rootdic[omtref] == 0:
                    self.rootdic[omtref]= transaction.getUdi()['omtroot'].gethashvalue()

            return True

    def ICalcintersection(self,seg1, seg2):
        s1x1,s1y1,s1x2,s1y2 = seg1
        s2x1, s2y1, s2x2, s2y2 = seg2
        seg2 = ((s2x1, s2y1), (s2x2, s2y2))
        seg1 = ((s1x1, s1y1), (s1x2, s1y2))
        try:
            p = seg1[0]
            r = (seg1[1][0] - seg1[0][0], seg1[1][1] - seg1[0][1])
            q = seg2[0]
            s = (seg2[1][0] - seg2[0][0], seg2[1][1] - seg2[0][1])
            denom = r[0] * s[1] - r[1] * s[0]
        except:
            return None
        if denom == 0:
            return None
        numer = float(q[0] - p[0]) * s[1] - (q[1] - p[1]) * s[0]
        t = numer / denom
        numer = float(q[0] - p[0]) * r[1] - (q[1] - p[1]) * r[0]
        u = numer / denom
        if (t < 0 or t > 1) or (u < 0 or u > 1):
            return None
        x = p[0] + t * r[0]
        y = p[1] + t * r[1]
        return (round(x,2), round(y,2))

    def IverifyKeyEnclosedby(self,k,nk,searchkey):
        return True

    def IpointOnline(self):
        return True

    def IverifySignature(self,message,signature,secret):
        '''return true if a signature belongs to an authenticated party.'''
        return True
    def IcomputeMerkleRoot(self,comphashes,collection_item):
        '''computes a possible merkle root given complementary hashes: comphashes, and collection-item'''
        proof_root = None
        return proof_root
    def IverifyMerkleRoots(self,referenced_omttree,proof_root):
        return self.rootdic[referenced_omttree] == proof_root

class Certificate(object):

    certtypes={'nodeverify':0, 'equivalentcert':1}
    def __init__(self,certtype,*args):
        self.certtype = certtype

    def __str__(self):
        s = self.certtype
        return s

class BlockChain(object):
    def __init__(self):
        self.blocklist =[]


#Driver classe for Omt.
class UDI:

    def __init__(self):
        self.operationcode={}

    def __setitem__(self, key, item):
        self.operationcode[key] = item

    def __getitem__(self, key):
        return self.operationcode[key]

    def addoperands(self,key,value):
        self.operationcode[key] = value

    def __str__(self):
        s = ''
        for k,v in self.operationcode.items():
            s += k +":"
            if isinstance(v,S):
                s += str(v)
            elif isinstance(v,list):
                if v:
                    for item in v:
                        if isinstance(item,tuple):
                            s += str(item[0])+','+str(item[1])+'\n'
            else:
                s += str(v)
            s +='\n'
        return s

    def toJson(self,outfilename):
        import json
        with open(outfilename, 'w') as outfile:
            json.dump(self.operationcode, outfile)
        #json.dumps(self.operationcode)
