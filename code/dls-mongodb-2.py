'''
Author: Naresh Adhikari and Dr. Ramkumar Mahalingam. (advisor)
Net-ID: Na542
Date: Oct 19,  2018

422101088
'''
import sys
import hashlib
if sys.version_info < (3, 6):
    import sha3
import datastructure2 as DS
reload(DS)

GD = {}  # giant dictionary
import sys
class DLS:
    def __init__(self):
        self.debug = True
        global GD
        GD = {}
        self.maincolid = self.ID(1, 0, 0)  # id of main collection with parent 0 and key 0
        self.depths = {}
        self.nhashes =0
        self.nleafs = 0
        self.xminkey = sys.maxint
        self.xmaxkey = -sys.maxint
        self.yminkey = sys.maxint
        self.ymaxkey = -sys.maxint
        pass
    def dlsstats(self):
        s =''
        s +='main-collection-id:' +self.maincolid+'\n'
        s +='hash-tree-root:' +self.GetMainCollectionCommitment()+'\n'
        s += 'hash-tree-depth:' +str(self.depths[self.maincolid])+'\n'
        s += 'number-of-hashes:' + str(self.nhashes)+'\n'
        s += 'number-of-leaves:' + str(self.nleafs)+'\n'
        s += 'minimum-key,maximum-key:' +str(self.xminkey)+','+str(self.xmaxkey)+'\n'
        return s

    def put(self, key, record):
        # a record can be [pcid,pkey] or structure S or a hash, pcid,pkey = parent collection id, parent key
        if type(record) is dict:
            GD[key] = record  # if mongodb is used, get a record with column "key"=key
        elif isinstance(record, DS.St):
            GD[key] = record.castToDict()
        else:
            GD[key] = record
        # NOTE:
        ''' for a python dictionary, does not matter if record exists for (key,record).
        for mongodb, if (key,record) exists: do-update, else do-insert
        '''
    def delete(self,key):
        try:
            del GD[key]
        except Exception, e:
            print(str(e))

    def get(self, key):
        # will get the record of 3 differnt types.
        return GD[key]

    def ID(self, *args):
        s = hashlib.sha3_224()
        for p in args:
            hasher =hashlib.sha3_224()
            hasher.update(str(p))
            s.update(hasher.hexdigest()[0:12])
        return s.hexdigest()[0:12]

    def h(self, *args):
        s = hashlib.sha3_224()
        for p in args:
            #hasher =hashlib.sha3_224()
            #hasher.update(str(p))
            s.update(str(p))#s.update(hasher.hexdigest()[0:12])
        return s.hexdigest()[0:12]

    def Init(self):
        cid = self.ID(1, 0, 0)  # id of main collection with parent 0 and key 0
        self.put(cid, [0, 0])  # even main collection has an entry
        self.put(self.ID(cid, 0, 0), 0)  # initialize root of main collection to zero
        self.depths[self.maincolid] = 0
        self.nhashes =0
        return cid

    def InsertFirstKey(self, cid, xk, yk,type=1):  # cid,key,type; default is LUT; type=0 is dict.
        oldhash = self.get(self.ID(cid, 0, 0))
        if oldhash == 0:
            #newhash = self.h(xk, xk, yk,yk,0, t)
            S = DS.St(xk,xk,yk,yk,val={},type_=0,posx=0,posy=0,colid=cid)
            newhash = S.h()
            print("put"),self.ID(cid, S.xkey,S.xnkey,S.ykey)
            print cid, S.xkey,S.xnkey,S.ykey
            self.put(self.ID(cid, S.xkey,S.xnkey,S.ykey), S) #id is cid+h(xk)+h(xnk)+h(yk)
            self.put(self.ID(cid, 0, 0), newhash)
            [p_cid, p_key] = self.get(cid)  # see if cid has parents
            if p_cid != 0:  # has parent
                self.Recursive(p_cid, p_key, oldhash, newhash)
            #stat
            self.nhashes +=1
            self.nleafs +=1
            self.xminkey = min(self.xminkey,xk)
            self.xmaxkey = max(self.xmaxkey,xk)
            self.yminkey = min(self.yminkey, yk)
            self.ymaxkey = max(self.ymaxkey, yk)
            return S

    def CreateCollection(self, cid, k):  # parent cid, parent key
        # a new colleciton is created under any parent collection 'cid', leaf structure S with key 'k'
        Sdic = self.get(self.ID(cid, k))
        S = DS.St(dic=Sdic)  # coz get returns dic type,
        if S.val == 0:  # value os zero so we can initialize collection by adding the record
            childcid = self.ID(1, cid, k)
            self.put(childcid, [cid, k])  # new collection created
            self.put(self.ID(childcid, 0, 0), 0)  # initializing root of the childcid hashes.
            self.depths[childcid] = 0
            return childcid

    def InsertKey(self, collectid=None, structS=None,xnewkey=None,ynewkey=None,dimension='y'):  # sk is the key that is split to insert key k
        import copy
        #print(".."),collectid, structS.xkey, structS.xnkey,structS.ykey,self.ID(collectid, structS.xkey, structS.xnkey,structS.ykey)
        Sdic = self.get(self.ID(collectid, structS.xkey, structS.xnkey,structS.ykey))  # cast to S obj
        S = DS.St(dic=Sdic)
        SN = copy.deepcopy(S);  # new item
        if SN.type_ == 0:
            SN.val = {}
        ##
        if dimension == 'y':
            if ynewkey is None:
                return

            if not ((S.ykey < ynewkey < S.ynkey) or (ynewkey > S.ykey >= S.ynkey) or (S.ykey >= S.ynkey > ynewkey)):
                return  None,None# new key k is not enclosed by sk, sk ......(k < skp <= sk) or (sk <= skp < k):
            oldyk = S.ykey #back-up
            S.ynkey = ynewkey;
            SN.ykey = ynewkey

            ox = S.x;
            oy = S.y  # old x, y values of item sk, it's going to be pulled down.
            oldhash_sk = self.get(self.ID(collectid, ox, oy))
            newhash_sk = S.h() #hash of split-key
            hash_k = SN.h() #self.h(SN.key, SN.nkey, SN.val, SN.type_)  # hash of new item
            # pulling down both items
            S.x = 2 * ox;
            SN.x = 2 * ox + 1;
            S.y = SN.y = oy + 1
            #write back items

            self.put(self.ID(collectid, S.xkey,S.xnkey,S.ykey), S)  # cast S and SN objects,rep
            self.put(self.ID(collectid, SN.xkey,SN.xnkey,SN.ykey), SN)  # new-1

            #delete previous S
            #self.delete(self.ID(collectid,S.xkey,S.xnkey,oldyk))
            # store new hashes
            self.put(self.ID(collectid, S.x, S.y), newhash_sk)
            self.put(self.ID(collectid, SN.x, SN.y), hash_k)

            # compute parent
            newparent = self.h(newhash_sk, hash_k)
            self.put(self.ID(collectid, ox, oy), newparent)  # oldhash_sk was the original value at ox, oy
            [oldroot, newroot] = self.UpdateHashes(collectid, ox, oy, oldhash_sk, newparent)
            [p_cid, p_key] = self.get(collectid)
            if p_cid != 0:
                self.Recursive(p_cid, p_key, oldroot, newroot)
        else:
            if xnewkey is None:
                return
            #S.xnkey next to the split key
            if not ((S.xkey < xnewkey < S.xnkey) or (xnewkey > S.xkey >= S.xnkey) or (
                    S.xkey >= S.xnkey > xnewkey)):
                return None, None  # new key k is not enclosed by sk, sk ......(k < skp <= sk) or (sk <= skp < k):
            old_xnk = S.xnkey #back-up
            S.xnkey = xnewkey;
            SN.xkey = xnewkey
            ox = S.x;
            oy = S.y  # old x, y values of item sk, it's going to be pulled down.
            oldhash_sk = self.get(self.ID(collectid, ox, oy))
            newhash_sk = S.h()  # hash of split-key
            hash_k = SN.h()  # self.h(S SN.key, SN.nkey, SN.val,N.type_)  # hash of new item
            #xpulling down both items
            S.x = 2 * ox;
            SN.x = 2 * ox + 1;
            S.y = SN.y = oy + 1

            # write back items
            self.put(self.ID(collectid, S.xkey, S.xnkey ,S.ykey), S)  # cast S and SN objects,rep
            self.put(self.ID(collectid, SN.xkey, SN.xnkey,SN.ykey), SN)  # new-1
            self.delete(self.ID(collectid,S.xkey,old_xnk,S.ykey))
            # store new hashes
            self.put(self.ID(collectid, S.x, S.y), newhash_sk)
            self.put(self.ID(collectid, SN.x, SN.y), hash_k)

            # compute parent
            newparent = self.h(newhash_sk, hash_k)
            self.put(self.ID(collectid, ox, oy), newparent)  # oldhash_sk was the original value at ox, oy
            [oldroot, newroot] = self.UpdateHashes(collectid, ox, oy, oldhash_sk, newparent)
            [p_cid, p_key] = self.get(collectid)
            if p_cid != 0:
                self.Recursive(p_cid, p_key, oldroot, newroot)

            self.xminkey = min(self.xminkey,xnewkey)
            self.xmaxkey = max(self.xmaxkey,xnewkey)
        #end-else

        # stats
        self.depths[self.maincolid] = max(S.y, self.depths[self.maincolid])
        self.nhashes += 2
        self.nleafs += 1
        return S, SN

    def UpdateHashes(self, cid, x, y, old, new):
        #if self.debug: print("updating hashes."), cid, x, y, old, new

        # in collection cid, hash x, y is modified from old to new
        # this functions updates all parent hashes accordingly
        oldroot = self.get(self.ID(cid, 0, 0))  # this is the root of the collection that will be
        # changed by this function.#along with all hashes on the way to the root

        #done in insertKey fun.
        #self.put(self.ID(cid, x, y), new)  # first we will set the new leaf hash;this is duplicate fun.

        xc = x;  # current values of x and y as we move up the tree.
        yc = y;
        ch = new  # current hash. we will keep hash extending with complementary nodes
        while yc > 0:
            if (xc % 2 == 1):  # xc is odd, the sibiling hash is to the left)
                sh = self.get(self.ID(cid, xc - 1, yc));  # sh is siblihg hash
                ch = self.h(sh, ch)  # current hash ch, left or right information
            else:
                # sc is even so sibling is to the right
                sh = self.get(self.ID(cid, xc + 1, yc))
                ch = self.h(ch, sh)
            # end el
            yc = yc - 1  # climb up
            xc = xc >> 1  # dividing by 2, move left/right
            # @TODO: collect all (sibling hashes) complementary hashes (comp.nodes) for creating certs
            ##also need to add things to collect all complementary hashes to help create certificates

            self.put(self.ID(cid, xc, yc), ch)  # update, this is the new current hash   a3w
        # endwhile
        # when wile loop is done, the root should be set
        return oldroot, ch  # ch is the new root.

    def Recursive(self, cid, key, old, new):  # old hash, new hash
        print("recursive")
        while cid != 0:
            if self.debug: print("recursive"), self.ID(cid, key), new
            Sdic = self.get(self.ID(cid, key))
            S = DS.St(dic=Sdic)  # bcoz get returns dic type,
            oldhash = self.get(self.ID(cid, S.x, S.y))
            S.val = new;
            newhash = self.ID(S.key, S.nkey, S.val, S.type_)
            self.put(self.ID(cid, S.x, S.y), newhash)
            self.put(self.ID(cid, key), S)
            self.UpdateHashes(cid, S.x, S.y, oldhash, newhash)
            # [cid,key] = self.get(self.ID(1,cid,key))
            [cid, key] = self.get(cid)
    def InsertKeyBalanced(cid, sk, k, kb):
        import copy
        # ks is the key that is split
        # k is the inserted key
        # kb is the key pulled down to become sibling of k
        Ssk = self.get(self.ID(cid, sk)) #Structure for split key sk.
        Sk = copy.deepcopy(Ssk)
        Skb = self.get(h(cid, kb))

        # First check k enclosed by sk and sk
        skp = Ssk.nkey
        if not (sk < k < skp) or (k < skp <= sk) or (sk <= skp < k):
            return

        Sk.key = k;
        Ssk.nkey = k;
        if (Sk.type_ == 0):
            Sk.val = 0;

        oldhash_ks = self.get(h(cid, Ssk.x, Ssk.y));
        newhash_ks = self.ID(Ssk.key, Ssk.nkey, Ssk.val, Ssk.type_);
        self.put(h(cid, Ssk.x, Ssk.y), newhash_ks)

        [oldroot, newroot] = self.UpdateHashes(cid, Ssk.x, Ssk.y, oldhash_ks, newhash_ks)
        #fix hash
        ox = Skb.x;
        oy = Skb.y;
        oldhash_kb = self.get(h(cid, ox, oy))

        Skb.x = Skb.x << 1; #double x
        Sk.y = Skb.y = Skb.y + 1;
        Sk.x = Skb.x + 1;  # make Skb and Sb siblings
        self.put(self.ID(cid,Sk.k),Sk)

        #put Sk
        hash_k = self.ID(Sk, key.Sk.nkey, Sk.val, Sk.type);
        self.put(self.ID(cid,Sk.x,Sk.y), hash_k)

        newparent = self.ID(oldhash_kb, hash_k)  # this will replace oldhash_
        [newroot, newnewroot] = self.UpdateHashes(cid, ox, oy, oldhash_kb, newparent)
        [p_cid, p_key] = self.get(cid);
        if p_cid != 0:
            self.Recursive(p_cid, p_key, old, newnew);
    def BalanceTree(cid, k1, k2, k3):

        # k2 and k3 are siblings. we want to move k2 up one level and make k3 a sibling of k1
        # (by pulling k1 down one level)
        S1 = self.get(self.ID(cid, k1));
        S2 = self.get(self.ID(cid, k2));
        S3 = get(self.ID(cid, k3));
        if (S2.y != S3.y) or (S3.x - S2.x != 1):
            return;
        oy = S2.y - 1;
        ox = S2.x << 1  # ox,oy will be new position of S2
        hashS2S3 = self.get(self.ID(cid, ox, oy));  # get parent hash of S2S3 will need to change it to hash of S2
        hashS2 = self.get(self.ID(cid, S2.x, S2.y));  # value to be assigned to ox,oy

        self.put(self.ID(cid, ox, oy), hashS2);
        self.put(self.ID(cid, k2), S2);  # and done

        self.delete(self.ID(cid, ox >> 1, oy + 1));
        self.delete(self.ID(cid, ox >> 1 + 1, oy + 1));  # old S2 and S3
        [old, new] = self.UpdateHashes(cid, ox, oy, hashS2S3, hashS2);

        hashS3 = self.get(self.ID(cid, S3.x, S3.y));
        hashS1 = self.get(self.ID(cid, S1.x, S1.y));
        hashS1S3 = self.h(hashS1, hashS3)

        ox = S1.x;
        oy = S1.y  # this is the position that will change from hashS1 -> hashS1S3

        S1.x = ox >> 1;
        S3.x = ox >> 1 + 1;
        S1.y = S3.y = oy + 1  # //change x and y in records

        self.put(h(cid, k1), S1);
        self.put(h(cid, k3), S3);  # //write back modified records

        # now write back new hashes
        self.put(h(cid, S1.x, S1.y), hashS1);
        self.put(h(cid, S3.x, S3.y), hashS3);
        self.put(h(cid, oox, ooy), hashS1S3);

        [new, newnew] = self.UpdateHashes(cid, ox, oy, hashS1, hashS1S3);
        [p_cid, p_key] = self.get(cid);
        if p_cid != 0:
            self.Recursive(p_cid, p_key, old, newnew)
    def SwapLeaves(cid, k1, k2):  # Swap k1 and k2

        S1 = self.get(self.ID(cid, k1));
        S2 = self.get(self.ID(cid, k2));
        h1 = self.get(self.ID(cid, S1.x, S1.y));
        h2 = self.get(self.ID(cid, S2.x, S2.y));
        y1 = S1.y;
        x1 = S1.x;

        S1.x = S2.x;
        S1.y = S2.y;
        S2.x = x1;
        S2.y = y1;  # swap positions
        self.put(self.ID(cid, k1), S1);
        self.put(self.ID(cid, k2), S2);
        self.put(self.ID(cid, S1.x, S1.y), h1);

        [old, new] = self.UpdateHashes(cid, S1.x, S1.y, h2, h1);
        self.put(self.ID(cid, S2.x, S2.y), h2);
        [new, newnew] = self.UpdateHashes(cid, S2.x, S2.y, h1, h2);
        [p_cid, p_key] = self.get(cid);
        if p_cid != 0:
            self.Recursive(p_cid, p_key, old, newnew);
    def GetMainCollectionId(self):
        return self.ID(1, 0, 0)
    def GetMainCollectionCommitmentKey(self):
        return self.ID(self.GetMainCollectionId(), 0, 0)
    def GetMainCollectionCommitment(self):
        return GD[self.GetMainCollectionCommitmentKey()]
    def Get2LevelCollectionId(self, pkey=None):
        # handle if pkey is none: as:
        # if self.debug: print("parent key value for parameter 'pkey' is required.")
        return self.ID(1, self.GetMainCollectionId(), pkey)
    def Get2LevelCollectionCommitmentKey(self, pkey=None):
        return self.ID(self.Get2LevelCollectionId(pkey), 0, 0)
    def Get2LevelCollectionCommitment(self, pkey=None):
        return GD[self.Get2LevelCollectionCommitmentKey(pkey)]
    def GetSizes(self):
        cntint = cntdic = cntlev = 0
        for k, value in GDL.items():
            if type(value) == list:
                cnt += 1
            elif type(value) == int:
                cntint += 1
            elif type(value) == dict:
                cntdic += 1
        return cntint, cntdic, cntlev
    def GetSizeGDL(self):
        return len(GDL)
    def GetHashes(self):
        hashes = []
        for k, item in GD.items():
            if type(item) not in [list, dict]:
                hashes += [item]
        return hashes
    def GetStructures(self):
        sts = []
        for k, item in GD.items():
            if type(item) is dict:
                sts += [item]
        return sts
    def GetTopLevelStructures(self): #leaf nodes
        all_sts = self.GetStructures()
        tlsts = []
        tlcolid = self.GetMainCollectionId()  # top level collection ids
        for st in all_sts:
            if st["colid"] == tlcolid:
                S = DS.St(dic=st)
                tlsts += [S]
        return tlsts
    def Get2LevelStructures(self):
        tlsts = self.GetTopLevelStructures()
        l2sts = {}
        for st in tlsts:
            l2colid = self.Get2LevelCollectionId(st["key"])
            for ast in self.GetStructures():
                if ast["colid"] == l2colid:
                    try:
                        l2sts[st["key"]] += [ast]
                    except:
                        l2sts[st["key"]] = [ast]
        return l2sts
    def GetHashId(self, colid, x, y):
        '''returns a hashid for a structure in collection id @colid, and position @x,y'''
        return self.ID(colid, x, y)

    def GetTopLevelLeafHashes(self):
        '''Top level leaf hashes are hashes of the top-level structures.'''
        leafhashes = []
        tlsts = self.GetTopLevelStructures()
        topcolid = self.GetMainCollectionId()  # i.e. top level collection's id.
        for st in tlsts:
            S = DS.St(dic=st)
            hashid = self.GetHashId(topcolid, S.x, S.y)
            leafhashes += [GD[hashid]]
        return leafhashes

    def Get2LevelLeafHashes(self):
        '''Top level leaf hashes are hashes of the 2-level structures.'''
        leafhashes = []
        level2ls = self.Get2LevelStructures()
        for key, l2sts in level2ls.items():
            lev2colid = self.Get2LevelCollectionId(key)
            for sts in l2sts:
                S = DS.St(dic=sts)
                hashid = self.GetHashId(lev2colid, S.x, S.y)
                print("l2,parent-key,sts"), "l2", key, sts
                leafhashes += [GD[hashid]]  # or [self.h(S.key,S.nkey,S.val,S.type_)]
        return leafhashes

    def keyExists(self, key):
        try:
            return GD[key]
        except:
            return False
    def GetHashTree(self, colid, x, y,):
        result = self.keyExists(self.GetHashId(colid, x, y))
        if not result:
            return
        else:
            hlist = [(x,y,result)]
            left = self.GetHashTree(colid, 2 * x, y + 1)
            if left is not None:
                hlist.extend((2*x,y+1,left))
            right = self.GetHashTree(colid, 2 * x + 1, y + 1)
            if right is not None:
                hlist.extend((2*x,y,right))
        return hlist
    def GetHashTree2(self, colid, x, y,lii):
        li =[]
        if self.keyExists(self.GetHashId(colid,x,y)):
            #print x,y,GD[self.GetHashId(colid, x, y)]
            li =[(x,y,GD[self.GetHashId(colid, x, y)])]
            lii +=li
            self.GetHashTree2(colid, 2 * x, y + 1,lii)
            self.GetHashTree2(colid, 2 * x + 1, y + 1,lii)
    def GetTopLevelHashTree(self):
        li =[]
        self.GetHashTree2(self.GetMainCollectionId(), 0, 0,li)
        return li
    def Get2LevelHashes(self):
        return dlsObj.GetHashTree(self.Get2LevelCollectionId(1), 0, 0)
    def GetCollectionIdentifiers(self):
        cids = {}
        for k, item in GD.items():
            if type(item) is list:
                cids[k] = item
        return cids
    def GetComplementaryNodesForS(self,leafNode):
        '''for an entry structure S with x and y position values,
        it returns complementary nodes of a hash-node @ position x, y in in hash-tree at '''
        compnodes = []
        StructureObj =leafNode
        cury = StructureObj.y  # current y
        curx = StructureObj.x
        print cury,curx
        while cury > 0:
            if curx % 2 == 1:
                # left side is sibling
                compnodes += [(curx-1,cury,self.get(self.ID(self.GetMainCollectionId(), curx - 1, cury)))]
            else:
                compnodes += [(curx+1,cury,self.get(self.ID(self.GetMainCollectionId(), curx + 1, cury)))]
            cury = cury - 1
            curx = curx / 2

        return compnodes

    def GetCommitmentFromCompNodes(self,compnodes,leaf_node):
        '''Returns a commitment for a leaf-node given compplementary nodes compnodes.'''
        #leaf_hash = self.h(leaf_node.xkey, leaf_node.xnkey, leaf_node.val, leaf_node.type_);
        leaf_hash = leaf_node.h()
        stack  = compnodes
        result = leaf_node.h()
        for i in range(0,len(stack)):
            if stack[i][0] % 2 == 1:
                result=self.h(result,stack[i][2])
            else:
                result=self.h(stack[i][2],result)
        return result

    def segment_rec_maping(self, segment):
        '''given a segment <(x1,y1,x2,y2),Ra || Rb>, map to a rectangle (x1,x2,y1,y2,v=0) in DLS dictionary to assign
        the value v={segKye:segValue}'''
        # find a rectangle (x1,x2,y1,y2,v=0) such that a segment is a diagonal of the rec.
        # set v=segment.label
        # update hash of the rectangle entry. Update root of the tree.
        key,matchedS = self.get_container_rec(segment)

        if not matchedS:
            return False
        print("--matched with:"),matchedS
        #@TO DO: create key for the segment and assign appropriate value of region code above/below the segment.
        if type(segment) == list or type(segment) is tuple:
            seg_key = Hasher().hash_segment(segment)
            seg_value= segment[4:]
            segment={seg_key:seg_value}
        try:
            matchedS.val.update(segment)
        except:
            matchedS.val = {}#
            matchedS,val.update(segment) #"Ra+Rb" #here insert segment as an entry in a dictionary value.

        #put matchedS after updating value
        self.put(self.ID(matchedS.colid,matchedS.xkey,matchedS.xnkey,matchedS.ykey), matchedS)
        #update hashes.
        ox = matchedS.x;
        oy = matchedS.y  # old x, y values of rectangle matchedS whose value is updated.
        oldhash = self.get(self.ID(matchedS.colid, ox, oy))
        newhash= matchedS.h()  # hash of updated entry
        self.UpdateHashes(matchedS.colid,ox,oy,oldhash,newhash)
        return True

    def get_container_rec(self,segment):
        '''return a leaf-level entry S that is rectanlge that contains segment as one of its diagonal.'''
        S = DS.St()
        for k, val in GD.items():
            if type(val) is dict:
                S = DS.St(dic=val)
                if self.isdiag(segment, S):
                    return k,S
        return None,None
    def isdiag(self, segment, Sobj):
        self.toforward_mode(segment)

        r1 = (segment[0] == Sobj.xkey and segment[2] == Sobj.xnkey) or (segment[0] == segment[2] == Sobj.xkey)
        if segment[1] < segment[3]:
            r2 = (segment[1] >= Sobj.ykey and segment[3] <= Sobj.ynkey)
        else:
            r2 = (segment[3] >= Sobj.ykey and segment[1] <= Sobj.ynkey)
        return r1 and r2

    def toforward_mode(self, segment):
        '''a segment (x1,y1,x2,y2) is received and returned (x1,x2,y1,y2):x1<x2'''
        x1, y1, x2, y2 = segment[0], segment[1], segment[2], segment[3]
        if x1 < x2:
            return True
        else:
            segment[0], segment[1], segment[2], segment[3] = x2, y2, x1, y1
        return False

    def plot_hash_tree(self, hash_tree):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        max_depth = self.depths[self.maincolid]
        pts={}
        for node in hash_tree:
            try:
                pts[node[1]] +=[node]
            except:
                pts[node[1]] =[]
                pts[node[1]] += [node]

        #fill up last level.
        l = [(i,max_depth,None,5*i,1) for i in range(0,2**max_depth)]
        for node in pts[max_depth]:
            idx = node[0]
            x,y = l[idx][-2],l[idx][-1]
            l[idx] = node + (x,y);
        pts[max_depth] = l
        ###
        #fill up upper levels.
        for d in range(max_depth-1,-1,-1):
            l = [(i,d,None) for i in range(0,2**d)]
            for i in range(0,len(l)):#node in pts[d]
                node = l[i]
                idx = node[0]*2
                xp = (pts[d+1][idx][-2] + pts[d+1][idx+1][-2])/float(2)
                yp = max_depth - node[1] +1
                l[i] = node + (xp,yp)

            for node in pts[d]:
                idx = node[0]
                n = l[idx]
                l[idx] = (n[0],n[1],node[2],n[3],n[4])
            pts[d] = l

        # graphing
        pts_dict_by_levels = pts
        fig = plt.figure(1, figsize=(5, 5), dpi=100)
        ax = fig.add_subplot(111)
        plt.gca()
        marks = ['ROOT'] +['#', '$', '@', '^', '0', '#', '^']
        marks = marks + marks[1:]
        maxx, maxy = -sys.maxint, -sys.maxint
        for k, v in pts_dict_by_levels.items():
            for node in v:
                if node[2] is None:
                    continue
                plt.text(node[-2], node[-1], marks[k])  # (node[-2],node[-1]))
                maxx, maxy = max(maxx, node[-2]), max(maxy, node[-1])
                # print ("\t"),node
        ax.set_xlim([0, maxx + 5])
        ax.set_ylim([0, maxy + 1])
        ax.set_title('split points')
        ax.grid(color='grey', linestyle='--', linewidth=.5)
        plt.show()
        return pts

import hashlib
from Hasher import Hasher
class MapProcess:

    def getHasher(self):
        return Hasher().getHasher()
    #functions
    def hash_segment(self,segment):
        s = self.getHasher()
        for e in segment:
            h = self.getHasher()
            h.update(str(e))
            s.update(h.hexdigest()[0:12])
        return s.hexdigest()[0:12]

    def shoelace_area(self,polygon_segs):
        '''polygon-segs is a sequence of connecting point in a polygon'''
        if len(polygon_segs)<3:
            return 0
        if polygon_segs[0] != polygon_segs[-1]:
            polygon_segs +=[polygon_segs[0]] #the first point must be last point in the series.

        lra = 0 #left-to-right product
        rla = 0
        for i in range(len(polygon_segs)-1):
            lra += polygon_segs[i][0]* polygon_segs[i+1][1]
            rla += polygon_segs[i][1]* polygon_segs[i+1][0]
        return abs(lra-rla)/float(2)


    def yblocks(self,xun=[], splits=[]):
        '''xun is list of unique x values; splits are split segments (x1,y1,x2,y2) in a polygon
        returns: {x:[(y1,y2),...]}'''
        xdict = {}
        for i in range(len(xun) - 1):
            xl = xun[i]
            xdict[xl] = []
        xdict[xun[-1]] = []  # add last ite,
        for i in range(len(splits)):
            x1 = splits[i][0]
            y1 = splits[i][1]
            x2 = splits[i][2]
            y2 = splits[i][3]
            try:
                slope= (y2-y1)/float(x2-x1)
            except:
                slope = 'inf'
            if y1 < y2:
                xdict[x1].append((y1, y2,slope))
            else:
                xdict[x1].append((y2, y1,slope))
        return xdict

    def polygon_dict(self,segment_tuples):
        '''n-tuple segment in segment_tuples represent a line segment.
        for ex. (x1,y1,x2,y2,prop1,prop2..) is a line segment with two end-points.
        The segment's inner boundary is denoted by inboundary_tag in the list inboundary_tags.
        This function returns a dictionary whose keys are 4-tuple segments, and values are the area code
        for above and below the segment.'''
        seg_dict = {}
        for i in range(0,len(segment_tuples)):
            if len(segment_tuples[i]) <4:
                raise "Error: Size of segment tuple must be 4."
            segKey =self.hash_segment(segment_tuples[i][0:4])
            seg_dict.update({segKey:segment_tuples[4:]})

        return seg_dict

    def load_segmentsfromtxt(self,infile):
        '''read infile txt file to read segments of a closed boundary line by line.'''
        segs = []
        #LOAD TXT FILE USING  NUMPY
        import numpy as np
        arr = np.loadtxt(infile)
        segsarr = arr[:,1:].astype(int) #ignore first column.
        return segsarr

    def segments_to_polyseq(self,segments):
        'segments of a polygon in either clockwise order or anticlockwise ordering.'
        import numpy
        polypts =[]
        if type(segments) == numpy.ndarray:
            segments = segments[0:,0:5].tolist()

        n = len(segments)
        for idx in xrange(0,5):
            seg = segments[idx]
            x1,y1,x2,y2,isswp = seg[0:5]
            #swapped.
            if isswp:
                x1,y1,x2,y2 = x2,y2,x1,y1
            polypts +=[(x1,y1)]
        #append start pt, append start pt as end pt to complete loop.
        polypts += [polypts[0]]
        return polypts

    def polysegments(self,points_list):
        pp = points_list
        segs = []
        n = len(pp) - 1  # the last point is also the first point.
        for i in range(n):  # make lines from points
            x1 = pp[i][0]
            x2 = pp[i + 1][0]
            y1 = pp[i][1]
            y2 = pp[i + 1][1]
            # collect all x-s in a dictionary
            # undict[x1]=1
            # undict[x2]=1
            # arrange lines so first x-coordinate is less that second
            if x1 <= x2:
                tt = (x1, y1, x2, y2,0) #0 means line from A-B is original
            else:
                tt = (x2, y2, x1, y1,1) #1 means direction change. B to A is original
            if tt not in segs:
                segs.append(tt)
        return segs
    def splitsegments(self,segs, xun):
        '''splits segments in segs where vertical lines pass through x-values in xun.
        that is, at points {(x,0): x in xun}'''
        ttt = []  # to hold all split lines
        for i in range(len(segs)):  # split lines
            x1 = segs[i][0]
            y1 = segs[i][1]
            x2 = segs[i][2]
            y2 = segs[i][3]
            li = xun.index(x1)
            hi = xun.index(x2)
            nn = hi - li
            if nn < 2:  # difference 0/1 means no split
                ttt.append((x1, y1, x2, y2)+segs[i][4:])
            else:  # chopping necessary
                yc = y1
                for j in range(li, hi):
                    xc = xun[j]
                    xn = xun[j + 1]
                    if xn < x2:
                        y = float(y2 - y1) / (x2 - x1) * (xn - x1) + y1
                        yn =  int(round(y))
                        ttt.append((xc, yc, xn, yn)+segs[i][4:])
                        yc = yn
                    elif xn == x2:  # last one
                        ttt.append((xc, yc, xn, y2)+segs[i][4:])
                # end for
        # end for
        return ttt

    def test_insert_sequence(self,dlsObj,cid,pkey,xun):

        dlsObj.InsertFirstKey(cid,pkey1,1)
        splitkey = pkey1

        for idx in range(1,len(xun),1):#exclude first, already inserted.
            print("\t K:"),xun[idx]
            try:
                dlsObj.InsertKey(cid,splitkey,xun[idx])
                splitkey = xun[idx]
            except Exception, e:
                print("Exception:"),str(e),xun[idx]

        for k in GD:
            if type(GD[k]) is dict:
                print GD[k]
            if type(GD[k]) is list:
                print("\t"),GD[k]
            if type(GD[k]) is str:
                print("\t\t"),GD[k]

    def enclosedby(self,sk,skp,k):
        #key,nextkey,new-key
        if (sk < k < skp) or (k > sk >= skp) or (sk >= skp > k):
            return True
        return False

    def getsplitkey(self,newkey,refxkey=None,dimension='y'):
        #search for keys to insert newkey in x-dimension or y-dimensin.
        #To Do: all leaves can be stored in different dictionary so that we won't have to search big dictionary GD
        split_keys = []
        if dimension == 'y':
            #return ykeys as split keys for the xkey in x dimension.
            if refxkey == None:
                print("xkey is missing.")
                return None
            for item, value in GD.items():
                if type(value) == dict:
                    S = DS.St(dic = value)
                    if S.xkey == refxkey and self.enclosedby(S.ykey,S.ynkey,newkey):
                        split_keys +=[S]
            return split_keys
        else:
            for id,value in GD.items():
                if type(value) == dict:
                    S = DS.St(dic=value)
                    if self.enclosedby(S.xkey, S.xnkey, newkey):
                        split_keys +=[S]#+=[id]
            return split_keys
        return None

    def create_dls_for_slabbed_polygon(self,dlsObj,cid,xun,yblocks_onx):
        #Random way of entering keys in the Dictionary
        import random,copy
        backxun = copy.deepcopy(xun) #back up.
        ikx = xun.pop(random.randint(0,len(xun)-1)) #random x-value
        Si= dlsObj.InsertFirstKey(cid,ikx,0,type=1) #first key cid,iks,ykey=0
        '''Here a random x-key->xnew is chosen from list of unique-x values. 
        It is inserted by splitting existing keys stored by St objects in listS in the GD.'''

        while len(xun) > 0:
            xnew = xun.pop(random.randint(0,len(xun)-1))
            listS= self.getsplitkey(xnew,dimension='x') #list of Struct St's #y-split-key,y-split-key-next
            if not listS:
                continue
            for splitS in listS:
                #@TO DO: do not split wrap-up entries.
                newS,newSN=dlsObj.InsertKey(collectid = cid, structS=splitS,xnewkey = xnew, ynewkey = None, dimension = 'x')
            #endfor
        #end while
        '''For each x-keys insert y-values.'''

        while len(backxun) > 0:
            ikx = backxun.pop(random.randint(0,len(backxun)-1))
            yuniqkx = yblocks_onx[ikx][1]  # unique-y values on x-ikx; the second list.
            while len(yuniqkx) > 0:
                ynew = yuniqkx.pop(random.randint(0,len(yuniqkx)-1))
                S= self.getsplitkey(ynew,refxkey=ikx, dimension='y') #y-split-key,y-split-key-next key
                if not S: continue
                for stS in S:
                    S,SN=dlsObj.InsertKey(collectid = cid,structS=stS,xnewkey = None, ynewkey = ynew, dimension = 'y')
                    #split the entry structS=S and add new ykey


    def  non_overlaping_yspans(self,ytuples):
        '''given list of (yi,yj,li) that means line li has y-span yi-yj, construct a coverup y-span
        such that the y-span's form a complete lut-entries. For instance:
        (1,3,a) and (2,4,b) are overlapping y-spans; it must result a tuple (1,4,{a,b}) and (4,1,{})
        tuples.'''
        #sort each tuple in ytuples.ytuples[jdx-1][1]

        for i in range(len(ytuples)):
            if ytuples[i][1] < ytuples[i][0]:
                ytuples[i] = (ytuples[i][1],ytuples[i][0]) + ytuples[i][2:]

        #sort list of tuples by a tuple's first element.
        from operator import itemgetter
        #slower version:
        #ytuples = sorted(ytuples, key=lambda tuple: tuple[0])
        #faster version:
        ytuples = sorted(ytuples,key=itemgetter(0))
        recs ={}
        #now merge the overlapping y-spans
        idx = 0
        jdx = idx
        while idx < len(ytuples):
            y0,y1,line = ytuples[idx]
            try:
                recs[y0][0].append(ytuples[idx])
            except:
                recs[y0] = []
                #recs[y0] += [ytuples[idx]]
                recs[y0].append([ytuples[idx]])
                bot,top = y0,y1 #bottom and top of y-value.
                recs[y0].append((bot,top))
            jdx = idx+1
            while jdx < len(ytuples) and ytuples[jdx][0] < top: #ytuples[jdx-1][1]:
                recs[y0][0].append(ytuples[jdx])
                #recs[y0][1]=(min(recs[y0][1][0],ytuples[jdx][0]),max(recs[y0][1][1],ytuples[jdx][1]))
                bot,top = min(bot,ytuples[jdx][0]),max(top,ytuples[jdx][1])
                recs[y0][1] = bot,top
                jdx +=1
            #endwhile
            idx = jdx
        #endwhile
        from collections import OrderedDict
        #ordred = OrderedDict(sorted(recs.items(),key=lambda t:t[0]))
        yvalues=[]
        #extract bot/top values of each rectangles.
        for k,v in recs.items():
            ybot,ytop = v[1]
            yvalues +=[ybot,ytop]
        yvalues = sorted(list(set(yvalues)))
        print("\t"),yvalues
        return recs,yvalues

    def graph(self,pp,splitsegs,xun,yblocks_onx=None):
        import numpy as np, sys
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from Hasher import Helper
        from collections import OrderedDict

        fig = plt.figure(1, figsize=(5, 5), dpi=100)
        ax = fig.add_subplot(111)
        #original polygon
        #ax.plot([t[0] for t in pp], [t[1] for t in pp], color='#6699cc', alpha=1,linewidth=1.5, solid_capstyle='round', zorder=2)
        #splits
        maxy = -sys.maxint
        for split in splitsegs:
            maxy = max(maxy,max([split[1], split[3]]))
            ax.plot([round(split[0]), round(split[2])], [split[1], split[3]], '-',linewidth=0.5, c='red')
            ax.plot([split[0], split[2]], [split[1], split[3]], 'x', c='blue') #split-points
            #ax.plot([split[0], split[0]], [0, 14], '--', c='red', linewidth=0.2)
            pass
        ax.set_title('split points')
        major_xticks = np.arange(-1, max(xun)+1, 1)
        major_yticks = np.arange(-1, maxy+1, 1)
        ax.set_xticks(major_xticks)
        ax.set_yticks(major_yticks)
        ax.yaxis.grid(which="major", color='blue', linestyle='--', linewidth=0.2)

        #ax.set_xticks(minor_yticks, minor=True)
        #ax.grid(color='gray', which="major",linestyle='--', linewidth=0.5)
        #if yblocks_onx plot vertical bars for vertial slabs; use y-spans to draw recs on the slabs.
        orderedyblocks_onx = OrderedDict(sorted(yblocks_onx.items(),key=lambda t:t[0]))
        Helper.addItem(orderedyblocks_onx)
        for xkey,yblocks in orderedyblocks_onx.items():
            #plot vertical bar.
            ax.plot([xkey,xkey],[0,maxy], '-',linewidth=0.75, c='green')
            yspans = yblocks[1]
            for yvalue in yspans:
                if Helper.getNextKey(xkey) !=None:
                    ax.plot([xkey, Helper.getNextKey(xkey)], [yvalue, yvalue], '-', linewidth=0.5,c='black')
                    pass
            pass
        plt.show()

    def test_create_dls_for_slabbed_polygon(self):

        import datetime, copy, numpy as np, hashlib, sys, collections, random

        from_segfile = True
        home = "D:/workspace/sqdm-repo/sqdm/out/tmp"
        infile1 = home + "/usa.prj.lbl.txn.int.txt"
        infile2 = home + "/states.prj.lbl.txn.int.txt"

        poly1 = [(5, 6), (2, 5), (4, 3), (5, 2), (6, 3), (7, 1), (9, 2), (11, 3), (13, 5), (14, 6), (13, 7),
                 (12, 8), (11, 9), (9, 12), (8, 13), (5, 12), (4, 10), (3, 9), (5, 6)]
        poly2 = [(11,9),(10,8.5), (9.5,8),(9,7),(8,7.5),(6,8),(4,9),(7,6), (6,5),(7,4),(9,5),(8,3),(12,1),(14,2),(15,1),(16,2),(17,4),(15,6),
                 (14,4),(11,5),(13,7),(11,9)]
        poly3 = [(11,15),(11,16),(7,16),(5,14),(7,15),(5,11),(7,12),(5,9 ),(9,7),(11,9),(12,7),(10,5),(12,4),(14,6),
                 (18,4),(20,7),(23,6),(24,9),(22,11),(20,10),(21,12),(19,10),(19,13),(18,12),(17,14),(15,12),(18,9),(16,8),
                 (14,11),(15,14),(13,15),(11,12),(9,14),(11,15)]

        MP = MapProcess()
        segs = []
        if from_segfile:
            segs = MP.load_segmentsfromtxt(infile1)
            segs =  tuple(map(tuple, segs))
            pp = MP.segments_to_polyseq(segs)
            print segs[0:10]
            print pp[0:10]
        else:
            pp = poly1
            pp = [(int(t[0]), int(t[1])) for t in pp]  # convert to integers.
            segs = MP.polysegments(pp)  # collection of tuple (x1,y1,x2,y2)
            #tag each segments:
            segs = [t+('A','B','s') for t in segs]

        xlist = [[t[0], t[2]] for t in segs]
        xlist = [x for sublist in xlist for x in sublist]
        xun = sorted(list(set(xlist)))

        splitsegs = MP.splitsegments(segs, xun)  # to hold all split lines

        yblocks_onx = MP.yblocks(xun=xun, splits=splitsegs) #To DO: remove slope from yblocks.
        print("yblocks_onx"),yblocks_onx

        for x,yblocklst in yblocks_onx.items():
            #sort all y-values in the x-interval.
            print("key x:"),x, yblocklst
            nonoverlaping_yspans,uniqyvalues = self.non_overlaping_yspans(yblocklst)
            yblocks_onx[x] = [nonoverlaping_yspans,uniqyvalues]

        dlsObj = DLS()
        cid = dlsObj.Init()
        print("collected unique x-values in map."),xun
        # graphing
        MP.create_dls_for_slabbed_polygon(dlsObj, cid,copy.deepcopy(xun),copy.deepcopy(yblocks_onx))
        MP.graph(pp,splitsegs,xun,yblocks_onx)

        print("dls stats:"),dlsObj.dlsstats()
        leaf_nodes = dlsObj.GetTopLevelStructures()
        print(len(leaf_nodes))
        for lnode in leaf_nodes:
            print dlsObj.ID(dlsObj.maincolid,lnode.xkey,lnode.ykey),lnode
            pass
        print("hash tree"), len(dlsObj.GetTopLevelHashTree())

        for hash in dlsObj.GetTopLevelHashTree():
            print hash
            pass

        lidx = random.randint(0,len(leaf_nodes)-1)
        print "leaf",lidx,leaf_nodes[lidx]
        compnodes= dlsObj.GetComplementaryNodesForS(leaf_nodes[lidx])
        print("complementary-nodes"),compnodes
        print("root"),dlsObj.GetCommitmentFromCompNodes(compnodes,leaf_nodes[lidx])==dlsObj.GetMainCollectionCommitment()
        assert dlsObj.GetMainCollectionCommitment() == dlsObj.GetCommitmentFromCompNodes(compnodes,leaf_nodes[lidx]),"root did not matched!!"
        #dlsObj.plot_hash_tree(dlsObj.GetTopLevelHashTree())

        #save set of segments as as a dictionary
        #save GD as a dictionary
        gdfile = "GD"
        self.save(GD, gdfile)
        print("GD saved in"), gdfile
        self.save(splitsegs,'segments')
        print("splitsegs saved in "),'segments'

    def load_data(self,infile):
        import json
        data = None
        infile = infile.replace(".json",'')
        infile = infile.replace(".p",'')
        infile = infile.replace(".pickle",'')
        for ext in ['.json','.pickle']:
            try:
                with open(infile + ext, 'r') as fp:
                    data = json.load(fp)
            except:
                continue
        return data

    def save(self,data,outfile):
        import json,sys,os.path,random
        #delete if exists
        if os.path.exists(outfile):
            if not os.path.isdir(outfile): os.remove(outfile)
            else:outfile += '/data-'
        outfile = outfile.replace(".json",'')
        outfile = outfile.replace(".p",'')
        outfile = outfile.replace(".pickle",'')
        with open(outfile+'.json', 'w') as fp:
            json.dump(data, fp,sort_keys=True,indent=4)

        print("json data saved in "),outfile

    def swap_endpoints(self,segments):
        '''for each (x1,y1,x2,y2) in segments, construct (x2,y2,x1,y1) if x2<x1.'''
        for i in range(len(segments)):
            x1,y1,x2,y2 = segment[i][0],segment[i][1],segment[i][2],segment[i][3]
            if x2<x1:
                segments[i] = (x2,y2,x1,y1,x2<x1) #isswapped end points A and B.
            else:
                segments[i] = (x1,y1,x2,y2,x2<x1)
        return segments

    def region_labeling(self,segments):
        '''for each segment (x1,y1,x2,y2, isswapped) that are end points for a
        segment AB such that x1<x2; construct  a tuple (x1,y1,x2,y2,Ra,Rb)
        where Ra and Rb are labels for regions above and below the segment A(x1,y1)--B(x2,y2)'''
        import osgeo, ogr, shapely


class TestMapProcess:
    mp  = MapProcess()
    def test_hash_segment(self):
        segment = (3341,41,50,13)
        print(self.mp.hash_segment(segment))
    def test_polygon_dict(self):
        segments =[(1,2,34,5),(1,24,65)]
        print self.mp.polygon_dict(segments)
    def test_shoelace_area(self):
        polysegs = [(2,4),(3,-8),(1,2)]
        poly1 = [(5, 6), (2, 5), (4, 3), (5, 2), (6, 3), (7, 1), (9, 2), (11, 3), (13, 5), (14, 6), (13, 7),
                 (12, 8), (11, 9), (9, 12), (8, 13), (5, 12), (4, 10), (3, 9), (5, 6)]
        poly2 = [(11,9),(10,8.5), (9.5,8),(9,7),(8,7.5),(6,8),(4,9),(7,6), (6,5),(7,4),(9,5),(8,3),(12,1),(14,2),(15,1),(16,2),(17,4),(15,6),
                 (14,4),(11,5),(13,7),(11,9)]
        print("area by shoelace"),self.mp.shoelace_area(poly2)

    def test_non_overlaping_yspans(self):

        yspans = [(11,14,'h'),(1,3,'a'),(3,4,'b'),(4,6,'c'),(5,6,'d'),(5,8,'e'),(7,8,'f'),(7,10,'g')]
        yspans = [(8, 9, -0.5), (9, 12, 1.5), (11, 15, 2.0), (11, 12, 0.5), (14, 16, 1.0), (14, 15, 0.5)]
        print("initial yspans"),yspans
        yspans,yvalues=self.mp.non_overlaping_yspans(yspans)
        print("after sorting each tuples:")
        print(yspans)
        print("yspans-keys,yvalues"),yspans.keys(),yvalues
        ##
        #sort tuples by first y-value.

    def test_segment_mapping(self):
        global GD
        dlsObj = DLS()
        GD = self.mp.load_data("GD")
        segments = self.mp.load_data("segments")
        missedmapping=[]
        #horizsegs = [t for t in segments if t[1]==t[3]]
        vertsegs = [t for t in segments if t[0]==t[2]]
        segs = segments
        for segment in segs:
            print("--"),segment
            ismapped=dlsObj.segment_rec_maping(segment)
            if not ismapped:
                missedmapping+=[segment]

        print("Missing mapping:"), 'cause: horizontal and vertical lines. '
        for seg in missedmapping:
            print seg
        #save GD
        self.mp.save(GD,"GD1")
        print("Total Segments"),len(segments)
        print("Segments Missing Matching"), len(missedmapping)
        print("Segments Matched"), len(segments)-len(missedmapping)

        ##plot mapping.

        leaf_nodes = dlsObj.GetTopLevelStructures()

    def test_segments_to_polyseq(self):
        poly2 = [(11,9),(10,8.5), (9.5,8),(9,7),(8,7.5),(6,8),(4,9),(7,6), (6,5),(7,4),(9,5),(8,3),(12,1),(14,2),(15,1),(16,2),(17,4),(15,6),
                 (14,4),(11,5),(13,7),(11,9)]

        segs = self.mp.polysegments(poly2)
        segs = self.test_load_segmentstxt()
        print segs[0:4,0:]
        print("-----"),type(segs)
        pts = self.mp.segments_to_polyseq(segs)
        print pts
        print len(pts)


    def test_load_segmentstxt(self):
        home ="D:/workspace/sqdm-repo/sqdm/out/tmp"
        infile1 = home+"/usa.prj.lbl.txn.int.txt"
        infile2 =home +"/states.prj.lbl.txn.int.txt"
        segsarr = self.mp.load_segmentsfromtxt(infile1)
        segs =segsarr[:,1:-3]/10000
        segs = segs.tolist()[0:1000]

        #plot
        #self.plot_boundary(segs)
        return segsarr

    def plot_boundary(self,segments):
        import numpy as np, sys
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        fig = plt.figure(1, figsize=(20,20), dpi=100)
        ax = fig.add_subplot(111)
        #splits
        maxy = maxx=-sys.maxint
        for split in segments:
            x1,x2,y1,y2 = split[0], split[2],split[1], split[3]
            maxy = max(maxy,max([y1,y2]))
            maxx = max(maxx, max([x1,x2]))
            #ax.plot([x1,x2], [y1,y2], '-',linewidth=0.5, c='red')
            #ax.plot([split[0], split[2]], [split[1], split[3]], 'x', c='blue') #split-points
            #ax.plot([split[0], split[0]], [0, 14], '--', c='red', linewidth=0.2)
            pass
        ax.set_title('boundary points')
        print(maxx,maxy)
        #major_xticks = np.arange(-1, maxx+1, 10000)
        #major_yticks = np.arange(-1, maxy+1, 10000)
        #ax.set_xticks(major_xticks)
        #ax.set_yticks(major_yticks)
        ax.yaxis.grid(which="major", color='blue', linestyle='--', linewidth=0.2)

        #ax.set_xticks(minor_yticks, minor=True)
        #ax.grid(color='gray', which="major",linestyle='--', linewidth=0.5)
        #if yblocks_onx plot vertical bars for vertial slabs; use y-spans to draw recs on the slabs.

        plt.show()

    def run_tests(self):
        import traceback
        from colorama import Fore, init
        init()
        try:
            self.test_hash_segment()
            print("Testing hash_segment,..test passed.")
        except Exception, e:
            print(str(e))
        try:
            self.test_polygon_dict()
            print("Testing polygon_dict,...test passed.")
        except Exception, e:
            print(str(e))
        try:
            self.test_shoelace_area()
            print("Testing shoelace_area,...test passed.")
        except Exception, e:
            print(str(e))
        try:
            self.test_non_overlaping_yspans()
        except Exception, e:
            print "Fetal:In test_non_overlapping_yspans",str(e)

        self.test_segment_mapping()
        try:
            1==1
        except Exception, e:
            print Fore.BLUE+"Fetal: In test_segment_mapping--"+str(e)+Fore.RESET
        #self.test_load_segmentstxt()
        self.test_segments_to_polyseq()

if __name__ == "__main__":
    MapProcess().test_create_dls_for_slabbed_polygon()
    #tmp = TestMapProcess()
    #tmp.run_tests()

    '''
    1) Given a series of points in a polygon, develop an algorithm to compute area of 
    polygon.
    2) Given a point, return 10 closest points/attributes to the point.
    3) Delegate (P, K1 -->K2)
    4) Merge the areas belonging to same K1. (P1 =~ P2). Remove line shared by P1 and P2. needs history
    of the ownerships.
    5) Add structures like Road(speedlimit,divided,one way)
    /utility/ lines for creating Auxilary maps.
    6) Find queries like 5 nearest gas stations for a point.
    7) Point Location problem. Given a point,which is the area that contains the point.'''



