import copy
import sys
import hashlib
if sys.version_info < (3, 6):
    import sha3

from datatypes import *

class TestOmt(object):

    def testInsertSNode(self,tobalancetree=False):
        collection_dict = {22: 0, 3: 0, 24:2,45:6,8:9,100:99,89:90, 22:333}
        print collection_dict
        mainOmt = OMT(collection=[]);
        mainOmt.bulkinsert(collection_dict, balancetree=tobalancetree)
        print mainOmt
        mainOmt.stats()

        #Test Existence
        #TODO: To test the existence of a key 'k', i) search a snode (k,nk,v) in the collection.
        #TODO: ii) derive complementary hashes for this snode.
        searchsnode = S(k=3,val=2)
        proofsnode, f = mainOmt.serchsnode(searchsnode.k) #found snode
        print proofsnode
        commlementaries = mainOmt.complementaries(proofsnode)
        comproot = mainOmt.rootofcomphashes(commlementaries, proofsnode)
        print str(searchsnode), " exists in mainOMT?", mainOmt.verifywithomtroot(comproot) and searchsnode == proofsnode

        udi = UDI()
        udi.addoperands('searchkey',searchsnode.k)
        #udi.addoperands('searchval',searchsnode.val)
        #udi.addoperands('proofsnode',proofsnode)
        #udi.addoperands('proofvos',commlementaries)
        #udi.addoperands('omtref',mainOmt.omtid)
        udi.addoperands('omtroot',mainOmt.omtroot().hashv)
        udi.toJson("../out/testjson.json")
        print udi
        # Test of non-existence
        #TODO: To test non-existence of a key k, i) search a snode (k',nk',v) in the collection such that 'k < k < nk'.
        #TODO: ii) derive complementary hashes for this snode.
        print("Test of non-existence")
        searchsnode = S(k=22,val=0)
        proofsnode, f = mainOmt.serchsnode(searchsnode.k) #found snode
        print proofsnode
        commlementaries = mainOmt.complementaries(proofsnode)
        comproot = mainOmt.rootofcomphashes(commlementaries, proofsnode)
        print str(searchsnode), " exists in mainOMT?", mainOmt.verifywithomtroot(comproot) and searchsnode == proofsnode

        #Test of min-max key
        print mainOmt.extractminmax()

    def testPolygonSegmentsOmt(self,tobalancetree=False):
        tcb = TCB()
        seg_dict =[]
        from simple_polygon.simplepoly import Polygon, SlSegment
        poly6 = [(4, 3), (2, 1), (4, 0), (5, 2), (6.5, 1), (6, 3), (8, 3.5), (5, 4),(8, 5), (1, 6), (2, 4), (4, 3)]

        '''
        #do these when using a conneted polygon.
        poly = Polygon(poly6)
        poly.xorderednamedsegs()   
        seg_dict = poly.tosegsdict(segment_id_type='linehash')        
        '''

        segments = [(1,2, 5,5),(2,3,10,1),(3,2,8,4),(4,3,6,1),(7,1,9,2)]#unique x-values x-ordinates
        segments = [(1,1,18,11),(3,7,5,3),(2,10,16,5),(10,6,13,4),(11,12,17,9)]
        segments = [(11, 12, 17, 9), (1, 1, 18, 11), (2, 10, 16, 5), (10, 6, 13, 4), (3, 7, 5, 3)]  # uniq x and y

        segObjlist=[]
        index =0
        for tup in segments:
            sego = SlSegment()
            sego.tupleToSeg(index,tup)
            seg_dict.append(sego.dictentry_linehashkey())
            index +=1

        print segObjlist
        print("segment-dict"),seg_dict

        segsomt = OMT(omtref='segomt', collection=[])
        segsomt.bulkinsertsegs(seg_dict, balancetree=tobalancetree)

        #initialize TCB or broadcaset commitment to this segment's set.
        udiomtinit = UDI()
        udiomtinit['omtref'] = segsomt.getomtref()
        udiomtinit['omtroot'] = segsomt.omtroot()
        udiomtinit['vos'] = []
        trans = Transaction('omtinit',udiomtinit)
        certificate = tcb.processTransac(trans)

        #for a valid certificate, do start setting new values for snodes.
        if certificate:
            #Update Segment OMT
            for seg_entry in seg_dict:
                print("--"),seg_entry
                searchkey,newval = seg_entry
                updatesnode = segsomt.serchsnodebykey(searchkey)#going to update this snode
                udinodeupdate = UDI()
                udinodeupdate['updatesnode'] = updatesnode
                udinodeupdate['updatekey'] = searchkey
                udinodeupdate['updatevalue'] = newval
                udinodeupdate['omtref'] = segsomt.getomtref()
                udinodeupdate['vos'] = segsomt.complementaries(updatesnode)
                trans = Transaction('nodeupdate', udinodeupdate)
                certificate =tcb.processTransac(trans) #now verify inputs,update omt root etc.
                #for valid certificate update the hashes.
                if certificate:
                    segsomt.update_snode(updatesnode, newval) #then update hashes locally.
                print("Updated segments' OMT.")
        del certificate
        del seg_dict
        print segsomt
        return segsomt,tcb

    def createEventOmt(self,segsomt,tcb):
        print("\n\n Creating Event Tree")
        #create an event tree.

        eventomt = OMT('eventomt',collection=[])
        for segsnode in segsomt.collection:
            x1,y1,x2,y2= segsnode.val[0], segsnode.val[1], segsnode.val[2], segsnode.val[3]

            starteventvalue = endeventvalue = 0 #segsnode.val when updating this omt.
            if eventomt.isempty():
                eventomt.insertfirstkey(x1, starteventvalue, collid=1, stype=STYPE['LUT'])  # inserted node returned
                # get enclosing snode for k,v
                splitsnode = eventomt.getsplitsnode(x2)
                balancingnode = eventomt.getbalancingnode(splitsnode)
                eventomt.insertkeybalanced(splitsnode, balancingnode, x2, endeventvalue)

            else:
                #create start event
                # get enclosing snode for k,v
                splitsnode = eventomt.getsplitsnode(x1)
                if not splitsnode:
                    print("S-node with this key may already exists in the OMT")
                    continue
                balancingnode = eventomt.getbalancingnode(splitsnode)
                eventomt.insertkeybalanced(splitsnode, balancingnode, x1, starteventvalue)

                #create end event
                # get enclosing snode for k,v
                splitsnode = eventomt.getsplitsnode(x2)
                if not splitsnode:
                    print("S-node with this key may already exists in the OMT")
                    continue
                balancingnode = eventomt.getbalancingnode(splitsnode)
                eventomt.insertkeybalanced(splitsnode, balancingnode, x2, endeventvalue)

        #initialize TCB or broadcaset commitment to this segment's set.
        udiomtinit = UDI()
        udiomtinit['omtref'] = eventomt.getomtref()
        udiomtinit['omtroot'] = eventomt.omtroot()
        trans = Transaction('omtinit', udiomtinit)
        certificate = tcb.processTransac(trans)
        tcb.ICurrentState()
        del trans
        del certificate
        #Update event-Snode to put value as (x1,y1,x2,y2,START/END).
        for segsnode in segsomt.collection:
            x1,y1,x2,y2= segsnode.val[0], segsnode.val[1], segsnode.val[2], segsnode.val[3]
            searchkey, eventype,newval = x1,EVENT_TYPE['START'],segsnode.val + (EVENT_TYPE['START'],) #adding item in tuple.
            updatesnode = eventomt.serchsnodebykey(searchkey)  # going to update this snode
            udiCreateEvent = UDI()
            udiCreateEvent['updatesnode'] = updatesnode
            udiCreateEvent['updatekey'] = searchkey
            udiCreateEvent['updatevalue'] = newval
            udiCreateEvent['eventtype'] = eventype
            udiCreateEvent['voseventomt'] = eventomt.complementaries(updatesnode)
            udiCreateEvent['eventomtref'] = eventomt.getomtref()
            udiCreateEvent['segmentsnode'] = segsnode
            udiCreateEvent['segomtref'] = segsomt.getomtref()
            udiCreateEvent['vossegomt'] = segsomt.complementaries(segsnode)
            trans = Transaction('eventnodeupdate', udiCreateEvent)
            certificate=tcb.processTransac(trans) #calling tcb to process this transaction
            if certificate:
                eventomt.update_snode(updatesnode, newval)  # then update hashes locally.

            #update event for Stop end point.
            searchkey, eventype,newval = x2,EVENT_TYPE['STOP'],segsnode.val+ (EVENT_TYPE['STOP'],) #adding item in tuple.
            updatesnode = eventomt.serchsnodebykey(searchkey)  # going to update this snode
            udiCreateEvent = UDI()
            udiCreateEvent['updatesnode'] = updatesnode
            udiCreateEvent['updatekey'] = searchkey
            udiCreateEvent['updatevalue'] = newval
            udiCreateEvent['eventtype'] = eventype
            udiCreateEvent['voseventomt'] = eventomt.complementaries(updatesnode)
            udiCreateEvent['eventomtref'] = eventomt.getomtref()
            udiCreateEvent['segmentsnode'] = segsnode
            udiCreateEvent['segomtref'] = segsomt.getomtref()
            udiCreateEvent['vossegomt'] = segsomt.complementaries(segsnode)
            trans = Transaction('eventnodeupdate', udiCreateEvent)
            certificate=tcb.processTransac(trans)  #calling tcb to process this transaction
            if certificate:
                eventomt.update_snode(updatesnode, newval)  # then update hashes locally
        print("Completed updating eventomt")
        #print eventomt
        tcb.ICurrentState()
        return segsomt,tcb, eventomt

    def processEvents(self,segsomt,tcb,eventomt):
        '''

        :param segsomt: stores segments committed for finding intersections
        :param tcb: is a TCB that monitors the important omt root values, executes some secure operations
        for verifications
        :param eventomt: stores events ordered by x-value of the segment's end-points.
        :return: Writes interesection points in BlockChain which is not modified by any external entities.
        Q. Is root a commitment to the value for a key in a collection?
        '''
        activeomt = OMT('activeomt',collection=[])

        #first event.
        initeventsnode = eventomt.serchsnodebykey(eventomt.omtstats['minkey'])
        print initeventsnode
        udiProcessEvent = UDI()
        udiProcessEvent['event'] = initeventsnode
        udiProcessEvent['eventvos'] = eventomt.complementaries(initeventsnode)
        udiProcessEvent['activeomtref'] = activeomt.getomtref()
        udiProcessEvent['eventomtref'] = eventomt.getomtref()

        trans = Transaction("initprocessevent",udiProcessEvent)
        certificate = tcb.processTransac(trans)
        if certificate:
            k,v = initeventsnode.k,initeventsnode.val
            print k, v[1], v[-1]
            activeomt.insertfirstkey(v[1],v,collid=1,stype=STYPE['LUT'])
        print activeomt
        print tcb.ICurrentState()


        #second event and other events
        currenteventnode = initeventsnode
        while currenteventnode.nk > currenteventnode.k:
            nxteventsnode = eventomt.serchsnodebykey(currenteventnode.nk)
            k, y1 = nxteventsnode.k, currenteventnode.val[1]

            splitactivesnode = activeomt.getsplitsnode(y1,searchkind=STYPE['LUT'])#split node in activeomt
            activebalancingnode = activeomt.getbalancingnode(splitactivesnode)
            activeomt.insertkeybalanced(splitactivesnode, activebalancingnode, y1, currenteventnode.val)

            activesegabove = segsomt.serchsnodebykey(nxteventsnode.nk)
            activesegbelow = segsomt.serchsnodebykey(nxteventsnode.k,findprev=True)

            '''
            udiProcessEvent = UDI()
            udiProcessEvent['event'] = nxteventsnode
            udiProcessEvent['eventvos'] = eventomt.complementaries(currenteventnode)
            udiProcessEvent['activeomtref'] = activeomt.getomtref()
            udiProcessEvent['eventomtref'] = eventomt.getomtref()
            udiProcessEvent['segabove'] = activesegabove
            udiProcessEvent['segabovevos'] = eventomt.complementaries(activesegabove)
            udiProcessEvent['segbelow'] = activesegbelow
            udiProcessEvent['segbelowvos'] = eventomt.complementaries(activesegbelow)

            trans = Transaction("furhterprocessevent",udiProcessEvent)
            certificate = tcb.processTransac(trans)
            '''
            currenteventnode = nxteventsnode



print("Completed...............")
testOmt = TestOmt()
segsomt,tcb =testOmt.testPolygonSegmentsOmt(True)
segsomt.stats()
print
segsomt,tcb,eventomt = testOmt.createEventOmt(segsomt,tcb)
print
eventomt.stats()
print
testOmt.processEvents(segsomt,tcb,eventomt)

