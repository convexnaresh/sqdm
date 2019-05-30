
import math
import ogr
PI = math.pi
import numpy as np
import random
import copy
class Metri:

    def __init__(self):
         pass

    @classmethod
    def upolsby_popper(self,area, perimeter):
        '''
        unit ratio value
        :param area: area value
        :param perimeter: perimieter value
        :return: area/perimeter-squaredd normalized
        '''
        try:
            return 4 * PI * area / perimeter ** 2
        except ZeroDivisionError:
            return math.nan

    @classmethod
    def inv_root_upolsby_popper(self,area, perimeter):
        '''
        unit ratio value
        :param area: area value
        :param perimeter: perimieter value
        :return: area/perimeter-squaredd normalized
        '''
        try:
            return pow(4 * PI * area / float(perimeter ** 2), 0.5)
        except ZeroDivisionError:
            return math.nan

    @classmethod
    def dpolsby_popper(self, dist_prop):
        '''
        A(D)/PERI(D)^2
        :param dist_prop:
        :return:
        '''
        for distkey in dist_prop.keys():
            dist_prop[distkey]["ISO_PERIMETRIC_IDX"] = self.upolsby_popper(dist_prop["AREA"][distkey],
                                                                         dist_prop["PERI"][distkey])

    @classmethod
    def dpolsby_popper2(self, distkey,dist_prop):
        '''
        A(D)/PERI(D)^2
        :param dist_prop:
        :return:
        '''
        dist_prop[distkey]["ISO_PERIMETRIC_IDX"]= self.upolsby_popper(dist_prop[distkey]["AREA"],
                            dist_prop[distkey]["PERI"])

    @classmethod
    def moment_area(self,dist_blocks_keys,block_prop):
        Ia = 0
        A = 0
        #each blockid bid
        for bid in dist_blocks_keys:
            Ia += block_prop[bid]["AREA"]*block_prop[bid]["BLKCTR_TO_DIST_CTR"]
            A += block_prop[bid]["AREA"]
            #block to district length
        return 2*PI*Ia/float(A*A)

    @classmethod
    def moment_popu(self,dist_blocks_keys,blocks_prop):
        Ip = 0
        A = 0
        P = 0
        for bid in dist_blocks_keys:
            Ip += blocks_prop[bid]["POPULATION"]*blocks_prop[bid]["BLKCTR_TO_DIST_CTR"]#block to district length
            A += blocks_prop[bid]["AREA"]
            P += blocks_prop[bid]["POPULATION"]

        #normalize
        return 2*PI*Ip/float(A*P)

    @classmethod
    def distance_sq(self,x1,y1,x2=None,y2=None):
        if x2 and y2:
            return (x1-x2)**2 + (y1-y2)**2
        else:
            #treat x1 and y1 as 2-tuple (x1,y1) for point
            if (type(x1) in [list,tuple]) and (type(y1) in [list,tuple]):
                x1, y1, x2, y2 = x1[0], x1[1], y1[0], y1[1]
                return (x1 - x2) ** 2 + (y1 - y2) ** 2
            else:
                return None

    @classmethod
    def distance_tpts_sq(self,tuple_pt1,tuple_pt2):
        x1,y1,x2,y2 = tuple_pt1[0],tuple_pt1[1],tuple_pt2[0],tuple_pt2[1]
        return (x1-x2)**2 + (y1-y2)**2

    @classmethod
    def dist_to_point(self,blocks_prop,other_point):
        '''
        Computes distance between two points for each features.
        :param blocks_prop:
        :param other_point:
        :return:
        '''

        return Metri.distance_sq(blocks_prop["CENTROID"][bid][0],
                                       blocks_prop["CENTROID"][bid][1],
                                       other_point[0],
                                       other_point[1])



    @classmethod
    def get_convexhull(self,polypts):
        poly = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in polypts:
            x, y = float(x), float(y)
            ring.AddPoint(x, y)
        poly.AddGeometry(ring)

        return poly.ConvexHull()

    @classmethod
    def convex_hull_area_ratio(self,distkey,dist_prop):
        return dist_prop[distkey]["AREA"]/float(dist_prop[distkey]["HULL_AREA"])

    @classmethod
    def comp_equal_area_circle(self,area,center,spatialref):
        from osgeo import ogr
        R = math.sqrt(area/math.pi)
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(center[0],center[1])
        bufferDistance = R
        circle = pt.Buffer(bufferDistance)
        #print "%s buffered by %d is %s" % (pt.ExportToWkt(), bufferDistance,'')
        return circle

    @classmethod
    def overlaping_area(self, eqarea_circle, polypts):
        from simple_polygon import util

        home = "D:/workspace/sqdm-repo/sqdm/out/tmp/redist/"

        poly = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in polypts:
            x, y = float(x), float(y)
            ring.AddPoint(x, y)
        poly.AddGeometry(ring)

        I = eqarea_circle.Intersection(poly)
        i = 0
        lA = []
        for p in I:
            #util.poly_ptstoshp(p.GetBoundary().GetPoints(), home + "cong_2011/"+str(i), save_as_multipt=False)
            i +=1
            lA +=[p.GetArea()]

        AI = sum(lA)
        return AI

    @classmethod
    def bctr_to_dctr(self,st_centroid,block_keys,block_prop):
        '''
        Computes distance between two points for each features.
        :param blocks_prop:
        :param other_point:
        :return:
        '''
        for block_id in block_keys:

            d=Metri.distance_sq(st_centroid,block_prop[block_id]["CENTROID"])
            if d == None:
                print type(block_prop[block_id]["CENTROID"])
                print block_prop[block_id]["CENTROID"]
            block_prop[block_id]["BLKCTR_TO_DIST_CTR"] = d

    @classmethod
    def blkctr_to_dist_peri(self,district_poly_pts,building_block_keys,refblock_prop):
        # distances-to-perimeter
        for blockid in building_block_keys:
            len_blk_to_peri = 0

            x2, y2 = refblock_prop[blockid]["CENTROID"]
            for boundar_pt in district_poly_pts[0:]:
                x1, y1 = boundar_pt
                dp = Metri.distance_sq(x1, y1, x2, y2)
                if len_blk_to_peri > dp:
                    len_blk_to_peri = dp

            refblock_prop[blockid]["BLKCTR_TO_DIST_PERI"] = len_blk_to_peri


    @classmethod
    def blkctr_to_dist_peri_np(self,district_poly_pts,building_block_keys,refblock_prop):
        # distances-to-perimeter
        N = 10*len(district_poly_pts)/100 #10 percent of data pts

        dp = copy.deepcopy(district_poly_pts)
        random.shuffle(dp)
        for blockid in building_block_keys:
            len_blk_to_peri = 0
            arb = np.asarray(dp[0:N])
            x2, y2 = refblock_prop[blockid]["CENTROID"]
            arc = np.asarray([(x2,y2)]* arb.shape[0])
            del2 = (arc - arb)**2  #sq diff
            del2 = del2.sum(axis=0) #sum
            len_blk_to_peri = min(del2)
            refblock_prop[blockid]["BLKCTR_TO_DIST_PERI"] = len_blk_to_peri
            del arb
            del arc
            del del2


    @classmethod
    def dist_population(self,building_block_keys,refblock_prop):
        # distances-to-perimeter
        dpop =0
        for blockid in building_block_keys:
            dpop +=refblock_prop[blockid]["POPULATION"]

        return dpop


    @classmethod
    def mean_radius_to_distctr(self,district_blockkeys,block_prop):
        dpda = 0 #prod of blockctr-dist-to-ctr and block's area
        A = 0 #total area
        for blockid in district_blockkeys:
            dpda += block_prop[blockid]["AREA"]\
                    * math.sqrt(block_prop[blockid]["BLKCTR_TO_DIST_CTR"])
            A +=block_prop[blockid]["AREA"]

        return (2*math.sqrt(A/float(PI))/3)/(dpda/float(A))

    @classmethod
    def dynamic_radius_to_distctr(self,district_blockkeys,block_prop):
        dpda = 0 #prod of blockctr-dist-to-ctr and block's area
        A = 0 #total area
        for blockid in district_blockkeys:
            dpda += block_prop[blockid]["AREA"] * block_prop[blockid]["BLKCTR_TO_DIST_CTR"]
            A +=block_prop[blockid]["AREA"]

        return math.sqrt((A/(2 * PI)) / (dpda/A));

    @classmethod
    def harmonic_radius_to_distctr(self,district_blockkeys,block_prop):
        dpda = 0 #prod of blockctr-dist-to-ctr and block's area
        A = 0 #total area
        for blockid in district_blockkeys:
            dpda += block_prop[blockid]["AREA"]/float(math.sqrt(block_prop[blockid]["BLKCTR_TO_DIST_CTR"]))
            A +=block_prop[blockid]["AREA"]

        return math.sqrt(A/PI)/2 / (A/dpda);

    @classmethod
    def rohrbach_index(self,district_blockkeys,block_prop):
        '''
        for each block b, area * min(distance-to-perimeter of district.)
        distance-to-perimeter = DIST(centroid-of-b, point-on-boundary-of-b)
        :param blocks_prop:
        :return:
        '''
        dpda = 0 #prod of dist-to-perimeter and block's area
        A = 0 #total area
        for blockid in district_blockkeys:
            dpda += block_prop[blockid]["AREA"] * block_prop[blockid]["BLKCTR_TO_DIST_PERI"]
            A +=block_prop[blockid]["AREA"]

        try:
            R3 = pow(math.sqrt(A/float(PI)),3) #radius whose area is A
        except Exception, e:
            print(""),str(e)
        return dpda/float(PI * R3 * 3)

    @classmethod
    def exchange_index(self,distkey,dist_prop):
        '''
        A(D Intersection EQAC)/A(EQAC), EQAC=Equal Area circle centered in CENTROID(D)
        :param dist_prop:
        :return:
        '''
        return dist_prop[distkey]["OVERLAP_AREA_EQCIRCLE"]/float(dist_prop[distkey]["AREA"])


    @classmethod
    def interpersonal_distance(self,dist_blockkeys,block_prop):
        tdist_pop = 0
        A = 0
        N = 0
        for this_blockid in dist_blockkeys:
            A += block_prop[this_blockid]["POPULATION"]
            for other_blockid in dist_blockkeys:
                if this_blockid != other_blockid:
                    pp = block_prop[this_blockid]["POPULATION"]* block_prop[other_blockid]["POPULATION"]
                    personal_dist = Metri.distance_sq(block_prop[this_blockid]["CENTROID"],
                                                      block_prop[other_blockid]["CENTROID"])
                    tdist_pop += pp*personal_dist

                    N +=1
        tdist_pop =tdist_pop/float(N) #avg dist
        # divide by 128R/45PI, R is equal area circle.
        return 45*PI*(tdist_pop/(128*math.sqrt(A/float(PI))))

    @classmethod
    def circum_circle(self, polygon):
        cradi, ccenter = 0, 0
        return cradi, ccenter

    @classmethod
    def reock(self,AD,mCSCA):
        '''
        A(D)/A(minCC), minCC=minimum Circum-Circle for D
        :param AD: area of district
        :param mCCA: area of minimum CircumScribing circle.
        :return:
        '''
        return AD/mCSCA

    @classmethod
    def darea_by_ac_eqpd(self, dist_prop):
        '''
        A(district)/A(circle-whose-perimeter-is-equal-to-district)
        :return:
        '''
        for distid in dist_prop['FEAT_IDS'].keys():
            dist_prop[distid]["DAREA_EQPC"]= 4 * PI * dist_prop["AREA"][distid] / float(dist_prop["PERI"][distid] ** 2)

    @classmethod
    def darea_by_ac_eqpd2(self, distkey,dist_prop):
        '''
        A(district)/A(circle-whose-perimeter-is-equal-to-district)
        :return:
        '''
        dist_prop[distkey]["DAREA_EQPC"]= 4 * PI * dist_prop[distkey]["AREA"] / float(dist_prop[distkey]["PERI"] ** 2)


    @classmethod
    def pd_to_peqac(self,dist_prop):
        '''
        Perimeter-of-District/Perimeter-of-Circle-Whose-Area-Eq-District
        Also called: Schwartzberg ratio
        Equal to Inverse-of-(pols-by-poler ratio)
        :param district_poly:
        :return:
        '''
        return {distid: self.inv_root_upolsby_popper(dist_prop["AREA"][distid],
                                          dist_prop["PERI"][distid])
                for distid in dist_prop['FEAT_IDS'].keys()}

    @classmethod
    def population_polygon(self, district_hull,dist_pop,district_blockkeys,blocks_prop):
        '''
        PoP(D)/PoP(district-Hull)
        :param district_hull: Convex Hull (polygon) for this district
        :param district_blockkeys: block ids composing all the districts
        :param district_pop: population of the district.
        :param blocks_prop:
        :return:
        '''
        popHull = 0
        #get total population inside the hull
        for blockid in district_blockkeys:
            center=blocks_prop[blockid]["CENTROID"]
            ctr = ogr.Geometry(ogr.wkbPoint)
            ctr.AddPoint(center[0], center[1])
            if ctr.Within(district_hull):
                popHull += blocks_prop[blockid]["POPULATION"]

        return dist_pop/float(popHull)
