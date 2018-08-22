import numpy as np
from  osgeo import ogr, osr


def workpolygon (polygonRef):
	#feature = layer.GetFeature(fid)
	#polygonRef = feature.GetGeometryRef()
	undict = {}
	hl = []
	for eachpoly in polygonRef: #polygonref can be list of refs.
		pp = eachpoly.GetBoundary().GetPoints()
		n=len(pp)-1#the last point is also the first
		for i in range(n): #make lines from points
			x1=pp[i][0]*1
			x2=pp[i+1][0]*1
			y1=pp[i][1]*1
			y2=pp[i+1][1]*1
			#collect all x-s in a dictionary
			#undict[x1]=1
			#undict[x2]=1
			#arrange lines so first x-coordinate is less that second
			if x1<=x2:
				tt=(x1,y1,x2,y2)
			else:
				tt=(x2,y2,x1,y1)
			if tt not in hl:
				hl.append(tt)

	arr = np.asarray(hl).astype(np.uint32)
	minx, miny = min(arr[:,0].min(), arr[:,2].min()), min(arr[:,1].min(), arr[:,3].min())
	if min(minx,miny) < 0:
		TM = np.asarray([minx,miny,minx,miny])
		arr = arr + TM
	arr = np.round(arr*100,0)
	hl = arr
	#find unique x-coords and split lines
	xun = sorted(list(set((list(arr[:,0])+list(arr[:,2])))))
	iv = [1 for i in range(len(xun))]
	undict = dict(zip(xun,iv))
	ttt=[]# to hold all split lines
	for i in range(len(hl)):#split lines
		x1=hl[i][0]
		y1=hl[i][1]
		x2=hl[i][2]
		y2=hl[i][3]
		li=xun.index(x1)
		hi=xun.index(x2)
		nn=hi-li
		if nn<2 : #difference 0/1 means no split
			ttt.append([x1,y1,x2,y2])
		else:  #chopping necessary
			yc=y1
			for j in range(li,hi):
				xc=xun[j]
				xn=xun[j+1]
				if xn < x2:
					y=float(y2-y1)/(x2-x1)*(xn-x1)+y1
					yn=round(y)
					ttt.append([xc,yc,xn,yn])
					yc=yn
				elif xn==x2: #last one
					ttt.append([xc,yc,xn,y2])
	xdict={}##now add all split lines to a dictionary
	for i in range(len(xun)-1):
		xl=xun[i]
		xh=xun[i+1]
		xdict[xl]=[]
	xdict[xh]=[] #add last ite,
	for i in range(len(ttt)):
		x1=ttt[i][0]
		y1=ttt[i][1]
		x2=ttt[i][2]
		y2=ttt[i][3]
		if y1<y2:
			xdict[x1].append((y1,y2))
		else:
			xdict[x1].append((y2,y1))
	n = len(hl)
	return xdict,n,len(ttt)  #n is original lines, len(ttt) is total splitlines, n=unique ls, len(ttt)=total chops

def ycoordtotypes(y): #convert y coordinates for a key to rectangle types
	y=sorted(y)
	yy=[] #we are going to convert this into a list of chars (for rect type)
	if len(y) == 0:
		pass
	elif len(y)==1:
		if y[0][0]==y[0][1]: #horizontal line
			yy.append('w')
		else:
			yy.append('g')# something like (3,5) -- need to add white like (5,3) for
			              #wrap-around
			yy.append('w')
		pass
	elif len(y)>1:
		for j in range(len(y)-1):
			if (y[j][0]<y[j][1]) and (y[j][1]<y[j+1][0]): #eg. (3,5),(7,?)
				yy.append('g') 
				yy.append('w')        #need to add (5,7) inbetween
			if (y[j][0]<y[j][1]) and (y[j][1]==y[j+1][0]): #eg (3,5),(5,?)
				yy.append('g')
			if (y[j][0]<y[j][1]) and (y[j][1]==y[j+1][1]) and (y[j][1]>y[j+1][0]): #eg. (3,5),(4,5) 
				yy.append('b') #acute angle double line
			if (y[j][0]==y[j][1]) and (y[j][1]<y[j+1][0]):#eg. (3,3),(5,?)
				yy.append('w')
		if len(y)>0:
			k=len(y)-1#index of last tuple
			if y[k][1]==y[k][0]: #...(8,8)
				yy.append('w') #wrapped around
			if (y[k][1]>y[k][0]): 
				yy.append('g') 
				yy.append('w')
			#if (y[k][1]>y[k][0]) and (y[k][0]>= y[k-1][1]):#...(?,6),(8,10)
				#yy.append('g') 
				#yy.append('w') #wrapped around
			#if (y[k][1]>y[k][0]) and (y[k][1]== y[k-1][1]):
				#yy.append('w') #wrapped around after a previous blue
	return yy
			
import DataSources as mydatasource
import sys
reload(mydatasource)
datatype = 'usstate'
polygons, inputsspatialref,feature_names = mydatasource.getPolygons(datatype=datatype,dofilter=True)

print("len of features_names"), len(polygons)
summ=[]
attr = ['data','pid','inputlines','splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes']
print attr
domapall= True #divides a complete map of polygons.
if domapall:
   howmany=1
else:howmany = len(polygons)
datarows = []
for fid in range(0, howmany):
	summ = []
        #feature = layer.GetFeature(fid)
        #polygonRef = feature.GetGeometryRef()
	if domapall:
		polygon = polygons #we are passing list of polys
		fidx = 'map'
	else:
		polygon = [polygons[fid]]
		fidx = fid
	#print fid,polygon.GetBoundary().GetPoints()[0:3]
	xdict,n,tn=workpolygon(polygon)
	#xdict is the dictionary of chopped lines
	d2f={}#dictionary indicating type 'w','g','b' etc
	      #'w' is clear rect. 'g' is for diagonals,'b' acute angles
	cb=0
	cw=0
	cg=0
	for key in xdict:
		y=xdict[key]
		yy=ycoordtotypes(y)
		d2f[key]=yy
		cb=cb+yy.count('b')
		cw=cw+yy.count('w')
		cg=cg+yy.count('g')
	#print "blue, green, white", cb,cg,cw
	summ.append((fid,n,tn,cb,cg,cw))
	summar = np.asarray(summ)
	datarows +=["{0},{1},{2},{3},{4},{5},{6},{7}\n".format(datatype,feature_names[fid],sum(summar[:,1]), sum(summar[:,2]),sum(summar[:,3]),sum(summar[:,4]), sum(summar[:,5]), sys.getsizeof(xdict))]
outfile = './results/' + datatype + '.csv'

with open(outfile,'w') as f:
    for line in datarows:
	f.write(line)
print("data saved in "), outfile


