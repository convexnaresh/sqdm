Required python packages:
import os
from  osgeo import ogr, osr
import pickle
from shapely.geometry import LinearRing,Point,mapping
import numpy as np
import time
import ntpath
import math


1) Merge counties to construct state's boundary. Save the states.shp file
2) Merge states to construct country boundary. Save the usa.shp file
3) Take states.shp or usa.shp file to :
    --2) transform the co-ordinate system to projected coordinate system
    --1) transform the projected co-ordinate system to universal co-ordinate system which is worldmercater.
    --3) convert the projected co-ordinate to an integer value.
    --4) now label the segments in the boundary by labels for top/bottom region codes
    --5) save the file as *.prj.lbl.txn.int000.shp; where *=states or usa for states or usa labeled boundary segments
    --6) Plot the *.int000.shp file to see if the boundary shape is preserved. must be similar to the original .shp file.
4) 