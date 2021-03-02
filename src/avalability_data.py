#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import sys
import re

#===================
data=sys.argv[1]
west=float(sys.argv[2])
south=float(sys.argv[3])

east=west+10.0
north=south+10.0
#--
flag=0
#--
fname="./inp/"+data+"Station_list.txt"
with open(fname) as f:
    lines = f.readlines()
    for line in lines[1::]:
        line = filter(None, re.split(" ",line))
        lon  = float(line[5])
        lat  = float(line[6])
        #print lon, lat
        if lon >= west and lon <= east and lat >= south and lat <= north:
            flag=1

#----
print (flag)