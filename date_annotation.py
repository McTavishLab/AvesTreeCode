import dendropy
import sys
import os
import csv
import copy
import json
import random
from opentree import OT, annotations
from helpers import crosswalk_to_dict


node_ages_json = sys.argv[1] 
outfile = open(sys.argv[2],'w')

ages = json.load(open(node_ages_json))

template ="""DATASET_SYMBOL
# See https://itol.embl.de/help.cgi#symbols
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL,number of dates

#dataset color (can be changed later)
COLOR,#ffff00


#largest symbol will be displayed with this size, others will be proportionally smaller.
MAXIMUM_SIZE,50

DATA
"""

outfile.write(template)

max_count = 0
for node in ages['node_ages']:
    if len(ages['node_ages'][node]) > max_count:
        max_count = len(ages['node_ages'][node])


for node in ages['node_ages']:
    count = len(ages['node_ages'][node])
    outfile.write("".join([node,",2,",str((count/max_count)*10), ",#0000ff,1,0\n"]))
    #100379,3,5,#0000ff,0,0

