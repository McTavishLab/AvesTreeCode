import dendropy
import sys
import os
import csv
import copy
import json
import random
from opentree import OT, annotations
from helpers import crosswalk_to_dict


support_annot = sys.argv[1] 
outfile = open(sys.argv[2],'w')

support = json.load(open(support_annot))

template ="""DATASET_PIECHART
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL,conflict and support

#dataset color (can be changed later)
COLOR,#ff0000

#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)
FIELD_COLORS,#00FF00,#82A67D,#FF0000


#field labels
FIELD_LABELS,strict support,traversed by,conflict

#maximum pie chart radius will be displayed with this size, others will be proportionally smaller. This setting applies to internal pie charts only.
MAXIMUM_SIZE,50

#pie chart height factor; For external pie charts, default radius will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
#HEIGHT_FACTOR,1

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the pie chart segments 
#BORDER_WIDTH,0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR,#0000ff

#display data as polar area diagrams (equal arcs with varying radius), instead of standard pie charts (varying arcs of equal radius)
#POLAR_AREA_DIAGRAM,0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
#the following fields are required for each node:
#ID,position,radius,value1,value2,value3...
#position defines the position of the pie chart on the tree:
#  -1 = external pie chart
#  a number between 0 and 1 = internal pie chart positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)

DATA
"""

outfile.write(template)

studies = 0
for node in support['nodes']:
    ss = len(support['nodes'][node].get("supported_by",[]))
    con = len(support['nodes'][node].get("conflicts_with",[]))
    trav = len(support['nodes'][node].get("partial_path_of",[]))
    total = ss+con+trav
    if con>ss:
        outfile.write("".join([node,",0,",str(total)+",",str(ss)+",", str(trav)+",",str(con)+"\n"]))

