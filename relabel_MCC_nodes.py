#!/usr/bin/env python
import sys
import os
import csv
import json
import dendropy
import subprocess
from opentree import OT
from chronosynth import chronogram
from helpers import crosswalk_to_dict

"""
Pass in set of trees with CLO labels
"""


input_tree_file = sys.argv[1] 
labelled_tree_file = sys.argv[2] ## phylo only
taxonomy_crosswalk = sys.argv[3]
outputfile = sys.argv[4]


## Generates a dictionary to get Clements names from OTT_ids
clements_name_map = crosswalk_to_dict(taxonomy_crosswalk)
ott_id_map = {v: k for k, v in clements_name_map.items()}



mcc = dendropy.Tree.get(path=input_tree_file, schema = 'Nexus')
mcc.is_rooted = True

labelled_tree = dendropy.Tree.get(path=labelled_tree_file, schema = 'Newick')

for tax in mcc.taxon_namespace:
     tax.label = ott_id_map.get(tax.label, tax.label)

#

label_map = {}
for node in labelled_tree:
    lab = None
    if node.label:
        lab = node.label
    if node.taxon:
        lab = node.taxon.label
    label_map[lab] = [leaf.taxon.label for leaf in node.leaf_iter()]


for lab in label_map:
    mccnode = mcc.mrca(taxon_labels = label_map[lab])
    mccnode.label = lab

for node in mcc:
    if node.label:
        if node.label.startswith("ott") or node.label.startswith("mrca"):
            pass
        else:
            node.label = "NA"
    else:
        node.label = "NA"

mcc.write(path=outputfile, schema="nexus")
