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


custom_synth_dir = os.path.abspath(sys.argv[1])
taxonomy_crosswalk = sys.argv[2] 
input_trees_file = sys.argv[3] 
filename = sys.argv[4]

## Generates a dictionary to get Clements names from OTT_ids
clements_name_map = crosswalk_to_dict(taxonomy_crosswalk)
ott_id_map = {v: k for k, v in clements_name_map.items()}


custom_synth = dendropy.TreeList.get(path=input_trees_file, schema = 'newick')
## Estimate dates


#label as ott_ids where possible
for tax in custom_synth.taxon_namespace:
    tax.label = ott_id_map.get(tax.label, tax.label)



dates_dir = "{}/dates_add_taxa".format(custom_synth_dir)
if not os.path.exists(dates_dir):
    os.mkdir(dates_dir)


all_chronograms = chronogram.find_trees(search_property="ot:branchLengthMode", value="ot:time")
all_birds = chronogram.find_trees(search_property="ot:ottId", value="81461")
bird_chrono = set(all_chronograms).intersection(set(all_birds))
problematic_chrono = set(["ot_2008@tree5", "pg_2829@tree6579"])
bird_chrono = bird_chrono - problematic_chrono



base_tree = dendropy.Tree.get_from_path(custom_synth_dir+"/phylo_only.tre", schema="Newick")
custom_str = chronogram.conflict_tree_str(base_tree)
all_dates_file = "{}/all_node_ages.json".format(dates_dir)
#if os.path.exists(all_dates_file):
#    all_dates = json.load(open(all_dates_file))
#else:

internal_label_map = {}

base_tree_leaves = set(tip.taxon.label for tip in base_tree.leaf_node_iter())

for node in base_tree:
    internal_label_map[node.label] = set(tip.taxon.label for tip in node.leaf_iter())

if os.path.exists(all_dates_file):
    all_dates = json.load(open(all_dates_file))
else:
    all_dates = chronogram.combine_ages_from_sources(list(bird_chrono),
                                                     json_out = all_dates_file,
                                                     compare_to = base_tree)


tree_iter = 0

for tree in custom_synth:
    internal_label_map_new ={}
    tree_iter+=1
    ##relabel
    leaves = [tip.taxon.label for tip in tree.leaf_node_iter()]
    root_node = "ott81461"
    tree.seed_node.label = root_node
#     node_iter = 0
#     for node in tree:
#         if node.is_leaf() == False:
#             node_iter += 1
#             if node.taxon:
#                 node_label = node.taxon.label
#             else:
#                 node_label = node.label
#             if node_label == "NA" or None:
#                 phylo_tips = base_tree_leaves.intersection(set(tip.taxon.label for tip in node.leaf_iter()))
#                 if len(phylo_tips) < 2:
#                     node.label = "node{}".format(node_iter)
#                 else:
#                     mrca = base_tree.mrca(taxon_labels=phylo_tips)
#                     node.label = mrca.label
#             assert base_tree_leaves.intersection(set(tip.taxon.label for tip in node.leaf_iter())) == internal_label_map.get(node.label, set())

#             internal_label_map_new[node.label] = set(tip.taxon.label for tip in node.leaf_iter())
#     assert set(internal_label_map.keys()).issubset(set(internal_label_map_new.keys()))
    max_age_est = 130
    print("dating full tree, all dates, mean")
    treesfile, sources = chronogram.date_tree(tree,
                                                all_dates,
                                                root_node,
                                                max_age_est,
                                                method='bladj',
                                                output_dir="{}/dates_all_mean_{}_{}".format(dates_dir, filename, tree_iter),
                                                select = "mean",
                                                reps = 1
                                                )
    dated_phylo = dendropy.Tree.get_from_path(treesfile, schema = "newick")
    dated_phylo.write(path="{}/dated_mean_all_dates_ott_labels_{}_tree{}.tre".format(dates_dir, filename, tree_iter), schema="newick")
    for tax in dated_phylo.taxon_namespace:
        tax.label = clements_name_map.get(tax.label, tax.label)
    dated_phylo.write(path="{}/dated_mean_all_dates_clements_labels_{}_tree{}.tre".format(dates_dir, filename, tree_iter), schema="newick")


#--------------------Generating citations -------------------------------------


node_date_count = 0
matched_date_studies = set()
for node in all_dates['node_ages']:
    node_date_count += 1
    for source in all_dates['node_ages'][node]:
        matched_date_studies.add(source['source_id'].split('@')[0])


dates_cite_file = open("{}/dates/full_dates_citations.txt".format(custom_synth_dir), "w")
cites = OT.get_citations(matched_date_studies)
dates_cite_file.write(cites)
dates_cite_file.close()


# #--------------------------------------------------------------------------------------------------

print("""date information for {ld} nodes in the tree
          was summarized from {lds} published studies""".format(ld=node_date_count,
                                                                 lds=len(matched_date_studies)))

for tree in custom_synth:
    tree_iter+=1
    ##relabel
    
    leaves = [tip.taxon.label for tip in tree.leaf_node_iter()]
    root_node = "ott81461"
    max_age_est = 130
    print("dating full tree, all dates, rand")
    for i in range(10):
        treesfile, sources = chronogram.date_tree(tree,
                                          all_dates,
                                          root_node,
                                          max_age_est,
                                          method='bladj',
                                          output_dir="{}/dates_all_rand_sample{}_{}_{}".format(dates_dir, filename, tree_iter, i),
                                          select = "random",
                                          reps = 1
                                          )
        dated_phylo = dendropy.Tree.get_from_path(treesfile, schema = "newick")
#        dated_phylo.write(path="{}/dated_mean_all_dates_ott_labels_{}_tree{}.tre".format(dates_dir, filename, tree_iter), schema="newick")
        for tax in dated_phylo.taxon_namespace:
            tax.label = clements_name_map.get(tax.label, tax.label)
        dated_phylo.write(path="{}/dated_rand_all_dates_clements_labels_{}_tree{}_{}.tre".format(dates_dir, filename, tree_iter, i), schema="newick")



