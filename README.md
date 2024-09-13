# AvesTreeCode

This repo stores scripts used to process the outputs of an OpenTree custom synthesis run.


### Step 1:
Download the raw custom synth data and taxonomy crosswalk data  using

```git clone https://github.com/McTavishLab/AvesData```

If you want to run these steps on your own tree, run custom sythnesis with Aves as root. 
You can use any tree collection as long as it has birds in it!

Download the custom synth directory to your computer, and unpack it

e.g. for v1.2
```wget  https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmpl8uyqu1p.tar.gz```

```tar -xzvf https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmpl8uyqu1p.tar.gz```

You should now have a folder with all the synthetic tree artifacts in it. 

This folder OpenTree synth tree folder download is packaged with each tree version in the AvesData repo
https://github.com/McTavishLab/AvesData


e.g. https://github.com/McTavishLab/AvesData/tree/main/Tree_versions/Aves_1.2/OpenTreeSynth
is the unpacked snacktavish_aves_81461_tmpl8uyqu1p.tar.gz


For later steps I will use the paths to data in the AvesData repo as examples.

### Step 2: Takes custom synth output and prunes to taxa in Clements taxonony,  relabels tips to Clements names.

estimates (very rough) dates for internal nodes, and writes iTOL annotation files.

Requires as inputs the mapping from OTT_ids to Clements names available at: 
https://github.com/McTavishLab/OpenTreeCLO/blob/main/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv

And the custom synth directory, untarred.
To perfom custom synth see:
https://opentreeoflife.github.io/CustomSynthesis/

Example:
process_custom_synth_birds.py ../AvesData/Tree_versions/Aves_1.2/OpenTreeSynth ../AvesData/Taxonomy_versions/Clements2023/OTT_crosswalk_2023.csv output_example


Output_example will contain:
* Treefiles:  
    - pruned.tre -> the synth tree pruned to taxa found in the taxonomy crosswalk  
    - phylo_only.tre   -> the tree pruned to taxa found in the taxonomy corsswalk AND in phylogenetic inputs  
* Dates information in subdirectory dates/: (rough ultrametricization using only Kimball et al 2019 and Oliveros et al. 2019)  
    - dates/custom_node_ages.json  -> the date input data. These dates are based only on   
    - dates/phylo_only_select_dates_mean_clements_labels.tre  -> the phylogeny only tree, dated, with Clements taxonomy labels  
    - dates/phylo_only_select_dates_mean_ott_labels.tre  -> the phylogeny only tree, dated, with OpenTree taxonomy labels  
    - dates/select_dates_citations.txt <- citations for the studies used to date  
    - dates/dates_select_phylo_only  <- directory of intermediate dating files  
* Metdata:  
    - citation_node_counts.tsv  <- Table of citations and how many nodes are supported by each  
    - tips_without_phylo.txt <- taxa in crosswalk with no phylogenetic information  
    - studies_per_tip.txt <- Mapping of what input studies inform placement of each tip  
* Annotation files for viewing on itol:  
    - ottlabel.txt <- Annotation file to re-label tips from ott ids to Clements labels  
    - jetz_support.txt   <- agreement with Jetz et al. 2012  
    - jetz_conflict.txt  <- disagreement with Jetz et al. 2012  
    - conflict_12.txt   <- Conflict scaled to max color at 12       
    - conflict_3.txt   <- Conflict scaled to max color at 3              
    - support_10.txt <- Support scaled to max color at 10   
    - support_20.txt <- Support scaled to max color at 20   


### Step 3: Add taxa not in phylogenies


### Step 4: Dates