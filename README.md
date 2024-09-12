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

??estimates dates for internal nodes, 
??and writes iTOL annotation files.

Requires as inputs the mapping from OTT_ids to Clements names available at: 
https://github.com/McTavishLab/OpenTreeCLO/blob/main/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv

And the custom synth directory, untarred.
To perfom custom synth see:
https://opentreeoflife.github.io/CustomSynthesis/

Example:
process_custom_synth_birds.py custom_synth_runs/snacktavish_aves_81461_tmpw3m8cs_b taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv 

