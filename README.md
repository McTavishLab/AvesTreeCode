# AvesTreeCode

This repo stores scripts used to process the outputs of an OpenTree custom synthesis run.

Steps to reproduce the results found in [McTavish et al. 2024](https://www.biorxiv.org/content/10.1101/2024.05.20.595017v1) are shown here.


### Step 1:
Download the raw custom synth data and taxonomy crosswalk data using

```git clone https://github.com/McTavishLab/AvesData```

If you want to run these steps on your own tree, run custom sythnesis with Aves as root. 
You can use any tree collection as long as it has birds in it!

Download the custom synth directory to your computer, and unpack it

e.g. for v1.2
```wget  https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmpl8uyqu1p.tar.gz```

```tar -xzvf https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmpl8uyqu1p.tar.gz```

You should now have a folder with all the synthetic tree artifacts in it. 

This OpenTree synth folder is packaged with each tree version in the AvesData repo
https://github.com/McTavishLab/AvesData


e.g. https://github.com/McTavishLab/AvesData/tree/main/Tree_versions/Aves_1.2/OpenTreeSynth
is the unpacked snacktavish_aves_81461_tmpl8uyqu1p.tar.gz


For later steps I will use the paths to data in the AvesData repo as examples.

### Python requirements
Python requirements are listed in requirements.txt

I recommend using a python virtual environment to avoid any conflicts, and installing the requrements there.

    $ virtualenv -p python3 venv-aves
    $ source venv-aves/bin/activate
    $ pip install -r requirements.txt

### Step 2: Takes custom synth output and prunes to taxa in Clements taxonony, relabels tips to Clements names.

Estimates rough dates for internal nodes, and writes iTOL annotation files.

Requires as inputs the mapping from OTT_ids to Clements names available at: 
https://github.com/McTavishLab/AvesData/blob/main/Taxonomy_versions/Clements2021/OTT_crosswalk_2021.csv

And the custom synth directory, untarred.
To perfom custom synth see:
https://opentreeoflife.github.io/CustomSynthesis/

```python process_custom_synth_birds.py \<Open Tree Synth directory\> \<Taxonomic crosswalk\> \<Name of output\>```


Example:
```python process_custom_synth_birds.py ../AvesData/Tree_versions/Aves_1.2/OpenTreeSynth ../AvesData/Taxonomy_versions/Clements2023/OTT_crosswalk_2023.csv output_example```


Output_example will contain:
* Treefiles:  
    - pruned.tre -> the synth tree pruned to taxa found in the taxonomy crosswalk  
    - phylo_only.tre   -> the tree pruned to taxa found in the taxonomy crosswalk AND in phylogenetic inputs 
    - phylo_only_clements_labels.tre -> phylo only tree with Clements taxon labels 
    - phylo_only_ott_labels_ultrametric.tre -> phylo only tree scaled to time
    - phylo_only_clements_labels_ultrametric.tre -> phylo only tree scaled to time with Clements taxon labels
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


run addTaxa_pipeline.R, in R

Takes as input:
* an ultrametric tree with Clements taxonomy labels, (e.g. phylo_only_clements_labels_ultrametric.tre)
* the eBird taxonomy for the year (e.g. AvesData/Taxonomy_versions/Clements2021/eBird_Taxonomy_v2021.csv)
* a set of taxon addition statements (e.g. AvesData/Taxonomy_versions/Clements2021/taxonAddition_2021taxonomy_v1_4.csv)


This script randomly resolves polytomies 100 times, and adds any missing taxa to each of those trees.
As the polytomy resolution and the taxon addition process are both stochastic, this generates a distribution of trees.

Outputs:
* taxon_addition_treeset.tre <- a file containing 100 complete trees
* taxonAdditionPlots.pdf <- a figure with for each family where species you added taxonomically are colored in red

### Step 4: Dates
Date estimates are translated from the dates mapped to the phylogeny only custom synthesis input tree onto each of the complete trees generated by the taxon addition algorithim.

The branch smoothing algorithim, BLADJ (Webb et al. 2008) is applied to the whole tree.
This is a simple algorithim which places undated nodes evernly between nodes with date estimates.

Takes as input:
* the phylogeny only tree from step 1.
* the taxon addition treeset from step 3.
* the taxonomy crosswalk table from the AvesData repo.

This dating aproach attempts to account for uncertainty in two ways -
- We date each of the trees generated by 100 iterations of the stochastic taxon addition process
- For each tree we estimated a full dated tree in two ways - 
    * node mean date calibration: where the avarage age of evey dated node is used to calibrate
    * random node date sampling: For each node with dates we randomly select one of the available dates to cailbrate. If there is only one date for that node, it is used with probability 0.5.

More information on Chronosynth at the [wiki](https://github.com/OpenTreeOfLife/chronosynth/wiki/Chronosynth-methods-overview)


```python date_complete_treeset.py <Complete treeset file from step 3> <Labelled phylo only tree from step 2>  <Taxonomic cross walk from AvesData>  <output directory>```

```python date_complete_treeset.py ../AvesData/Tree_versions/Aves_1.2/Clements2021/taxon_addition_treeset.tre ../AvesData/Tree_versions/Aves_1.2/Clements2021/phylo_only.tre  ../AvesData/Taxonomy_versions/Clements2021/OTT_crosswalk_2021.csv  dated_treeset_2021

```

#### Outputs:
The output directory with contain:

    * full_dates_citations.txt <- a file containing the citations for all the studies used in dating
    * dates_add_taxa/ <- a folder
        - all_nodes.json <- a json file containing all the node dates for the tree  
        - dated_all_mean_dates_clements\<i\>.tre <- Dated tree using mean node age for each dated node. Topologies are numbered 1-100.  Labels are Clements names. (a set with ott labels is created as well)  
        - dated_all_rand_dates_clements\<i\>.tre <- Dated tree using a randomly sampled node age for each dated node. Toplogies are numbered 1-100. Labels are Clements names. (a set with ott labels is created as well)  
        - A folder for each run, containing the files to run bladj, and the bladj output.  


I summarized these trees using RevBayes:


```
$ echo "trees" > dated_treeset_2021/dated_rand_sample_clements.tre ##Rev Bayes expects a header line, or you lose your first tree
$ cat dated_treeset_2021/dates_add_taxa/dated_rand_all_dates_clements_labels_tree*.tre >> dated_treeset_2021/dated_rand_sample_clements.tre
$ rb
> tt = readTreeTrace("dated_rand_sample_clements.tre", "clock", burnin=0)
> mcc_tree = mccTree(trace=tt, file="mcc_dated_clements.nex", positiveBranchLengths=TRUE)
```

This summary tree does not include the internal node labels from the phylogeny that we use to track citations.

We re-apply internal node labels using a script that takes them from the phylogeny tree and applies them to the summary tree:

```
python relabel_MCC_nodes.py mcc_dated_clements.nex phylo_only.tre ../../../Taxonomy_versions/Clements2021/OTT_crosswalk_2021.csv  summary_dated_clements.nex
```