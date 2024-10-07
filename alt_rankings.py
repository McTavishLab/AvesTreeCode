import dendropy

jetz = "jetz_output/phylo_only_ott_labels_ultrametric.tre"
wu = "wu_output/phylo_only_ott_labels_ultrametric.tre"
prum = "prum_output/phylo_only_ott_labels_ultrametric.tre"
stiller = "output_example_2021/dates/phylo_only_ott_labels_ultrametric.tre"


aves_tree = dendropy.Tree.get(path=stiller, schema="newick")
jetz_tree = dendropy.Tree.get(path=jetz, schema="newick", taxon_namespace = aves_tree.taxon_namespace)
wu_tree = dendropy.Tree.get(path=wu, schema="newick", taxon_namespace = aves_tree.taxon_namespace)
prum_tree = dendropy.Tree.get(path=prum, schema="newick", taxon_namespace = aves_tree.taxon_namespace)


dendropy.calculate.treecompare.symmetric_difference(aves_tree, prum_tree, is_bipartitions_updated=False)

dendropy.calculate.treecompare.symmetric_difference(aves_tree, wu_tree, is_bipartitions_updated=False)

dendropy.calculate.treecompare.symmetric_difference(aves_tree, jetz_tree, is_bipartitions_updated=False)