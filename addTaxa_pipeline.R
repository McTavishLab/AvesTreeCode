## Developed by Eliot Miller
#start your libraries
install.packages("devtools")
install.packages("RRphylo")
library(devtools)
# AddTaxa dependecy Laser was recently removed from cran and needs to be installed from github
install_github("cran/laser")
install_github("eliotmiller/addTaxa")
library(addTaxa)
library(RRphylo)

#load the tree from the Aves Data Repo

setwd("AvesTreeCode")
myTree <- read.tree("../AvesData/Tree_versions/Aves_1.3/Clements2021/phylo_only_clements_labels.tre")

#load the grouping file and taxonomy
groupings <- read.csv("../AvesData/Taxonomy_versions/Clements2021/TaxonomytaxonAddition_2021taxonomy_v1_4.csv")
tax <- read.csv("../AvesData/Taxonomy_versions/Clements2021/eBird_Taxonomy_v2021.csv") 


#force the trees to be ultrametric
#this comes from a jonathan chang's blog
#https://www.r-bloggers.com/2021/07/three-ways-to-check-and-fix-ultrametric-phylogenies/

#convenience variables
N <- Ntip(myTree)
myTree$edge.length <- rep(.1, N) # add edge length
myTree <- reorder(myTree, "postorder")
e1 <- myTree$edge[, 1] # parent node
e2 <- myTree$edge[, 2] # child node
EL <- myTree$edge.length
ages <- numeric(N + myTree$Nnode)

for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) {
    ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
  } else {
    recorded_age <- ages[e1[ii]]
    new_age <- ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) {
      cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
      EL[ii] <- recorded_age - ages[e2[ii]]
    }
  }
}

ultraTree <- myTree
ultraTree$edge.length <- EL
is.ultrametric(ultraTree)
is.rooted(ultraTree)
is.binary(ultraTree)


#see whether every species that is missing has an addition statement
tax <- tax[tax$CATEGORY == "species",]
tax$underscores <- sub(" ", "_", tax$SCI_NAME)
missing <- setdiff(tax$underscores, ultraTree$tip.label)
hasStatement <- groupings$species[groupings$add.to != ""]
setdiff(missing, hasStatement) #yes

#now randomly resolve the polytomies 100 times, and generate 100 trees primed
#for taxon addition
resolved <- list()

for(i in 1:100)
{
  resolved[[i]] <- multi2di(ultraTree, random=TRUE)
}

#bumping into an issue where taxa in add.to and do.not.break statements no longer exist.
#let's try and identify issues and fix them. families screw things up here with spaces,
#but this still mostly works. can work through these manually
allAdd <- unique(unlist(strsplit(groupings$add.to, ", |; ")))
allDoNotBreak <- unique(unlist(strsplit(groupings$do.not.break, ", |; ")))
allValid <- unique(c(groupings$family, groupings$genus, groupings$species))
setdiff(allAdd, allValid)
setdiff(allDoNotBreak, allValid)

#loop over missing and try addition tests with each and confirm there is a valid
#position for every one (i assume it's possible to create situations where after
#one addition it becomes impossible to add another, but haven't encountered yet)
randomized <- sample(missing)
for(i in 1:length(randomized))
{
  print(randomized[i])
  additionPrep(resolved[[1]], groupings, randomized[i])
}

#if you want to find polytomies
desc <- plyr::count(ultraTree$edge[,1])
poly <- desc[desc$freq > 2,]

#create the actual complete trees. you previously threw an error in additionPrep
#if the tree was not binary. commented that check out. going for it here, not sure
#what outcome will be, but appears to run. takes 433s for one tree, so should
#take 12h for 100 trees
completeSet <- list()

for(i in 1:length(resolved))
{
  completeSet[[i]] <- customAdd(resolved[[i]], addition.statements=groupings, no.trees=1)
}

#save these out
toSave <- list()

for(i in 1:length(completeSet))
{
  toSave[[i]] <- completeSet[[i]]$trees[[1]]
}

class(toSave) <- "multiPhylo"

#somehow we get polytomies again. resolve again
toSave <- lapply(toSave, multi2di, random=TRUE)

#write.tree(toSave, "data/v1_3_tree_v1_4_add_2021_taxonomy_resolved.tre")

#try plotting a tree for every family where species you added taxonomically are colored
#in red
families <- unique(tax$FAMILY)

#drop NAs
families <- families[!is.na(families)]

for(i in 1:length(families))
{
  spp <- tax$underscores[tax$FAMILY==families[i]]
  
  if(length(spp) <= 1)
  {
    next()
  }
  
  subtree <- ladderize(drop.tip(toSave[[1]], setdiff(toSave[[1]]$tip.label, spp)))
  
  noTaxa <- length(subtree$tip.label)
  
  #figure out which taxa were added taxonomically  
  tipColors <- rep("black", length(subtree$tip.label))
  tipColors[subtree$tip.label %in% missing] <- "red"
  
  if(noTaxa > 1 & noTaxa < 10)
  {
    textSize <- 1
  }
  else if(noTaxa >= 10 & noTaxa < 50)
  {
    textSize <- 0.7
  }
  else if(noTaxa >= 50 & noTaxa < 100)
  {
    textSize <- 0.5
  }
  else if(noTaxa >= 100 & noTaxa < 150)
  {
    textSize <- 0.4
  }
  else if(noTaxa >= 150 & noTaxa < 300)
  {
    textSize <- 0.2
  }
  else if(noTaxa >= 300 & noTaxa < 500)
  {
    textSize <- 0.1
  }
  else if(noTaxa >= 500 & noTaxa < 1000)
  {
    textSize <- 0.8
  }
  else if(noTaxa >= 1000 & noTaxa < 3000)
  {
    textSize <- 0.05
  }
  else
  {
    textSize <- 0.02
  }
  
  pdf(file=paste("taxonAdditionPlots/",families[i], ".pdf", sep=""), height=10, width=5)
  plot(subtree, cex=textSize, tip.color=tipColors)
  dev.off()
}
