library(phylotools)
library(ggtree)
library(googlesheets4)
library(phytools)
library(tidyverse)
library(TreeTools)

# Read in our focal taxa:
focalTaxa <- read_sheet("https://docs.google.com/spreadsheets/d/166mvT3i85SJO25s155Dm_SBkzlRPCRHsQDl4Sj9SOLA/edit?usp=sharing")
focalTaxa$speciesUnderscore <- gsub(pattern = " ",
                                    replacement = "_",
                                    focalTaxa$species)

# Read in phylogeny from Kawahara et al. 2023:
# This particular file is the one they use in their HISSE/BISSE analyses:
butterflyTree <- read.tree("./kawahara2023/Data.S21.MLtrees/tree7_AA154partitions_renamed.tre")

# Get the tips in the phylogeny that match our focal species:
filterTips <- function(speciesName) {
  match <- butterflyTree[["tip.label"]] %>% 
    str_subset(pattern = speciesName)
  
  match <- str_flatten(match,
                       collapse = ",")
  
  match <- ifelse(match == "",
                  "noMatch",
                  match)
  
  result <- c(speciesName,
              match)
  
  return(result)
}

matches <- purrr::map(focalTaxa$speciesUnderscore,
                      filterTips)
matches <- as.data.frame(do.call(rbind, matches)) 

# Some species don't have an exact match; get a congener for these. 
speciesWithoutExactMatch <- filter(matches, 
                                   V2 == "noMatch")
speciesWithoutExactMatch$genus <- str_split_i(string = speciesWithoutExactMatch$V1,
                                              pattern = "_",
                                              i = 1)
genusMatches <- purrr::map(speciesWithoutExactMatch$genus,
                           filterTips)
genusMatches <- as.data.frame(do.call(rbind, genusMatches)) 
genusMatches <- full_join(genusMatches,
                          speciesWithoutExactMatch,
                          by = c("V1" = "genus")) %>%
  distinct()

# Two genera both (a) don't have exact species matches and (b) have >1 species. For these, pick a different congener for each species:
genusMatches <- genusMatches %>%
  group_by(V1) %>% 
  mutate(index = row_number())

speciesSelection <- purrr::map2(genusMatches$V2.x, 
                                genusMatches$index, 
                                ~ str_split_i(string = .x, 
                                              pattern = ",", 
                                              i = .y))
speciesSelection <- as.data.frame(do.call(rbind, speciesSelection)) 
genusMatches$species <- speciesSelection$V1
genusMatches <- genusMatches %>% 
  ungroup %>%
  select(c(V1.y,
           species))

# Did these still fail to get a match for any species? (e.g. they have no congener in the tree?)
noMatches <- filter(genusMatches,
                    species == "noMatch")

# We end up with no match for Fabriciana_adippe,
# But per Moya et al. 2017, it is quite closely related to Argynnis_paphia, 
# So I'll pick that for its match. 
genusMatches <- rbind(genusMatches, 
                      setNames(c("Fabriciana_adippe", 
                                 "BN001485_MS2016005_Nymphalidae_Heliconiinae_Argynnini_Argynnis_paphia"), 
                               names(genusMatches)))

# Combine all of those:
colnames(matches) <- c("originalSpecies",
                       "speciesInTree")
colnames(genusMatches) <- c("originalSpecies",
                            "speciesInTree")

allMatches <- rbind(matches,
                    genusMatches) %>%
  filter(speciesInTree != "noMatch") %>%
  group_by(speciesInTree) %>% 
  mutate(index = row_number())

# Prune the Kawahara phylogeny:
prunedTree <- keep.tip(butterflyTree,
                       tip = allMatches$speciesInTree)
ggtree(prunedTree) +
  geom_tiplab()

# Sometimes two species in my focal taxa matched to the same species in the tree, so I need to split those tips:
tipsToSplit <- allMatches[duplicated(allMatches$speciesInTree),]
tipsToSplit <- tipsToSplit$speciesInTree

# Now split the tips:
for (i in tipsToSplit) {
  indexToSplit <- grep(pattern = i,
                       x = prunedTree[["tip.label"]])
  label <- paste(i,
                 "_2",
                 sep = "")
  
  prunedTree <- AddTip(prunedTree,
                       where = indexToSplit,
                       label = label,
                       edgeLength = 0.5 * prunedTree$edge.length[which(prunedTree$edge[,2] == indexToSplit)])
  
} 

ggtree(prunedTree) +
  geom_tiplab() +
  xlim(c(0, 0.5))

# Now replace the names so they match our focal taxa:
# Make a column in the dataframe with names that match those in the pruned tree:
allMatches <- allMatches %>%
  group_by(speciesInTree) %>% 
  mutate(index = paste(row_number()))
allMatches$index <- paste("_",
                          allMatches$index,
                          sep = "")
allMatches$index <- gsub(pattern = "_1",
                         replacement = "",
                         allMatches$index)
allMatches$speciesNameInPrunedTree <- paste(allMatches$speciesInTree,
                                            allMatches$index,
                                            sep = "")
# Do the renaming:
prunedTree <- treeio::rename_taxa(tree = prunedTree,
                            data = allMatches,
                            key = speciesNameInPrunedTree,
                            value = originalSpecies)
ggtree(prunedTree) +
  geom_tiplab() +
  xlim(c(0, 0.5))

# Check that all species are present and correct:
prunedTree[["tip.label"]] %in% focalTaxa$speciesUnderscore

# Save the tree in Newick format:
dir.create("./speciesTree/")
ape::write.tree(prunedTree, 
                file = "./speciesTree/focalSpeciesTree.txt")
