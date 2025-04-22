library(ggimage)
library(ggtree)

readTree <- function(file) {
  # Read in a single tree and get the tip labels of that tree:
  text <- readr::read_file(file)
  
  # In case hyphy has already run on this, add back in the semicolon that ape requires:
  text <- gsub("\n",
               "",
               text,
               perl = TRUE)
  text <- gsub("$(?<!;)",
               ";",
               text,
               perl = TRUE)
  tree <- ggtree::read.tree(text = text)
  return(tree)
}

imageDirectory <- "./speciesImages"

testTree <- readTree(file = "./06_labelledPhylogenies/tails/labelled_HOG0011193_tree.txt")

ggtree(testTree) + 
  geom_tiplab(aes(image = paste0(imageDirectory, 
                               '/', 
                               str_split_i(label,
                                           pattern = "_",
                                           i = 1), 
                               '.jpg')), 
              geom = "image", 
              offset = 0.05, 
              align = 2, 
              size = 0.08) + 
  geom_tiplab(geom = 'label',
              offset = 0.1) + 
  xlim(0, 1) 





