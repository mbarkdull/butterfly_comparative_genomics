library(rBLAST)
library(tidyverse)
library(googlesheets4)
library(rentrez)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                        "/usr/local/ncbi/blast/bin/", 
                        sep = .Platform$path.sep))

### Download the candidate genes from NCBI based on my google sheet: ####
# Read in the information from my public google sheet, getting NCBI gene IDs:
candidates <- read_sheet("https://docs.google.com/spreadsheets/d/1D_cnkNskGndVqMIfkKU3dpjHX3_6OD4Cv78YIiRiGco/edit?usp=drive_link")
candidates$`NCBI locus ID` <- as.character(candidates$`NCBI locus ID`)
ncbiIDs <- candidates$`NCBI locus ID`

# Get the list of already-downloaded proteins:
geneFiles <- dir("candidate_genes", 
                 recursive = TRUE, 
                 full.names = TRUE)

# List the ones that are empty:
emptyGeneFiles <- geneFiles[file.info(geneFiles)[["size"]]==0]

# Remove empty files
unlink(emptyGeneFiles, 
       recursive = TRUE, 
       force = FALSE)

# List the downloaded genes:
alreadyDownloaded <- list.files(path = "./candidate_genes/",
                                pattern = "*fasta") %>% 
  str_split_i(pattern = "\\.",
              i = 1)

# Get a list of only proteins that still need to be downloaded:
ncbiIDs <- base::setdiff(ncbiIDs,
                         alreadyDownloaded)

# Create an output directory to which the files can be downloaded:
dir.create("./candidate_genes/")

# Write a function to download a single protein sequence:
gettingProteinSequence <- function(ID) {
  # Link the gene ID to a protein ID:
  geneToProtein <- entrez_link(dbfrom = 'gene', 
                               id = ID, 
                               db = 'protein')
  proteinID <- geneToProtein[["links"]][["gene_protein_refseq"]][[1]]
  # Get that protein sequence as a fasta:
  entrez_fetch(db = "protein", 
               id = proteinID, 
               rettype = "fasta") %>%
    write(file = paste("./candidate_genes/",
                       ID,
                       ".fasta",
                       sep = ""))
}

# Make a safe version with possibly:
possiblyGettingProteinSequence <- purrr::possibly(gettingProteinSequence,
                                                  otherwise = "Error")
# Download the list of protein sequences:
purrr::map(ncbiIDs,
           possiblyGettingProteinSequence)

#### Combine all the proteomes against which we'll BLAST ####
# Use an if-else statement so this step doesn't have to be repeated:
if (file.exists("allProteomes.fasta")) {
  print("The file exists!")
} else {
  print("The file does not exist.")
  # Construct a list of the proteomes:
  allProteomes <- list.files("03_OrthoFinder/fasta",
                             full.names = TRUE) %>% 
    str_subset(pattern = "_aa.fasta") %>% 
    str_subset(pattern = "bomo",
               negate = TRUE)
  
  # Read in all the proteomes and combine them:
  proteomes <- purrr::map(allProteomes,
                          phylotools::read.fasta)
  proteomes <- as.data.frame(do.call(rbind, proteomes)) 
  phylotools::dat2fasta(proteomes,
                        outfile = "allProteomes.fasta")
  rm(proteomes)
}

#### BLAST the candidates against my proteomes ####
# Construct the database:
makeblastdb(file = './allProteomes.fasta', dbtype = "prot")

# Process the orthogroups so they can be identified:
# Read in the orthogroups file:
orthogroupMembers <- read_delim(file = "./03_OrthoFinder/fasta/OrthoFinder/Results_Apr15/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", 
                                col_names = TRUE,
                                delim = "\t") %>%
  dplyr::select(-c("OG",
                   "Gene Tree Parent Clade")) 

orthogroupMembers <- orthogroupMembers %>%
  tidyr::pivot_longer(cols = -c(HOG),
                      names_to = "junk") %>%
  dplyr::select(-junk) %>%
  dplyr::distinct() %>%
  filter(!is.na(value))

orthogroupMembers$HOG <- str_split_i(orthogroupMembers$HOG,
                                     pattern = "\\.",
                                     i = 2)

# Write a function to blast candidate genes against my genomes and identify the matching orthogroup:
fetchingOrthogroups <- function(i) {
  # Read in my query sequence:
  querySequence <- readAAStringSet(i, 
                                   format = 'fasta')
  
  # Prep the search:
  blastSearch <- blast(db = './allProteomes.fasta', 
                       type = 'blastp')
  
  # Run the search:
  searchResults <- predict(blastSearch, 
                           querySequence)
  
  # Make a column that is the product of length and identity for selecting the best match:
  searchResults$bestMatch <- searchResults$pident*searchResults$length
  
  # Filter the search to include only matches with 100% identity:
  searchResults <- dplyr::filter(searchResults, 
                                 bestMatch == max(bestMatch))
  searchResults <- head(searchResults, 
                        n = 1)
  
  # Extract the percent identity of the match:
  percentID <- max(searchResults$pident)
  
  # Get out the matching sequence name:
  matchingSequence <- searchResults$sseqid
  
  # Get the orthogroup that contains the sequence from our BLAST:
  orthogroup <- filter(orthogroupMembers, 
                       grepl(pattern = matchingSequence,
                             x = value))
  if (length(orthogroup$HOG) == 1) {
    orthogroup <- orthogroup$HOG
    print(paste("Found a match:",
                matchingSequence))
  } else {
    orthogroup <- "noMatch"
  }
  
  #Construct an object that contains that orthogroup, and the file it matches to, so you know which gene of interest goes with each orthogroup:
  results <- c(orthogroup, i, percentID)
  return(results)
}

# Make a safe version with possibly:
possiblyFetchingOrthogroups <- possibly(fetchingOrthogroups,
                                        otherwise = "Error")

# List all of the genes to potentially BLAST for:
allGenes <- list.files("./candidate_genes",
                       full.names = TRUE)

# Read in any past results, so as not to duplicate analyses:
if(file.exists("./candidate_genes/blastSummary.csv") == TRUE){
  pastResults <- read_delim(file = "./candidate_genes/blastSummary.csv")
  pastResults <- pastResults$V2
  
  # Get a list of genes to search:
  genesToSearch <- base::setdiff(allGenes, pastResults)
  
} else {
  genesToSearch <- allGenes
}

# BLAST for the new candidates:
allGenesMatches <- purrr::map(genesToSearch,
                              possiblyFetchingOrthogroups)

allGenesMatches <- as.data.frame(do.call(rbind, allGenesMatches)) 

copy <- allGenesMatches

allGenesMatches <- allGenesMatches %>%
  filter(V1 != "noMatch") %>%
  filter(V1 != "Error")

# Combine new results with old results, if any exist:
if(file.exists("./candidate_genes/blastSummary.csv") == TRUE) {
  pastResults <- read_delim(file = "./candidate_genes/blastSummary.csv") %>%
    distinct()
  pastResults$V3 <- as.character(pastResults$V3)
  
  allResults <- dplyr::bind_rows(allGenesMatches,
                                 pastResults) %>%
    distinct()
  
  # Export all results:
  write_delim(allResults,
              file = "./candidate_genes/blastSummary.csv",
              delim = ",")
} else {
  # Export all results:
  write_delim(allGenesMatches,
              file = "./candidate_genes/blastSummary.csv",
              delim = ",")
} 

# Read in all results:
allResults <- read_delim(file = "./candidate_genes/blastSummary.csv") %>%
  distinct()

# Now we have the list of orthogroups that we want to analyze. 
#### Need to get alignments with MAFFT: ####
dir.create("./04_aligned_hogs")

alignHOGS <- function(i) {
  pathToHOGs <- "03_OrthoFinder/fasta/OrthoFinder/Results_Apr15/N1/HOG_Sequences/"
  HOGfile <- paste(pathToHOGs,
                   "N1.",
                   i,
                   ".fa",
                   sep = "")
  outputFile <- paste("./04_aligned_hogs/",
                      i,
                      "_alignment.fasta",
                      sep = "")
  
  if (file.exists(outputFile)) {
    print("Already aligned")
  } else {
    mafftCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/mafft --localpair --maxiterate 1000 --anysymbol ",
                          HOGfile,
                          " > ",
                          outputFile,
                          sep = "")
    system(mafftCommand)
  }
  
}

possiblyAlignHOGS <- purrr::possibly(alignHOGS,
                                     otherwise = "Error")
purrr::map(allResults$V1,
           possiblyAlignHOGS)


#### And gene trees with FastTree: ####
dir.create("./05_hog_trees/")

hogTrees <- function(i) {
  alignedHOGfile <- paste("./04_aligned_hogs/",
                          i,
                          "_alignment.fasta",
                          sep = "")
  outputFile <- paste("./05_hog_trees/",
                      i,
                      "_tree.txt",
                      sep = "")
  
  if (file.exists(outputFile)) {
    print("Already aligned")
  } else {
    fasttreeCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/fasttree ",
                             alignedHOGfile,
                             " > ",
                             outputFile, 
                             sep = "")
    system(fasttreeCommand)
  }
}

possiblyHogTrees <- purrr::possibly(hogTrees,
                                     otherwise = "Error")
purrr::map(allResults$V1,
           possiblyHogTrees)

#### Now we need to label the gene trees for each analysis ####
phenotypes <- read_sheet("https://docs.google.com/spreadsheets/d/166mvT3i85SJO25s155Dm_SBkzlRPCRHsQDl4Sj9SOLA/edit?usp=sharing") %>%
  select(c("abbreviation",
           "border_oceli_spots",
           "eyespots",
           "non_border_oceli_eyespots",
           "yellow_color",
           "green_color",
           "blue_color",
           "tails",
           "venous_striping"))

# Create an output directory:
dir.create("./06_labelledPhylogenies/")

# List all of the tree files:
orthofinderTreeFiles <- list.files(path = "./05_hog_trees",
                                   full.names = TRUE)

# Label trees for each trait:
labellingForSpecificTrait <- function(trait) {
  trait <- sym(trait)
  # Create an output folder for this trait:
  dir.create(paste("./06_labelledPhylogenies/",
                   trait,
                   "/",
                   sep = ""))
  # Get a list of the foreground species for this trait:
  speciesToLabel <- phenotypes %>%
    filter(!!trait == "yes")
  speciesToLabel <- speciesToLabel$abbreviation
  
  # Write a function to label just the tips on a single tree with Hyphy:
  labellingTrees <- function(i) {
    filename <- str_split_i(i,
                            pattern = "/",
                            i = 3)
    orthogroup <- str_split_i(filename,
                              pattern = "_",
                              i = 1)
    # Copy the unlabelled tree to the output folder; we'll be editing this tree three times. 
    treePaste <- paste("./06_labelledPhylogenies/",
                       trait,
                       "/labelled_",
                       orthogroup,
                       "_tree.txt",
                       sep = "")
    file.copy(from = i,
              to = treePaste)
    
    # Read in a single tree and get the tip labels of that tree:
    text <- readr::read_file(paste("./06_labelledPhylogenies/",
                                   trait,
                                   "/labelled_",
                                   orthogroup,
                                   "_tree.txt",
                                   sep = ""))
    
    # In case hyphy has already run on this, add back in the semicolon that ape requires:
    text <- gsub("\n",
                 "",
                 text,
                 perl = TRUE)
    text <- gsub("$(?<!;)",
                 ";",
                 text,
                 perl = TRUE)
    tree <- ape::read.tree(text = text)
    
    #plot(tree)
    tips <- as.data.frame(tree[["tip.label"]])
    tips 
    
    # Get the list of tip labels with matches to the species of interest:
    tipsToLabel <- filter(tips,
                          str_split_i(tips$`tree[["tip.label"]]`, pattern = "_", i = 1) %in% speciesToLabel)
    tipsToLabel <- tipsToLabel$`tree[["tip.label"]]`
    tipsToLabel
    
    write(tipsToLabel, 
          paste("tipsToLabel_",
                filename,
                sep = ""))
    
    # Run the Hyphy labelling script:
    hyphyCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/hyphy ../hyphy-analyses/LabelTrees/label-tree.bf --tree ./06_labelledPhylogenies/",
                          trait,
                          "/labelled_",
                          orthogroup,
                          "_tree.txt --list tipsToLabel_",
                          orthogroup,
                          "_tree.txt --output ./06_labelledPhylogenies/",
                          trait,
                          "/labelled_",
                          orthogroup,
                          "_tree.txt --internal-nodes \"All descendants\"",
                          sep = "")
    cat(hyphyCommand)
    system(hyphyCommand)
    file.remove(paste("tipsToLabel_",
                      filename,
                      sep = ""))
  }
  
  # Apply it over all trees:
  possiblyLabellingTrees <- purrr::possibly(labellingTrees,
                                            otherwise = "Error")
  purrr::map(orthofinderTreeFiles,
             possiblyLabellingTrees)
}

possiblyLabellingForSpecificTrait <- possibly(labellingForSpecificTrait,
                                              otherwise = "error")

traits <- c("border_oceli_spots",
            "eyespots",
            "non_border_oceli_eyespots",
            "yellow_color",
            "green_color",
            "blue_color",
            "tails",
            "venous_striping")

purrr::map(traits,
           possiblyLabellingForSpecificTrait)



