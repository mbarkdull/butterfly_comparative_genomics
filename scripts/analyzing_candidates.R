library(rBLAST)
library(tidyverse)
library(googlesheets4)
library(rentrez)
library(splitstackshape) 
library(furrr)

future::plan(multisession)
options(future.globals.maxSize= +Inf)

Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                        "/usr/local/ncbi/blast/bin/", 
                        sep = .Platform$path.sep))

### Download the candidate genes from NCBI based on my google sheet: ####
# Read in the information from my public google sheet, getting NCBI gene IDs:
candidates <- read_sheet("https://docs.google.com/spreadsheets/d/1D_cnkNskGndVqMIfkKU3dpjHX3_6OD4Cv78YIiRiGco/edit?usp=drive_link")
candidates$`NCBI locus ID` <- as.character(candidates$`NCBI locus ID`)
ncbiIDs <- candidates$`NCBI locus ID`

# Create an output directory to which the files can be downloaded:
dir.create("./04_candidate_genes/")

# Get the list of already-downloaded proteins:
geneFiles <- dir("04_candidate_genes", 
                 recursive = TRUE, 
                 full.names = TRUE)

# Remove empty files
unlink(geneFiles[file.info(geneFiles)[["size"]]==0], 
       recursive = TRUE, 
       force = FALSE)

# List the downloaded genes:
alreadyDownloaded <- list.files(path = "./04_candidate_genes/",
                                pattern = "*fasta") %>% 
  str_split_i(pattern = "\\.",
              i = 1)

# Get a list of only proteins that still need to be downloaded:
ncbiIDs <- base::setdiff(ncbiIDs,
                         alreadyDownloaded)

# Write a function to download a single protein sequence:
gettingGeneSequence <- function(ID) {
  # Link the gene ID to a nucleotide sequence ID:
  geneToProt <- entrez_link(dbfrom = 'gene', 
                            id = ID, 
                            db = 'protein')
  protID <- geneToProt[["links"]][["gene_protein_refseq"]][[1]]
  
  
  # Get the gene sequence as a fasta:
  entrez_fetch(db = "protein", 
               id = protID, 
               rettype = "fasta") %>%
    write(file = paste("./04_candidate_genes/",
                       ID,
                       ".fasta",
                       sep = ""))
}

# Make a safe version with possibly:
possiblygettingGeneSequence <- purrr::possibly(gettingGeneSequence,
                                               otherwise = "Error")
# Download the list of protein sequences:
purrr::map(ncbiIDs,
           possiblygettingGeneSequence)

#### Combine all the proteomes against which we'll BLAST ####
# Use an if-else statement so this step doesn't have to be repeated:
if (file.exists("allProteomes.fasta")) {
  print("The file exists!")
} else {
  print("The file does not exist.")
  # Construct a list of the proteomes:
  allProteomes <- list.files("03_OrthoFinder_AA/fasta",
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
makeblastdb(file = './allProteomes.fasta', 
            dbtype = "prot")

# Process the orthogroups so they can be identified:
# Read in the orthogroups file:
orthogroupMembers <- read_delim(file = "./03_OrthoFinder_AA/fasta/OrthoFinder/Results_Apr15/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", 
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
allGenes <- list.files("./04_candidate_genes",
                       full.names = TRUE)

# Read in any past results, so as not to duplicate analyses:
if(file.exists("./04_candidate_genes/blastSummary.csv") == TRUE){
  pastResults <- read_delim(file = "./04_candidate_genes/blastSummary.csv")
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
if(file.exists("./04_candidate_genes/blastSummary.csv") == TRUE) {
  pastResults <- read_delim(file = "./04_candidate_genes/blastSummary.csv") %>%
    distinct()
  pastResults$V3 <- as.character(pastResults$V3)
  
  allResults <- dplyr::bind_rows(allGenesMatches,
                                 pastResults) %>%
    distinct()
  
  # Export all results:
  write_delim(allResults,
              file = "./04_candidate_genes/blastSummary.csv",
              delim = ",")
} else {
  # Export all results:
  write_delim(allGenesMatches,
              file = "./04_candidate_genes/blastSummary.csv",
              delim = ",")
} 

# Read in all results:
allResults <- read_delim(file = "./04_candidate_genes/blastSummary.csv") %>%
  distinct()

# Now we have the list of orthogroups that we want to analyze. 
#### Get the corresponding nucleotide sequences ####
# Read in all of the BLAST results:
allResults <- read_delim(file = "./04_candidate_genes/blastSummary.csv") %>%
  distinct()

# Read in all of the orthogroups:
readOrthogroup <- function(orthogroupNumber) {
  orthogroup <- phylotools::read.fasta(paste("./03_OrthoFinder_AA/fasta/OrthoFinder/Results_Apr15/N1/HOG_Sequences/N1.",
                                             orthogroupNumber,
                                             ".fa",
                                             sep = "")) 
  orthogroup$orthogroup <- orthogroupNumber
  orthogroup$species <- str_split_i(orthogroup$seq.name,
                                    pattern = "_",
                                    i = 1)
  orthogroup$originalName <- gsub(pattern = "^[^_]*_", 
                                  replacement = "", 
                                  x = orthogroup$seq.name)
  orthogroup$originalName <- gsub(pattern = "\\.1", 
                                  replacement = "", 
                                  x = orthogroup$originalName)
  return(orthogroup)
}
possiblyReadOrthogroup <- purrr::possibly(readOrthogroup,
                                          otherwise = "error")
allOrthogroups <- purrr::map(allResults$V1,
                             possiblyReadOrthogroup)
allOrthogroups <- as.data.frame(do.call(rbind, allOrthogroups)) 

speciesList <- list.files("./01_RawCDSFiles/") %>%
  str_split_i(pattern = "_",
              i = 1)

findAllProteins <- function(speciesName) {
  # Filter the orthogroups to just the one species:
  speciesOrthogroup <- filter(allOrthogroups,
                              species == speciesName)
  
  # Read in one annotation:
  annotations <- list.files(pattern = "_annotation")
  annotation <- annotations[grepl(speciesName, annotations)]
  annotation <- ape::read.gff(gzfile(annotation))
  annotation$species <- speciesName
  
  #Read in the gene sequences:
  geneSequences <- phylotools::read.fasta(paste("./01_RawCDSFiles/",
                                                speciesName,
                                                "_cds.fasta",
                                                sep = ""))
  geneSequences$geneName <- str_split_i(geneSequences$seq.name,
                                        pattern = " ",
                                        i = 1)
  
  # This function will run over speciesOrthogroup$originalName to find genes:
  findProteinMatches <- function(geneName,
                                 orthogroupNumber) {
    # Look for the protein names in the annotation:
    print("finding protein in annotation")
    proteinMatches <- filter(annotation,
                             grepl(geneName, attributes)) %>%
      distinct()
    proteinMatches <- cSplit(proteinMatches, 
                             'attributes', 
                             ';') 
    proteinMatches <- pivot_longer(proteinMatches, 
                                   cols = starts_with("attributes"),
                                   names_to = NULL) %>%
      filter(!is.na(value)) %>%
      filter(!grepl(pattern = "gbkey=|Note=|product=|version=|start_range=|_eAED=|_QI=|_AED=|transl_except=",
                    x = value))
    
    if (speciesName == "baph") {
      proteinMatches <- proteinMatches %>%
        filter(!grepl(pattern = "gene=",
                      x = value))
    }
    
    proteinMatches$nameInAnnotation <- str_split_i(proteinMatches$value,
                                                   pattern = "=",
                                                   i = 2) 
    
    if(!unique(proteinMatches$species) %in% c('juco', 'pael')) {
      proteinMatches$nameInAnnotation <- gsub(pattern = "^[^:]*:", 
                                              replacement = "", 
                                              x = proteinMatches$nameInAnnotation)
    }
    
    
    featuresToKeep <- c("mRNA",
                        "CDS")
    proteinMatches <- proteinMatches %>%
      filter(type %in% featuresToKeep) %>% 
      dplyr::filter(nchar(nameInAnnotation) > 2)
    
    
    
    # Now use that to find the coding sequences:
    uniqueSequencesToSearch <- unique(proteinMatches$nameInAnnotation)
    findCodingSequences <- function(oneSearch) {
      if (unique(proteinMatches$species) %in% c('hede', 'pael')) {
        print(paste("Finding",
                    oneSearch))
        matchingCodingSequences <- filter(geneSequences,
                                          geneName == oneSearch) %>%
          distinct() %>%
          filter(seq.name != "error")
        print(paste("&&&&&&&&&&&&&&& The gene",
                    oneSearch,
                    "matched:",
                    matchingCodingSequences$seq.name))
        
        matchingCodingSequences$searchSequence <- geneName
        matchingCodingSequences$orthogroup <- orthogroupNumber
        return(matchingCodingSequences)
      } else {
        print(paste("Finding",
                    oneSearch))
        matchingCodingSequences <- filter(geneSequences,
                                          grepl(oneSearch, geneName))
        matchingCodingSequences$searchSequence <- geneName
        matchingCodingSequences$orthogroup <- orthogroupNumber
        matchingCodingSequences<- matchingCodingSequences %>%
          distinct() %>%
          filter(seq.name != "error")
        print(paste("&&&&&&&&&&&&&&& The gene",
                    oneSearch,
                    "matched:",
                    matchingCodingSequences$seq.name)
        )
        return(matchingCodingSequences)
      } 
    }
    
    possiblyFindCodingSequences <- purrr::possibly(findCodingSequences,
                                                   otherwise = "error")
    
    matchingCodingSequences <- purrr::map(uniqueSequencesToSearch,
                                          possiblyFindCodingSequences) 
    matchingCodingSequences <- as.data.frame(do.call(rbind, matchingCodingSequences)) 
    matchingCodingSequences <- matchingCodingSequences 
    
  }
  possiblyFindProteinMatches <- purrr::possibly(findProteinMatches,
                                                otherwise = "error")
  matches <- purrr::map2(speciesOrthogroup$originalName,
                         speciesOrthogroup$orthogroup,
                         possiblyFindProteinMatches) 
  matches <- as.data.frame(do.call(rbind, matches)) %>%
    select(c(searchSequence,
             seq.text,
             orthogroup)) %>%
    filter(orthogroup != "error") %>%
    distinct()
  matches$species <- speciesName
  matches <- full_join(matches, 
                       speciesOrthogroup,
                       by = c("searchSequence" = "originalName",
                              "orthogroup" = "orthogroup")) %>%
    distinct()
  return(matches)
}
possiblyFindAllProteins <- purrr::possibly(findAllProteins,
                                           otherwise = "error")

allProteins <- purrr::map(speciesList,
                          possiblyFindAllProteins)
allProteins <- as.data.frame(do.call(rbind, allProteins)) 
write.csv(allProteins,
          "./allNucleotideOrthogroups.csv")

#### Now get codon-aware alignment (MAFFT then PAL2NAL) ####
# Read in all of the BLAST results:
allResults <- read_delim(file = "./04_candidate_genes/blastSummary.csv") %>%
  distinct()

# Create output directories:
dir.create("./06_alignments")
dir.create("./06_alignments/01_aligned_AAs")

alignHOGS <- function(i) {
  pathToHOGs <- "./03_OrthoFinder_AA/fasta/OrthoFinder/Results_Apr15/N1/HOG_Sequences/"
  HOGfile <- paste(pathToHOGs,
                   "N1.",
                   i,
                   ".fa",
                   sep = "")
  
  sequence <- phylotools::read.fasta(file = HOGfile) 
  sequence$species <- str_split_i(sequence$seq.name,
                                  pattern = "_",
                                  i = 1)
  sequence <- sequence %>%
    filter(species != "cace") %>%
    filter(species != "hela") %>%
    select(-c("species"))
  
  phylotools::dat2fasta(dat = sequence,
                        outfile = "temp.fasta")
  
  outputFile <- paste("./06_alignments/01_aligned_AAs/",
                      i,
                      "_alignment.fasta",
                      sep = "")
  
  if (file.exists(outputFile)) {
    print("Already aligned")
    file.remove("temp.fasta")
  } else {
    mafftCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/mafft --localpair --maxiterate 1000 --anysymbol temp.fasta > ",
                          outputFile,
                          sep = "")
    system(mafftCommand)
    
    file.remove("temp.fasta")
  }
  
}

possiblyAlignHOGS <- purrr::possibly(alignHOGS,
                                     otherwise = "Error")
purrr::map(allResults$V1,
           possiblyAlignHOGS)

# Then generate codon-aware CDS alignments with pal2nal: 
dir.create("./06_alignments/02_individual_amino_acid_sequences")
dir.create("./06_alignments/03_individual_cds_sequences")
dir.create("./06_alignments/04_pal2nal_alignments")

# First we need to get each individual amino acid sequence as a fasta file:
aminoAcidAlignments <- list.files("./06_alignments/01_aligned_AAs",
                                  full.names = TRUE)

possiblyRead.fasta <- possibly(phylotools::read.fasta,
                               otherwise = "Error")
allaminoAcidAlignments <- purrr::map(aminoAcidAlignments,
                                     possiblyRead.fasta)
allaminoAcidAlignments <- as.data.frame(do.call(rbind, allaminoAcidAlignments))

# Export each amino acid sequence as a separate file:
exportIndividualSequences <- function(sequenceName,
                                      sequencesToFilter,
                                      outputDirectory) {
  individualSequence <- filter(sequencesToFilter,
                               seq.name == sequenceName)
  
  outputFile <- paste(outputDirectory,
                      sequenceName,
                      ".fasta",
                      sep = "")
  
  phylotools::dat2fasta(dat = individualSequence,
                        outfile = outputFile)
}


possiblyExportIndividualSequences <- purrr::possibly(exportIndividualSequences,
                                                     otherwise = "error")
purrr::map(allaminoAcidAlignments$seq.name,
           ~possiblyExportIndividualSequences(.x, 
                                              sequencesToFilter = allaminoAcidAlignments,
                                              outputDirectory = "./06_alignments/02_individual_amino_acid_sequences/"))

# Get every coding sequence as a separate fasta:
allCdsFiles <- read.csv("./allNucleotideOrthogroups.csv") %>%
  select(c("seq.name",
           "seq.text.x"))

purrr::map(allCdsFiles$seq.name,
           ~possiblyExportIndividualSequences(.x, 
                                              sequencesToFilter = allCdsFiles,
                                              outputDirectory = "./06_alignments/03_individual_cds_sequences/"))

# Now run pal2nal individually on pairs of AA/CDS files:
cdsToAlign <- list.files("./06_alignments/03_individual_cds_sequences/") 
cdsToAlign <- gsub(pattern = ".fasta",
                   replacement = "",
                   cdsToAlign)

runPAL2NAL <- function(gene) {
  alignedProteinsFile <- paste("./06_alignments/02_individual_amino_acid_sequences/", 
                               gene, 
                               ".fasta", 
                               sep = "")
  cdsSequencesFile <- paste("./06_alignments/03_individual_cds_sequences/", 
                            gene, 
                            ".fasta", 
                            sep = "")
  outputFileName <- paste("./06_alignments/04_pal2nal_alignments/", 
                          gene, 
                          "_nucleotideAlignments.fasta", 
                          sep = "")
  
  if (file.exists(outputFileName)) {
    print("Already aligned.")
  } else {
    # Construct a system command to run pal2nal
    pal2nalCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/pal2nal.pl  ",
                            alignedProteinsFile,
                            "  ",
                            cdsSequencesFile,
                            "  -output fasta > ",
                            outputFileName,
                            sep = "")
    cat(pal2nalCommand)
    system(pal2nalCommand)
  }
}

possiblyrunPAL2NAL <- possibly(runPAL2NAL,
                               otherwise = "Error")
# Run the alignment:
purrr::map(cdsToAlign,
           possiblyrunPAL2NAL)

# Now recombine the individual cds sequences into orthogroups:
dir.create("./06_alignments/05_aligned_orthogroup_cds_sequences")

aminoAcidAlignments <- list.files("06_alignments/01_aligned_AAs/",
                                  full.names = TRUE)

recombiningNucleotideOrthogroups <- function(i) {
  # Read in one orthogroup:
  aminoAcidOrthogroup <- phylotools::read.fasta(i)
  orthogroupNumber <- str_split_i(i,
                                  pattern = "/",
                                  i = 4) %>%
    str_split_i(pattern = "_",
                i = 1)
  # Get a list of the gene names in the orthogroup:
  geneNames <- aminoAcidOrthogroup$seq.name
  
  geneNames <- paste("./06_alignments/04_pal2nal_alignments/",
                     geneNames,
                     "_nucleotideAlignments.fasta",
                     sep = "")
  possiblyRead.fasta <- possibly(phylotools::read.fasta,
                                 otherwise = "Error")
  allGenesInOrthogroup <- purrr::map(geneNames,
                                     possiblyRead.fasta)
  allGenesInOrthogroup <- as.data.frame(do.call(rbind, allGenesInOrthogroup)) %>%
    filter(seq.name != "Error") %>%
    filter(seq.text != "NA")
  
  phylotools::dat2fasta(allGenesInOrthogroup,
                        outfile = paste("./06_alignments/05_aligned_orthogroup_cds_sequences/",
                                        orthogroupNumber,
                                        ".fasta",
                                        sep = ""))
}

possiblyrecombiningNucleotideOrthogroups <- possibly(recombiningNucleotideOrthogroups,
                                                     otherwise = "Error")
# Run the alignment in parallel, using furrr:
purrr::map(aminoAcidAlignments,
           possiblyrecombiningNucleotideOrthogroups)

#### Get gene trees with FastTree: ####
dir.create("./07_hog_trees/")

hogTrees <- function(i) {
  alignedHOGfile <- paste("./06_alignments/01_aligned_AAs/",
                          i,
                          "_alignment.fasta",
                          sep = "")
  outputFile <- paste("./07_hog_trees/",
                      i,
                      "_tree.txt",
                      sep = "")
  
  if (file.exists(outputFile)) {
    print("Already aligned")
  } else {
    nucleotideFile <- phylotools::read.fasta(file = paste("./06_alignments/05_aligned_orthogroup_cds_sequences/",
                                                          i,
                                                          ".fasta",
                                                          sep = ""))
    aaHogFile <- phylotools::read.fasta(file = alignedHOGfile) %>%
      filter(seq.name %in% nucleotideFile$seq.name)
    
    phylotools::dat2fasta(aaHogFile,
                          outfile = alignedHOGfile)
    
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
dir.create("./08_labelledPhylogenies/")

# List all of the tree files:
orthofinderTreeFiles <- list.files(path = "./07_hog_trees",
                                   full.names = TRUE)

# Label trees for each trait:
labellingForSpecificTrait <- function(trait) {
  trait <- sym(trait)
  # Create an output folder for this trait:
  dir.create(paste("./08_labelledPhylogenies/",
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
    # Copy the unlabelled tree to the output folder: 
    treePaste <- paste("./08_labelledPhylogenies/",
                       trait,
                       "/labelled_",
                       orthogroup,
                       "_tree.txt",
                       sep = "")
    file.copy(from = i,
              to = treePaste)
    
    # Read in a single tree and get the tip labels of that tree:
    text <- readr::read_file(paste("./08_labelledPhylogenies/",
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
    hyphyCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/hyphy ../hyphy-analyses/LabelTrees/label-tree.bf --tree ./08_labelledPhylogenies/",
                          trait,
                          "/labelled_",
                          orthogroup,
                          "_tree.txt --list tipsToLabel_",
                          orthogroup,
                          "_tree.txt --output ./08_labelledPhylogenies/",
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
    
    # Read in the labelled tree and replace any periods with underscores:
    text <- readr::read_file(paste("./08_labelledPhylogenies/",
                                   trait,
                                   "/labelled_",
                                   orthogroup,
                                   "_tree.txt",
                                   sep = ""))
    text <- gsub("\n",
                 "",
                 text,
                 perl = TRUE)
    text <- gsub("$(?<!;)",
                 ";",
                 text,
                 perl = TRUE)
    tree <- ape::read.tree(text = text)
    tree[["tip.label"]] <- gsub(pattern = "\\.",
                                replacement = "_",
                                x = tree[["tip.label"]])
    tree[["tip.label"]] <- gsub(pattern = "-",
                                replacement = "_",
                                x = tree[["tip.label"]])
    write.tree(phy = tree,
               file = paste("./08_labelledPhylogenies/",
                            trait,
                            "/labelled_",
                            orthogroup,
                            "_tree.txt",
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


