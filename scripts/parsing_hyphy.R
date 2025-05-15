library(tidyverse)
library(googlesheets4)
library(splitstackshape) 
library(furrr)

#### Read in BUSTED-PH results and do GO enrichment ####
# List all results file with size greater than zero:
bustedphResults <- list.files(path = "./09_bustedph",
                              pattern = "*.txt",
                              full.names = TRUE,
                              recursive = TRUE)
bustedphResults <- bustedphResults[sapply(bustedphResults, file.size) > 0]

# Write a function to parse one file:
parsingBustedph <- function(i) {
  singleBustedphResult <- RJSONIO::fromJSON(content = i)
  
  correctedNames <- gsub(pattern = "_",
                         replacement = "-",
                         names(singleBustedphResult[["branch attributes"]][["0"]]))
  
  inputGenes <- paste(unique(correctedNames),
                      collapse = "|")
  
  singleBustedphResultSimple <- c(singleBustedphResult[["input"]][["file name"]],
                                  singleBustedphResult[["test results"]][["p-value"]],
                                  singleBustedphResult[["test results background"]][["p-value"]],
                                  singleBustedphResult[["test results shared distributions"]][["p-value"]],
                                  inputGenes,
                                  i)
  singleBustedphResultSimple <- as.data.frame(t(singleBustedphResultSimple))
}

# Apply that to all the results:
possiblyparsingBustedph <- purrr::possibly(parsingBustedph,
                                           otherwise = "Error")
allBustedphResults <- purrr::map(bustedphResults,
                                 possiblyparsingBustedph)
allBustedphResults <- as.data.frame(do.call(rbind, allBustedphResults))   

# Give columns meaningful names:
colnames(allBustedphResults) <- c("inputFile",
                                  "testPvalue",
                                  "backgroundPvalue",
                                  "differencePvalue",
                                  "HOGmembers",
                                  "resultFile")

# Get a column with just the orthogroup ID:
allBustedphResults$orthogroup <- stringr::str_split_i(allBustedphResults$inputFile, 
                                                      pattern = "/",
                                                      i = 9) %>%
  stringr::str_split_i(pattern = "\\.", 
                       i = 1)

# Get a column with just the trait:
allBustedphResults$trait <- stringr::str_split_i(allBustedphResults$resultFile, 
                                                      pattern = "/",
                                                      i = 3) 

# Make sure columns are numeric:
allBustedphResults$testPvalue <- as.numeric(as.character(allBustedphResults$testPvalue)) 
allBustedphResults$backgroundPvalue <- as.numeric(as.character(allBustedphResults$backgroundPvalue)) 
allBustedphResults$differencePvalue <- as.numeric(as.character(allBustedphResults$differencePvalue)) 

# Do FDR corrections on p-values:
allBustedphResults$testPvalueFDR <- p.adjust(allBustedphResults$testPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$backgroundPvalueFDR <- p.adjust(allBustedphResults$backgroundPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$differencePvalueFDR <- p.adjust(allBustedphResults$differencePvalue, method='BH') %>%
  as.numeric(as.character()) 

# Get a descriptive column for the gene's selective regime:
allBustedphResults <- allBustedphResults %>% mutate(selectionOn =
                                                      case_when(as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "ForegroundOnly",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) <= 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "BackgroundOnly",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) <= 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) <= 0.05 ~ "NoEvidenceOfSelection",
                                                                
                                                                as.numeric(as.character(testPvalueFDR)) > 0.05 & 
                                                                  as.numeric(as.character(backgroundPvalueFDR)) > 0.05 &
                                                                  as.numeric(as.character(differencePvalueFDR)) > 0.05 ~ "NoEvidenceOfSelection"))

# Add info about the original candidate gene:
candidateGeneInfo <- read.csv("./04_candidate_genes/blastSummary.csv")
allBustedphResults <- left_join(allBustedphResults,
                                candidateGeneInfo,
                                by = c("orthogroup" = "V1"))
allBustedphResults$candidate <- str_split_i(allBustedphResults$V2,
                                            pattern = "/",
                                            i = 3) %>%
  str_split_i(pattern = "\\.",
              i = 1)
candidateGeneInfo <- read_sheet("https://docs.google.com/spreadsheets/d/1D_cnkNskGndVqMIfkKU3dpjHX3_6OD4Cv78YIiRiGco/edit?usp=drive_link")
candidateGeneInfo$`NCBI locus ID` <- as.character(candidateGeneInfo$`NCBI locus ID`)
allBustedphResults <- left_join(allBustedphResults,
                                candidateGeneInfo,
                                by = c("candidate" = "NCBI locus ID"))

# Export the results to a csv:
readr::write_delim(allBustedphResults,
                   file = "./allBustedPHResults.csv",
                   delim = ",",
                   quote = "none")
