library(tidyverse)
library(googlesheets4)
library(splitstackshape) 
library(furrr)


future::plan(multisession)
options(future.globals.maxSize= +Inf)

#### Run BUSTED-PH ####
dir.create("./09_bustedph/")
dir.create("./06_alignments/06_aligned_orthogroup_cds_sequences_no_stops")

# List all of the tree files:
orthofinderTreeFiles <- list.files(path = "./07_hog_trees",
                                   full.names = TRUE)

# Write a function to run BUSTED-PH on a single trait:
bustedPhSingleTrait <- function(trait) {
  dir.create(paste("./09_bustedph/",
                   trait,
                   sep = ""))
  
  trait <- sym(trait)
  
  # Write a function to analyze a single orthogroup:
  bustedPh <- function(i) {
    filename <- str_split_i(i,
                            pattern = "/",
                            i = 3)
    orthogroup <- str_split_i(filename,
                              pattern = "_",
                              i = 1)
    outputFile <- paste("./09_bustedph/",
                        trait,
                        "/",
                        trait,
                        "_",
                        orthogroup,
                        ".txt",
                        sep = "")
    
    
    if (file.exists(outputFile)) {
      print("Already analyzed.")
    } else {
      # First remove any stop codons:
      removeStopsCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/hyphy /Users/mbarkdull/Projects/hyphy-develop/res/TemplateBatchFiles/CleanStopCodons.bf Universal ./06_alignments/05_aligned_orthogroup_cds_sequences/",
                                  orthogroup,
                                  ".fasta No/Yes ./06_alignments/06_aligned_orthogroup_cds_sequences_no_stops/",
                                  orthogroup,
                                  ".fasta ",
                                  sep = "")
      cat(removeStopsCommand)
      system(removeStopsCommand)
      
      
      # Run the BUSTED-PH analysis:
      hyphyCommand <- paste("/opt/anaconda3/envs/bioinformatics/bin/hyphy ../hyphy-analyses/BUSTED-PH/BUSTED-PH.bf --alignment ./06_alignments/06_aligned_orthogroup_cds_sequences_no_stops/",
                            orthogroup,
                            ".fasta --tree ./08_labelledPhylogenies/",
                            trait,
                            "/labelled_",
                            orthogroup,
                            "_tree.txt --srv Yes --branches Foreground --output ",
                            outputFile,
                            sep = "")
      cat(hyphyCommand)
      write(paste("starting", trait, "and", orthogroup),
            file = "bustedph_progress.txt",
            append = TRUE)
      system(hyphyCommand)
      write(paste("*****Finished", trait, "and", orthogroup),
            file = "bustedph_progress.txt",
            append = TRUE)
    } 
  }
  
  # Apply it over all orthogroups for this trait:
  possiblyBustedPh <- purrr::possibly(bustedPh,
                                      otherwise = "Error")
  purrr::map(orthofinderTreeFiles,
             possiblyBustedPh)
}

possiblyBustedPhSingleTrait <- possibly(bustedPhSingleTrait,
                                        otherwise = "error")

traits <- c("border_oceli_spots",
            "eyespots",
            "non_border_oceli_eyespots",
            "yellow_color",
            "green_color",
            "blue_color",
            "tails",
            "venous_striping")

furrr::future_map(traits,
                  possiblyBustedPhSingleTrait)

