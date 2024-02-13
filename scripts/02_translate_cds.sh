#!/bin/bash

# I want to run Transdecoder on all of the downloaded transcript files.
# My current problem is that not all of the CDS sequences are appropriate for running through TransDecoder; some of them are based on gene predictions, so they have the 3' UTR removed already, meaning TransDecoder can't read them and they should just be run through my simpler DataTranslating script. I need to write some sort of if-then statement to handle these two cases.
  # I could do this by including another column in the input url file that says whether the CDS are based on gene predictions; however, this puts additional work up-front on the researcher.
  # I think that I could alternatively do this by simply running TransDecoder on all of the files, then noting which species don't have the appropriate output and re-running DataTranslating on just those files.

mkdir 02_translated_cds_files

#make a list of cds file names 
while read -r line;
do
	export speciesCode=`echo "$line" | awk -F',' '{print $2}'`
	export cdsFileName=$speciesCode'_cds.fasta'
	echo $cdsFileName >> cdsFileList.txt
done < $1

# Split the input file list into a user-specified number of chunks:
export chunkNumber="l/"$2
split --number=$chunkNumber --additional-suffix=.txt -d cdsFileList.txt "cdsFileList"
# Remove any chunks that are empty:
find . -name 'cdsFileList*' -type f -empty -delete
# Create a file listing those chunks:
ls "cdsFileList"* > "cdsChunkList.txt"

# While the list of chunks is not empty,
while [ -s "cdsChunkList.txt" ]
do
  # Run transdecoder, meaning:
  # Create a holder for the chunks:
  export batchSize=$2
  export currentBatch=0
  export batchFileNames=()
  # Then, while reading each line in our list of chunked files,
  while read -r line;
  do
    export batchFile=$line
    batchFileNames+=($batchFile)
    if [ ${#batchFileNames[@]} -eq $batchSize ]; then
      for batchFile in ${batchFileNames[@]} ; do
        #sleep 10 &
        ./Scripts/02_01_single_transdecoder_run.sh $batchFile  &
      done
      wait
      batchFileNames=()
    fi
  done < "cdsChunkList.txt"

 done
