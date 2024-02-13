#!/bin/bash

# I want to run Transdecoder on all of the downloaded transcript files.
# My current problem is that not all of the CDS sequences are appropriate for running through TransDecoder; some of them are based on gene predictions, so they have the 3' UTR removed already, meaning TransDecoder can't read them and they should just be run through my simpler DataTranslating script. I need to write some sort of if-then statement to handle these two cases.
  # I could do this by including another column in the input url file that says whether the CDS are based on gene predictions; however, this puts additional work up-front on the researcher.
  # I think that I could alternatively do this by simply running TransDecoder on all of the files, then noting which species don't have the appropriate output and re-running DataTranslating on just those files.

mkdir 02_translated_cds_files

while read -r line;
do
	export speciesCode=`echo "$line" | awk -F',' '{print $2}'`
	export cdsFileName=$speciesCode'_cds.fasta'
	cp ./01_RawCDSFiles/$cdsFileName 02_translated_cds_files
	cd 02_translated_cds_files
	# Now we can run Transdecoder on the cleaned file:
    echo "First, attempting TransDecoder run on $cdsFileName"
    /programs/TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $cdsFileName
    /programs/TransDecoder-v5.5.0/TransDecoder.Predict -t $cdsFileName --single_best_only
done < $1