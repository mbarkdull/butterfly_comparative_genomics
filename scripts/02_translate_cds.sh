#!/bin/bash

# I want to run Transdecoder on all of the downloaded transcript files.
# My current problem is that not all of the CDS sequences are appropriate for running through TransDecoder; some of them are based on gene predictions, so they have the 3' UTR removed already, meaning TransDecoder can't read them and they should just be run through my simpler DataTranslating script. I need to write some sort of if-then statement to handle these two cases.
  # I could do this by including another column in the input url file that says whether the CDS are based on gene predictions; however, this puts additional work up-front on the researcher.
  # I think that I could alternatively do this by simply running TransDecoder on all of the files, then noting which species don't have the appropriate output and re-running DataTranslating on just those files.

mkdir 02_translated_cds_files

while read -r line;
do
	export speciesCode=`echo "$line" | awk -F',' '{print $2}'`
	export cdsFile=
	# Now we can run Transdecoder on the cleaned file:
    echo "First, attempting TransDecoder run on $cleanName"
    /programs/TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $cleanName
    /programs/TransDecoder-v5.5.0/TransDecoder.Predict -t $cleanName --single_best_only




done < $1


# The command to run this is `./DataTransdecoder [input file]`. This will read in each line of the input file; create a variable that is the url from which to download the data, based on the first field of the line of the input file; it will create a temporary file that holds the individual tree for that gene; and will then send those pieces of information to HyPhy to be run.
  while read -r line;
  do
    # This uses the second field on each line of the input urls file to get the file name where the filtered, longest transcripts data was sent, ex. nful_filteredTranscripts.fasta:
    export abbrev=`echo "$line" | awk -F',' '{print $4}'`
    export rawName=$abbrev'_filteredTranscripts.fasta'
    # Now we construct the cleaned file name, which is the one that we want to run Transdecoder on, ex. cleanedacept_transcripts:
    export cleanName=cleaned$rawName

    echo "Translating the nucleotide sequences in $cleanName"

    # Now we can make a directory for translating the inputs and copy the cleaned files there:
    mkdir ./4_1_TranslatedData
    cd ./4_1_TranslatedData
    cp ../3_CleanedData/$cleanName ./$cleanName
    cp ../scripts/TranscriptFilesTranslateScript.py ./

    # Now we can run Transdecoder on the cleaned file:
    echo "First, attempting TransDecoder run on $cleanName"
    /programs/TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $cleanName
    /programs/TransDecoder-v5.5.0/TransDecoder.Predict -t $cleanName --single_best_only

    # Now we need to run the TranscriptFilesTranslateScript.py script just on those files that didn't successfully go through TransDecoder.
    mkdir ./OutputFiles
    export transdecoderFile=$cleanName".transdecoder.pep"
    echo "Checking for $transdecoderFile"

    #Check if the Transdecoder file exists; if it does, move the amino acid and nucleotide files to appropriate directories.
    if [ -f "$transdecoderFile" ]; then
        echo "TransDecoder run successful; $transdecoderFile exists"
        # Copy the TranDecoder file to something with a more readable name, ex. translatedacep_transcripts
        cp $transdecoderFile ./OutputFiles/translated$rawName
        echo "Copying $transdecoderFile to ./OutputFiles/translated$rawName"
        cat ./OutputFiles/translated$rawName | tr  "."  "_" > temp.fasta
        mv temp.fasta ./OutputFiles/translated$rawName
        rm temp.fasta
        sed -i'.original' -e "s|\.|_|g" ./OutputFiles/translated$rawName

        # copy the created coding sequence file to a coding sequences folder:
        export cdsFile=$cleanName".transdecoder.cds"
        mkdir ../4_2_TransdecoderCodingSequences
        cp $cdsFile ../4_2_TransdecoderCodingSequences
        mv ../4_2_TransdecoderCodingSequences/$cdsFile ../4_2_TransdecoderCodingSequences/cds_$rawName
        cat ../4_2_TransdecoderCodingSequences/cds_$rawName | tr  "."  "_" > temp.fasta
        mv temp.fasta ../4_2_TransdecoderCodingSequences/cds_$rawName
        rm temp.fasta
        sed -i'.original' -e "s|\.|_|g" ../4_2_TransdecoderCodingSequences/cds_$rawName

    # If the Transdecoder file does not exist, run my translating script on the file instead.
    else
        echo "$transdecoderFile was not produced because the input file already consisted of coding sequences (CDS). Running TranscriptFilesTranslateScript.py on $cleanName to produce amino acid sequences"
        # This will make a list of all of the cleaned names and pass them into a text file; then that text file is the input for the python script.
        echo "$cleanName" >> CleanedNucleotideSequences.txt

        echo __________________________________________________

        python ./TranscriptFilesTranslateScript.py CleanedNucleotideSequences.txt
        mv ./"translated"$rawName ./OutputFiles/"translated"$rawName
        sed -i'.original' -e "s|\.|_|g" ./OutputFiles/"translated"$rawName

        mkdir ../4_2_TransdecoderCodingSequences
        cp ../3_CleanedData/$cleanName ../4_2_TransdecoderCodingSequences
        mv ../4_2_TransdecoderCodingSequences/$cleanName ../4_2_TransdecoderCodingSequences/cds_$rawName
        sed -i'.original' -e "s|\.|_|g" ../4_2_TransdecoderCodingSequences/cds_$rawName
    fi

    cd ../

  done < $1

  # The "raw name" from the input file is XXX_transcripts.fasta
  # The python script produces files named like: translatedXXX_transcripts.fasta
  # TransDecoder produces files named like: cleanedXXX_transcripts.fasta.transdecoder.pep
  # I then copy those result files to files named like: translatedXXX_transcripts.fasta
    # in the folder ./TransDecoderResults