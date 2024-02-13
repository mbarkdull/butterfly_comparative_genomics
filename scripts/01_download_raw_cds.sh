#!/bin/bash

# I am going to have an input file that is a tab-delimited text file with data urls and desired file speciesCodes to send the download.

# The command to run this script is `./DataDownload [input file]`. This will read in the input urls file and download the transcripts, proteins, and GFF files for each species.
  mkdir ./01_RawCDSFiles

  while read -r line;
  do
    # This creates a variable, url, that holds the information about the download url on this line of the input file:
    export cdsUrl=`echo "$line" | awk -F',' '{print $3}'`
    export speciesCode=`echo "$line" | awk -F',' '{print $2}'`

    # This uses the second field on that line to get the file speciesCode where that data should ultimately be sent:
    export cdsFileName=$speciesCode'_cds.fasta'

    cd ./01_RawCDSFiles

    # This tells us what the zipped, downloaded file should be sent to:
    #if [[ $cdsUrl == *.zip ]]
      #then
        #export zipFn=$speciesCode.zip
      #else [[ $cdsUrl == *.gz ]]
        #export zipFn=$speciesCode.gz
    #fi
    export zipFn=$speciesCode.gz

    # These lines print the file speciesCode variable information to the console so you can see it:
    echo __________________________________________________
    #echo Dowloading data from the url: $cdsUrl
    #echo Downloading the zipped data to $zipFn
    echo wget $cdsUrl 

    # This sends the download url to wget, unzips the downloaded file, and respeciesCodes it to the speciesCode you indicated in the input file.
    wget $cdsUrl 
    #-O $zipFn

    #if [[ $zipFn == *.zip ]]
      #then
      #  unzip -p $zipFn > $cdsFileName
      #  echo Unzipping $zipFn to $cdsFileName
      #else [[ $zipFn == *.gz ]]
      #  gunzip -c $zipFn > $cdsFileName
      #  echo Unzipping $zipFn to $cdsFileName
    #fi

   # rm $zipFn
    echo ________________________Done with cds download__________________________

    cd ../

  done < $1

  # The command to run this script is `./DataDownload [input file]`. This will read in the input urls file and download the transcripts, proteins, and GFF files for each species.
