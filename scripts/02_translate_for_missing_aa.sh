# Script to run Transdecoder for any genomes for which we could not download a proteins file.

# List all of the cds files:
ls ./01_RawCDSFiles/ >> cds_files.txt

# list all of the aa files:
ls ./02_AAFiles/ >> aa_files.txt

# Make a directory in which to run Transdecoder:
mkdir ./02_transdecoder_run

while read -r line;
do
	export speciesCode=`echo "$line" | awk -F'_' '{print $1}'`
  export outputPepFile="./02_AAFiles/"$speciesCode"_aa.fasta"
  if [ -f "$outputPepFile" ]; then
    echo "*************************************************************"
    echo "$outputPepFile exists; an amino acid sequences file was downloaded,
    or TransDecoder has already been run on this species."
    echo "*************************************************************"
  else
    export cdsFileName=$line
	  # Now we can run Transdecoder on the nucleotide file:
    echo "First, attempting TransDecoder run on $cdsFileName"
    cd ./02_transdecoder_run
    ../TransDecoder/TransDecoder.LongOrfs -t ../01_RawCDSFiles/$cdsFileName
    ../TransDecoder/TransDecoder.Predict -t ../01_RawCDSFiles/$cdsFileName --single_best_only
    cp "./"$speciesCode"_cds.fasta.transdecoder.pep" ../02_AAFiles/$speciesCode"_aa.fasta"
    cd ../
  fi
done < cds_files.txt

rm -r ./02_transdecoder_run
rm cds_files.txt
rm aa_files.txt