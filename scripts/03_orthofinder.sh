#!/bin/bash

# Script to run OrthoFinder. Takes two command line options:
  # $1: how many threads to use for highly parallelizable steps of OrthoFinder
  # Path to the reference species' transcripts
mkdir ./03_OrthoFinder/
mkdir ./03_OrthoFinder/fasta

for directory in ./02_translated_cds_files/*/
do
  # Remove the trailing "/"
  export directory=${directory%*/}
  # Extract the species code by printing everything after the final "/"
  export sampleCode="${directory##*/}"
  echo This is the sample code: $sampleCode
  # Copy the proteins file into the orthofinder directory:
  cp $directory"/"$sampleCode"_cds.fasta.transdecoder.pep" ./03_OrthoFinder/fasta/
  mv ./03_OrthoFinder/fasta/$sampleCode"_cds.fasta.transdecoder.pep" ./03_OrthoFinder/fasta/$sampleCode"_proteins.fasta"
  # Append a species code to each gene name:
  export samplePrefix=$sampleCode"_"
  sed -i "s/>/>$samplePrefix/g" "./03_OrthoFinder/fasta/"$sampleCode"_proteins.fasta"
done

mv ./llun_proteins.fasta ./03_OrthoFinder/fasta/llun_proteins.fasta

mkdir -p /workdir/$USER/tmp
source /programs/miniconda3/bin/activate orthofinder-2.5.4

# Options:
  # -S: sequence search option; here using diamond.
  # -t: number of threads for upstream processes
  # -f: directory with input data
  # -d: run on nucleotide sequences rather than amino acid sequences
  # -M: infer multiple sequence aligments and gene trees.
orthofinder -S diamond -t $1 -f ./03_OrthoFinder/fasta/ -M msa

# Now we need to use an orthofinder utility to extract just the orthogroups we are interested in (the N1 clade)
cd ./03_OrthoFinder/
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder_source.tar.gz
tar -xzvf OrthoFinder_source.tar.gz
chmod u+rwx OrthoFinder_source/tools/create_files_for_hogs.py
export resultsDate="$(date +'%b%d')"
OrthoFinder_source/tools/create_files_for_hogs.py ./fasta/OrthoFinder/Results_$resultsDate/ ./fasta/OrthoFinder/Results_$resultsDate/ N1