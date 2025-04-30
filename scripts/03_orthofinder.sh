#!/bin/bash

# Script to run OrthoFinder. Takes two command line options:
  # $1: how many threads to use for highly parallelizable steps of OrthoFinder
  # Path to the reference species' transcripts
mkdir ./03_OrthoFinder/
mkdir ./03_OrthoFinder/fasta

for file in ./02_AAFiles/*
do
  # Extract the species code:
  export speciesCode=`echo "$file" | awk -F'_' '{print $2}' | awk -F'/' '{print $2}'`
  echo This is the species code: $speciesCode
  # Copy the proteins file into the orthofinder directory:
  cp ./02_AAFiles/$speciesCode"_aa.fasta" ./03_OrthoFinder/fasta/
  # Append a species code to each gene name:
  export samplePrefix=$speciesCode"_"
  sed -i .backup "s/>/>$samplePrefix/g" "./03_OrthoFinder/fasta/"$speciesCode"_aa.fasta"
  rm ./03_OrthoFinder/fasta/*.backup
  # One of the species also has periods (.) for some reason in the gene sequences, so get rid of these:
  sed -i .backup '/^>/! s/\./-/g' "./03_OrthoFinder/fasta/"$speciesCode"_aa.fasta"
  rm ./03_OrthoFinder/fasta/*.backup
done

# Download the Bombyx mori proteins as an outgroup:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/269/925/GCF_030269925.1_ASM3026992v2/GCF_030269925.1_ASM3026992v2_protein.faa.gz -O bomo.gz
gunzip -c bomo.gz > ./03_OrthoFinder/fasta/bomo_aa.fasta

mkdir -p ./tmp
#source /programs/anaconda3/bin/activate orthofinder-3.0.1b1

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