#!/bin/bash 

#Prepare a probe set file (Phthiraptera-2.8Kv1.fasta) in the initial working folder
#modify all genome assemblies endding with .fa ("species_names.fa") and copy them to working_folder/genomes/
#Tools MAFFT, seqkit, trimAl, faToTwoBit and FASconCAT are used in this script and the former three have been installed in environmental paths

#Define variables
DIR_fa2bit="/usr/local/bin"
DIR_TRIMAL="/home/zf/install/trimal-1.4.1"
DIR_FASconCAT="/home/zf/install/FASconCAT-G"
THREADS="8"
GROUP="Phthiraptera"
PROBE="Phthiraptera-2.8Kv1.fasta"

#Initiate the phyluce environment
source activate phyluce

##Change the first letter of the assembly name in the folder genomes/ into upper case
cd genomes/
a=`ls |xargs -n1`
echo "$a" |while read line
do
  b=`echo $line | cut -c1`
  c=`echo $line |cut -c2-`
  d=`echo $b |grep -c '[a-z]'`
  if [ $d -eq 1 ]
    then
     e=`echo $b | tr 'a-z' 'A-Z'`
     mv $b$c $e$c
  fi
done
cd ..

#prepare a species list
echo "prepare a species list..."
ls genomes/ | sed "s/'.fa'//g" > species.list

SPECIES_NAME=$(cat species.list)

#Generate a python script simplify_headers.py to simplify the sequence head
echo "simplify the sequence head..."
echo '#!/usr/bin/python' >> simplify_headers.py
echo 'from Bio import SeqIO' >> simplify_headers.py
for SPECIES in $SPECIES_NAME
  do
    echo -e '\n' >> simplify_headers.py
    echo "# $SPECIES" >> ../simplify_headers.py
    echo "with open(\"$SPECIES.fa\", "'"rU"'") as infile:" >> simplify_headers.py
    echo "  with open(\"$SPECIES.fasta\", \"w\") as outf:" >> simplify_headers.py
    echo "    for seq in SeqIO.parse(infile, 'fasta'):" >> simplify_headers.py
    echo "      seq.name = \"\"" >> simplify_headers.py
    echo "      seq.description = \"\"" >> simplify_headers.py
    echo "      outf.write(seq.format('fasta'))" >> simplify_headers.py
  done

#Simplify the FASTA heads
cd genomes/
python ../simplify_headers.py
rm *.fa

#Move each genome assembly to its own directory
for SPECIES in *; do mkdir ${SPECIES%.*}; mv $SPECIES ${SPECIES%.*}; done

#Convert genomes to 2bit format
echo "Convert genomes to 2bit format..."
for SPECIES in *; do $DIR_fa2bit/faToTwoBit $SPECIES/$SPECIES.fasta $SPECIES/${SPECIES%.*}.2bit; done

#Generate the configure file for the genome data
echo '[scaffolds]' >> ../$GROUP-genome.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:../genomes/$SPECIES/$SPECIES.2bit" >> ../$GROUP-genome.conf; done

#Align bait set to the extant genome sequences
echo "Align bait set to the extant genome sequences..."
cd .. && mkdir -p uces/$GROUP-genome-lastz && cd uces/
phyluce_probe_run_multiple_lastzs_sqlite --db uces.sqlite --output $GROUP-genome-lastz --probefile ../$PROBE --scaffoldlist $(cat ../species.list) --genome-base-path ../genomes --identity 50 --cores $THREADS

#Extracting FASTA sequence matching UCE loci from genome sequences
echo "Extracting FASTA sequence matching UCE loci from genome sequences..."
phyluce_probe_slice_sequence_from_genomes --conf ../$GROUP-genome.conf --lastz $GROUP-genome-lastz --output $GROUP-genome-fasta --flank 400 --name-pattern ""$PROBE"_v_{}.lastz.clean"

#prepare the configure file for species list
echo '[all]' >> ../taxon-sets.conf
cat ../species.list >> ../taxon-sets.conf

#Match contigs to baits
mkdir log
phyluce_assembly_match_contigs_to_probes --contigs $GROUP-genome-fasta --probes ../$PROBE --output in-silico-lastz --min-coverage 67 --log-path log

#Get match counts
mkdir -p taxon-sets/insilico-incomplete
phyluce_assembly_get_match_counts --locus-db in-silico-lastz/probe.matches.sqlite --taxon-list-config ../taxon-sets.conf --taxon-group 'all' --output taxon-sets/insilico-incomplete/insilico-incomplete.conf --log-path log --incomplete-matrix

#Change the first letter of the species name in the folder ../../$GROUP-genome-fasta into upper case
cd $GROUP-genome-fasta/
a=`ls |xargs -n1`
echo "$a" |while read line
do
  b=`echo $line | cut -c1`
  c=`echo $line |cut -c2-`
  d=`echo $b |grep -c '[a-z]'`
  if [ $d -eq 1 ]
    then
     e=`echo $b | tr 'a-z' 'A-Z'`
     mv $b$c $e$c
  fi
done

#Extract FASTA information of UCEs
echo "Extract final UCEs..."
cd ../taxon-sets/insilico-incomplete/
phyluce_assembly_get_fastas_from_match_counts --contigs ../../$GROUP-genome-fasta --locus-db ../../in-silico-lastz/probe.matches.sqlite --match-count-output insilico-incomplete.conf --output insilico-incomplete.fasta --incomplete-matrix insilico-incomplete.incomplete --log-path ../../log

source deactivate

#Generate loci list
cd .. && mkdir -p concat/raw_loci && cd concat/
sed -n '/uce-/p' ../insilico-incomplete/insilico-incomplete.conf > loci.list

#Generate the sequence file for each locus
LOCI_NAME=$(cat loci.list)
for LOCI in $LOCI_NAME; do cat ../insilico-incomplete/insilico-incomplete.fasta | seqkit grep -r -p "$LOCI"_ > raw_loci/$LOCI.fa; done

#Simplify the sequence heads
for LOCI in $LOCI_NAME; do sed -i "s/"$LOCI"_\| |"$LOCI"//g" raw_loci/$LOCI.fa; done

#Summarize the loci number for each species
SPECIES_NAME=$(cat ../../../species.list)
for SPECIES in $SPECIES_NAME
do
  for LOCI in $LOCI_NAME
  do
  if [ -z $(grep -o $SPECIES raw_loci/$LOCI.fa) ]
    then
      echo "$LOCI in $SPECIES does not exist" | tee -a raw_loci/log.txt
  fi
  done
done

#Filter loci having too few taxa (less than three)
mkdir loci_filter
TOTAL_TAXA=$(cat ../../../species.list | wc -l)
for LOCI in $LOCI_NAME
do
  if [ $(grep -o $LOCI raw_loci/log.txt | wc -l) -lt `expr $TOTAL_TAXA - 2` ]
    then
        cp raw_loci/$LOCI.fa loci_filter
        echo "$LOCI" >> loci_filter/loci_name_filter.log
  fi 
done

#Align sequences for each locus
mkdir loci_align
LOCI_NAME_FILTER=$(cat loci_filter/loci_name_filter.log)
for LOCI_FILTER in $LOCI_NAME_FILTER; do linsi loci_filter/$LOCI_FILTER.fa > loci_align/$LOCI_FILTER.fa; done

#Trim the alignments
mkdir loci_align_trim
for LOCI_FILTER in $LOCI_NAME_FILTER
  do
    echo "trim nucleotide sequence of loci $LOCI_FILTER ......"
    trimal -in loci_align/$LOCI_FILTER.fa -out loci_align_trim/$LOCI_FILTER.fas -automated1 -matrix $DIR_TRIMAL/dataset/matrix.degenerated_dna
  done

#Concatenate all the nucleotide/protein alignments in phylip format and generate partition files
mkdir -p matrices/all
cd loci_align_trim/
perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
mv FcC* ../matrices/all
cd ../matrices/

#Generate the final matrices
for percent in 50 60 70 80 90 100
do
  echo "generating $percent% matrix ......"
  mkdir matrix$percent
  for LOCI in $LOCI_NAME_FILTER
  do
    num1=$(echo "scale=1;$TOTAL_TAXA*(1-$percent/100)"|bc)
    num2=$(grep -o $LOCI ../raw_loci/log.txt | wc -l)
    c=$(echo "$num1 >= $num2" | bc)
    if [ $c -eq 1 ]
      then
        cp ../loci_align_trim/$LOCI.fas matrix$percent
    fi 
  done
  cd matrix$percent/
  perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
  cd ..
done

#Summarize the matrices
echo "Infmormation of matrices is saved in matrices/extraction.log" | tee -a extraction.log
echo -e '\n'
echo "nucleotide matrices" | tee -a extraction.log
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
done
echo -e '\n'
echo "All the supermatrices and partition files are kept in /uces/taxon-sets/concat/matrices"
