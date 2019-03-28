#!/bin/bash 
#2019.03.28 by ZF


#Type 'sh script4_UCE_extraction.sh assemblies_folder probe_sequence', e.g. sh script3_UCE_probe_design.sh assemblies probe.fasta
#Prepare a probe set file (e.g. Phthiraptera-2.8Kv1.fasta)
#Modify all genome assemblies endding with .fa ("SPECIES_NAME.fa") and copy them to the same folder
#Tools faToTwoBit and FASconCAT, as well as PHYLUCE conda enviorment, are used and will be automatically checked prior to formal analyses in this script. MAFFT and trimal are not necessary to install because they have been included in the PHYLUCE conda enviorment.


##check the packages
echo "Checking the package dependency......" | tee -a log.txt
echo -e "\n" >> log.txt

#Check the phyluce conda environment
read -p "Please input the name of PHYLUCE conda environment ( e.g. phyluce):      " ENV
source activate $ENV
if [ $(which phyluce_probe_strip_masked_loci_from_set) ]
    then
      echo "PHYLUCE conda environment ...... OK" | tee -a log.txt
    else
      until [ $(which phyluce_probe_strip_masked_loci_from_set) ]
        do
          read -p "PHYLUCE conda environment is not found. Please input the correct conda environment name ( e.g. phyluce):      " ENV
          source activate $ENV
        done
      echo "PHYLUCE conda environment ...... OK" | tee -a log.txt
fi

#check faToTwoBit
if [ $(which faToTwoBit) ]
    then
      echo "faToTwoBit ...... OK" | tee -a log.txt
      EXE_FATOTWOBIT=$(which faToTwoBit)
      DIR_FATOTWOBIT=${EXE_FATOTWOBIT%/*}
    else
      until [ -x $DIR_FATOTWOBIT/faToTwoBit ]
        do
          read -p "faToTwoBit is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " DIR_FATOTWOBIT
        done
      echo "faToTwoBit ...... OK" | tee -a log.txt
fi

#check FASconCAT
until [ -s $DIR_FASconCAT/FASconCAT*.pl ]
    do
      read -p "DIR_FASconCAT is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/FASconCAT-G):      " DIR_FASconCAT
    done
echo "FASconCAT ...... OK"

#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done
echo "THREADS=$THREADS" >> parameters.cfg

#Check the name of the lineage
read -p "Please input the name of the lineage for the prefix of UCE probe (e.g. Collembola, Diptera):      " GROUP
echo "GROUP=$GROUP" >> parameters.cfg


#copy probe set to the working folder
echo -e "\n" | tee -a log.txt
echo "copy probe set to the working folder..." | tee -a log.txt
cat $2 > probe.fasta

#copy the assemblies to 0-genomes/
echo -e "\n" | tee -a log.txt
echo "copy the assemblies to 0-genomes/..." | tee -a log.txt
cp -r $1 0-genomes


#Change the first letter of the assembly name in the folder 0-genomes/ into upper case
cd 0-genomes/
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
echo "prepare a species list..." | tee -a log.txt
ls 0-genomes/ | sed "s/.fa//g" > species.list

SPECIES_NAME=$(cat species.list)

#Generate a python script simplify_headers.py to simplify the sequence head
echo "simplify the sequence head..." | tee -a log.txt
echo '#!/usr/bin/python' >> simplify_headers.py
echo 'from Bio import SeqIO' >> simplify_headers.py
for SPECIES in $SPECIES_NAME
  do
    echo -e '\n' >> simplify_headers.py
    echo "# $SPECIES" >> simplify_headers.py
    echo "with open(\"$SPECIES.fa\", "'"rU"'") as infile:" >> simplify_headers.py
    echo "  with open(\"$SPECIES.fasta\", \"w\") as outf:" >> simplify_headers.py
    echo "    for seq in SeqIO.parse(infile, 'fasta'):" >> simplify_headers.py
    echo "      seq.name = \"\"" >> simplify_headers.py
    echo "      seq.description = \"\"" >> simplify_headers.py
    echo "      outf.write(seq.format('fasta'))" >> simplify_headers.py
  done

#Simplify the FASTA heads
cd 0-genomes/
python ../simplify_headers.py
rm *.fa

#Move each genome assembly to its own directory
for SPECIES in *; do mkdir ${SPECIES%.*}; mv $SPECIES ${SPECIES%.*}; done

#Convert genomes to 2bit format
echo "Convert genomes to 2bit format..." | tee -a ../log.txt
for SPECIES in *; do $DIR_FATOTWOBIT/faToTwoBit $SPECIES/$SPECIES.fasta $SPECIES/${SPECIES%.*}.2bit; done

#Generate the configure file for the genome data
echo "Generate the configure file for the genome data..." | tee -a ../log.txt
echo '[scaffolds]' >> ../$GROUP-genome.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:../0-genomes/$SPECIES/$SPECIES.2bit" >> ../$GROUP-genome.conf; done

#Align bait set to the extant genome sequences
echo -e "\n" | tee -a ../log.txt
echo "Align bait set to the extant genome sequences..." | tee -a ../log.txt
cd .. && mkdir 1-uce_align && cd 1-uce_align/
phyluce_probe_run_multiple_lastzs_sqlite --db uces.sqlite --output $GROUP-genome-lastz --probefile ../probe.fasta --scaffoldlist $(cat ../species.list) --genome-base-path /home/zf/Desktop/scripts/test_UCE_extraction/test/0-genomes --identity 50 --cores $THREADS

#Extracting FASTA sequence matching UCE loci from genome sequences
echo -e "\n" | tee -a ../log.txt
echo "Extracting FASTA sequence matching UCE loci from genome sequences..." | tee -a ../log.txt
phyluce_probe_slice_sequence_from_genomes --conf ../$GROUP-genome.conf --lastz $GROUP-genome-lastz --output $GROUP-genome-fasta --flank 400 --name-pattern "probe.fasta_v_{}.lastz.clean" | tee -a ../log.txt

#prepare the configure file for species list
echo '[all]' >> ../taxon-sets.conf
cat ../species.list >> ../taxon-sets.conf

#Match contigs to baits
echo -e "\n" | tee -a ../log.txt
phyluce_assembly_match_contigs_to_probes --contigs $GROUP-genome-fasta --probes ../probe.fasta --output in-silico-lastz --min-coverage 67 | tee -a ../log.txt

#Get match counts
echo -e "\n" | tee -a ../log.txt
cd .. && mkdir 2-uces && cd 2-uces/
phyluce_assembly_get_match_counts --locus-db ../1-uce_align/in-silico-lastz/probe.matches.sqlite --taxon-list-config ../taxon-sets.conf --taxon-group 'all' --output insilico-incomplete.conf --incomplete-matrix | tee -a ../log.txt

#Change the first letter of the species name in the folder ../../$GROUP-genome-fasta into upper case
cd ../1-uce_align/$GROUP-genome-fasta/
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
cd ../../2-uces/
echo -e "\n" | tee -a ../log.txt
echo "Extract final UCEs..." | tee -a ../log.txt
phyluce_assembly_get_fastas_from_match_counts --contigs ../1-uce_align/$GROUP-genome-fasta --locus-db ../1-uce_align/in-silico-lastz/probe.matches.sqlite --match-count-output insilico-incomplete.conf --output insilico-incomplete.fasta --incomplete-matrix insilico-incomplete.incomplete | tee -a ../log.txt


#Generate loci list
cd .. && mkdir 3-raw_loci
sed -n '/uce-/p' 2-uces/insilico-incomplete.conf > loci.list

#Generate the sequence file for each locus
LOCI_NAME=$(cat loci.list)
cd 3-raw_loci/
for LOCI in $LOCI_NAME
  do
    for SPECIES in $(cat ../species.list)
      do
        echo ""$LOCI"_"$SPECIES" |"$LOCI"" >> $LOCI.seq_name
      done
  done
for LOCI in $LOCI_NAME
  do
    awk 'FNR==NR{a[$0];next} /^>/{val=$0;sub(/^>/,"",val);flag=val in a?1:0} flag' $LOCI.seq_name ../2-uces/insilico-incomplete.fasta > $LOCI.fa
  done
rm *name

#Simplify the sequence heads
for LOCI in $LOCI_NAME; do sed -i "s/"$LOCI"_\| |"$LOCI"//g" $LOCI.fa; done

#Summarize the loci number for each species
echo -e "\n" | tee -a ../log.txt
SPECIES_NAME=$(cat ../species.list)
for SPECIES in $SPECIES_NAME
do
  for LOCI in $LOCI_NAME
  do
  if [ -z $(grep -o $SPECIES $LOCI.fa) ]
    then
      echo "$LOCI in $SPECIES does not exist" | tee -a log.txt
  fi
  done
done

#Filter loci having too few taxa (less than three)
cd .. && mkdir 4-loci_filter
TOTAL_TAXA=$(cat species.list | wc -l)
for LOCI in $LOCI_NAME
do
  if [ $(grep -o $LOCI 3-raw_loci/log.txt | wc -l) -lt `expr $TOTAL_TAXA - 2` ]
    then
        cp 3-raw_loci/$LOCI.fa 4-loci_filter
        echo "$LOCI" >> 4-loci_filter/loci_name_filter.log
  fi 
done

#Align sequences for each locus
mkdir 5-loci_align
LOCI_NAME_FILTER=$(cat 4-loci_filter/loci_name_filter.log)
for LOCI_FILTER in $LOCI_NAME_FILTER
  do
    linsi --thread $THREADS 4-loci_filter/$LOCI_FILTER.fa > 5-loci_align/$LOCI_FILTER.fa
    test -s 5-loci_align/$LOCI_FILTER.fa && (echo "loci $LOCI_FILTER has been aligned" | tee -a log.txt) || mafft --thread $THREADS 4-loci_filter/$LOCI_FILTER.fa > 5-loci_align/$LOCI_FILTER.fa
  done

#Trim the alignments
mkdir 6-loci_trim
for LOCI_FILTER in $LOCI_NAME_FILTER
  do
    echo "trim nucleotide sequence of loci $LOCI_FILTER ......"
    trimal -in 5-loci_align/$LOCI_FILTER.fa -out 6-loci_trim/$LOCI_FILTER.fas -automated1
  done

#Concatenate all the nucleotide/protein alignments in phylip format and generate partition files
mkdir -p 7-matrices/all
cd 6-loci_trim/
perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
mv FcC* ../7-matrices/all
cd ../7-matrices/

#Generate the final matrices
for percent in 50 60 70 80 90 100
do
  echo "generating $percent% matrix ......"
  mkdir matrix$percent
  for LOCI in $LOCI_NAME_FILTER
  do
    num1=$(echo "scale=1;$TOTAL_TAXA*(1-$percent/100)"|bc)
    num2=$(grep -o $LOCI ../3-raw_loci/log.txt | wc -l)
    c=$(echo "$num1 >= $num2" | bc)
    if [ $c -eq 1 ]
      then
        cp ../6-loci_trim/$LOCI.fas matrix$percent
    fi 
  done
  cd matrix$percent/
  perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
  cd ..
done

#Summarize the matrices
echo "Infmormation of matrices is saved in 7-matrices/summary.extraction" | tee -a summary.extraction
echo -e '\n'
echo "nucleotide matrices" | tee -a summary.extraction
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
done
echo -e '\n'
echo "All the supermatrices and partition files are kept in 7-matrices"

source deactivate
