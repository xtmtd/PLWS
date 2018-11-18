#!/bin/bash

#Tools MAFFT, trimAl and FASconCAT-G are used in this script and the former two have been installed in environmental paths

#mkdir raw_busco

#Copy all the BUSCO results (run_* folders) into the same subfolder ./raw_busco

#Define variables
DIR_TRIMAL="/home/zf/trimal-1.4.1"
DIR_FASconCAT="/home/zf/FASconCAT-G"
THREADS="48"

#Generate a file containing species list
ls raw_busco/ | sed 's/run_//g' > species.list

#Generate a file of loci list
SPECIES_NAME=$(cat species.list)
for SPECIES in $SPECIES_NAME
do
cd raw_busco/run_$SPECIES/single_copy_busco_sequences
ls | sed 's/.fna\|.faa//g' >> ../../../temp.list; cd ../../..
done
sort -n temp.list | uniq > loci.list
rm temp.list

#Modify the head name of the fasta files for each locus
for SPECIES in $SPECIES_NAME
	do
	  cd raw_busco/run_$SPECIES/single_copy_busco_sequences/
	  sed -i "1c >$SPECIES" EOG*
	  cd ../../..
	done

#Merge sequences of the same locus into the fasta files
mkdir -p extraction/raw_loci/fna extraction/raw_loci/faa
LOCI_NAME=$(cat loci.list)
for LOCI in $LOCI_NAME; do touch ./extraction/raw_loci/fna/$LOCI.fna && touch ./extraction/raw_loci/faa/$LOCI.faa; done
for LOCI in $LOCI_NAME
do
    for SPECIES in $SPECIES_NAME
    do
	  if [ -f raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna ]
          then
            cat raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna >> ./extraction/raw_loci/fna/$LOCI.fna 
            cat raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.faa >> ./extraction/raw_loci/faa/$LOCI.faa
          else
            echo "$LOCI in $SPECIES does not exist" | tee -a ./extraction/raw_loci/log.txt
       fi
    done
done

#Filter loci having too few taxa (less than three)
mkdir -p extraction/loci_filter/fna extraction/loci_filter/faa
TOTAL_TAXA=$(cat species.list | wc -l)
for LOCI in $LOCI_NAME
do
    echo -e $LOCI'\t'$(grep -o $LOCI extraction/raw_loci/log.txt | wc -l) | tee -a extraction/loci_filter/sequence_number.log
    if [ $(grep -o $LOCI extraction/raw_loci/log.txt | wc -l) -lt `expr $TOTAL_TAXA - 2` ]
      then
        cp extraction/raw_loci/fna/$LOCI.fna extraction/loci_filter/fna
        cp extraction/raw_loci/faa/$LOCI.faa extraction/loci_filter/faa
        echo "$LOCI" >> extraction/loci_filter/loci_name_filter.log
    fi
done

#Align all the faa and fna files
mkdir -p extraction/loci_filter/fna_align extraction/loci_filter/faa_align
LOCI_NAME_FILTER=$(cat extraction/loci_filter/loci_name_filter.log)
for LOCI_FILTER in $LOCI_NAME_FILTER
do
  linsi --thread $THREADS extraction/loci_filter/fna/$LOCI_FILTER.fna > extraction/loci_filter/fna_align/$LOCI_FILTER.fna
  linsi --thread $THREADS extraction/loci_filter/faa/$LOCI_FILTER.faa > extraction/loci_filter/faa_align/$LOCI_FILTER.faa
done

#Trim the alignments
mkdir -p extraction/loci_trim/fna extraction/loci_trim/faa
for LOCI_FILTER in $LOCI_NAME_FILTER
do
  trimal -in extraction/loci_filter/fna_align/$LOCI_FILTER.fna -out extraction/loci_trim/fna/$LOCI_FILTER.nuc.fas -automated1 -matrix $DIR_TRIMAL/dataset/matrix.degenerated_dna
    trimal -in extraction/loci_filter/faa_align/$LOCI_FILTER.faa -out extraction/loci_trim/faa/$LOCI_FILTER.aa.fas -automated1 
done

#Concatenate all the nucleotide/protein alignments in phylip format and generate partition files
mkdir -p extraction/loci_concat/all/fna extraction/loci_concat/all/faa
cd extraction/loci_trim/fna
perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
mv FcC* ../../loci_concat/all/fna/
cd ../faa
perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
mv FcC* ../../loci_concat/all/faa/
cd ../../loci_concat/

#Generate the final matrices
for percent in 50 60 70 80 90 100
do
    mkdir -p matrix$percent/fna matrix$percent/faa
  for LOCI in $LOCI_NAME_FILTER
  do
    num1=$(echo "scale=1;$TOTAL_TAXA*(1-$percent/100)"|bc)
    num2=$(grep -o $LOCI ../raw_loci/log.txt | wc -l)
    c=$(echo "$num1 >= $num2" | bc)
    if [ $c -eq 1 ]
      then
        cp ../loci_trim/fna/$LOCI.nuc.fas matrix$percent/fna
        cp ../loci_trim/faa/$LOCI.aa.fas matrix$percent/faa
        echo $LOCI >> matrix$percent/busco$percent.loci.list
    fi 
  done
  cd matrix$percent/fna/
  perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
  cd ../faa
  perl $DIR_FASconCAT/FASconCAT-G_v1.04.pl -a -p -p -s -l
  cd ../../
done

#Summerize the matrices
echo "Infmormation of matrices is saved in /loci_cat/extraction.log" | tee -a extraction.log
echo -e '\n'
echo "protein BUSCO matrices" | tee -a extraction.log
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/faa/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/faa/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/faa/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/faa/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
done
echo -e '\n'
echo "nucleotide BUSCO matrices" | tee -a extraction.log
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/fna/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/fna/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/fna/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/fna/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a extraction.log
done
echo -e '\n'
echo "All the supermatrices and partition files are kept in extraction/loci_concat"

